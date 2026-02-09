# Script for wine-type data anlaysis (Section 7)
library(data.table)
library(ROCR)

set.seed(2026)  
your_path <- getwd()
# read the wine-type datasets
df1 <- read.csv(
  paste0(your_path,"/wine+quality/winequality-red.csv"),
  sep = ";"
)
df1$outcome <- 1 # code red wine with 1

df2 <- read.csv(
  paste0(your_path,"/wine+quality/winequality-white.csv"),
  sep = ";"
)
df2$outcome <- 0 # code white wine with 0

dat_ppi <- rbind(df1,df2) # combine the datasets

# ramdomly sample 2000 observations as the training set
idx <- sample(nrow(dat_ppi), 2000) 
dat_train <- dat_ppi[idx, ]

# the remaining observation as the evaluation set
dat_test  <- dat_ppi[-idx, ]

# fit a logistic regression model f within the training set, using all the features
fit <- glm(outcome ~ ., 
           data = dat_train,
           family = binomial())

# obtain the prefictions of f on the evaluation set
prob_test <- predict(fit, dat_test, type = "response")

# compute the AUC of the logistic regression model f on the evaluation set for sanity check
pred_obj <- prediction(prob_test, dat_test$outcome)
perf_obj <- performance(pred_obj, "auc")
auc_f <- as.numeric(perf_obj@y.values)
auc_f

# record the predictions 
dat_test$f <- prob_test

# the prespecified \pi(X) using the quality feature
dat_test$labprob <- ifelse(dat_test$quality <= 6, 0.2, 0.3)

# generate the indicator variables
dat_test$labind <- rbinom(nrow(dat_test), size = 1, prob = dat_test$labprob)
dat_lab <- dat_test[dat_test$labind == 1, ]
dat_unlab <- dat_test[dat_test$labind == 0, ]

nrow(dat_lab) # check the size of the labeled data
nrow(dat_unlab) # check the size of the unlabeled data

## tpr
# Functions for estimators and weighted influence functions computation for TPR
tpr.true.yn <- function(dat){
  return(sum((dat$R>alpha)*dat$Y)/sum(dat$Y))
}
tpr.est.yn <- function(dat){
  return(sum((dat$R>alpha)*dat$Y*dat$labind/dat$labprob)/sum(dat$Y*dat$labind/dat$labprob))
}
tpr.est.fxall <- function(dat){
  return(sum((dat$R>alpha)*dat$f)/sum(dat$f))
}
tpr.est.fxn <- function(dat){
  return(sum((dat$R>alpha)*dat$f*dat$labind/dat$labprob)/sum(dat$f*dat$labind/dat$labprob))
}
tpr.phi.yn <- function(dat){
  return((1*(dat$R>alpha)*dat$Y/mean(dat$Y) - mean((dat$R>alpha)*dat$Y)*dat$Y/(mean(dat$Y))^2)*dat$labind/dat$labprob)
}

tpr.phi.fxn <- function(dat){
  return((1*(dat$R>alpha)*dat$f/mean(dat$f) - mean((dat$R>alpha)*dat$f)*dat$f/(mean(dat$f))^2)*(dat$labind/dat$labprob-1))
}

## fpr
# Functions for estimators and weighted influence functions computation for FPR
fpr.true.yn <- function(dat){
  return(sum((dat$R>alpha)*(1-dat$Y))/sum(1-dat$Y))
}
fpr.est.yn <- function(dat){
  return(sum((dat$R>alpha)*(1-dat$Y)*dat$labind/dat$labprob)/sum((1-dat$Y)*dat$labind/dat$labprob))
}
fpr.est.fxall <- function(dat){
  return(sum((dat$R>alpha)*(1-dat$f))/sum(1-dat$f))
}
fpr.est.fxn <- function(dat){
  return(sum((dat$R>alpha)*(1-dat$f)*dat$labind/dat$labprob)/sum((1-dat$f)*dat$labind/dat$labprob))
}
fpr.phi.yn <- function(dat){
  return((1*(dat$R>alpha)*(1-dat$Y)/mean(1-dat$Y) - mean((dat$R>alpha)*(1-dat$Y))*(1-dat$Y)/(mean(1-dat$Y))^2)*dat$labind/dat$labprob)
}
fpr.phi.fxn <- function(dat){
  return((1*(dat$R>alpha)*(1-dat$f)/mean(1-dat$f) - mean((dat$R>alpha)*(1-dat$f))*(1-dat$f)/(mean(1-dat$f))^2)*(dat$labind/dat$labprob-1))
}

## AUC
# Functions for estimators and weighted influence functions computation for AUC
auc.true.yn <- function(dat){
  
  score.order <- order(dat$score, decreasing=TRUE)
  score.sorted <- dat$score[score.order]
  
  dups <- rev(duplicated(rev(score.sorted)))
  
  ## TPR
  TPR <- cumsum(dat$Y[score.order])/sum(dat$Y[score.order])
  TPR <- c(0,TPR[!dups])
  
  ## FPR
  FPR <- cumsum((1-dat$Y[score.order]))/sum((1-dat$Y[score.order]))
  FPR <- c(0,FPR[!dups])
  
  auc <- 0
  for (i in 2:length(FPR)) {
    auc <- auc + 0.5 * (FPR[i] - FPR[i-1]) * (TPR[i] + TPR[i-1])
  }
  return(auc)
}
auc.est.yn <- function(dat){
  
  score.order <- order(dat$score, decreasing=TRUE)
  score.sorted <- dat$score[score.order]
  
  dups <- rev(duplicated(rev(score.sorted)))
  
  dat$w <- dat$labind / dat$labprob # IPW weights
  ## TPR
  TPR <- cumsum(dat$Y[score.order] * dat$w[score.order])/sum(dat$Y[score.order] * dat$w[score.order])
  
  TPR <- c(0,TPR[!dups])
  
  ## FPR
  FPR <- cumsum((1-dat$Y[score.order])* dat$w[score.order])/sum((1-dat$Y[score.order])* dat$w[score.order])
  FPR <- c(0,FPR[!dups])
  
  auc <- 0
  for (i in 2:length(FPR)) {
    auc <- auc + 0.5 * (FPR[i] - FPR[i-1]) * (TPR[i] + TPR[i-1])
  }
  return(auc)
}
auc.est.fxall <- function(dat){
  
  score.order <- order(dat$score, decreasing=TRUE)
  score.sorted <- dat$score[score.order]
  
  dups <- rev(duplicated(rev(score.sorted)))
  
  ## TPR
  TPR <- cumsum(dat$f[score.order])/sum(dat$f[score.order])
  TPR <- c(0,TPR[!dups])
  
  ## FPR
  FPR <- cumsum((1-dat$f[score.order]))/sum((1-dat$f[score.order]))
  FPR <- c(0,FPR[!dups])
  
  auc <- 0
  for (i in 2:length(FPR)) {
    auc <- auc + 0.5 * (FPR[i] - FPR[i-1]) * (TPR[i] + TPR[i-1])
  }
  return(auc)
}

auc.est.fxn <- function(dat){
  
  score.order <- order(dat$score, decreasing=TRUE)
  score.sorted <- dat$score[score.order]
  
  dups <- rev(duplicated(rev(score.sorted)))
  
  dat$w <- dat$labind / dat$labprob # IPW weights
  ## TPR
  TPR <- cumsum(dat$f[score.order] * dat$w[score.order])/sum(dat$f[score.order] * dat$w[score.order])
  
  TPR <- c(0,TPR[!dups])
  
  ## FPR
  FPR <- cumsum((1-dat$f[score.order])* dat$w[score.order])/sum((1-dat$f[score.order])* dat$w[score.order])
  FPR <- c(0,FPR[!dups])
  
  auc <- 0
  for (i in 2:length(FPR)) {
    auc <- auc + 0.5 * (FPR[i] - FPR[i-1]) * (TPR[i] + TPR[i-1])
  }
  return(auc)
}

## The influence functions computation for AUC mostly follows from the cvAUC package 
## with some small adjustments to accommodate our setups
auc.phi.yn <- function(dat){
  library(cvAUC)
  predictions <- dat$score  # Length-V list of predicted values
  labels <- dat$Y  # Length-V list of true labels
  n_obs <- nrow(dat)  # Number of observations
  
  # Inverse probability weights across entire data set
  w1 <- 1/(sum(labels==1)/n_obs)  # Inverse weights for positive class 1/E[Y]
  w0 <- 1/(sum(labels==0)/n_obs)  # Inverse weights for negative class 1/E[1-Y]
  
  # This is required to cleanly get past R CMD CHECK
  # https://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check
  pred <- label <- NULL
  fracNegLabelsWithSmallerPreds <- fracPosLabelsWithLargerPreds <- icVal <- NULL  
  
  .IC <- function(fold_preds, fold_labels, pos, neg, w1, w0) { 
    # Applied to a single fold's (preds, labels)
    n_rows <- length(fold_labels)
    n_pos <- sum(fold_labels == pos)
    n_neg <- n_rows - n_pos
    auc <- AUC(fold_preds, fold_labels)
    DT <- data.table(pred = fold_preds, label = fold_labels)
    DT <- DT[order(pred, -xtfrm(label))]  #Sort by asc(pred), desc(label)
    DT[, fracNegLabelsWithSmallerPreds := cumsum(label == neg)/n_neg] #fracNegLabelSWithSmallerPreds: E[I(R>\alpha)(1-Y)]/E[1-Y] at all R_i
    DT <- DT[order(-pred, label)] 
    DT[, fracPosLabelsWithLargerPreds := cumsum(label == pos)/n_pos] #fracPosLabelsWithLargerPreds: (1-E[I(R<=\alpha)Y]/E[Y]) at all R_i
    DT[, icVal := ifelse(label == pos, w1 * (fracNegLabelsWithSmallerPreds - auc),
                         w0 * (fracPosLabelsWithLargerPreds - auc))]
    #DT[, icVal :=  label*w1 * (fracNegLabelsWithSmallerPreds - auc)+(1-label)*w0 * (fracPosLabelsWithLargerPreds - auc)]
    #return(mean(DT$icVal))
    return(DT$icVal)
  }
  return(.IC(fold_preds=predictions, fold_labels=labels, pos=1, neg=0, w1=w1, w0=w0)*(dat$labind / dat$labprob)) # weighted
}

auc.phi.fxn <- function(dat){
  predictions <- dat$score  # Length-V list of predicted values
  labels <- dat$f  # Length-V list of f
  n_obs <- nrow(dat)  # Number of observations
  
  # Inverse probability weights across entire data set
  w1 <- 1/(mean(labels))  # Inverse weights for f
  w0 <- 1/(mean(1-labels))  # Inverse weights for 1-f
  
  # This is required to cleanly get past R CMD CHECK
  # https://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check
  pred <- label <- NULL
  fracNegLabelsWithSmallerPreds <- fracPosLabelsWithLargerPreds <- icVal <- NULL  
  
  .IC <- function(fold_preds, fold_labels, w1, w0) { 
    # Applied to a single fold's (preds, labels)
    
    n_pos <- sum(fold_labels)
    n_neg <- sum(1-fold_labels)
    auc_f <- auc.est.fxn(dat)
    DT <- data.table(pred = fold_preds, label = fold_labels)
    DT <- DT[order(pred, -xtfrm(label))]  #Sort by asc(pred), desc(label)
    DT[, fracNegLabelsWithSmallerPreds := cumsum(1-label)/n_neg]
    DT <- DT[order(-pred, label)] 
    DT[, fracPosLabelsWithLargerPreds := cumsum(label)/n_pos]
    DT[, icVal :=  label*w1 * (fracNegLabelsWithSmallerPreds - auc_f)+(1-label)*w0 * (fracPosLabelsWithLargerPreds - auc_f)]
    #return(mean(DT$icVal))
    return(DT$icVal)
  }
  return(.IC(fold_preds=predictions, fold_labels=labels, w1=w1, w0=w0)*(dat$labind / dat$labprob-1)) # weighted
}


inf_dat <- dat_test

inf_dat$Y <- inf_dat$outcome

# Use density as the biomarker to be evaluated
inf_dat$score <- inf_dat$density
inf_dat$R <- inf_dat$density

alpha <- 0.998
size_all <- nrow(inf_dat)

dat.i <- inf_dat

tpr_phi_n.i <- tpr.phi.yn(dat.i)
fpr_phi_n.i <- fpr.phi.yn(dat.i)
auc_phi_n.i <- auc.phi.yn(dat.i)

# compute the estimators that do not utilize the predictions of the unlabeled data
tpr_est_n.i <- tpr.est.yn(dat.i)
fpr_est_n.i <- fpr.est.yn(dat.i)
auc_est_n.i <- auc.est.yn(dat.i)

# compute the se estimates of the estimators that do not utilize the predictions of the unlabeled data
tpr_sd_n.i <- sqrt(var(tpr_phi_n.i))/sqrt(size_all)
fpr_sd_n.i <- sqrt(var(fpr_phi_n.i))/sqrt(size_all)
auc_sd_n.i <- sqrt(var(auc_phi_n.i))/sqrt(size_all)

tpr_phi_f_all.i <- tpr.phi.fxn(dat.i)
fpr_phi_f_all.i <- fpr.phi.fxn(dat.i)
auc_phi_f_all.i <- auc.phi.fxn(dat.i)

# compute \omega
tpr_omega.i <- cov(tpr_phi_n.i,tpr_phi_f_all.i)/(var(tpr_phi_f_all.i))
fpr_omega.i <- cov(fpr_phi_n.i,fpr_phi_f_all.i)/(var(fpr_phi_f_all.i))
auc_omega.i <- cov(auc_phi_n.i,auc_phi_f_all.i)/(var(auc_phi_f_all.i))

tpr_est_f_all.i <- tpr.est.fxall(dat.i)
tpr_est_f_n.i <- tpr.est.fxn(dat.i)
fpr_est_f_all.i <- fpr.est.fxall(dat.i)
fpr_est_f_n.i <- fpr.est.fxn(dat.i)
auc_est_f_all.i <- auc.est.fxall(dat.i)
auc_est_f_n.i <- auc.est.fxn(dat.i)

# construct the PPI estimators
tpr_est_ppi.i <- tpr_est_n.i + tpr_omega.i * (tpr_est_f_all.i - tpr_est_f_n.i)
fpr_est_ppi.i <- fpr_est_n.i + fpr_omega.i * (fpr_est_f_all.i - fpr_est_f_n.i)
auc_est_ppi.i <- auc_est_n.i + auc_omega.i * (auc_est_f_all.i - auc_est_f_n.i)

# compute the se estimates of the PPI estimators
tpr_sd_ppi.i <- sqrt((var(tpr_phi_n.i)-cov(tpr_phi_f_all.i,tpr_phi_n.i)^2/(var(tpr_phi_f_all.i)))/(size_all))
fpr_sd_ppi.i <- sqrt((var(fpr_phi_n.i)-cov(fpr_phi_f_all.i,fpr_phi_n.i)^2/(var(fpr_phi_f_all.i)))/(size_all))
auc_sd_ppi.i <- sqrt((var(auc_phi_n.i)-cov(auc_phi_f_all.i,auc_phi_n.i)^2/(var(auc_phi_f_all.i)))/(size_all))



# computing the ideal tpr, fpr, and AUC estimates when all test data labeled 
all_dat <- inf_dat
all_dat$labind <- 1 # all the Y are observed
all_dat$labprob <- 1 # no unlabeled data
tpr_phi_all.i <- tpr.phi.yn(all_dat)
fpr_phi_all.i <- fpr.phi.yn(all_dat)
auc_phi_all.i <- auc.phi.yn(all_dat)

# record the ideal tpr, fpr, and AUC estimates when all test data labeled 
tpr_est_all.i <- tpr.est.yn(all_dat)
fpr_est_all.i <- fpr.est.yn(all_dat)
auc_est_all.i <- auc.est.yn(all_dat)

# compute the se estimates for the ideal tpr, fpr, and AUC estimates when all test data labeled 
tpr_sd_all.i <- sqrt(var(tpr_phi_all.i))/sqrt(size_all)
fpr_sd_all.i <- sqrt(var(fpr_phi_all.i))/sqrt(size_all)
auc_sd_all.i <- sqrt(var(auc_phi_all.i))/sqrt(size_all)

# record the results
res.i <- list(
  tpr_est_n_i=tpr_est_n.i,
  tpr_est_ppi_i=tpr_est_ppi.i,
  tpr_est_all_i=tpr_est_all.i,
  
  tpr_sd_n_i=tpr_sd_n.i,
  tpr_sd_ppi_i=tpr_sd_ppi.i,
  tpr_sd_all_i=tpr_sd_all.i,
  
  fpr_est_n_i=fpr_est_n.i,
  fpr_est_ppi_i=fpr_est_ppi.i,
  fpr_est_all_i=fpr_est_all.i,
  
  fpr_sd_n_i=fpr_sd_n.i,
  fpr_sd_ppi_i=fpr_sd_ppi.i,
  fpr_sd_all_i=fpr_sd_all.i,
  
  auc_est_n_i=auc_est_n.i,
  auc_est_ppi_i=auc_est_ppi.i,
  auc_est_all_i=auc_est_all.i,
  
  
  auc_sd_n_i=auc_sd_n.i,
  auc_sd_ppi_i=auc_sd_ppi.i,
  auc_sd_all_i=auc_sd_all.i
)


# Create the table
tbl <- data.frame(
  Metric = c("TPR", "FPR", "AUC"),
  
  `Point (i)`   = c(res.i$tpr_est_n_i,
                    res.i$fpr_est_n_i,
                    res.i$auc_est_n_i),
  
  `Point (ii)`  = c(res.i$tpr_est_ppi_i,
                    res.i$fpr_est_ppi_i,
                    res.i$auc_est_ppi_i),
  
  `Point (iii)` = c(res.i$tpr_est_all_i,
                    res.i$fpr_est_all_i,
                    res.i$auc_est_all_i),
  
  `SE (i)`   = c(res.i$tpr_sd_n_i,
                 res.i$fpr_sd_n_i,
                 res.i$auc_sd_n_i),
  
  `SE (ii)`  = c(res.i$tpr_sd_ppi_i,
                 res.i$fpr_sd_ppi_i,
                 res.i$auc_sd_ppi_i),
  
  `SE (iii)` = c(res.i$tpr_sd_all_i,
                 res.i$fpr_sd_all_i,
                 res.i$auc_sd_all_i)
)


tbl_fmt <- tbl

point_cols <- grep("^Point", names(tbl_fmt), value = TRUE)
se_cols    <- grep("^SE", names(tbl_fmt), value = TRUE)

tbl_fmt[point_cols] <- lapply(tbl_fmt[point_cols], function(x)
  sprintf("%.4f", x)
)

tbl_fmt[se_cols] <- lapply(tbl_fmt[se_cols], function(x)
  sprintf("%.5f", x)
)

library(knitr)

caption_text <- paste(
  "Estimates of TPR, FPR, and AUC with estimated standard errors.",
  "(i) Estimator that does not use predictions in the unlabeled data;",
  "(ii) the PPI estimator;",
  "(iii) the ideal estimator using labels for the unlabeled data."
)

kable(tbl_fmt, caption = caption_text)

