# simulation code for Section 6.2 
# PPI estimators under covariate distribution shift in the context of Section 3.1

expit <- function(a){ return(exp(a)/(1+exp(a)))}

## mean
# compute the estimators and the weighted influence functions
mean.true.yn <- function(dat){
  return(mean(dat$Y))
}
mean.est.yn <- function(dat){
  return(sum(dat$Y*dat$labind/dat$labprob)/sum(dat$labind/dat$labprob))
}
mean.est.fxall <- function(dat){
  return(mean(dat$f))
}
mean.est.fxn <- function(dat){
  return(sum(dat$f*dat$labind/dat$labprob)/sum(dat$labind/dat$labprob))
}
mean.phi.yn <- function(dat){
  observed <- dat$labind == 1
  ipw <- dat$labind/dat$labprob
  mean.y <- sum(dat$Y[observed]*ipw[observed])/sum(ipw[observed])
  phi <- numeric(nrow(dat))
  phi[observed] <- (dat$Y[observed] - mean.y)*ipw[observed]
  return(phi)
}
mean.phi.fxn <- function(dat){
  return((dat$f-mean(dat$f))*(dat$labind/dat$labprob-1))
}


## tpr
# compute the estimators and the weighted influence functions
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
  observed <- dat$labind == 1
  ipw <- dat$labind/dat$labprob
  mean.y <- sum(dat$Y[observed]*ipw[observed])/nrow(dat)
  mean.ry <- sum((dat$R[observed]>alpha)*dat$Y[observed]*ipw[observed])/nrow(dat)
  phi <- numeric(nrow(dat))
  phi[observed] <- ((dat$R[observed]>alpha)*dat$Y[observed]/mean.y -
                      mean.ry*dat$Y[observed]/(mean.y^2))*ipw[observed]
  return(phi)
}
tpr.phi.fxn <- function(dat){
  return((1*(dat$R>alpha)*dat$f/mean(dat$f) - mean((dat$R>alpha)*dat$f)*dat$f/(mean(dat$f))^2)*(dat$labind/dat$labprob-1))
}

## fpr
# compute the estimators and the weighted influence functions
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
  observed <- dat$labind == 1
  ipw <- dat$labind/dat$labprob
  mean.one.minus.y <- sum((1-dat$Y[observed])*ipw[observed])/nrow(dat)
  mean.r.one.minus.y <- sum((dat$R[observed]>alpha)*(1-dat$Y[observed])*ipw[observed])/nrow(dat)
  phi <- numeric(nrow(dat))
  phi[observed] <- ((dat$R[observed]>alpha)*(1-dat$Y[observed])/mean.one.minus.y -
                      mean.r.one.minus.y*(1-dat$Y[observed])/(mean.one.minus.y^2))*ipw[observed]
  return(phi)
}
fpr.phi.fxn <- function(dat){
  return((1*(dat$R>alpha)*(1-dat$f)/mean(1-dat$f) - mean((dat$R>alpha)*(1-dat$f))*(1-dat$f)/(mean(1-dat$f))^2)*(dat$labind/dat$labprob-1))
}

## AUC
# compute the estimators and the weighted influence functions
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
  observed <- dat$labind == 1
  predictions <- dat$score[observed]
  labels <- dat$Y[observed]
  case.weights <- (dat$labind/dat$labprob)[observed]
  
  # Inverse weighted probabilities for the positive and negative classes
  w1 <- 1 / mean(dat$labind * dat$Y / dat$labprob)
  w0 <- 1 / mean(dat$labind * (1-dat$Y) / dat$labprob)
  auc_y <- auc.est.yn(dat)
  
  pred <- label <- obs.weight <- row.id <- NULL
  fracNegLabelsWithSmallerPreds <- fracPosLabelsWithLargerPreds <- icVal <- NULL
  
  .IC <- function(fold_preds, fold_labels, fold_weights, pos, neg, w1, w0, auc) {
    n_pos <- sum(fold_weights*(fold_labels == pos))
    n_neg <- sum(fold_weights*(fold_labels == neg))
    DT <- data.table(row.id = seq_along(fold_labels), pred = fold_preds,
                     label = fold_labels, obs.weight = fold_weights)
    DT <- DT[order(pred, -xtfrm(label))]
    DT[, fracNegLabelsWithSmallerPreds :=
         cumsum(obs.weight*(label == neg))/n_neg]
    DT <- DT[order(-pred, label)]
    DT[, fracPosLabelsWithLargerPreds :=
         cumsum(obs.weight*(label == pos))/n_pos]
    DT[, icVal := ifelse(label == pos,
                         w1*(fracNegLabelsWithSmallerPreds - auc),
                         w0*(fracPosLabelsWithLargerPreds - auc))]
    DT <- DT[order(row.id)]
    return(DT$icVal)
  }
  
  phi <- numeric(nrow(dat))
  phi[observed] <- .IC(fold_preds = predictions, fold_labels = labels,
                       fold_weights = case.weights, pos = 1, neg = 0,
                       w1 = w1, w0 = w0, auc = auc_y)*case.weights
  return(phi)
}

auc.phi.fxn <- function(dat){
  predictions <- dat$score
  labels <- dat$f
  n_obs <- nrow(dat)
  
  # Inverse probabilities for the soft positive and negative classes
  w1 <- 1/(mean(labels))
  w0 <- 1/(mean(1-labels))
  auc_f <- auc.est.fxall(dat)
  
  pred <- label <- row.id <- NULL
  fracNegLabelsWithSmallerPreds <- fracPosLabelsWithLargerPreds <- icVal <- NULL
  
  .IC <- function(fold_preds, fold_labels, w1, w0, auc) {
    n_pos <- sum(fold_labels)
    n_neg <- sum(1-fold_labels)
    DT <- data.table(row.id = seq_along(fold_labels), pred = fold_preds,
                     label = fold_labels)
    DT <- DT[order(pred, -xtfrm(label))]
    DT[, fracNegLabelsWithSmallerPreds := cumsum(1-label)/n_neg]
    DT <- DT[order(-pred, label)]
    DT[, fracPosLabelsWithLargerPreds := cumsum(label)/n_pos]
    DT[, icVal := label*w1*(fracNegLabelsWithSmallerPreds - auc) +
         (1-label)*w0*(fracPosLabelsWithLargerPreds - auc)]
    DT <- DT[order(row.id)]
    return(DT$icVal)
  }
  
  return(.IC(fold_preds = predictions, fold_labels = labels,
             w1 = w1, w0 = w0, auc = auc_f)*
           (dat$labind/dat$labprob-1))
}


library(data.table)
library(ROCR)
library(MASS)

f_types <- c("all","noise","bad","ideal","pure_noise")# different prediction models f
# "all" is a random forest model predicting Y from (R,X_1,...,X_5)
# "bad" is a random forest model predicting Y from (R,X_1,X_4,X_4)
# "noise" is a random forest model predicting Y from (X_4,X_5), which are independent of Y
# "ideal" represents E[Y|R,X_1,...,X_5]
# "pure_noise" is Unif[0.01,0.99], uniformly distributed random noise

alpha <- 0.6

results_dir <- "results/covshift_ipw"
set.seed(2025)

# The covariance matrix of (R,X_1,...,X_5)
Var_mat <- matrix(c(
  1.0, 0.5, 0.3, 0.2, 0, 0,
  0.5, 1.0, 0.4, 0.3, 0, 0,
  0.3, 0.4, 1.0, 0.2, 0, 0,
  0.2, 0.3, 0.2, 1.0, 0, 0,
  0,   0,    0,   0,  1, 0,
  0,  0,    0,   0,   0, 2
), nrow=6, byrow=TRUE)

# prespecified E[Y|R,X_1,...,X_5]
beta_intcpt <- -0.5
beta_R  <-  1.0
beta_X1q <- -0.9   # X_1^2
beta_X2a <-  0.6   # |X_2|
beta_X3c <-  0.5   # X_3^3
beta_Rx3 <-  1.5   # R * X_3
beta_X1x2 <- -0.7  # X_1 * X_2

compute_eta_p <- function(df){
  eta <- with(df,
              beta_intcpt +
                beta_R  * R +
                beta_X1q * X1^2 +
                beta_X2a * abs(X2) +
                beta_X3c * X3^3 +
                beta_Rx3 * (R * X3) +
                beta_X1x2 * (X1 * X2)
  )
  
  # transform to probability
  plogis(eta)
}


m <- 100000 # size of the training data for estimating the random forest prediction models

# generate (R,X_1,...,X_5) for the training data
train_dat <- as.data.frame(mvrnorm(m, mu = c(0,0,0,0,3,20), Sigma = Var_mat)) 

names(train_dat) <- c("R","X1","X2","X3","X4","X5")

# compute E[Y|R,X_1,...,X_5]
train_dat$ideal <- compute_eta_p(train_dat)

# generate Y for the training data, following the prespecified model E[Y|R,X_1,...,X_5]
train_dat$Y <- rbinom(m, size=1, prob=train_dat$ideal)
train_dat$Y <- factor(train_dat$Y, levels=0:1)  # for classification

# prespecified function for \pi(X)
lp <- with(train_dat, 
           R 
           - 0.9 * X1 
           + 0.7 * X2 * X3 - 0.5
)

# raw \pi(X)
p_raw <- expit(lp)
eps <- 0.2 # set \epsilon s.t. \pi(X) \in [\epsilon, 1−\epsilon]

# rescaled \pi(X) (\pi(X) \in [\epsilon, 1−\epsilon])
train_dat$labprob <- 0.2 + (1 - 2*0.2) * p_raw
train_dat$labind <- rbinom(m,1,train_dat$labprob) # generate the indicator variable based on \pi(X)

library(randomForest)

# A random forest model predicting Y from (R,X_1,...,X_5)
rf_all <- randomForest(Y ~ R + X1 + X2 + X3 + X4 + X5,
                       data = train_dat)

# A random forest model predicting Y from (X_4,X_5), which are independent of Y
rf_noise <- randomForest(Y ~ X4 + X5,
                         data = train_dat)

# 3) A random forest model predicting Y from (R,X_1,X_4,X_4)
rf_bad <- randomForest(Y ~ R + X1 + X4 + X5,
                       data = train_dat)

# data generation function
dat_gen_full <- function(size,f_type){
  
  # generate the data in the context of Section 3.1
  inf_dat <- as.data.frame(mvrnorm((size), mu=c(0,0,0,0,3,20), Sigma=Var_mat))
  names(inf_dat) <- c("R","X1","X2","X3","X4","X5")
  
  inf_dat$ideal <- compute_eta_p(inf_dat) # the prespecified ideal f: E[Y|R,X_1,...,X_5]
  
  inf_dat$Y <- rbinom(size, size=1, prob=inf_dat$ideal) # generating Y following the prespecified model
  
  # obtain predictions using the "all" random forest model
  inf_dat$all <- predict(rf_all, newdata = inf_dat,type = "prob")[ ,"1"]
  # obtain predictions using the "noise" random forest model
  inf_dat$noise <- predict(rf_noise, newdata = inf_dat,type = "prob")[ ,"1"]
  # obtain predictions using the "pure_noise" f
  inf_dat$pure_noise <- runif(size,min=0.01, max=0.99)
  # obtain predictions using the "bad" random forest model
  inf_dat$bad <- predict(rf_bad, newdata = inf_dat,type = "prob")[ ,"1"]
  
  # the f used for constructing PPI estimators
  inf_dat$f <- inf_dat[,f_type]
  
  inf_dat$score <- inf_dat$R # this is for facilitating the AUC computation
  
  # prespecified function for \pi(X)
  lp <- with(inf_dat, 
             R 
             - 0.9 * X1 
             + 0.7 * X2 * X3 - 0.5
  )
  
  # raw \pi(X)
  p_raw <- expit(lp)
  
  # rescaled \pi(X) (\pi(X) \in [\epsilon, 1−\epsilon])
  inf_dat$labprobtrue <- 0.2 + (1 - 2*0.2) * p_raw
  inf_dat$labprob <- inf_dat$labprobtrue
  inf_dat$labind <- rbinom(size,1,inf_dat$labprobtrue)
  
  return(inf_dat)
}


# simulation function
ppi.sim <- function(nrep,dat_gen,size_all){
  
  ## 10000 replications using fully labeled data to compute the targets
  target_rep <- 10000
  mean_target_0_vec <- tpr_target_0_vec <- fpr_target_0_vec <- auc_target_0_vec <- vector("numeric",length = target_rep)
  for(index in 1:target_rep){
    
    set.seed(index+13145)
    
    inf_dat.index <- as.data.frame(mvrnorm(size_all, mu=c(0,0,0,0,3,20), Sigma=Var_mat)) # generating (R,X_1,...,X_5)
    names(inf_dat.index) <- c("R","X1","X2","X3","X4","X5")
    
    inf_dat.index$ideal <- compute_eta_p(inf_dat.index)
    inf_dat.index$Y <- rbinom(size_all, size=1, prob=inf_dat.index$ideal) # generating Y
    inf_dat.index$score <- inf_dat.index$R # facilitating the computation of AUC
    
    mean_target_0_vec[index] <- mean.true.yn(inf_dat.index)
    tpr_target_0_vec[index] <- tpr.true.yn(inf_dat.index)
    fpr_target_0_vec[index] <- fpr.true.yn(inf_dat.index)
    auc_target_0_vec[index] <- auc.true.yn(inf_dat.index)
  }
  # compute the mean of the 10000 replications as the targets
  mean_target_0 <- mean(mean_target_0_vec)
  tpr_target_0 <- mean(tpr_target_0_vec)
  fpr_target_0 <- mean(fpr_target_0_vec)
  auc_target_0 <- mean(auc_target_0_vec)
  
  
  for(k in 1:length(f_types)){
    # for saving the estimators that do not utilize the predictions of the unlabeled data
    mean_est_n_vec <- tpr_est_n_vec <- fpr_est_n_vec <- auc_est_n_vec <- vector("numeric",length = nrep)
    # for saving the PPI estimators
    mean_est_ppi_vec <- tpr_est_ppi_vec <- fpr_est_ppi_vec <- auc_est_ppi_vec <- vector("numeric",length = nrep)
    # for saving the \omega
    mean_omega_vec <- tpr_omega_vec <- fpr_omega_vec <- auc_omega_vec  <- vector("numeric",length = nrep)
    # for saving the se estimates of the estimators that do not utilize the predictions of the unlabeled data
    mean_sd_n_vec <- tpr_sd_n_vec <- fpr_sd_n_vec <- auc_sd_n_vec <- vector("numeric",length=nrep)
    # for saving the se estimates of the PPI estimators
    mean_sd_ppi_vec <- tpr_sd_ppi_vec <- fpr_sd_ppi_vec <- auc_sd_ppi_vec <-  vector("numeric",length=nrep)
    # for saving the coverage of the estimators that do not utilize the predictions of the unlabeled data
    mean_cov_n_vec <- tpr_cov_n_vec <- fpr_cov_n_vec <- auc_cov_n_vec <- vector("numeric",length=nrep)
    # for saving the coverage of the PPI estimators
    mean_cov_ppi_vec <- tpr_cov_ppi_vec <- fpr_cov_ppi_vec <- auc_cov_ppi_vec <-  vector("numeric",length=nrep)
    
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    
    current_path <- file.path(
      results_dir,
      sprintf("f_%s_size_%d", f_types[k], size_all)
    )
    
    dir.create(
      current_path,
      recursive = TRUE,
      showWarnings = FALSE
    )

    
    # simulations -------------------------------------------------------------
    
    for( i in 1:nrep){
      
      set.seed(i)
      dat.i <- dat_gen(size_all,f_types[k])
      
      mean_phi_n.i <- mean.phi.yn(dat.i)
      tpr_phi_n.i <- tpr.phi.yn(dat.i)
      fpr_phi_n.i <- fpr.phi.yn(dat.i)
      auc_phi_n.i <- auc.phi.yn(dat.i)
      
      mean_est_n.i <- mean.est.yn(dat.i)
      tpr_est_n.i <- tpr.est.yn(dat.i)
      fpr_est_n.i <- fpr.est.yn(dat.i)
      auc_est_n.i <- auc.est.yn(dat.i)
      
      # record the estimators that do not utilize the predictions of the unlabeled data
      mean_est_n_vec[i] <- mean_est_n.i
      tpr_est_n_vec[i] <- tpr_est_n.i
      fpr_est_n_vec[i] <- fpr_est_n.i
      auc_est_n_vec[i] <- auc_est_n.i
      
      # record th se estimates of the estimators that do not utilize the predictions of the unlabeled data
      mean_sd_n_vec[i] <- sqrt(var(mean_phi_n.i))/sqrt(size_all)
      tpr_sd_n_vec[i] <- sqrt(var(tpr_phi_n.i))/sqrt(size_all)
      fpr_sd_n_vec[i] <- sqrt(var(fpr_phi_n.i))/sqrt(size_all)
      auc_sd_n_vec[i] <- sqrt(var(auc_phi_n.i))/sqrt(size_all)
      
      # record the coverage of the estimators that do not utilize the predictions of the unlabeled data
      if(mean_est_n.i+1.96*mean_sd_n_vec[i] >= mean_target_0 & mean_est_n.i-1.96*mean_sd_n_vec[i] <= mean_target_0){
        mean_cov_n_vec[i] <- 1
      }else{
        mean_cov_n_vec[i] <- 0
      }
      if(tpr_est_n.i+1.96*tpr_sd_n_vec[i] >= tpr_target_0 & tpr_est_n.i-1.96*tpr_sd_n_vec[i] <= tpr_target_0){
        tpr_cov_n_vec[i] <- 1
      }else{
        tpr_cov_n_vec[i] <- 0
      }
      if(fpr_est_n.i+1.96*fpr_sd_n_vec[i] >= fpr_target_0 & fpr_est_n.i-1.96*fpr_sd_n_vec[i] <= fpr_target_0){
        fpr_cov_n_vec[i] <- 1
      }else{
        fpr_cov_n_vec[i] <- 0
      }
      if(auc_est_n.i+1.96*auc_sd_n_vec[i] >= auc_target_0 & auc_est_n.i-1.96*auc_sd_n_vec[i] <= auc_target_0){
        auc_cov_n_vec[i] <- 1
      }else{
        auc_cov_n_vec[i] <- 0
      }
      
      mean_phi_f_all.i <- mean.phi.fxn(dat.i)
      tpr_phi_f_all.i <- tpr.phi.fxn(dat.i)
      fpr_phi_f_all.i <- fpr.phi.fxn(dat.i)
      auc_phi_f_all.i <- auc.phi.fxn(dat.i)
      
      # compute \omega
      mean_omega_vec[i] <- cov(mean_phi_n.i,mean_phi_f_all.i)/(var(mean_phi_f_all.i))
      tpr_omega_vec[i] <- cov(tpr_phi_n.i,tpr_phi_f_all.i)/(var(tpr_phi_f_all.i))
      fpr_omega_vec[i] <- cov(fpr_phi_n.i,fpr_phi_f_all.i)/(var(fpr_phi_f_all.i))
      auc_omega_vec[i] <- cov(auc_phi_n.i,auc_phi_f_all.i)/(var(auc_phi_f_all.i))
      
      mean_est_f_all.i <- mean.est.fxall(dat.i)
      mean_est_f_n.i <- mean.est.fxn(dat.i)
      tpr_est_f_all.i <- tpr.est.fxall(dat.i)
      tpr_est_f_n.i <- tpr.est.fxn(dat.i)
      fpr_est_f_all.i <- fpr.est.fxall(dat.i)
      fpr_est_f_n.i <- fpr.est.fxn(dat.i)
      auc_est_f_all.i <- auc.est.fxall(dat.i)
      auc_est_f_n.i <- auc.est.fxn(dat.i)
      
      # construct the PPI estimator
      mean_est_ppi_vec[i] <- mean_est_n.i + mean_omega_vec[i]*(mean_est_f_all.i - mean_est_f_n.i)
      tpr_est_ppi_vec[i] <- tpr_est_n.i + tpr_omega_vec[i]*(tpr_est_f_all.i - tpr_est_f_n.i)
      fpr_est_ppi_vec[i] <- fpr_est_n.i + fpr_omega_vec[i]*(fpr_est_f_all.i - fpr_est_f_n.i)
      auc_est_ppi_vec[i] <- auc_est_n.i + auc_omega_vec[i]*(auc_est_f_all.i - auc_est_f_n.i)
      
      # record the se estimates of the PPI estimators
      mean_sd_ppi_vec[i] <- sqrt((var(mean_phi_n.i)-cov(mean_phi_f_all.i,mean_phi_n.i)^2/(var(mean_phi_f_all.i)))/(size_all))
      tpr_sd_ppi_vec[i] <- sqrt((var(tpr_phi_n.i)-cov(tpr_phi_f_all.i,tpr_phi_n.i)^2/(var(tpr_phi_f_all.i)))/(size_all))
      fpr_sd_ppi_vec[i] <- sqrt((var(fpr_phi_n.i)-cov(fpr_phi_f_all.i,fpr_phi_n.i)^2/(var(fpr_phi_f_all.i)))/(size_all))
      auc_sd_ppi_vec[i] <- sqrt((var(auc_phi_n.i)-cov(auc_phi_f_all.i,auc_phi_n.i)^2/(var(auc_phi_f_all.i)))/(size_all))
      
      if(mean_sd_ppi_vec[i]=='NaN'|tpr_sd_ppi_vec[i]=='NaN'|fpr_sd_ppi_vec[i]=='NaN'|auc_sd_ppi_vec[i]=='NaN'){
        print(i)
        stop()
      }
      
      # record the coverage of the PPI estimators
      if(mean_est_ppi_vec[i]+1.96*mean_sd_ppi_vec[i] >= mean_target_0 & mean_est_ppi_vec[i]-1.96*mean_sd_ppi_vec[i] <= mean_target_0){
        mean_cov_ppi_vec[i] <- 1
      }else{
        mean_cov_ppi_vec[i] <- 0
      }
      if(tpr_est_ppi_vec[i]+1.96*tpr_sd_ppi_vec[i] >= tpr_target_0 & tpr_est_ppi_vec[i]-1.96*tpr_sd_ppi_vec[i] <= tpr_target_0){
        tpr_cov_ppi_vec[i] <- 1
      }else{
        tpr_cov_ppi_vec[i] <- 0
      }
      if(fpr_est_ppi_vec[i]+1.96*fpr_sd_ppi_vec[i] >= fpr_target_0 & fpr_est_ppi_vec[i]-1.96*fpr_sd_ppi_vec[i] <= fpr_target_0){
        fpr_cov_ppi_vec[i] <- 1
      }else{
        fpr_cov_ppi_vec[i] <- 0
      }
      if(auc_est_ppi_vec[i]+1.96*auc_sd_ppi_vec[i] >= auc_target_0 & auc_est_ppi_vec[i]-1.96*auc_sd_ppi_vec[i] <= auc_target_0){
        auc_cov_ppi_vec[i] <- 1
      }else{
        auc_cov_ppi_vec[i] <- 0
      }
      
      # save the i_th simulation results
      res.i <- list(mean_est_n_i=mean_est_n_vec[i],
                    mean_est_ppi_i=mean_est_ppi_vec[i],
                    
                    mean_sd_n_i=mean_sd_n_vec[i],
                    mean_sd_ppi_i=mean_sd_ppi_vec[i],
                    
                    tpr_est_n_i=tpr_est_n_vec[i],
                    tpr_est_ppi_i=tpr_est_ppi_vec[i],
                    
                    tpr_sd_n_i=tpr_sd_n_vec[i],
                    tpr_sd_ppi_i=tpr_sd_ppi_vec[i],
                    
                    fpr_est_n_i=fpr_est_n_vec[i],
                    fpr_est_ppi_i=fpr_est_ppi_vec[i],
                    
                    fpr_sd_n_i=fpr_sd_n_vec[i],
                    fpr_sd_ppi_i=fpr_sd_ppi_vec[i],
                    
                    auc_est_n_i=auc_est_n_vec[i],
                    auc_est_ppi_i=auc_est_ppi_vec[i],
                    
                    
                    auc_sd_n_i=auc_sd_n_vec[i],
                    auc_sd_ppi_i=auc_sd_ppi_vec[i]
      )
      save(res.i, file = file.path(current_path, paste0("res.", i, ".Rdata")))
    }
   
    # save the summary results for one simulation scenario
    res <- list(mean_target_0=mean_target_0,
                
                mean_est_n=mean(mean_est_n_vec),
                mean_est_ppi=mean(mean_est_ppi_vec),
                
                mean_sd_n=mean(mean_sd_n_vec),
                mean_sd_ppi=mean(mean_sd_ppi_vec),
                
                mean_cov_n= mean(mean_cov_n_vec),
                mean_cov_ppi= mean(mean_cov_ppi_vec),
                
                
                tpr_target_0=tpr_target_0,
                
                tpr_est_n=mean(tpr_est_n_vec),
                tpr_est_ppi=mean(tpr_est_ppi_vec),
                
                tpr_sd_n=mean(tpr_sd_n_vec),
                tpr_sd_ppi=mean(tpr_sd_ppi_vec),
                
                tpr_cov_n= mean(tpr_cov_n_vec),
                tpr_cov_ppi= mean(tpr_cov_ppi_vec),
                
                
                fpr_target_0=fpr_target_0,
                
                fpr_est_n=mean(fpr_est_n_vec),
                fpr_est_ppi=mean(fpr_est_ppi_vec),
                
                fpr_sd_n=mean(fpr_sd_n_vec),
                fpr_sd_ppi=mean(fpr_sd_ppi_vec),
            
                fpr_cov_n= mean(fpr_cov_n_vec),
                fpr_cov_ppi= mean(fpr_cov_ppi_vec),
                
                
                auc_target_0=auc_target_0,
                
                auc_est_n=mean(auc_est_n_vec),
                auc_est_ppi=mean(auc_est_ppi_vec),
                
                auc_sd_n=mean(auc_sd_n_vec),
                auc_sd_ppi=mean(auc_sd_ppi_vec),
                
                auc_cov_n= mean(auc_cov_n_vec),
                auc_cov_ppi= mean(auc_cov_ppi_vec)
                
    )
    save(res, file = file.path(current_path, "res.summary.Rdata"))
  }
  
  
}

# Simulations for Section 6.2
ppi.sim(nrep=2500,dat_gen=dat_gen_full,size_all = 50000)



# Code for generating Figure 3

library(tidyverse)
library(ggplot2)

models      <- c("all", "ideal", "noise", "bad","pure_noise")
metrics     <- c("mean", "tpr", "fpr", "auc")

ratio_list <- list()

for (model in models) {
  folder_path <- file.path(
    results_dir,
    sprintf("f_%s_size_50000", model)
  )
  
  # find all simulation files
  sim_files <- list.files(
    path        = folder_path,
    pattern     = "^res\\.[0-9]+\\.(RData|rda)$",
    full.names  = TRUE,
    ignore.case = TRUE
  )
  if (length(sim_files) == 0) {
    stop("No simulation files found in: ", folder_path)
  }

  # load simulations 
  n_rep <- length(sim_files)
  sims  <- vector("list", n_rep)
  for (i in seq_len(n_rep)) {
    tmp_env  <- new.env()
    obj_name <- load(sim_files[i], envir = tmp_env)
    sims[[i]] <- tmp_env[[obj_name]]
  }
  
  # compute mean(sd_ppi/sd_n) for each metric
  for (metric in metrics) {
    est_n   <- sapply(sims, function(x) x[[paste0(metric, "_sd_n_i")]])
    est_ppi <- sapply(sims, function(x) x[[paste0(metric, "_sd_ppi_i")]])
    ratio_vec   <- est_ppi / est_n
    mean_ratio  <- mean(ratio_vec, na.rm = TRUE)
    
    ratio_list[[length(ratio_list) + 1]] <-
      data.frame(
        model      = model,
        metric     = metric,
        mean_ratio = mean_ratio,
        stringsAsFactors = FALSE
      )
  }
}

df_ratio <- bind_rows(ratio_list)
df_ratio$metric <- recode(df_ratio$metric,
                          mean = "Mean",
                          tpr  = "TPR",
                          fpr  = "FPR",
                          auc  = "AUC")
df_ratio$metric <- factor(df_ratio$metric, levels = c("Mean", "TPR", "FPR", "AUC"))
df_ratio$model[df_ratio$model=="ideal"] <- "aaideal"

ratio_sec3 <- ggplot(df_ratio,
                     aes(x     = model,
                         y     = mean_ratio,
                         color = model,
                         shape = model,
                         group = model)) +
  geom_point(size = 3) +
  facet_wrap(~ metric, nrow = 1) +
  labs(
    x     = "Prediction model",
    y     = expression(mean(hat(sigma)[PPI] / hat(sigma)[lab])),
    color = "Prediction model",
    shape = "Prediction model",
    title = ""
  ) +
  scale_color_manual(
    values = c("all"="salmon","bad"="chartreuse4","aaideal"="turquoise3","noise"="orchid","pure_noise"="orange2"),
    labels = c(
      "all"= expression(RF*"("*Y* plain("~")*bold(X)*")"),
      "bad"= expression(RF*"("*Y* plain("~")*R+X[1]+X[4]+X[5]*")"),
      "noise"= expression(RF*"("*Y* plain("~")*X[4]+X[5]*")"),
      "aaideal"= expression(Pr(Y==1~"|"~ bold(X))),
      "pure_noise"="Unif[0.01,0.99]"
    )
  )+
  scale_shape_manual(
    values = c("all"=16,"bad"=17,"aaideal"=15,"noise"=18,"pure_noise"=4),
    labels = c(
      "all"= expression(RF*"("*Y* plain("~")*bold(X)*")"),
      "bad"= expression(RF*"("*Y* plain("~")*R+X[1]+X[4]+X[5]*")"),
      "noise"= expression(RF*"("*Y* plain("~")*X[4]+X[5]*")"),
      "aaideal"= expression(Pr(Y==1~"|"~ bold(X))),
      "pure_noise"="Unif[0.01,0.99]"
    )
  )+
  scale_x_discrete(
    breaks = c("aaideal","all","bad","noise", "pure_noise"),        
    labels = c(
      "all"= expression(RF*"("*Y* plain("~")*bold(X)*")"),
      "bad"= expression(RF*"("*Y* plain("~")*R+X[1]+X[4]+X[5]*")"),
      "noise"= expression(RF*"("*Y* plain("~")*X[4]+X[5]*")"),
      "aaideal"= expression(Pr(Y==1~"|"~ bold(X))),
      "pure_noise"="Unif[0.01,0.99]"
    )
  ) +
  theme_minimal() +
  theme_bw()+
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14),
    panel.spacing = unit(1, "lines"),
    strip.text    = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )
ggsave(
  filename = file.path(results_dir, "ratio_cov_shift.png"),
  plot = ratio_sec3,
  width = 8,       
  height = 5,      
  units = "in",    
  dpi = 300        
)


# Code for generating Figure 4
models     <- c("all", "ideal", "noise", "bad","pure_noise")
alpha_lvls <- seq(0.01, 0.99, by = 0.05) # nominal coverage levels
metrics    <- c("mean", "tpr", "fpr", "auc")
methods    <- c("n", "ppi")

cov_results <- map_dfr(models, function(model) {
  # build the folder path
  folder_name <- sprintf("f_%s_size_50000", model)
  folder_path <- file.path(results_dir, sprintf("f_%s_size_50000", model))
  
  if (!dir.exists(folder_path)) {
    stop("Folder not found: ", folder_path)
  }
  
  # load the targets
  sum_env <- new.env()
  load(file.path(folder_path, "res.summary.Rdata"), envir = sum_env)
  targets <- sum_env$res
  
  # find all sim files 
  sim_files <- list.files(
    path        = folder_path,
    pattern     = "^res\\.[0-9]+\\.(RData|rda)$",
    full.names  = TRUE,
    ignore.case = TRUE
  )
  if (length(sim_files) == 0) {
    stop("No simulation files found in: ", folder_path)
  }
  cat("  Loading simulation files\n")
  
  # load each one into its own env
  sims <- map(sim_files, function(f) {
    e <- new.env()
    load(f, envir = e)
    # assume the first object loaded is your list
    e[[ls(e)[1]]]
  })
  cat("  Loaded all simulations\n")
  
  cat("  Computing coverage...\n")
  # 2c) compute coverage
  expand.grid(alpha  = alpha_lvls,
              metric = metrics,
              method = methods,
              stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    mutate(
      Z = qnorm(1 - alpha / 2),
      coverage = pmap_dbl(
        list(metric, method, Z),
        function(metric, method, Z) {
          true_val <- targets[[paste0(metric, "_target_0")]]
          ests     <- map_dbl(sims, ~ .x[[paste0(metric, "_est_",  method, "_i")]])
          sds      <- map_dbl(sims, ~ .x[[paste0(metric, "_sd_",   method, "_i")]])
          mean((true_val >= ests - Z * sds) & (true_val <= ests + Z * sds),
               na.rm = TRUE)
        }
      ),
      model = model
    ) %>%
    select(model, alpha, method, metric, coverage)
})

df <- cov_results
df_n_collapsed <- df %>%
  filter(method == "n") %>%
  group_by(alpha,metric) %>%
  summarise(
    coverage = first(coverage),
    method = "n",
    model = "zonlylabest",
    .groups = "drop"
  )
df_others <- df %>%
  filter(!(method == "n"))

df_final <- bind_rows(df_others, df_n_collapsed)
df_final$metric <- recode(df_final$metric,
                          mean = "Mean",
                          tpr  = "TPR",
                          fpr  = "FPR",
                          auc  = "AUC")
df_final$metric <- factor(df_final$metric, levels = c("Mean", "TPR", "FPR", "AUC"))

cov_sec3 <- ggplot(df_final,
                   aes(x = 1 - alpha, y = coverage,
                       color = model, linetype = model)) +
  geom_point(aes(shape = model), size = 1.5) +
  facet_wrap(~ metric, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_color_manual(
    values = c(all="salmon", bad="chartreuse4", ideal="turquoise3",
               noise="orchid", zonlylabest="orange2", pure_noise="navy"),
    labels = c(
      all = expression(RF*"("*Y* plain("~")*R+X[1]+...+X[5]*")"),
      bad = expression(RF*"("*Y* plain("~")*R+X[1]+X[4]+X[5]*")"),
      noise = expression(RF*"("*Y* plain("~")*X[4]+X[5]*")"),
      ideal = expression(Pr(Y==1~"|"~ R,X[1],...,X[5])),
      pure_noise = "Unif[0.01,0.99]",
      zonlylabest = "Only using labeled data"
    )
  ) +
  scale_linetype_manual(
    values = c(all="solid", bad="dashed", ideal="dotted",
               noise="dotdash", zonlylabest="longdash", pure_noise="twodash"),
    labels = c(
      all = expression(RF*"("*Y* plain("~")*R+X[1]+...+X[5]*")"),
      bad = expression(RF*"("*Y* plain("~")*R+X[1]+X[4]+X[5]*")"),
      noise = expression(RF*"("*Y* plain("~")*X[4]+X[5]*")"),
      ideal = expression(Pr(Y==1~"|"~ R,X[1],...,X[5])),
      pure_noise = "Unif[0.01,0.99]",
      zonlylabest = "Only using labeled data"
    )
  ) +
  scale_shape_manual(
    values = c(all=16, bad=17, ideal=15, noise=18, zonlylabest=4, pure_noise=6),
    labels = c(
      all = expression(RF*"("*Y* plain("~")*R+X[1]+...+X[5]*")"),
      bad = expression(RF*"("*Y* plain("~")*R+X[1]+X[4]+X[5]*")"),
      noise = expression(RF*"("*Y* plain("~")*X[4]+X[5]*")"),
      ideal = expression(Pr(Y==1~"|"~ R,X[1],...,X[5])),
      pure_noise = "Unif[0.01,0.99]",
      zonlylabest = "Only using labeled data"
    )
  ) +
  labs(color = NULL, linetype = NULL, shape = NULL,
       x = "Nominal coverage", y = "Empirical coverage", title = "") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "red", linewidth = 0.6) +
  theme_minimal() +
  theme_bw()+
  theme(panel.spacing = unit(2, "lines"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(
  filename =  file.path(results_dir, "cov_sec3.png"),
  plot = cov_sec3,
  width = 8,      
  height = 4,     
  units = "in",    
  dpi = 300       
)

# Code for generating Figure S.7
models  <- c("all", "ideal", "noise", "bad", "pure_noise")
metrics <- c("mean", "tpr", "fpr", "auc")

base_dir <- results_dir

diff_list <- list()

for (model in models) {
  print(model)
  
  folder_path <- file.path(base_dir, sprintf("f_%s_size_50000", model))
  
  if (!dir.exists(folder_path)) {
    stop("Folder not found: ", folder_path)
  }
  
  # load res.summary.RData
  sum_env <- new.env()
  load(file.path(folder_path, "res.summary.Rdata"), envir = sum_env)
  res <- sum_env$res
  
  # compute mean(est_ppi - target_0)
  for (metric in metrics) {
    
    est_name    <- paste0(metric, "_est_ppi")
    target_name <- paste0(metric, "_target_0")
    
    if (!(est_name %in% names(res))) {
      stop("Missing ", est_name, " in ", folder_path)
    }
    if (!(target_name %in% names(res))) {
      stop("Missing ", target_name, " in ", folder_path)
    }
    
    diff_vec  <- res[[est_name]] - res[[target_name]]
    mean_diff <- mean(diff_vec, na.rm = TRUE)
    
    diff_list[[length(diff_list) + 1]] <-
      data.frame(
        model     = model,
        metric    = metric,
        mean_diff = mean_diff,
        stringsAsFactors = FALSE
      )
  }
}

df_diff <- bind_rows(diff_list)
df_diff$metric <- recode(df_diff$metric,
                         mean = "Mean",
                         tpr  = "TPR",
                         fpr  = "FPR",
                         auc  = "AUC")
df_diff$metric <- factor(df_diff$metric,
                         levels = c("Mean", "TPR", "FPR", "AUC"))

df_diff$model[df_diff$model == "ideal"] <- "aaideal"

diff_sec3 <- ggplot(df_diff,
                    aes(x     = model,
                        y     = mean_diff,
                        color = model,
                        shape = model,
                        group = model)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  facet_wrap(~ metric, nrow = 1) +
  labs(
    x     = "Prediction model",
    y     = expression(mean(hat(theta)[PPI] - theta[0])),
    color = "Prediction model",
    shape = "Prediction model",
    title = ""
  ) +
  scale_color_manual(
    values = c("all"="salmon",
               "bad"="chartreuse4",
               "aaideal"="turquoise3",
               "noise"="orchid",
               "pure_noise"="orange2"),
    labels = c(
      "all"= expression(RF*"("*Y* plain("~")*bold(X)*")"),
      "bad"= expression(RF*"("*Y* plain("~")*R+X[1]+X[4]+X[5]*")"),
      "noise"= expression(RF*"("*Y* plain("~")*X[4]+X[5]*")"),
      "aaideal"= expression(Pr(Y==1~"|"~ bold(X))),
      "pure_noise"="Unif[0.01,0.99]"
    )
  ) +
  scale_shape_manual(
    values = c("all"=16,"bad"=17,"aaideal"=15,"noise"=18,"pure_noise"=4),
    labels = c(
      "all"= expression(RF*"("*Y* plain("~")*bold(X)*")"),
      "bad"= expression(RF*"("*Y* plain("~")*R+X[1]+X[4]+X[5]*")"),
      "noise"= expression(RF*"("*Y* plain("~")*X[4]+X[5]*")"),
      "aaideal"= expression(Pr(Y==1~"|"~ bold(X))),
      "pure_noise"="Unif[0.01,0.99]"
    )
  ) +
  scale_x_discrete(
    breaks = c("aaideal","all","bad","noise","pure_noise"),
    labels = c(
      "all"= expression(RF*"("*Y* plain("~")*bold(X)*")"),
      "bad"= expression(RF*"("*Y* plain("~")*R+X[1]+X[4]+X[5]*")"),
      "noise"= expression(RF*"("*Y* plain("~")*X[4]+X[5]*")"),
      "aaideal"= expression(Pr(Y==1~"|"~ bold(X))),
      "pure_noise"="Unif[0.01,0.99]"
    )
  ) +
  theme_minimal() +
  theme_bw() +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    axis.text     = element_text(size = 10),
    axis.title    = element_text(size = 14),
    panel.spacing = unit(1, "lines"),
    strip.text    = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(results_dir, "bias_sec3.png"),
  plot     = diff_sec3,
  width    = 8,
  height   = 5,
  units    = "in",
  dpi      = 300
)
