# simulation code for Section B.2.2 (illustrating Section 5.3)
# PPI estimators for TPR, FPR, and AUC when f(X) itself is the biomarker
expit <- function(a){ return(exp(a)/(1+exp(a)))}

## tpr
# point estimation
tpr.est.yn <- function(dat){
  return(sum((dat$score>alpha)*dat$Y)/sum(dat$Y))
}
tpr.est.fxn <- function(dat){
  return(sum((dat$score>alpha)*dat$f)/sum(dat$f))
}
# influence function computation
tpr.phi.yn <- function(dat){
  return(1*(dat$score>alpha)*dat$Y/mean(dat$Y) - mean((dat$score>alpha)*dat$Y)*dat$Y/(mean(dat$Y))^2)
}
tpr.phi.fxn <- function(dat){
  return(1*(dat$score>alpha)*dat$f/mean(dat$f) - mean((dat$score>alpha)*dat$f)*dat$f/(mean(dat$f))^2)
}

## fpr
# point estimation
fpr.est.yn <- function(dat){
  return(sum((dat$score>alpha)*(1-dat$Y))/sum(1-dat$Y))
}
fpr.est.fxn <- function(dat){
  return(sum((dat$score>alpha)*(1-dat$f))/sum(1-dat$f))
}
# influence function computation
fpr.phi.yn <- function(dat){
  return(1*(dat$score>alpha)*(1-dat$Y)/mean(1-dat$Y) - mean((dat$score>alpha)*(1-dat$Y))*(1-dat$Y)/(mean(1-dat$Y))^2)
}
fpr.phi.fxn <- function(dat){
  return(1*(dat$score>alpha)*(1-dat$f)/mean(1-dat$f) - mean((dat$score>alpha)*(1-dat$f))*(1-dat$f)/(mean(1-dat$f))^2)
}

## AUC
# point estimation
auc.est.yn <- function(dat){
  
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
auc.est.fxn <- function(dat){
  
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
  return(.IC(fold_preds=predictions, fold_labels=labels, pos=1, neg=0, w1=w1, w0=w0))
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
  return(.IC(fold_preds=predictions, fold_labels=labels, w1=w1, w0=w0))
}



library(data.table)
library(ROCR)
library(MASS)

nvec <- c(1000) # n=1000
lambdavec <- c(0.1) # \lambda=0.1
f_types <- c("all","ideal")
# different prediction models f
# "all" is a random forest model predicting Y from (R,X_1,...,X_5)
# "ideal" represents E[Y|R,X_1,...,X_5]

alpha <- 0.6
results_dir <- paste0("your result path") # path to save the simulation results

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

m <- 100000 # size of the training data for estimating the random forest prediction model

# generate (R,X_1,...,X_5) for the training data, with the same data generation process in Section 6.1
train_dat <- as.data.frame(mvrnorm(m, mu = c(0,0,0,0,3,20), Sigma = Var_mat)) 

names(train_dat) <- c("R","X1","X2","X3","X4","X5")

# compute E[Y|R,X_1,...,X_5]
train_dat$ideal <- compute_eta_p(train_dat)

# generate Y for the training data, following the prespecified model E[Y|R,X_1,...,X_5]
train_dat$Y <- rbinom(m, size=1, prob=train_dat$ideal)
train_dat$Y <- factor(train_dat$Y, levels=0:1)  # for classification

library(randomForest)

# A random forest model predicting Y from (R,X_1,...,X_5)
rf_all <- randomForest(Y ~ R + X1 + X2 + X3 + X4 + X5,
                       data = train_dat)

# data generation -------------------------------------------------------------------
dat_gen <- function(n,N,f_type){
  
  # generate the data to construct the labeled data estimators and PPI estimators
  # following the same data generation process described in Section 6.1
  inf_dat <- as.data.frame(mvrnorm((n+N), mu=c(0,0,0,0,3,20), Sigma=Var_mat))
  names(inf_dat) <- c("R","X1","X2","X3","X4","X5")
  
  inf_dat$ideal <- compute_eta_p(inf_dat) # the prespecified ideal f: E[Y|R,X_1,...,X_5]
  
  inf_dat$Y <- rbinom(n+N, size=1, prob=inf_dat$ideal)  # generating Y following the prespecified model
  
  # obtain predictions using the "all" random forest model
  inf_dat$all <- predict(rf_all, newdata = inf_dat,type = "prob")[ ,"1"]
  
  # the f used for constructing PPI estimators
  inf_dat$f <- inf_dat[,f_type]
  
  # set f itself to be the biomarker
  inf_dat$score <- inf_dat[,f_type]
  
  return(inf_dat)
}

# simulation function
ppi.sim <- function(nrep,dat_gen,n,lambda){
  
  N <- n/lambda # compute N, the unlabeled dataset size
  
  ## 10000 replications using fully labeled data to compute the target TPR, FPR, AUC
  target_rep <- 10000
  
  for(k in 1:length(f_types)){ # different f's corresponding to different biomarkers, and thus different targets

    f_type <- f_types[k]
    tpr_target_0_vec <- fpr_target_0_vec <- auc_target_0_vec <- vector("numeric",length = target_rep)
    for(index in 1:target_rep){
      
      set.seed(index+4564356)
      
      inf_dat.index <- as.data.frame(mvrnorm((n+N), mu=c(0,0,0,0,3,20), Sigma=Var_mat)) # generating (R,X_1,...,X_5)
      names(inf_dat.index) <- c("R","X1","X2","X3","X4","X5")
      
      # obtaining E[Y|R,X_1,...,X_5] as one of the biomarker candidates
      inf_dat.index$ideal <- compute_eta_p(inf_dat.index)
      inf_dat.index$Y <- rbinom((n+N), size=1, prob=inf_dat.index$ideal) # generating Y
      
      # obtaining f_all(R,X_1,...,X_5) as one of the biomarker candidates
      inf_dat.index$all <- predict(rf_all, newdata = inf_dat.index,type = "prob")[ ,"1"]
      
      inf_dat.index$score <- inf_dat.index[,f_type] # set the f itself as the biomarker
      
      tpr_target_0_vec[index] <- tpr.est.yn(inf_dat.index)
      fpr_target_0_vec[index] <- fpr.est.yn(inf_dat.index)
      auc_target_0_vec[index] <- auc.est.yn(inf_dat.index)
    }
    # compute the median of the 10000 replications as the targets
    tpr_target_0 <- median(tpr_target_0_vec)
    fpr_target_0 <- median(fpr_target_0_vec)
    auc_target_0 <- median(auc_target_0_vec)
    
    # for saving the labeled data estimators
    tpr_est_n_vec <- fpr_est_n_vec <- auc_est_n_vec <- vector("numeric",length = nrep)
    # for saving the PPI estimators
    tpr_est_ppi_vec <- fpr_est_ppi_vec <- auc_est_ppi_vec <- vector("numeric",length = nrep)
    # for saving the \omega
    tpr_omega_vec <- fpr_omega_vec <- auc_omega_vec  <- vector("numeric",length = nrep)
    # for saving the se estimates for the labeled data estimators
    tpr_sd_n_vec <- fpr_sd_n_vec <- auc_sd_n_vec <- vector("numeric",length=nrep)
    # for saving the se estimates for the PPI estimators
    tpr_sd_ppi_vec <- fpr_sd_ppi_vec <- auc_sd_ppi_vec <-  vector("numeric",length=nrep)
    # for saving the coverage for the labeled data estimators
    tpr_cov_n_vec <- fpr_cov_n_vec <- auc_cov_n_vec <- vector("numeric",length=nrep)
    # for saving the coverage for the PPI estimators
    tpr_cov_ppi_vec <- fpr_cov_ppi_vec <- auc_cov_ppi_vec <-  vector("numeric",length=nrep)
    
    current_path <- paste(results_dir,'f_',f_types[k],'_n_',n,'_lambda_',lambda,'/',sep='')
    dir.create(file.path(current_path),showWarnings = FALSE)
    
    # simulations -------------------------------------------------------------
    
    for( i in 1:nrep){
      
      print(i)
      set.seed(i)
      dat.i <- dat_gen(n,N,f_types[k])
      
      
      # Y is observed in n observations
      dat_n.i <- dat.i[1:n,]
      
      tpr_phi_n.i <- tpr.phi.yn(dat_n.i)
      fpr_phi_n.i <- fpr.phi.yn(dat_n.i)
      auc_phi_n.i <- auc.phi.yn(dat_n.i)
      
      tpr_phi_f_n.i <- tpr.phi.fxn(dat_n.i)
      fpr_phi_f_n.i <- fpr.phi.fxn(dat_n.i)
      auc_phi_f_n.i <- auc.phi.fxn(dat_n.i)
      
      tpr_est_n.i <- tpr.est.yn(dat_n.i)
      fpr_est_n.i <- fpr.est.yn(dat_n.i)
      auc_est_n.i <- auc.est.yn(dat_n.i)
      
      # record the labeled data estimators
      tpr_est_n_vec[i] <- tpr_est_n.i
      fpr_est_n_vec[i] <- fpr_est_n.i
      auc_est_n_vec[i] <- auc_est_n.i
      
      # record the se estimates of the labeled data estimator
      tpr_sd_n_vec[i] <- sqrt(var(tpr_phi_n.i))/sqrt(n)
      fpr_sd_n_vec[i] <- sqrt(var(fpr_phi_n.i))/sqrt(n)
      auc_sd_n_vec[i] <- sqrt(var(auc_phi_n.i))/sqrt(n)
      
      # record the coverage of the labeled data estimator
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
      
      tpr_phi_f_all.i <- tpr.phi.fxn(dat.i)
      fpr_phi_f_all.i <- fpr.phi.fxn(dat.i)
      auc_phi_f_all.i <- auc.phi.fxn(dat.i)
      
      # compute \omega
      tpr_omega_vec[i] <- cov(tpr_phi_n.i,tpr_phi_f_n.i)/(var(tpr_phi_f_all.i))
      fpr_omega_vec[i] <- cov(fpr_phi_n.i,fpr_phi_f_n.i)/(var(fpr_phi_f_all.i))
      auc_omega_vec[i] <- cov(auc_phi_n.i,auc_phi_f_n.i)/(var(auc_phi_f_all.i))
      
      tpr_est_f_all.i <- tpr.est.fxn(dat.i)
      tpr_est_f_n.i <- tpr.est.fxn(dat_n.i)
      fpr_est_f_all.i <- fpr.est.fxn(dat.i)
      fpr_est_f_n.i <- fpr.est.fxn(dat_n.i)
      auc_est_f_all.i <- auc.est.fxn(dat.i)
      auc_est_f_n.i <- auc.est.fxn(dat_n.i)
      
      # construct the PPI estimator
      tpr_est_ppi_vec[i] <- tpr_est_n.i + tpr_omega_vec[i]*(tpr_est_f_all.i - tpr_est_f_n.i)
      fpr_est_ppi_vec[i] <- fpr_est_n.i + fpr_omega_vec[i]*(fpr_est_f_all.i - fpr_est_f_n.i)
      auc_est_ppi_vec[i] <- auc_est_n.i + auc_omega_vec[i]*(auc_est_f_all.i - auc_est_f_n.i)
      
      # record the se estimates of the PPI estimators
      tpr_sd_ppi_vec[i] <- sqrt((var(tpr_phi_n.i)-cov(tpr_phi_f_n.i,tpr_phi_n.i)^2/(var(tpr_phi_f_all.i)+lambda*var(tpr_phi_f_all.i)))/n)
      fpr_sd_ppi_vec[i] <- sqrt((var(fpr_phi_n.i)-cov(fpr_phi_f_n.i,fpr_phi_n.i)^2/(var(fpr_phi_f_all.i)+lambda*var(fpr_phi_f_all.i)))/n)
      auc_sd_ppi_vec[i] <- sqrt((var(auc_phi_n.i)-cov(auc_phi_f_n.i,auc_phi_n.i)^2/(var(auc_phi_f_all.i)+lambda*var(auc_phi_f_all.i)))/n)
      
      if(tpr_sd_ppi_vec[i]=='NaN'|fpr_sd_ppi_vec[i]=='NaN'|auc_sd_ppi_vec[i]=='NaN'){
        print(i)
        stop()
      }
      
      # record the coverage of the PPI estimators
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
      res.i <- list(
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
      save(res.i, file = paste(current_path,"res.",i,'.Rdata',sep=''))
    }
    
    # save the simulation results for one f as the biomarker
    res <- list(
      tpr_target_0=tpr_target_0,
      
      tpr_est_n=median(tpr_est_n_vec),
      tpr_est_ppi=median(tpr_est_ppi_vec),
      
      tpr_sd_n=median(tpr_sd_n_vec),
      tpr_sd_ppi=median(tpr_sd_ppi_vec),
      
      tpr_cov_n= mean(tpr_cov_n_vec),
      tpr_cov_ppi= mean(tpr_cov_ppi_vec),

      
      fpr_target_0=fpr_target_0,
      
      fpr_est_n=median(fpr_est_n_vec),
      fpr_est_ppi=median(fpr_est_ppi_vec),
      
      fpr_sd_n=median(fpr_sd_n_vec),
      fpr_sd_ppi=median(fpr_sd_ppi_vec),
      
      fpr_cov_n= mean(fpr_cov_n_vec),
      fpr_cov_ppi= mean(fpr_cov_ppi_vec),
      
      
      auc_target_0=auc_target_0,
      
      auc_est_n=median(auc_est_n_vec),
      auc_est_ppi=median(auc_est_ppi_vec),
      
      auc_sd_n=median(auc_sd_n_vec),
      auc_sd_ppi=median(auc_sd_ppi_vec),
      
      auc_cov_n= mean(auc_cov_n_vec),
      auc_cov_ppi= mean(auc_cov_ppi_vec),
    )
    save(res,file=paste(current_path,'res.summary.Rdata',sep=''))
  }
  
}

# Simulation for Section B.2.4, for n=1000, \lambda=0.1, over 2500 replications
ppi.sim(nrep=2500,dat_gen,n=nvec[1],lambda = lambdavec[1])





