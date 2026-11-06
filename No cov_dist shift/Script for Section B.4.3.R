# simulation code using the Section 5.1 no-covariate-shift scenario

# For each target, let phi(X,Y) be the full-data influence function and let
# m(X) = E{phi(X,Y) | X}.  The efficient-rectifier-type estimator is
#
#   theta_n + mean_{n+N}{m_hat(X)} - mean_{n}{m_hat(X)}.

# m_hat is estimated by an additive natural-spline regression using (R,X1,...,X5,f(X)).  
# It is trained only on labeled observations and is cross-fitted.  
# The role of f is therefore only to help estimate m(X).

expit <- function(a){ return(exp(a)/(1+exp(a)))}

## TPR
# point estimation
tpr.est.yn <- function(dat){
  return(sum((dat$R>alpha)*dat$Y)/sum(dat$Y))
}

# influence function computation
tpr.phi.yn <- function(dat){
  mean.y <- mean(dat$Y)
  mean.ry <- mean((dat$R>alpha)*dat$Y)
  return((dat$R>alpha)*dat$Y/mean.y - mean.ry*dat$Y/mean.y^2)
}

## FPR
# point estimation
fpr.est.yn <- function(dat){
  return(sum((dat$R>alpha)*(1-dat$Y))/sum(1-dat$Y))
}

# influence function computation
fpr.phi.yn <- function(dat){
  mean.one.minus.y <- mean(1-dat$Y)
  mean.r.one.minus.y <- mean((dat$R>alpha)*(1-dat$Y))
  return((dat$R>alpha)*(1-dat$Y)/mean.one.minus.y -
           mean.r.one.minus.y*(1-dat$Y)/mean.one.minus.y^2)
}

## AUC point estimation
auc.est.yn <- function(dat){
  score.order <- order(dat$score, decreasing=TRUE)
  score.sorted <- dat$score[score.order]

  dups <- rev(duplicated(rev(score.sorted)))

  ## TPR
  TPR <- cumsum(dat$Y[score.order])/sum(dat$Y[score.order])
  TPR <- c(0,TPR[!dups])

  ## FPR
  FPR <- cumsum(1-dat$Y[score.order])/sum(1-dat$Y[score.order])
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
  
  predictions <- dat$score
  labels <- dat$Y
  n_obs <- nrow(dat)
  
  w1 <- 1 / (sum(labels == 1) / n_obs)
  w0 <- 1 / (sum(labels == 0) / n_obs)
  
  pred <- label <- row.id <- NULL
  fracNegLabelsWithSmallerPreds <- fracPosLabelsWithLargerPreds <- icVal <- NULL
  
  .IC <- function(fold_preds, fold_labels, pos, neg, w1, w0){
    
    n_rows <- length(fold_labels)
    n_pos <- sum(fold_labels == pos)
    n_neg <- n_rows - n_pos
    
    auc <- AUC(fold_preds, fold_labels)
    
    DT <- data.table::data.table(
      row.id = seq_along(fold_labels),
      pred = fold_preds,
      label = fold_labels
    )
    
    DT <- DT[order(pred, -xtfrm(label))]
    DT[, fracNegLabelsWithSmallerPreds := cumsum(label == neg) / n_neg]
    
    DT <- DT[order(-pred, label)]
    DT[, fracPosLabelsWithLargerPreds := cumsum(label == pos) / n_pos]
    
    DT[, icVal := ifelse(label == pos,
                         w1 * (fracNegLabelsWithSmallerPreds - auc),
                         w0 * (fracPosLabelsWithLargerPreds - auc))]
    
    DT <- DT[order(row.id)]
    return(DT$icVal)
  }
  return(.IC(fold_preds = predictions,
             fold_labels = labels, 
             pos = 1, neg = 0,w1 = w1, w0 = w0))
}
  
## Cross-fitted spline regression for m(X) = E{phi(X,Y) | X}
# The three influence functions share the same spline design, 
# so a multivariate linear model estimates them together within each training fold.
crossfit.conditional.if <- function(dat,n,nfolds=5,spline_df=4){
  N <- nrow(dat)-n
  labeled.index <- seq_len(n)
  unlabeled.index <- n+seq_len(N)

  # no outcome from a held-out labeled fold is used to fit the spline model that predicts that fold.
  labeled.fold <- sample(rep(seq_len(nfolds),length.out=n))
  unlabeled.fold <- sample(rep(seq_len(nfolds),length.out=N))

  mhat <- matrix(NA_real_,nrow=nrow(dat),ncol=3)
  colnames(mhat) <- c("tpr","fpr","auc")

  spline.terms <- paste0(
    "splines::ns(",c("R","X1","X2","X3","X4","X5","f"),
    ", df=",spline_df,")"
  )
  spline.formula <- as.formula(
    paste("cbind(phi_tpr,phi_fpr,phi_auc) ~",
          paste(spline.terms,collapse=" + "))
  )

  for(fold in seq_len(nfolds)){
    train.index <- labeled.index[labeled.fold!=fold]
    validation.index <- labeled.index[labeled.fold==fold]
    prediction.index <- c(validation.index,
                          unlabeled.index[unlabeled.fold==fold])

    spline.dat <- dat[train.index,]
    spline.dat$phi_tpr <- tpr.phi.yn(spline.dat)
    spline.dat$phi_fpr <- fpr.phi.yn(spline.dat)
    spline.dat$phi_auc <- auc.phi.yn(spline.dat)

    spline.fit <- lm(spline.formula,data=spline.dat)
    fold.prediction <- predict(spline.fit,newdata=dat[prediction.index,])
    fold.prediction <- as.matrix(fold.prediction)

    mhat[prediction.index,] <- fold.prediction
  }


  return(as.data.frame(mhat))
}


library(MASS)

nvec <- c(1000)
lambdavec <- c(0.8,0.5,0.25,0.1,0.01)

f_types <- c("all","bad","noise","ideal","pure_noise")
# "all" is a random forest model predicting Y from (R,X1,...,X5)
# "bad" is a random forest model predicting Y from (R,X1,X4,X5)
# "noise" is a random forest model predicting Y from (X4,X5)
# "ideal" is E[Y|R,X1,...,X5]
# "pure_noise" is Unif[0.01,0.99]

alpha <- 0.6
nfolds <- 5
spline_df <- 4
results_dir <- "results/noncovshift_tprfprauc_eff"
dir.create(results_dir,recursive=TRUE,showWarnings=FALSE)
set.seed(2025)

# The covariance matrix of (R,X1,...,X5)
Var_mat <- matrix(c(
  1.0, 0.5, 0.3, 0.2, 0, 0,
  0.5, 1.0, 0.4, 0.3, 0, 0,
  0.3, 0.4, 1.0, 0.2, 0, 0,
  0.2, 0.3, 0.2, 1.0, 0, 0,
  0,   0,   0,   0,   1, 0,
  0,   0,   0,   0,   0, 2
),nrow=6,byrow=TRUE)

# Prespecified E[Y|R,X1,...,X5]
beta_intcpt <- -0.5
beta_R <- 1.0
beta_X1q <- -0.9
beta_X2a <- 0.6
beta_X3c <- 0.5
beta_Rx3 <- 1.5
beta_X1x2 <- -0.7

compute_eta_p <- function(df){
  eta <- with(df,
              beta_intcpt+
                beta_R*R+
                beta_X1q*X1^2+
                beta_X2a*abs(X2)+
                beta_X3c*X3^3+
                beta_Rx3*(R*X3)+
                beta_X1x2*(X1*X2))
  plogis(eta)
}

m <- 100000

# Generate the independent training data for the five fixed prediction models f.
train_dat <- as.data.frame(mvrnorm(m,mu=c(0,0,0,0,3,20),Sigma=Var_mat))
names(train_dat) <- c("R","X1","X2","X3","X4","X5")
train_dat$ideal <- compute_eta_p(train_dat)
train_dat$Y <- rbinom(m,size=1,prob=train_dat$ideal)
train_dat$Y <- factor(train_dat$Y,levels=0:1)

library(randomForest)

rf_all <- randomForest(Y~R+X1+X2+X3+X4+X5,data=train_dat)
rf_noise <- randomForest(Y~X4+X5,data=train_dat)
rf_bad <- randomForest(Y~R+X1+X4+X5,data=train_dat)


# Data generation function: unchanged from the Section 5.1/6.1 scenario.
dat_gen <- function(n,N,f_type){
  inf_dat <- as.data.frame(
    mvrnorm(n+N,mu=c(0,0,0,0,3,20),Sigma=Var_mat)
  )
  names(inf_dat) <- c("R","X1","X2","X3","X4","X5")

  inf_dat$ideal <- compute_eta_p(inf_dat)
  inf_dat$Y <- rbinom(n+N,size=1,prob=inf_dat$ideal)

  inf_dat$all <- predict(rf_all,newdata=inf_dat,type="prob")[,"1"]
  inf_dat$noise <- predict(rf_noise,newdata=inf_dat,type="prob")[,"1"]
  inf_dat$pure_noise <- runif(n+N,min=0.01,max=0.99)
  inf_dat$bad <- predict(rf_bad,newdata=inf_dat,type="prob")[,"1"]

  inf_dat$f <- inf_dat[,f_type]
  inf_dat$score <- inf_dat$R

  return(inf_dat)
}


# Simulation function
efficient.rectifier.sim <- function(nrep,dat_gen,n,lambda,
                                    nfolds=5,spline_df=4,
                                    target_rep=10000){
  N <- as.integer(round(n/lambda))
  size_all <- n+N
  unlabeled.weight <- N/size_all

  ## Fully labeled Monte Carlo data are used only to compute fixed targets.
  tpr_target_0_vec <- fpr_target_0_vec <- auc_target_0_vec <-
    vector("numeric",length=target_rep)

  for(index in seq_len(target_rep)){
    set.seed(index+4564356)

    inf_dat.index <- as.data.frame(
      mvrnorm(size_all,mu=c(0,0,0,0,3,20),Sigma=Var_mat)
    )
    names(inf_dat.index) <- c("R","X1","X2","X3","X4","X5")
    inf_dat.index$ideal <- compute_eta_p(inf_dat.index)
    inf_dat.index$Y <- rbinom(size_all,size=1,prob=inf_dat.index$ideal)
    inf_dat.index$score <- inf_dat.index$R

    tpr_target_0_vec[index] <- tpr.est.yn(inf_dat.index)
    fpr_target_0_vec[index] <- fpr.est.yn(inf_dat.index)
    auc_target_0_vec[index] <- auc.est.yn(inf_dat.index)
  }

  
  tpr_target_0 <- mean(tpr_target_0_vec)
  fpr_target_0 <- mean(fpr_target_0_vec)
  auc_target_0 <- mean(auc_target_0_vec)

  for(k in seq_along(f_types)){
    # Labeled-data estimators
    tpr_est_n_vec <- fpr_est_n_vec <- auc_est_n_vec <-
      vector("numeric",length=nrep)
    # Efficient-rectifier-type estimators
    tpr_est_eff_vec <- fpr_est_eff_vec <- auc_est_eff_vec <-
      vector("numeric",length=nrep)
    # The fitted rectifier terms
    tpr_rectifier_vec <- fpr_rectifier_vec <- auc_rectifier_vec <-
      vector("numeric",length=nrep)
    # Estimated standard errors
    tpr_sd_n_vec <- fpr_sd_n_vec <- auc_sd_n_vec <-
      vector("numeric",length=nrep)
    tpr_sd_eff_vec <- fpr_sd_eff_vec <- auc_sd_eff_vec <-
      vector("numeric",length=nrep)
    # 95% confidence interval coverage
    tpr_cov_n_vec <- fpr_cov_n_vec <- auc_cov_n_vec <-
      vector("numeric",length=nrep)
    tpr_cov_eff_vec <- fpr_cov_eff_vec <- auc_cov_eff_vec <-
      vector("numeric",length=nrep)

    current_path <- file.path(
      results_dir,
      paste0("f_",f_types[k],"_n_",n,"_lambda_",lambda)
    )
    dir.create(current_path,recursive=TRUE,showWarnings=FALSE)

    for(i in seq_len(nrep)){
      set.seed(i)
      dat.i <- dat_gen(n,N,f_types[k])

      # Y is observed only in the first n observations.
      dat_n.i <- dat.i[seq_len(n),]
      unlabeled.index <- n+seq_len(N)

      # Cross-fitted estimates of E{phi(X,Y)|X}.  Only labeled outcomes enter the spline fits; the unlabeled observations are prediction-only rows.
      mhat.i <- crossfit.conditional.if(
        dat=dat.i,n=n,nfolds=nfolds,spline_df=spline_df
      )

      tpr_phi_n.i <- tpr.phi.yn(dat_n.i)
      fpr_phi_n.i <- fpr.phi.yn(dat_n.i)
      auc_phi_n.i <- auc.phi.yn(dat_n.i)

      tpr_est_n.i <- tpr.est.yn(dat_n.i)
      fpr_est_n.i <- fpr.est.yn(dat_n.i)
      auc_est_n.i <- auc.est.yn(dat_n.i)

      
      tpr_rectifier.i <- mean(mhat.i$tpr)-mean(mhat.i$tpr[seq_len(n)])
      fpr_rectifier.i <- mean(mhat.i$fpr)-mean(mhat.i$fpr[seq_len(n)])
      auc_rectifier.i <- mean(mhat.i$auc)-mean(mhat.i$auc[seq_len(n)])

      tpr_est_eff.i <- tpr_est_n.i+tpr_rectifier.i
      fpr_est_eff.i <- fpr_est_n.i+fpr_rectifier.i
      auc_est_eff.i <- auc_est_n.i+auc_rectifier.i

      
      tpr_sd_n.i <- sqrt(var(tpr_phi_n.i)/n)
      fpr_sd_n.i <- sqrt(var(fpr_phi_n.i)/n)
      auc_sd_n.i <- sqrt(var(auc_phi_n.i)/n)

      
      
      tpr_var_eff.i <- (
        mean(tpr_phi_n.i^2) -
          unlabeled.weight * var(mhat.i$tpr)
      ) / n
      
      fpr_var_eff.i <- (
        mean(fpr_phi_n.i^2) -
          unlabeled.weight * var(mhat.i$fpr)
      ) / n
      
      auc_var_eff.i <- (
        mean(auc_phi_n.i^2) -
          unlabeled.weight * var(mhat.i$auc)
      ) / n
      
      tpr_sd_eff.i <- sqrt(tpr_var_eff.i)
      fpr_sd_eff.i <- sqrt(fpr_var_eff.i)
      auc_sd_eff.i <- sqrt(auc_var_eff.i)


      tpr_est_n_vec[i] <- tpr_est_n.i
      fpr_est_n_vec[i] <- fpr_est_n.i
      auc_est_n_vec[i] <- auc_est_n.i

      tpr_est_eff_vec[i] <- tpr_est_eff.i
      fpr_est_eff_vec[i] <- fpr_est_eff.i
      auc_est_eff_vec[i] <- auc_est_eff.i

      tpr_rectifier_vec[i] <- tpr_rectifier.i
      fpr_rectifier_vec[i] <- fpr_rectifier.i
      auc_rectifier_vec[i] <- auc_rectifier.i

      tpr_sd_n_vec[i] <- tpr_sd_n.i
      fpr_sd_n_vec[i] <- fpr_sd_n.i
      auc_sd_n_vec[i] <- auc_sd_n.i

      tpr_sd_eff_vec[i] <- tpr_sd_eff.i
      fpr_sd_eff_vec[i] <- fpr_sd_eff.i
      auc_sd_eff_vec[i] <- auc_sd_eff.i

      tpr_cov_n_vec[i] <- as.numeric(
        tpr_est_n.i-1.96*tpr_sd_n.i<=tpr_target_0 &
          tpr_target_0<=tpr_est_n.i+1.96*tpr_sd_n.i
      )
      fpr_cov_n_vec[i] <- as.numeric(
        fpr_est_n.i-1.96*fpr_sd_n.i<=fpr_target_0 &
          fpr_target_0<=fpr_est_n.i+1.96*fpr_sd_n.i
      )
      auc_cov_n_vec[i] <- as.numeric(
        auc_est_n.i-1.96*auc_sd_n.i<=auc_target_0 &
          auc_target_0<=auc_est_n.i+1.96*auc_sd_n.i
      )

      tpr_cov_eff_vec[i] <- as.numeric(
        tpr_est_eff.i-1.96*tpr_sd_eff.i<=tpr_target_0 &
          tpr_target_0<=tpr_est_eff.i+1.96*tpr_sd_eff.i
      )
      fpr_cov_eff_vec[i] <- as.numeric(
        fpr_est_eff.i-1.96*fpr_sd_eff.i<=fpr_target_0 &
          fpr_target_0<=fpr_est_eff.i+1.96*fpr_sd_eff.i
      )
      auc_cov_eff_vec[i] <- as.numeric(
        auc_est_eff.i-1.96*auc_sd_eff.i<=auc_target_0 &
          auc_target_0<=auc_est_eff.i+1.96*auc_sd_eff.i
      )

      res.i <- list(
        tpr_est_n_i=tpr_est_n_vec[i],
        tpr_est_eff_i=tpr_est_eff_vec[i],
        tpr_rectifier_i=tpr_rectifier_vec[i],
        tpr_sd_n_i=tpr_sd_n_vec[i],
        tpr_sd_eff_i=tpr_sd_eff_vec[i],

        fpr_est_n_i=fpr_est_n_vec[i],
        fpr_est_eff_i=fpr_est_eff_vec[i],
        fpr_rectifier_i=fpr_rectifier_vec[i],
        fpr_sd_n_i=fpr_sd_n_vec[i],
        fpr_sd_eff_i=fpr_sd_eff_vec[i],

        auc_est_n_i=auc_est_n_vec[i],
        auc_est_eff_i=auc_est_eff_vec[i],
        auc_rectifier_i=auc_rectifier_vec[i],
        auc_sd_n_i=auc_sd_n_vec[i],
        auc_sd_eff_i=auc_sd_eff_vec[i]
      )
      save(res.i,file=file.path(current_path,paste0("res.",i,".Rdata")))
    }

    res <- list(
      tpr_target_0=tpr_target_0,
      tpr_est_n=mean(tpr_est_n_vec),
      tpr_est_eff=mean(tpr_est_eff_vec),
      tpr_sd_n=mean(tpr_sd_n_vec),
      tpr_sd_eff=mean(tpr_sd_eff_vec),
      tpr_empirical_sd_n=sd(tpr_est_n_vec),
      tpr_empirical_sd_eff=sd(tpr_est_eff_vec),
      tpr_cov_n=mean(tpr_cov_n_vec),
      tpr_cov_eff=mean(tpr_cov_eff_vec),

      fpr_target_0=fpr_target_0,
      fpr_est_n=mean(fpr_est_n_vec),
      fpr_est_eff=mean(fpr_est_eff_vec),
      fpr_sd_n=mean(fpr_sd_n_vec),
      fpr_sd_eff=mean(fpr_sd_eff_vec),
      fpr_empirical_sd_n=sd(fpr_est_n_vec),
      fpr_empirical_sd_eff=sd(fpr_est_eff_vec),
      fpr_cov_n=mean(fpr_cov_n_vec),
      fpr_cov_eff=mean(fpr_cov_eff_vec),

      auc_target_0=auc_target_0,
      auc_est_n=mean(auc_est_n_vec),
      auc_est_eff=mean(auc_est_eff_vec),
      auc_sd_n=mean(auc_sd_n_vec),
      auc_sd_eff=mean(auc_sd_eff_vec),
      auc_empirical_sd_n=sd(auc_est_n_vec),
      auc_empirical_sd_eff=sd(auc_est_eff_vec),
      auc_cov_n=mean(auc_cov_n_vec),
      auc_cov_eff=mean(auc_cov_eff_vec)
    )
    save(res,file=file.path(current_path,"res.summary.Rdata"))
  }
}


# Simulations for lambda in {0.8,0.5,0.25,0.1,0.01}

efficient.rectifier.sim(nrep=2500,dat_gen,n=nvec[1],lambda=lambdavec[1],nfolds=nfolds,spline_df=spline_df)

efficient.rectifier.sim(nrep=2500,dat_gen,n=nvec[1],lambda=lambdavec[2],nfolds=nfolds,spline_df=spline_df)

efficient.rectifier.sim(nrep=2500,dat_gen,n=nvec[1],lambda=lambdavec[3],nfolds=nfolds,spline_df=spline_df)

efficient.rectifier.sim(nrep=2500,dat_gen,n=nvec[1],lambda=lambdavec[4],nfolds=nfolds,spline_df=spline_df)



