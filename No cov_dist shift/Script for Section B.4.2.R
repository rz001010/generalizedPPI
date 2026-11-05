# simulation code for calibrated PPI-AIPW mean estimation 

expit <- function(a){ return(exp(a)/(1+exp(a)))}

## mean
# point estimation
mean.est.yn <- function(dat){
  return(mean(dat$Y))
}
mean.est.fxn <- function(dat){
  return(mean(dat$f))
}
# influence function computation
mean.phi.yn <- function(dat){
  return(dat$Y-mean(dat$Y))
}
mean.phi.fxn <- function(dat){
  return(dat$f-mean(dat$f))
}


library(MASS)
library(randomForest)
library(ppiAIPW)
library(ggplot2)

nvec <- c(1000)

lambdavec <- c(0.8,0.5,0.25,0.1,0.01) # different lambda=n/N values

f_types <- c("all","bad","noise","ideal","pure_noise") # different prediction models f
# "all" is a random forest model predicting Y from (R,X_1,...,X_5)
# "bad" is a random forest model predicting Y from (R,X_1,X_4,X_5)
# "noise" is a random forest model predicting Y from (X_4,X_5), which are independent of Y
# "ideal" represents E[Y|R,X_1,...,X_5]
# "pure_noise" is Unif[0.01,0.99], uniformly distributed random noise


m <- 100000

calibration_method <- "monotone_spline"

results_dir <- "results/noncovshift_mean"
dir.create(results_dir,recursive=TRUE,showWarnings=FALSE)

set.seed(2025)

# The covariance matrix of (R,X_1,...,X_5)
Var_mat <- matrix(c(
  1.0, 0.5, 0.3, 0.2, 0, 0,
  0.5, 1.0, 0.4, 0.3, 0, 0,
  0.3, 0.4, 1.0, 0.2, 0, 0,
  0.2, 0.3, 0.2, 1.0, 0, 0,
  0,   0,   0,   0,   1, 0,
  0,   0,   0,   0,   0, 2
),nrow=6,byrow=TRUE)

# prespecified E[Y|R,X_1,...,X_5]
beta_intcpt <- -0.5
beta_R <- 1.0
beta_X1q <- -0.9   # X_1^2
beta_X2a <- 0.6    # |X_2|
beta_X3c <- 0.5    # X_3^3
beta_Rx3 <- 1.5    # R * X_3
beta_X1x2 <- -0.7  # X_1 * X_2

compute_eta_p <- function(df){
  eta <- with(df,
              beta_intcpt+
                beta_R*R+
                beta_X1q*X1^2+
                beta_X2a*abs(X2)+
                beta_X3c*X3^3+
                beta_Rx3*(R*X3)+
                beta_X1x2*(X1*X2)
  )
  plogis(eta)
}


# generate the training data for the three random forest prediction models
train_dat <- as.data.frame(mvrnorm(m,mu=c(0,0,0,0,3,20),Sigma=Var_mat))
names(train_dat) <- c("R","X1","X2","X3","X4","X5")
train_dat$ideal <- compute_eta_p(train_dat)
train_dat$Y <- rbinom(m,size=1,prob=train_dat$ideal)
train_dat$Y <- factor(train_dat$Y,levels=0:1)

# A random forest model predicting Y from (R,X_1,...,X_5)
rf_all <- randomForest(Y~R+X1+X2+X3+X4+X5,data=train_dat)
# A random forest model predicting Y from (X_4,X_5)
rf_noise <- randomForest(Y~X4+X5,data=train_dat)
# A random forest model predicting Y from (R,X_1,X_4,X_5)
rf_bad <- randomForest(Y~R+X1+X4+X5,data=train_dat)


# data generation function
dat_gen <- function(n,N,f_type){
  inf_dat <- as.data.frame(mvrnorm((n+N),mu=c(0,0,0,0,3,20),Sigma=Var_mat))
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


# simulation function
ppi.sim <- function(nrep,dat_gen,n,lambda){
  
  N <- n/lambda # the unlabeled dataset size
  
  # Use the same mean-target calculation and seeds as Section 6.1.
  target_rep <- 10000
  mean_target_0_vec <- vector("numeric",length=target_rep)
  
  for(index in 1:target_rep){
    set.seed(index+4564356)
    
    inf_dat.index <- as.data.frame(
      mvrnorm((n+N),mu=c(0,0,0,0,3,20),Sigma=Var_mat)
    )
    names(inf_dat.index) <- c("R","X1","X2","X3","X4","X5")
    
    inf_dat.index$ideal <- compute_eta_p(inf_dat.index)
    inf_dat.index$Y <- rbinom(
      (n+N),size=1,prob=inf_dat.index$ideal
    )
    
    mean_target_0_vec[index] <- mean.est.yn(inf_dat.index)
  }
  
  mean_target_0 <- mean(mean_target_0_vec)
  
  for(k in 1:length(f_types)){
    
    
    mean_est_n_vec <- vector("numeric",length=nrep)
   
    mean_est_ppi_vec <- vector("numeric",length=nrep)
    
    mean_est_ppipp_vec <- vector("numeric",length=nrep)
    
    mean_est_calppi_vec <- vector("numeric",length=nrep)
    
    
    mean_omega_vec <- vector("numeric",length=nrep)
    
    
    mean_sd_n_vec <- vector("numeric",length=nrep)
    mean_sd_ppi_vec <- vector("numeric",length=nrep)
    mean_sd_ppipp_vec <- vector("numeric",length=nrep)
    mean_sd_calppi_vec <- vector("numeric",length=nrep)
    
    
    mean_cov_n_vec <- vector("numeric",length=nrep)
    mean_cov_ppi_vec <- vector("numeric",length=nrep)
    mean_cov_ppipp_vec <- vector("numeric",length=nrep)
    mean_cov_calppi_vec <- vector("numeric",length=nrep)
    
    current_path <- file.path(
      results_dir,
      paste0("f_",f_types[k],"_n_",n,"_lambda_",lambda)
    )
    dir.create(current_path,recursive=TRUE,showWarnings=FALSE)
    
    
    for(i in 1:nrep){
      
      set.seed(i)
      dat.i <- dat_gen(n,N,f_types[k])
      
      # Y is observed in the first n observations
      dat_n.i <- dat.i[1:n,]
      # Only f is used from the remaining N observations
      dat_N.i <- dat.i[(n+1):(n+N),]
      
      mean_phi_n.i <- mean.phi.yn(dat_n.i)
      mean_est_n.i <- mean.est.yn(dat_n.i)
      
      # record the labeled-data estimator and its estimated standard error
      mean_est_n_vec[i] <- mean_est_n.i
      mean_sd_n_vec[i] <- sqrt(var(mean_phi_n.i))/sqrt(n)
      
      mean_est_f_all.i <- mean.est.fxn(dat.i)
      mean_est_f_n.i <- mean.est.fxn(dat_n.i)
      
      # Original PPI
      original_ppi_fit <- ppiAIPW::mean_inference(
        Y=dat_n.i$Y,
        Yhat=dat_n.i$f,
        Yhat_unlabeled=dat_N.i$f,
        method="aipw",
        efficiency_maximization=FALSE,
        alpha=0.05,
        inference="wald"
      )
      mean_est_ppi_vec[i] <- mean_est_n.i+(mean_est_f_all.i-mean_est_f_n.i)
      mean_sd_ppi_vec[i] <- as.numeric(original_ppi_fit$se)
      
      
      # PPI++
      ppipp_fit <- ppiAIPW::mean_inference(
        Y=dat_n.i$Y,
        Yhat=dat_n.i$f,
        Yhat_unlabeled=dat_N.i$f,
        method="aipw",
        efficiency_maximization=TRUE,
        alpha=0.05,
        inference="wald"
      )
      mean_omega_vec[i] <- as.numeric(ppipp_fit$efficiency_lambda)
      mean_est_ppipp_vec[i] <- mean_est_n.i+
        mean_omega_vec[i]*(mean_est_f_all.i-mean_est_f_n.i)
      mean_sd_ppipp_vec[i] <- as.numeric(ppipp_fit$se)
      
      # Calibrated PPI-AIPW
      calibrated_fit <- ppiAIPW::mean_inference(
        Y=dat_n.i$Y,
        Yhat=dat_n.i$f,
        Yhat_unlabeled=dat_N.i$f,
        method=calibration_method,
        efficiency_maximization=FALSE,
        alpha=0.05,
        inference="wald"
      )
      calibrated_f_n.i <- as.numeric(predict(calibrated_fit$calibrator,dat_n.i$f))
      calibrated_f_N.i <- as.numeric(predict(calibrated_fit$calibrator,dat_N.i$f))
      calibrated_f_all_mean.i <- (
        sum(calibrated_f_n.i)+sum(calibrated_f_N.i)
      )/(n+N)
      
      mean_est_calppi_vec[i] <- mean_est_n.i+
        (calibrated_f_all_mean.i-mean(calibrated_f_n.i))
      mean_sd_calppi_vec[i] <- as.numeric(calibrated_fit$se)
      
      
      mean_cov_n_vec[i] <- as.numeric(
        abs(mean_est_n_vec[i]-mean_target_0)<=1.96*mean_sd_n_vec[i]
      )
      mean_cov_ppi_vec[i] <- as.numeric(
        abs(mean_est_ppi_vec[i]-mean_target_0)<=1.96*mean_sd_ppi_vec[i]
      )
      mean_cov_ppipp_vec[i] <- as.numeric(
        abs(mean_est_ppipp_vec[i]-mean_target_0)<=1.96*mean_sd_ppipp_vec[i]
      )
      mean_cov_calppi_vec[i] <- as.numeric(
        abs(mean_est_calppi_vec[i]-mean_target_0)<=1.96*mean_sd_calppi_vec[i]
      )
      
      res.i <- list(
        mean_est_n_i=mean_est_n_vec[i],
        mean_est_ppi_i=mean_est_ppi_vec[i],
        mean_est_ppipp_i=mean_est_ppipp_vec[i],
        mean_est_calppi_i=mean_est_calppi_vec[i],
        
        mean_sd_n_i=mean_sd_n_vec[i],
        mean_sd_ppi_i=mean_sd_ppi_vec[i],
        mean_sd_ppipp_i=mean_sd_ppipp_vec[i],
        mean_sd_calppi_i=mean_sd_calppi_vec[i],
        
        mean_omega_i=mean_omega_vec[i]
      )
      save(
        res.i,
        file=file.path(current_path,paste0("res.",i,".Rdata"))
      )
    }
    
    # save the simulation results for one simulation scenario
    res <- list(
      mean_target_0=mean_target_0,
      
      mean_est_n=mean(mean_est_n_vec),
      mean_est_ppi=mean(mean_est_ppi_vec),
      mean_est_ppipp=mean(mean_est_ppipp_vec),
      mean_est_calppi=mean(mean_est_calppi_vec),
      
      mean_sd_n=mean(mean_sd_n_vec),
      mean_sd_ppi=mean(mean_sd_ppi_vec),
      mean_sd_ppipp=mean(mean_sd_ppipp_vec),
      mean_sd_calppi=mean(mean_sd_calppi_vec),
      
      mean_emp_sd_n=sd(mean_est_n_vec),
      mean_emp_sd_ppi=sd(mean_est_ppi_vec),
      mean_emp_sd_ppipp=sd(mean_est_ppipp_vec),
      mean_emp_sd_calppi=sd(mean_est_calppi_vec),
      
      mean_bias_n=mean(mean_est_n_vec-mean_target_0),
      mean_bias_ppi=mean(mean_est_ppi_vec-mean_target_0),
      mean_bias_ppipp=mean(mean_est_ppipp_vec-mean_target_0),
      mean_bias_calppi=mean(mean_est_calppi_vec-mean_target_0),
      
      mean_cov_n=mean(mean_cov_n_vec),
      mean_cov_ppi=mean(mean_cov_ppi_vec),
      mean_cov_ppipp=mean(mean_cov_ppipp_vec),
      mean_cov_calppi=mean(mean_cov_calppi_vec),
      
      mean_omega=mean(mean_omega_vec),
      nrep=nrep,
      n=n,
      N=N,
      lambda=lambda,
      f_type=f_types[k],
      calibration_method=calibration_method
    )
    save(res,file=file.path(current_path,"res.summary.Rdata"))
    
    
    summary_table <- data.frame(
      method=c("n","ppi","ppipp","calppi"),
      point_estimate=c(
        res$mean_est_n,res$mean_est_ppi,res$mean_est_ppipp,res$mean_est_calppi
      ),
      estimated_sd=c(
        res$mean_sd_n,res$mean_sd_ppi,res$mean_sd_ppipp,res$mean_sd_calppi
      ),
      empirical_sd=c(
        res$mean_emp_sd_n,res$mean_emp_sd_ppi,res$mean_emp_sd_ppipp,
        res$mean_emp_sd_calppi
      ),
      bias=c(
        res$mean_bias_n,res$mean_bias_ppi,res$mean_bias_ppipp,res$mean_bias_calppi
      ),
      coverage_95=c(
        res$mean_cov_n,res$mean_cov_ppi,res$mean_cov_ppipp,res$mean_cov_calppi
      ),
      stringsAsFactors=FALSE
    )
    write.csv(
      summary_table,file=file.path(current_path,"mean_summary.csv"),
      row.names=FALSE
    )
  }
}


# Simulations for lambda in {0.8,0.5,0.25,0.1,0.01}
ppi.sim(nrep=2500,dat_gen,n=nvec[1],lambda=lambdavec[1])
ppi.sim(nrep=2500,dat_gen,n=nvec[1],lambda=lambdavec[2])
ppi.sim(nrep=2500,dat_gen,n=nvec[1],lambda=lambdavec[3])
ppi.sim(nrep=2500,dat_gen,n=nvec[1],lambda=lambdavec[4])
ppi.sim(nrep=2500,dat_gen,n=nvec[1],lambda=lambdavec[5])


