#eta1/eta_sd,eta2/eta_sd follows Beta(e,f) prior distribution
#one adjacent matrix are estimated from full Bayesian MCMC sampling of 
#sparse invser covariance and adjacent matrices assuming Gaussian Graphical model
#bayanal3_FBE will take input of hyperparameters for the parameters in model setup
#mu for gamma, beta0 for beta,alpha0 for alpha
#nu0,sigmasq0 for sigmasqbeta
#also input of the parameters-eta1, eta2 in variant MRF prior for gamma
#adjacent matrices-R and R2 from differnt sources of information
#typically one from correlation matrix of the data itself
#another from biological information such as KEGG pathway database
#MCMC update
bayanal3_FBE <- function(design,y,z_dat,R,seed,file_address) {
  #randomly generate eta1,eta2 through auxiliary variable-omega_eta1,omega_eta2
  library(Rcpp)
  library(BH)
  library(RcppArmadillo)
  #  library(RcppEigen)
  library(MCMCpack)

  
  
  
  p=ncol(design)
  #Sig=1.5*diag(p)
  R2=1.5*diag(p)
  
  sourceCpp("MollerRand2eta_random.cpp")
  ## Initialize
  nbeta<-dim(R)[1]
  ntau<-ncol(z_dat)
  omega <- exp(mu)/(1+exp(mu))
  gamma <- rbinom(nbeta, 1, omega)
  beta <- rep(0, nbeta)
  sigmasqbeta<-1/rgamma(1,nu0/2,sigmasq0*nu0/2)
  sigmabeta<-sqrt(sigmasqbeta)
  sigmaalpha<-sqrt(h*sigmasqbeta)
  sigmasqtau<-htau*sigmasqbeta
  sigmatau<-sqrt(sigmasqtau)
  tau<-rnorm(ntau,tau0,sigmatau)
  alpha<-rnorm(1,alpha0,sigmaalpha)
  eta1=eta_sd*rbeta(1,e,f)
  eta2=eta_sd*rbeta(1,e,f)
  beta[gamma==1] <- rnorm(sum(gamma), beta0,sigmabeta)
  #calculate the log posterior marginal probability at the initiazation of the parameters
  select <- which(gamma==1)
  select.zero<-which(gamma==0)
  nselect <- length(select)
  betasub <- beta[select]
  designsub <- design[, select, drop=F]
  z<-apply(designsub,1,function(x){return(x%*%betasub)})
  
  tau_zdat<-z_dat%*%tau
  tau_z_tot<-alpha+z+tau_zdat
  loglik <- sum(y*(tau_z_tot))-sum(log(1+exp(tau_z_tot)))
  logtau <- -sum((tau-tau0)^2)/(2*sigmasqtau) -0.5*ntau*log(2*pi*sigmasqtau)
  
  logalpha<--(alpha-alpha0)^2/(2*h*sigmasqbeta)-0.5*log(2*pi*h*sigmasqbeta)
  logbeta <- -sum((betasub-beta0)^2)/(2*sigmasqbeta) -0.5*nselect*log(2*pi*sigmasqbeta)
  loggamma <- c(mu*nselect + t(gamma)%*%(eta1*R+eta2*R2)%*%gamma)
  logsigmasqbeta<-log(dinvgamma(sigmasqbeta,nu0/2,sigmasq0*nu0/2))
  ## log posterior marginal joint distribution of alpha,beta,gamma
  logp=loglik + logbeta + loggamma+logalpha+logsigmasqbeta+logtau
  #set the parameter vectors for MCMC interations
  betaMC <- gammaMC <- matrix(NA, nrow=niter, ncol=nbeta)
  alphaMC<-matrix(NA,nrow=niter,ncol=1)
  tauMC<-matrix(NA,nrow=niter,ncol=ntau)
  
  sigmasqbetaMC<-matrix(NA,nrow=niter,ncol=1)
  for (i in 1:niter) {
  #read the MCMC samples of adjacent matrix generated from matlab code and make sure the file names are consistent
    R2.dat<-read.table(paste(file_address,'met_p_',as.character(seed),"_cov_",as.character(i),".txt",sep="")) 
    
    
    
    R2=as.matrix(R2.dat)
    
    diag(R2)<-0
    
    for (k in 1:3) {
      ## update main effects
      uid <- sample(1:nbeta, 1)
      gamma.prop <- gamma
      gamma.prop[uid] <- 1 - gamma.prop[uid]
      beta.prop <-  beta
      beta.prop[uid]<-ifelse(gamma.prop[uid]==0,0,rnorm(1, beta0, sigmabeta))
      nselect.prop<-ifelse(gamma.prop[uid]==1,nselect+1,nselect-1)
      #difference at beta
      logbeta.prop<-ifelse(gamma.prop[uid]==1,logbeta+log(dnorm(beta.prop[uid],beta0,sigmabeta)),
                           logbeta-log(dnorm(beta[uid],beta0,sigmabeta)))
      #differnece at likelihood
      if(gamma.prop[uid]==1){
        z.prop<-z+beta.prop[uid]*design[,uid]
      } else {
        z.prop<-z-beta[uid]*design[,uid]}
      
      
      tau_z_tot.prop<-alpha+z.prop+tau_zdat
      loglik.prop <- sum(y*(tau_z_tot.prop))-sum(log(1+exp(tau_z_tot.prop)))
      
      
      #difference at gamma
      loggamma.prop<-ifelse(gamma.prop[uid]==1,
                            loggamma+mu+sum(gamma.prop*(eta1*R[,uid]+eta2*R2[,uid]))+sum(gamma.prop*(eta1*R[uid,]+eta2*R2[uid,])),
                            loggamma-mu-sum(gamma*(eta1*R[,uid]+eta2*R2[,uid]))-sum(gamma*(eta1*R[uid,]+eta2*R2[uid,])))
      #difference at log-posterior marginal probability
      logp.prop<-logbeta.prop+ loglik.prop +loggamma.prop+logalpha+logsigmasqbeta+logtau
      AAA <- logp.prop - logp
      +ifelse(gamma.prop[uid]==0, log(dnorm(beta[uid], beta0,sigmabeta)), -log(dnorm(beta.prop[uid], beta0,sigmabeta)))
      if (log(runif(1)) < AAA) {
        gamma <- gamma.prop
        beta <- beta.prop
        logp <- logp.prop
        z<-z.prop
        tau_z_tot<- tau_z_tot.prop
        loglik<-loglik.prop
        logbeta<-logbeta.prop
        loggamma<-loggamma.prop
        nselect<-nselect.prop
      }
      ##################################################################################################
      include <- which(gamma==1)
      ninclude <- length(include)
      ## update regressional coefficients, update them one by one
      for (j in include) {
        beta.prop <- beta
        beta.prop[j] <- rnorm(1, beta[j], sigmabeta)
        #difference at beta
        logbeta.prop<-logbeta-(beta.prop[j]-beta0)^2/(2*sigmasqbeta)+(beta[j]-beta0)^2/(2*sigmasqbeta)
        #differnece at likelihood
        z.prop<-z+(beta.prop[j]-beta[j])*design[,j,drop=F]
        
        tau_z_tot.prop<-alpha+z.prop+tau_zdat
        loglik.prop <- sum(y*(tau_z_tot.prop))-sum(log(1+exp(tau_z_tot.prop)))
        
        #difference at log-posterior marginal probability
        logp.prop<-logbeta.prop+loglik.prop+loggamma+logalpha+logsigmasqbeta+logtau
        AAA <- logp.prop - logp
        if (log(runif(1)) < AAA) {
          beta <- beta.prop
          logp <- logp.prop
          logbeta<-logbeta.prop
          loglik<-loglik.prop
          z<-z.prop
          tau_z_tot<-tau_z_tot.prop
        }
        
      }
      
    }
    #update eta1, eta2
    list_test2<-moller_2eta_trnormal(R, R2,  mu, eta1, eta2,Tmax,
                                     mu_tilde,eta1_tilde,eta2_tilde,eta_sd,gamma,e,f)
    eta1<-unlist(list_test2["eta1"])
    eta2<-unlist(list_test2["eta2"])
    loggamma= c(mu*sum(gamma==1) + t(gamma)%*%(eta1*R+eta2*R2)%*%gamma)
    logp=loglik + logbeta + loggamma+logalpha+logsigmasqbeta+logtau
    
    
    ## update alpha
    alpha.prop <- rnorm(1,alpha,sigmaalpha)
    
    tau_z_tot.prop<-alpha.prop+z+tau_zdat
    loglik.prop <- sum(y*(tau_z_tot.prop))-sum(log(1+exp(tau_z_tot.prop)))
    
    
    
    
    logalpha.prop<-logalpha-(alpha.prop-alpha0)^2/(2*h*sigmasqbeta)+(alpha-alpha0)^2/(2*h*sigmasqbeta)
    logp.prop<-loglik.prop+logalpha.prop+loggamma+logbeta+logsigmasqbeta+logtau
    AAA <- logp.prop - logp
    if (log(runif(1)) < AAA) {
      alpha <- alpha.prop
      logp <- logp.prop
      loglik<-loglik.prop
      logalpha<-logalpha.prop
      tau_z_tot<-tau_z_tot.prop
      
    }
    
    ## update tau
    # alpha.prop <- rnorm(1,alpha,sigmaalpha)
    tau.prop<-rnorm(ntau,tau,sigmatau)
    tau_zdat.prop<-z_dat%*%tau.prop
    tau_z_tot.prop<-alpha+z+tau_zdat.prop
    loglik.prop <- sum(y*(tau_z_tot.prop))-sum(log(1+exp(tau_z_tot.prop)))
    #logalpha.prop<-logalpha-(alpha.prop-alpha0)^2/(2*h*sigmasqbeta)+(alpha-alpha0)^2/(2*h*sigmasqbeta)
    
    
    logtau.prop <- logtau+sum((tau-tau0)^2)/(2*sigmasqtau)-sum((tau.prop-tau0)^2)/(2*sigmasqtau) 
    
    
    
    logp.prop<-loglik.prop+logalpha+loggamma+logbeta+logsigmasqbeta+logtau.prop
    
    
    AAA <- logp.prop - logp
    if (log(runif(1)) < AAA) {
      # alpha <- alpha.prop
      tau<-tau.prop
      tau_zdat<-tau_zdat.prop
      tau_z_tot<-tau_z_tot.prop
      logp <- logp.prop
      loglik<-loglik.prop
      # logalpha<-logalpha.prop
      logtau<-logtau.prop
      
    }
    
    
    
    
    
    
    
    
    
    ##update sigmasqbeta
    sigmasqbeta.prop<-exp(-rnorm(1,-log(sigmasqbeta),1))
    logsigmasqbeta.prop=log(dinvgamma(sigmasqbeta.prop,nu0/2,sigmasq0*nu0/2))
    
    
    
    logalpha.prop<--(alpha-alpha0)^2/(2*h*sigmasqbeta.prop)-0.5*log(2*pi*h*sigmasqbeta.prop)
    
    sigmasqtau.prop<-htau*sigmasqbeta.prop
    sigmatau.prop<-sqrt(sigmasqtau.prop)
    
    
    logtau.prop<- -sum((tau-tau0)^2)/(2*sigmasqtau.prop) -0.5*ntau*log(2*pi*sigmasqtau.prop)
    
    
    
    betasub<- beta[which(gamma==1)]
    logbeta.prop <- -sum((betasub-beta0)^2)/(2*sigmasqbeta.prop) -0.5*nselect*log(2*pi*sigmasqbeta.prop)
    logp.prop<-logbeta.prop+loglik+loggamma+logalpha.prop+logsigmasqbeta.prop+logtau.prop
    AAA <- logp.prop - logp
    if (log(runif(1)) < AAA) {
      sigmasqbeta <- sigmasqbeta.prop
      sigmabeta<-sqrt(sigmasqbeta.prop)
      sigmaalpha<-sqrt(h*sigmasqbeta.prop)
      
      sigmasqtau<-sigmasqtau.prop
      sigmatau<-sigmatau.prop
      logp <- logp.prop
      logalpha<-logalpha.prop
      logbeta<-logbeta.prop
      logsigmasqbeta<- logsigmasqbeta.prop
      logtau<-logtau.prop
      
    }
    
    
    
    ## Store chains
    betaMC[i, ] <- beta
    gammaMC[i, ] <- gamma
    alphaMC[i,]<-alpha
    sigmasqbetaMC[i,]<-sigmasqbeta
    
  }
  parms <- list(design=design,y=y, mu=mu, nu0=nu0,sigmasq0=sigmasq0,e=e,f=f,eta_sd=eta_sd,alpha0=alpha0,tau0=tau0,beta0=beta0,R=R,h=h,htau=htau, niter=niter)
  list(parms=parms, betaMC=betaMC, alphaMC=alphaMC,tauMC=tauMC, gammaMC=gammaMC,sigmasqbetaMC=sigmasqbetaMC)
}
