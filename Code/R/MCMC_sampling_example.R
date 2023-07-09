#########################################################################################
#please update to your own working directory, set up the working directory
setwd("C:/Users/YYao/Downloads/Real Data")

set.seed(202307)
#set up the hyperparameters
nu0=2
sigmasq0=10
h=1.5
htau=2
mu=-log(1/0.5-1)
e<-2
f<-1
eta_sd=0.5
alpha0=0
tau0=0
beta0=0

Tmax=64
mu_tilde=-4.9
eta1_tilde=0.1
eta2_tilde=0.1
options(scipen=999)
niter=1000

file_address="C:/Users/YYao/Downloads/Real Data/"
#run MCMC sampling of our Bayesian variable selection method using the example data 
mc_test_p20_dat<-bayanal3_FBE(design=X,y=y,z_dat=matrix(Z_dat,ncol=1),seed=512,file_address=file_address,R=adj1_X)

#mean posterior inclusion probability of variable selection indicators of 20 variables for selection based on MCMC sampling
apply(mc_test_p20_dat[[5]][401:1000,],2, function(x){sum(x)/(niter-400)})
#the number of variables selected 
sum(apply(mc_test_p20_dat[[5]][401:1000,],2, function(x){sum(x)/(niter-400)})>0.5)
