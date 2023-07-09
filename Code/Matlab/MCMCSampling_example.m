%set up the random seed 
seed=512;
% load observation data for the predictors 
X1 = load('/MATLAB Drive/data/X_dat.txt');

%covariance for X1
S1=X1'*X1;
%covariance for X2
%S2=X2*X2';

%other hyperparameters setup in this simulation scenario with p=100

%% Fix some hyperparameters
  h = 50^2; v0 = 0.015^2; v1 = h*v0;
  p=20;
  n=100;
  lambda = 1; pii = 6/(p-1); 
  V0 = v0*ones(p); V1 = v1*ones(p); burnin = 100; nmc =1000;
  rng(seed);
 % output the posterior samples of adjacent matrix in txt format
 Z_save1 = BayesGGM_SSVS_FixedV0V1_mod2(S1,n,eye(p),V0,V1,lambda,pii,burnin,nmc,seed);
 
