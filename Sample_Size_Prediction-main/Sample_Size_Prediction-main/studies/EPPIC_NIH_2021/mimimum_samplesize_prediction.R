####---- Functions -------
#Functions for calculating mininum sample sizes for prediction. 
#Based on Riley's works

####################
# Function to calculate the minimum sample for binary response based on criterion:
# iii) precise estimation of the average outcome risk in the population, i.e., the prevalance
# Use the formula by given in Fig 1 of paper \cite{riley2020calculating}: 
# https://www.research.manchester.ac.uk/portal/files/161373531/bmj.m441.full.pdf
# Arguments: 
#      prevalence: Prevalence of outcome. 
#      delta: +/- proportion regarding the precise estimation.
# Example: 
# pmsampsize(type = "b", cstatistic =0.8, parameters = 30, prevalence=0.3);
# minSampleSize.aveRisk(prevalence=0.3, delta=0.05); #result is consistent with the above's "Criteria 3"
####################
minSampleSize.aveRisk = function(prevalence, delta=0.05) {
  zscore = qnorm(1-delta/2);  
  return(n = (zscore/delta)^2*prevalence*(1-prevalence));
}


####################
# Function to calculate the minimum sample for binary response based on criterion: 
# ii) small absolute difference in the model's apparent and adjusted Nagelkerke's R-squared value. 
# Based on 1) C-statistic or R2cs and maxR2cs; 
#          2) The formula by given in Fig 5 of paper \cite{riley2020calculating}: 
#              https://www.research.manchester.ac.uk/portal/files/161373531/bmj.m441.full.pdf
# Depend: approximate_R2 function for calculating R2cs if C-statistic is used. 
# Arguments: 
#     P: The number of parameters
# Cstat: C-statistic, i.e., the AUC. It is NULL if R2cs and maxR2cs are given. 
# prevalence: Prevalence of outcome. Required if calculation is based on Cstat 
# R2cs: Cox and Snell R2. It is NULL if Cstat is given. 
# maxR2cs: Maximum of Cox and Snell R2. It is NULL if Cstat is given. 
# delta: Small absolute difference in the model's apparent and adjusted Nagelkerke's R-squared value
# n.approx: sample sized used in approximating R2cs by Cstat
# Example: 
# P = 20;
# delta=0.05;
# R2cs = 0.2; #Cox and Snell R2;
# maxR2cs = 0.33; # maximum of Cox and Snell R2;
# minSampleSize.NagR2(P, R2cs=R2cs, maxR2cs=maxR2cs);
# 
# pmsampsize(type = "b", cstatistic =0.8, parameters = 30, prevalence=0.3, shrinkage=0.8, mmoe=1.2)
# minSampleSize.NagR2(P=30, Cstat=0.8, prevalence=0.3); #result is consistent with the above's "Criteria 2"
# minSampleSize.NagR2(P=30, Cstat=0.8, prevalence=0.3, delta=0.1); #Smaller n due to relaxed delta from 0.05 to 0.1. 
####################
minSampleSize.NagR2 = function(P, Cstat=NULL, prevalence=NULL, R2cs=NULL, maxR2cs=NULL, delta=0.05, n.approx=1000000) {
  if (!is.null(Cstat)){
    R2 = approximate_R2(auc = Cstat, prev = prevalence, n=n.approx);
    R2cs = R2$R2.coxsnell; 
    maxR2cs = R2cs/R2$R2.nagelkerke;
  }
  S = R2cs/(R2cs + delta*maxR2cs);
  return(n = P/((S-1)*log(1 - R2cs/S)));
}


####################
# Function for calculating Cox and Snell R2 based on C-statistic and prevalence. 
# Code copied from the paper: https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8806
# Depend: lrm function in the R library rms for the logit model fitting.  
# Arguments: auc: area under the ROC curve, i.e., the C-statistic
#            prev: prevalence of outcome in the population of the data.
# Example: 
# install.packages("rms");
# library(rms);
# set.seed(1234);
# approximate_R2(auc = 0.81, prev = 0.77, n=1000000);
# $R2.nagelkerke
# [1] 0.3183689
# $R2.coxsnell
# [1] 0.2100957
####################
approximate_R2 = function(auc, prev, n = 1000000){
  # define mu as a function of the C statistic
  mu = sqrt(2) * qnorm(auc);
  # simulate large sample linear prediction based on two normals for nonâ eventsN(0, 1), events and N(mu, 1)
  LP = c(rnorm(prev*n,  mean=0, sd=1), rnorm((1-prev)*n, mean=mu, sd=1));
  y = c(rep(0, prev*n), rep(1, (1-prev)*n));
  # Fit a logistic regression with LP as covariate;
  # this is essentially a calibration model, and the intercept and
  # slope estimate will ensure the outcome proportion is accounted
  # for, without changing C statistic
  fit = lrm(y ~ LP);
  max_R2 = function(prev){
    1-(prev^prev*(1-prev)^(1-prev))^2;
  }
  
  return(list(R2.nagelkerke = as.numeric(fit$stats['R2']),
              R2.coxsnell = as.numeric(fit$stats['R2']) * max_R2(prev)));
}
