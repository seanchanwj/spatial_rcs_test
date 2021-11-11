#####################################
## Set up and running Spatial NetRate
#####################################

library(VGAM)
library(tensorflow)
library(reticulate)
library(scales)

#set working directory to directory where the function and simulated data are stored
setwd( "INSERTPATH" )

# source function to generate RC
source("function_genRC.R")

# Set up miniconda environment and download tensorflow and tensorflow probability (old versions for compatibility with code)
# Note need to change paths 

reticulate::install_miniconda(path = "../AppData/Local/r-miniconda/", update = TRUE, force = FALSE)

Sys.setenv(RETICULATE_PYTHON="../AppData/Local/r-miniconda/envs/tensorflow")
tensorflow::install_tensorflow(version = "1.14",
                               extra_packages = c("tensorflow-probability==0.7" "scipy==1.5.0")  ,
                               conda_python_version = "3.6",
                               envname              = "tensorflow")
conda_list()
use_condaenv("tensorflow", required = TRUE)
tfp <- reticulate::import("tensorflow_probability",convert=FALSE)
tf_config() 

# Load simulated data

res<-read.csv('res.csv')

# make matrices and reformat for function  
dd<-list(n=nrow(res), I= c(rep(1,length.out = n_seed), rep(0, length.out = nrow(res)-n_seed)),t=res$inf_times, d=res$inf_dist)
tmat<-dmat<-matrix(0,nrow=dd$n,ncol=dd$n)

# shift for 15 day minimum serial interval 
for (i in 1:dd$n) { # for infectees
  for(j in 1:(dd$n)){ # for all cases/infectors
    tmat[i,j]=dd$t[i]-dd$t[j]-15 #otherwise calculate time difference		
    
  }
}


tmat[tmat<0]=0 # anything less than zero is in the past to set time difference to zero
tmat = tmat[-which(dd$I==1),] # imported cases cant be infectees and so remove these

# calculate distance in meters
for (i in 1:dd$n) { # for infectees
  for(j in 1:(dd$n)){ # for all cases/infectors
    dmat[i,j]= abs(dd$d[i]-dd$d[j]) #otherwise calculate time difference		
    
  }
}

dmat<-dmat/1000 # convert to kilometers
dmat = dmat[-which(dd$I==1),] # imported cases cant be infectees and so remove these

sd_mid<-spatialnetrate(tp=tmat, dp= dmat, fixed = "epsilon", alpha = c(0.002, 0.001), delta=c(0.01,0.001), SpatialKernel = "exponential", epsilon = 1e-20)

# inspect results 

# distribution of Rc
hist(sd_mid[[1]])

# delta parameter
sd_mid[[2]]

# alpha parameter 
hist(sd_mid[[3]])

#epsilon
sd_mid[[4]]
# AIC
sd_mid[[5]]

