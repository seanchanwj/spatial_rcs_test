#####################################
## Set up and running Spatial NetRate
#####################################

library(VGAM)
library(tensorflow)
library(reticulate)
library(scales)
library(tidyverse)
library(lubridate)

#set working directory to directory where the function and simulated data are stored
setwd( "/Users/seanchanwj/Library/CloudStorage/OneDrive-Personal/MSc_PhD Public Health/Research/Malaria/Spatial_NetRateMalBhu|_test" )

# source function to generate RC
source("function_genRC.R")

# Set up miniconda environment and download tensorflow and tensorflow probability (old versions for compatibility with code)
# Note need to change paths 

reticulate::install_miniconda(path = "../AppData/Local/r-miniconda/", update = TRUE, force = FALSE)

Sys.setenv(RETICULATE_PYTHON="../AppData/Local/r-miniconda/envs/tensorflow")
tensorflow::install_tensorflow(version = "1.14",
                               extra_packages = c("tensorflow-probability==0.7", "scipy==1.5.0")  ,
                               conda_python_version = "3.6",
                               envname              = "tensorflow")
conda_list()
use_condaenv("tensorflow", required = TRUE)
tfp <- reticulate::import("tensorflow_probability",convert=FALSE)

#For MacOS, try using the following installation: 

  #Install Miniconda: https://docs.conda.io/en/main/miniconda.html
  #Open Terminal
  #Create Conda env: conda create --name rstudio-tf-2.8
  #Activate env: conda activate rstudio-tf-2.8
  #Downgrade Python: conda install python=3.8
  #Install tensorflow-deps: conda install -c apple tensorflow-deps==2.8.0
  #Install tensorflow-macos: python -m pip install tensorflow-macos==2.8.0
  #Install tensorflow-metal: python -m pip install tensorflow-metal
  #Install tensorflow: conda install -c anaconda tensorflow
  
  ##(Depending on versions)
  

  #Edit ~/.Renviron: echo "RETICULATE_PYTHON=~/miniforge3/envs/rstudio-tf-2.8/bin/python" >> ~/.Renviron
  #Check installed packages: conda list | grep tensorflow

  #Go to R-Studio and restart R-Session Session -> Restart R
  #Init and test tensorflow with following R script:
  #reticulate::use_python("~/miniconda3/envs/rstudio-tf-2.8/bin/python")
  ##reticulate::use_condaenv("rstudio-tf-2.8", conda="~/miniconda3/bin/conda", required = TRUE)


reticulate::use_python("~/miniconda3/envs/rstudio-tf-2.8/bin/python")
reticulate::use_condaenv("rstudio-tf-2.8", conda="~/miniconda3/bin/conda", required = TRUE)

  #Run Script: 
library(tensorflow)
tf_config() 
tf_version() 
tf$config$list_logical_devices()
#If tf_config() crashes go back to Terminal and upgrade to python 3.9 conda install python=3.9

  #Check again if all deps are listed conda list | grep tensorflow, for me tensorflow-macos and metal were missing after the python upgrade, I had to repeat step 8 and 9
  #Restart 
  #If numpy has issues - to install numpy: conda install -c anaconda numpy


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

