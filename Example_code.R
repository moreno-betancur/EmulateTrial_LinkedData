#################################################################################
#                                                                               #
#  Example code for the paper                                                   #
#  "Emulating a trial of a nurse-home visiting program in the policy-relevant   #
#   target population using administrative linked data"                         #
#   by Moreno-Betancur et al. (submitted)                                       #
#                                                                               #
#  It shows how we ran a given method for a given outcome with a given approach #
#  to missingness on one (possibly imputed) dataset for the paper               #
#                                                                               #
#  Margarita Moreno-Betancur, 15 Dec 2020                                       #
#                                                                               #
#################################################################################

#### Load libraries and functions and set parameters ----

library(tmle)
library(SuperLearner)
library(boot)
source("funcs.R")

# Set outcome amonsgt "aedcdv1","aedcphy","aedcsoc","aedcemo","aedclan","aedccom", "napread","napnum","napwri","napspe","napgra"
  out<-"aedcdv1" 
# Set missingness approach amongst "CC", "MI", "MCIM", "CCov"
  miss<-"CC"     
# Set the imputed dataset number (if miss="MI")
  m<-100         
# Number of bootstrap runs
  R<-1000        
# Number ofcross-validation folds for SuperLearner in TMLE
  V<-10          

# Confounder names
MAIN<-c('yob','babysex','RA','irsad','matage','fatage','mothermarital','joblessfamily',
        'mothcob','babyoutc','smoking','antenataldic','preterm','hypertension','gestatdiabetes',
        'bwgaz','congenitalabnormBDR','DVeverCP',
        'everatsi','congenitalabnormFHV','drugalcoholFHV','mentalhealthFHV','dvFHV','poorattributionFHV',
        'cafhsarea','CP30days_indicator','housingSA_before1','max_educ_parent1','nbirths')

# Names of continuous covariate-missingness indicator terms for MCIM approach
MAINextra<-c("MI.fatage","MI.irsad")

# Algorithm library for SuperLearner for propensity score and missingness models
libq <-c("SL.mean","SL.glm","SL.step", "SL.glm.interaction",
         "SL.glmnet", "SL.xgboost",  "SL.bartMachine", "SL.randomForest",
         "SL.rpart","SL.stepAIC", "SL.bayesglm","SL.earth","SL.ipredbagg","SL.gbm")

# Algorithm library for SuperLearner for outcome model
libg <-libq[!libq%in%c("SL.randomForest","SL.bartMachine")]

RES<-data.frame()

outcome<-ifelse(out%in%c("aedcdv1","aedcphy","aedcsoc","aedcemo","aedclan","aedccom"), "AEDC", "NAPLAN")

#### Analyses according to outcome & missingness approach ----

if(miss=="CC")
{
  #Read in cleaned data - dataset not publicly available
  dat<-read.csv(paste("dat",outcome,"-clean.csv",sep=""))  
  
  #Remove any record with missingness in any covariate (anymissMAIN is an indicator of this) or outcome
  dat<-dat[!dat$anymissMAIN&!is.na(dat[,out]),]            
  
  Y<-dat[,out]
  A<-dat$trt
  W<-dat[,MAIN]
  
  # g-computation
  RES<-rbind(RES,runboot(Y,A,W,R=R,statgcomp))
  
  # IPW 
  RES<-rbind(RES,runboot(Y,A,W,R=R,statipw))
  
  # TMLE - complete outcome
  fit<-tmle(Y=Y,A=A,W=W,verbose=T,family="binomial",Q.SL.library =libq,g.SL.library =libg,V=V)
  RES<-rbind(RES, tmleres(fit))
  tmleother(fit,out, miss, m, R, V)
  
  methods<-c("gcomp","IPW","TMLE")
  
  
} else if(miss=="MI")
{  
  #Read in cleaned imputed data - dataset not publicly available
  dat<-read.csv(paste("MI datasets/dat",outcome,m,".csv",sep="")) 
  
  Y<-dat[,out]
  A<-dat$trt
  W<-dat[,MAIN]
  
  # g-computation
  RES<-rbind(RES,runboot(Y,A,W,R=R,statgcomp))
  
  # IPW 
  RES<-rbind(RES,runboot(Y,A,W,R=R,statipw))
  
  # TMLE - complete outcome
  fit<-tmle(Y=Y,A=A,W=W,verbose=T,family="binomial",Q.SL.library =libq,g.SL.library =libg,V=V)
  RES<-rbind(RES, tmleres(fit))  
  tmleother(fit,out, miss, m, R, V)
  
  methods<-c("gcomp","IPW","TMLE")
  
}  else if(miss=="MCIM")
{ 
  #Read in cleaned data prepared for MCIM - dataset not publicly available
  dat<-read.csv(paste("dat",outcome,"-MCIM.csv",sep=""))        
  
  Y<-dat[,out]
  Delta<-1*(!is.na(Y))
  A<-dat$trt
  W<-dat[,c(MAIN,MAINextra)]
  
  # TMLE - extended for missing outcome
  fit<-tmle(Y=Y,A=A,W=W,Delta=Delta,verbose=T,family="binomial",Q.SL.library =libq,g.SL.library =libg,V=V)
  RES<-rbind(RES, tmleres(fit))  
  tmleother(fit,out, miss, m, R, V)
  
  methods<-c("TMLE")
  
} else if(miss=="CCov")
{
  #Read in cleaned data - dataset not publicly available
  dat<-read.csv(paste("dat",outcome,"-clean.csv",sep=""))     
  
  #Remove any record with missingness in any covariate (anymissMAIN is an indicator of this) 
  dat<-dat[!dat$anymissMAIN,]                   
  
  Y<-dat[,out]
  Delta<-1*(!is.na(Y))
  A<-dat$trt
  W<-dat[,MAIN]
  
  # TMLE - extended for missing outcome
  fit<-tmle(Y=Y,A=A,W=W,Delta=Delta,verbose=T,family="binomial",Q.SL.library =libq,g.SL.library =libg,V=V)
  RES<-rbind(RES, tmleres(fit))  
  tmleother(fit,out, miss, m, R, V)
  
  methods<-c("TMLE")
  
}

RES<-cbind(rep(miss,nrow(RES)),rep(m,nrow(RES)),rep(R,nrow(RES)),rep(V,nrow(RES)),
           methods,RES)
names(RES)<-c("Missing","m","R","V","Method","ATT","SE","CIlow","CIupp","pvalue")

write.table(RES,paste("RES_",out,".csv"),append=T,row.names = F,col.names = F)



