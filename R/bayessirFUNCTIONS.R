library(mvtnorm)

enviforce=function(covs, obsDates, alphaPars){
  numCovs = dim(covs)[2]
  numDays = dim(covs)[1]
  
  firstday = obsDates[1]
  
  TempCovs=matrix(nrow=numDays,ncol=numCovs)
  
  for(i in 1:numCovs){
    TempCovs[,i]=covs[,i]*alphaPars[i+1]
  }
  
  allalphas<-c()
  
  for(j in 1:numDays){
    allalphas[j]=exp(alphaPars[1]+sum(TempCovs[j,]));
  }
  
  alphasWithBreaks<-list()
  for(j in 1:(length(obsDates)-1)){
    alphasWithBreaks[[j]] = cbind(obsDates[j]:(obsDates[j+1]-1),
                                  allalphas[(obsDates[j]-firstday+1):(obsDates[j+1]-1-firstday+1)])
  }
  return(alphasWithBreaks)
  
}

prior.fun<-function(param,priorMean,priorSD){
#   priorMean=priorMean[-c(which(is.na(param)))]
#   priorSD=priorSD[-c(which(is.na(param)))]
#   param=param[-c(which(is.na(param)))]
  
  priors<-rep(NA,length(param))
  
  for(i in 1:length(param)){
    priors[i]<-dnorm(param[i], mean=priorMean[i], sd=priorSD[i], log = T) 
  }
  return(sum(priors))
}

logit=function(p){log(p/(1-p))}           
expit=function(p){exp(p)/(1+exp(p))}

