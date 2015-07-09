#' @title PMMH algorithm for time-varying SIRS model
#' @param obscholLIST List containing the cholera counts for each phase of data 
#' @param obsdaysLIST List containing the observation days for each phase of data 
#' @param COVMATLIST List containing the matricies of daily covariates for each phase of data 
#' @param th Starting value for parameter vales
#' @param trans function to transform the parameter values 
#' @param invtrans inverse transformation function 
#' @param pmean Prior means for Normal prior distributions
#' @param psd Prior standard deviations for Normal prior distributions
#' @param betaIIndex Index of which \code{th} values correspond to \eqn{\beta_I} 
#' @param betaWIndex Index of which \code{th} values correspond to \eqn{\beta_W} 
#' @param gammaIndex Index of which \code{th} values correspond to \eqn{\gamma} 
#' @param kappaIndex Index of which \code{th} values correspond to \eqn{\kappa} 
#' @param etaIndex Index of which \code{th} values correspond to \eqn{\eta} 
#' @param muIndex Index of which \code{th} values correspond to \eqn{\mu} 
#' @param alphasIndex Index of which \code{th} values correspond to \eqn{\alpha} parameters
#' @param rhoIndex Vector of length equal to number of phases of data, Index of which \code{th} values correspond to \eqn{\rho}, the probability of infected individuals seeking treatment
#' @param startmeansIndex Index of which \code{th} values correspond to means of initial distributions for the numbers of susceptible and infected individuals respectively
#' @param nu1Index Index of which \code{th} value corresponds to \eqn{\nu_1} power
#' @param nu2Index Index of which \code{th} value corresponds to \eqn{\nu_2} power
#' @param nu3Index Index of which \code{th} value corresponds to \eqn{\nu_3} power
#' @param nu4Index Index of which \code{th} value corresponds to \eqn{\nu_4} power
#' @param burn Number of iterations for burn-in run
#' @param prelim Number of total iterations for preliminary run, preliminary run = burn-in run + secondary run
#' @param iters Number of iterations for final run
#' @param thin Amount to thin the chain; only every \code{thin}th iteration of the PMMH algorithm is saved
#' @param tune Tuning parameter for the covariance of the multivariate normal proposal distribution in the final run of the PMMH algorithm
#' @param ll Starting value for log-likelihood
#' @param psigma Standard deviation for independent normal proposal distribution in preliminary run
#' @param deltavalue Initial value to use for tau in tau-leaping algorithm
#' @param critical Critical size for modified tau-leaping algorithm; if the population of a compartment is lower than this number a single step algorithm is used until the population gets above the critical size.
#' @param PopSize Population size
#' @param theMU The rate at which immunity is lost, if setting this value
#' @param numParticles Number of particles 
#' @param resultspath File path for results
#' @param UseGill boolian; if 1, uses the gillespie algorithm. If 0, uses tau-leaping algorithm
#' @param UseSIWR boolian; if 1, uses the SIWR model. If 0, uses SIRS model
#' @param setBetaW boolian; if 1, sets BetaW parameter. If 0, estimates BetaW parameter.
#' @param setKappa boolian; if 1, sets Kappa parameter. If 0, estimates Kappa parameter.
#' @param setEta boolian; if 1, sets Eta parameter. If 0, estimates Eta parameter.
#' @param setRatio boolian; if 1, sets ratio of kappa and eta. 
#' @param setval 
#' @param setAlpha0 boolian; if 1, sets alpha0 parameter. If 0, estimates alpha0 parameter.
#' @param setAlpha0val
#' @param setstartmeansval
#' @param usetprior
#' @param alphadf
#' @param uselaplaceprior
#' @param maxWval Upper bound for W compartment. If 0, no bounding.
#' @return Posterior samples from the final run of the PMMH algorithm
#' 
#' Also, writes 4 files which are updated every 100th iteration: 
#' 
#' 1. prelimpmcmctimes.csv: times and acceptance ratios for preliminary PMMH run
#' 
#' 2. prelimthmat.csv: preliminary PMMH output
#' 
#' 3. FINALpmcmctimes: times and acceptance ratios for final PMMH run
#' 
#' 4. FINALthmat.csv: final PMMH output
#' @examples
#'
#' \dontrun{
#' 
#' library(bayessir)
#' ############################
#' ## simulate data
#' ############################
#' 
#' SimTimes=seq(0,365*4.5, by=14) 
#' 
#' ############################
#' # environmental force of infection
#' ############################
#' 
#' int<- -6
#' A<-2
#' sincovAmp<- c(2.1,1.8,2,2.2,2)
#' wave<-pi/(365/2)
#' t<-0:(max(SimTimes))
#' 
#' sincov<-sin(wave*t)
#' allsincov=matrix(NA,nrow=length(t),ncol=length(sincovAmp))
#' for(i in 1:length(sincovAmp)){
#'   allsincov[,i]<-sincovAmp[i]*sincov
#' }
#' 
#' sincov[1:365]=allsincov[1:365,1]
#' sincov[366:(365*2)]=allsincov[366:(365*2),2]
#' sincov[(365*2+1):(365*3)]=allsincov[(365*2+1):(365*3),3]
#' sincov[(365*3+1):(365*4)]=allsincov[(365*3+1):(365*4),4]
#' sincov[(365*4+1):length(sincov)]=allsincov[(365*4+1):length(sincov),5]
#' 
#' alpha<-exp(int+A*sincov)
#' 
#' ###########
#' pop=10000 #population size
#' 
#' phiS=2900
#' phiI=84
#' 
#' th1=.5/10000 #beta
#' th2=0.12 #gamma
#' th3=.0018 #mu
#' rho=90/10000 #reporting rate
#' 
#' nu1=1
#' nu2=1
#' nu3=0
#' nu4=0
#' 
#' set.seed(10)
#' sus0=rpois(1,phiS)
#' inf0=rpois(1,phiI)
#' shortstart<-as.matrix(c(sus0,inf0))
#' 
#' allcovs<-enviforce(as.matrix(c(sincov)),SimTimes,c(int,A))
#' 
#' simstates<-matrix(NA,nrow=(length(SimTimes)),ncol=2)
#' simstates[1,]<-shortstart
#' 
#' for (i in 2:length(SimTimes)){
#'   simstates[i,]<-inhomoSIRSGillespie(simstates[i-1,],pop,SimTimes[i-1],SimTimes[i]-SimTimes[i-1],
#'                                      c(th1,th2,th3,nu1,nu2,nu3,nu4),allcovs[[i-1]][,2],allcovs[[i-1]][,1])
#' }
#' set.seed(9)
#' SimData<-c()
#' for(i in 1:dim(simstates)[1]) SimData[i]<-rbinom(1,simstates[i,2],rho)
#' 
#' ##################################################
#' 
#' #####
#' # data for inference
#' #####
#' COVMATLIST=list(as.matrix(sincov))
#' obscholLIST=list(SimData)
#' obsdaysLIST=list(SimTimes)
#' 
#' numofcovs=1
#' #################################################
#' trans=function(p){
#'   c(log(p[1]),  #beta
#'     log(p[2]),  #gamma
#'     log(p[3]),  #mu
#'     p[4],       #alpha0
#'     p[5],       #alpha1
#'     logit(p[6]))#rho
#' }
#' 
#' invtrans=function(p){
#'   c(exp(p[1]),  #beta 
#'     exp(p[2]),  #gamma
#'     exp(p[3]),  #mu
#'     p[4],       #alpha0
#'     p[5],       #alpha1
#'     expit(p[6]))#rho
#' }
#' 
#' #prior means
#' pbetaI=log(1.25e-04) 
#' pgamma=log(.1)
#' pmu=log(.0009)
#' palpha0=-8
#' palphas=rep(0,numofcovs)
#' prho=logit(.03) 
#' 
#' pmean=c(pbetaI,pgamma,pmu,palpha0,palphas,prho)
#' 
#' #prior standard deviations
#' psd=c(5,  #beta
#'       .09,#gamma 
#'       .3, #mu
#'       5,  #alpha0
#'       5,  #alpha1
#'       2)  #rho
#' 
#' betaIIndex=1 #need one for each phase of data collection, we only simulated one phase
#' gammaIndex=2
#' muIndex=3
#' alphasIndex=4:5
#' rhoIndex=6
#' 
# We aren't estimating these in this example
#' startmeansIndex=nu1Index=nu2Index=nu3Index=nu4Index=NA
#' 
# These are parameters for the SIWR model
#' betaWIndex=kappaIndex=etaIndex=NA
#' 
#' # Iterations set small for example purposes; increase for applications
#' burn = 0
#' prelim = 10
#' iters =10
#' thin =1
#'  
#' tune=1
#' 
#' psigma<-diag(c(0.012, #beta
#'             0.012, #gamma
#'             0.180, #mu
#'             0.120, #alpha0
#'             0.120, #alpha1
#'             0.012)) #rho
#'             
#' 
#' #start values
#' #Names of th input are used for the column names in the matrix output
#' th=c(
#'   betaI=abs(rnorm(1,th1,th1/3)), 
#'   gamma=abs(rnorm(1,th2,th2/10)),
#'   mu=abs(rnorm(1,th3,th3/10)),
#'   alpha0=rnorm(1,int,1),
#'   alpha1=rnorm(1,A,1),
#'   rho=abs(rnorm(1,rho,rho)))
#'
#' 
#'  
#' resultspath<-getwd()
#'
#' deltavalue=1
#' critical=10
#' numParticles=100
#' UseGill=0
#' UseSIWR=0
#'
#' #Set the population size for inference
#' PopSize=10000
#' #Set the rate immunity is lost
#' theMU=NA #doesn't matter what this is since we are estimating mu in this example
#' 
#' 
#' setstartmeansval=list(c(10000*.21,10000*.0015)) 
#'   
#' setBetaW=setKappa=setEta=setRatio=setAlpha0=0 #don't want to set these right now
#' setval=setAlpha0val=NA #don't want to set these right now
#' 
#' 
#' uset=uselaplace=0 # not using shrinkage priors for alpha parameters
#' alphadf=5
#' 
#' maxW=50000
#' 
#' ll=-50000
#' 
#' 
#' bayessirOUT=bayessir(obscholLIST,obsdaysLIST,COVMATLIST,
#'                     th,trans,invtrans,pmean,psd,
#'                     betaIIndex,betaWIndex,gammaIndex,kappaIndex,etaIndex,muIndex,alphasIndex,rhoIndex,startmeansIndex,nu1Index,nu2Index,nu3Index,nu4Index,
#'                     burn,prelim,iters,thin,tune,ll,psigma,
#'                     deltavalue,critical,PopSize,theMU,numParticles,resultspath,UseGill,UseSIWR,setBetaW,setKappa,setEta,setRatio,setval,setAlpha0,
#'                     setAlpha0val,setstartmeansval,uset,alphadf,uselaplace,maxW)
#' 
#' #Output columns are posterior samples for parameters in th, in addition to the log-likelihood and accepted values of the hidden states susT and infT at the final observation time T
#' 
#' #Posterior histograms for parameter values
#' nvars=dim(bayessirOUT)[2]
#' par(mfrow=c(1,nvars-3))
#' for(i in 1:(nvars-3)) hist(bayessirOUT[,i],main="",xlab=colnames(bayessirOUT)[i])
#' 
#' #Trace plots for all output
#' par(mfrow=c(1,nvars))
#' for(i in 1:(nvars)) plot(ts(bayessirOUT[,i]),xlab=colnames(bayessirOUT)[i],ylab="")
#' }
#' @export
bayessir=function(obscholLIST,obsdaysLIST,COVMATLIST,
                  th,trans,invtrans,
                  pmean,psd,
                  betaIIndex,betaWIndex,gammaIndex,kappaIndex,etaIndex,muIndex,alphasIndex,rhoIndex,startmeansIndex,nu1Index,nu2Index,nu3Index,nu4Index,
                  burn,prelim,iters,thin,tune,
                  ll,psigma,
                  deltavalue,critical,PopSize,theMU,
                  numParticles,resultspath,UseGill,UseSIWR,setBetaW,setKappa,setEta,setRatio,setval,setAlpha0,setAlpha0val,setstartmeansval,
                  usetprior,alphadf,uselaplaceprior,maxWval){
  
  pacc=ptot=acc=tot=0
  
  if(!is.list(setstartmeansval)){
    message("setstartmeansval not list")
  }
  
  if(UseSIWR){
    numlatentstates=3
  }else{
    numlatentstates=2
  }  
  
  NumberOfPhases=length(COVMATLIST)
  
  lastx=rep(NA,numlatentstates*NumberOfPhases)

  transth<-trans(th)
  p=length(th)+numlatentstates*NumberOfPhases+1
  pthmat=matrix(0,nrow=prelim,ncol=p)
  
  statenames=c("susT","infT","WT")[1:numlatentstates]
  colnames(pthmat)=c(names(th),"ll",rep(statenames,NumberOfPhases))
  
  ptransthmat=matrix(0,nrow=prelim,ncol=p)
  colnames(ptransthmat)=c(names(th),"ll",rep(statenames,NumberOfPhases))
  
  thmat=matrix(0,nrow=iters,ncol=p)
  colnames(thmat)=c(names(th),"ll",rep(statenames,NumberOfPhases))
  
  transthmat=matrix(0,nrow=iters,ncol=p)
  colnames(transthmat)=c(names(th),"ll",rep(statenames,NumberOfPhases))
  
  maxdays=c()
  for(z in 1:NumberOfPhases){
    maxdays=c(maxdays,length(obsdaysLIST[[z]]))
  }
  allweights=matrix(NA,nrow=prelim,ncol=sum(maxdays)-NumberOfPhases)
  
  if((setBetaW+setKappa+setEta+setRatio > 0) && ((setBetaW+setKappa+setEta+setRatio) != length(setval))){
    stop("Error: Not enough values in for the fixed parameters")
  }
  
  ##############
  ##############
  
  begTime <- Sys.time()
  
  # prelim pMCMC loop
  for (i in 1:prelim) {
    if(i%%100==0){ 
      message(paste(i,""),appendLF=FALSE)
      write.table(matrix(c(i,difftime(Sys.time(),begTime,unit="hour"),pacc/ptot),nrow=1),
                  paste(resultspath,"/prelimpmcmctimes.csv",sep=""),
                  append=T,sep=",",col.names=F,row.names=F)
      write.csv(pthmat,paste(resultspath,"/prelimthmat.csv",sep=""))
      #write.csv(allweights,paste(resultspath,"/prelimweights.csv",sep="")) 
    }
    
    for (j in 1:thin) {    
      ptot<-ptot+1
      
      transthprop=rmvnorm(1,transth,psigma)
      thprop=invtrans(transthprop) 
      
      ##fix nu parameters
      
      if(is.na(nu1Index)){
        nu1prop=1
      }else{
        nu1prop=thprop[nu1Index]
      }
      
      if(is.na(nu2Index)){
        nu2prop=1
      }else{
        nu2prop=thprop[nu2Index]
      }

      if(is.na(nu3Index)){
        nu3prop=0
      }else{
        nu3prop=thprop[nu3Index]
      }
      
      if(is.na(nu4Index)){
        nu4prop=0
      }else{
        nu4prop=thprop[nu4Index]
      }
      
      if(any(is.na(startmeansIndex))){
        startmeansprop=setstartmeansval
      }else{
        startmeansprop=list()
        for(smI in 1:length(startmeansIndex)) startmeansprop[[smI]]=as.numeric(thprop[startmeansIndex[[smI]]])
      }

      
      if(is.na(muIndex)){
        MUprop=theMU
      }else{
        MUprop=thprop[muIndex]
      }
      
      if(i==1) print(c(nu1prop,nu2prop,nu3prop,nu4prop))

      llprops=rep(NA,NumberOfPhases)
      lastxprop=c()
      ESSweights=list()

for(z in 1:NumberOfPhases){
        
        
        if(setAlpha0){
          alphaBreakList<-enviforce(COVMATLIST[[z]],obsdaysLIST[[z]],c(setAlpha0val,thprop[alphasIndex]))
        }else{
          alphaBreakList<-enviforce(COVMATLIST[[z]],obsdaysLIST[[z]],thprop[alphasIndex])
        }
        
        if(UseSIWR){
          if(setBetaW && setEta && setKappa){
            theBetaW=setval[1]
            theKappa=setval[2]
            theEta=setval[3]
          }else if(setBetaW && setEta){
            theBetaW=setval[1]
            theKappa=thprop[kappaIndex]
            theEta=setval[2]
          }else if(setBetaW && setKappa){
            theBetaW=setval[1]
            theKappa=setval[2]
            theEta=thprop[etaIndex]
          }else if(setKappa && setEta){
            theBetaW=thprop[betaWIndex]
            theKappa=setval[1]
            theEta=setval[2]
          }else if(setBetaW){
            theBetaW=setval
            theKappa=thprop[kappaIndex]
            theEta=thprop[etaIndex]
          }else if(setKappa){
            theBetaW=thprop[betaWIndex]
            theKappa=setval
            theEta=thprop[etaIndex]
          }else if(setEta){
            theBetaW=thprop[betaWIndex]
            theKappa=thprop[kappaIndex]
            theEta=setval
          }else if(setRatio){
            theBetaW=thprop[betaWIndex]
            bdratio=setval
            theKappa=thprop[kappaIndex]
            theEta=theKappa/bdratio
          }else{
            theBetaW=thprop[betaWIndex]
            theKappa=thprop[kappaIndex]
            theEta=thprop[etaIndex]
          }           

          pf=SMCforSIWR(obscholLIST[[z]],obsdaysLIST[[z]],c(thprop[betaIIndex[z]],theBetaW,thprop[gammaIndex],MUprop,theKappa,theEta),alphaBreakList,
                                 startmeansprop[[z]],numParticles,PopSize,thprop[rhoIndex[z]],deltavalue,critical,UseGill,maxWval)    
        }else{
          pf=SMCwithModPossionTL(obscholLIST[[z]],obsdaysLIST[[z]],c(thprop[betaIIndex[z]],thprop[gammaIndex],MUprop,nu1prop,nu2prop,nu3prop,nu4prop),alphaBreakList,
                                 startmeansprop[[z]],numParticles,PopSize,thprop[rhoIndex[z]],deltavalue,critical,UseGill,UseSIWR)    
          
        }
   
        llprops[z]=pf$ll   
        lastxprop=c(lastxprop,pf$lastX)
        ESSweights[[z]]=pf$weights  

      }

llprop=sum(llprops)

      if (log(runif(1)) < llprop - ll + prior.fun(transthprop,pmean,psd,usetprior,uselaplaceprior,alphasIndex,alphadf) - prior.fun(transth,pmean,psd,usetprior,uselaplaceprior,alphasIndex,alphadf)){
#             +dmvnorm(transth,transthprop,psigma,log=T) - dmvnorm(transthprop,transth,psigma,log=T)) {
        transth=transthprop
        th=thprop
        ll=llprop
        lastx=lastxprop
        pacc<-pacc+1
      }
    }
    ptransthmat[i,]=c(transth,ll,lastx)
    pthmat[i,]=c(invtrans(transth),ll,lastx)
    allweights[i,1:length(unlist(ESSweights))]=unlist(ESSweights)


  }
  
  message("Done with prelim loop!")
  prelimtime<-Sys.time()-begTime
  message(Sys.time()-begTime)
  write.csv(pthmat,paste(resultspath,"/prelimthmat.csv",sep=""))
  
  ssigma<-cov(ptransthmat[burn:prelim,-c(which(colnames(ptransthmat) %in% c("ll",statenames)))])
  

allweights=matrix(NA,nrow=iters,ncol=sum(maxdays)-NumberOfPhases)

  ##final loop
  for (i in 1:iters) {
#     if(i%%10==0) message(paste(i,""),appendLF=FALSE)
    if(i%%100==0){ 
      message(paste(i,""),appendLF=FALSE)
      write.table(matrix(c(i,difftime(Sys.time(),begTime,unit="hour"),acc/tot),nrow=1),
                  paste(resultspath,"/FINALpmcmctimes.csv",sep=""),
                  append=T,sep=",",col.names=F,row.names=F)
      write.csv(thmat,paste(resultspath,"/FINALthmat.csv",sep=""))
      #write.csv(allweights,paste(resultspath,"/FINALweights.csv",sep=""))
      
    }
    
    for (j in 1:thin) {    
      tot<-tot+1
      
      transthprop=rmvnorm(1,transth,tune*ssigma,method="svd")    
      thprop=invtrans(transthprop) 
      
      if(is.na(nu1Index)){
        nu1prop=1
      }else{
        nu1prop=thprop[nu1Index]
      }
      
      if(is.na(nu2Index)){
        nu2prop=1
      }else{
        nu2prop=thprop[nu2Index]
      }
      
      if(is.na(nu3Index)){
        nu3prop=0
      }else{
        nu3prop=thprop[nu3Index]
      }
      
      if(is.na(nu4Index)){
        nu4prop=0
      }else{
        nu4prop=thprop[nu4Index]
      }
      
      if(any(is.na(startmeansIndex))){
        startmeansprop=setstartmeansval
      }else{
        startmeansprop=list()
        for(smI in 1:length(startmeansIndex)) startmeansprop[[smI]]=as.numeric(thprop[startmeansIndex[[smI]]])
      }      

      
      if(is.na(muIndex)){
        MUprop=theMU
      }else{
        MUprop=thprop[muIndex]
      }



      NumberOfPhases=length(COVMATLIST)
      llprops=rep(NA,NumberOfPhases)
      lastxprop=c()
      ESSweights=list()

      
      for(z in 1:NumberOfPhases){
        
        if(setAlpha0){
          alphaBreakList<-enviforce(COVMATLIST[[z]],obsdaysLIST[[z]],c(setAlpha0val,thprop[alphasIndex]))
        }else{
          alphaBreakList<-enviforce(COVMATLIST[[z]],obsdaysLIST[[z]],thprop[alphasIndex])
        }
        
        if(UseSIWR){
          
          if(setBetaW && setEta && setKappa){
            theBetaW=setval[1]
            theKappa=setval[2]
            theEta=setval[3]
          }else if(setBetaW && setEta){
            theBetaW=setval[1]
            theKappa=thprop[kappaIndex]
            theEta=setval[2]
          }else if(setBetaW && setKappa){
            theBetaW=setval[1]
            theKappa=setval[2]
            theEta=thprop[etaIndex]
          }else if(setKappa && setEta){
            theBetaW=thprop[betaWIndex]
            theKappa=setval[1]
            theEta=setval[2]
          }else if(setBetaW){
            theBetaW=setval
            theKappa=thprop[kappaIndex]
            theEta=thprop[etaIndex]
          }else if(setKappa){
            theBetaW=thprop[betaWIndex]
            theKappa=setval
            theEta=thprop[etaIndex]
          }else if(setEta){
            theBetaW=thprop[betaWIndex]
            theKappa=thprop[kappaIndex]
            theEta=setval
          }else if(setRatio){
            theBetaW=thprop[betaWIndex]
            bdratio=setval
            theKappa=thprop[kappaIndex]
            theEta=theKappa/bdratio
          }else{
            theBetaW=thprop[betaWIndex]
            theKappa=thprop[kappaIndex]
            theEta=thprop[etaIndex]
          }           
          
          
          
          pf=SMCforSIWR(obscholLIST[[z]],obsdaysLIST[[z]],c(thprop[betaIIndex[z]],theBetaW,thprop[gammaIndex],MUprop,theKappa,theEta),alphaBreakList,
                        startmeansprop[[z]],numParticles,PopSize,thprop[rhoIndex[z]],deltavalue,critical,UseGill,maxWval)    
        }else{
          pf=SMCwithModPossionTL(obscholLIST[[z]],obsdaysLIST[[z]],c(thprop[betaIIndex[z]],thprop[gammaIndex],MUprop,nu1prop,nu2prop,nu3prop,nu4prop),alphaBreakList,
                                 startmeansprop[[z]],numParticles,PopSize,thprop[rhoIndex[z]],deltavalue,critical,UseGill,UseSIWR)    
        }
        
        llprops[z]=pf$ll 
        lastxprop=c(lastxprop,pf$lastX)
        ESSweights[[z]]=pf$weights  
        
      }
      
      llprop=sum(llprops)

      
      if (log(runif(1)) < llprop - ll + prior.fun(transthprop,pmean,psd,usetprior,uselaplaceprior,alphasIndex,alphadf) - prior.fun(transth,pmean,psd,usetprior,uselaplaceprior,alphasIndex,alphadf)) {
        transth=transthprop
        th=thprop
        ll=llprop
        lastx=lastxprop
        acc<-acc+1
      }
    }
    transthmat[i,]=c(transth,ll,lastx)
    thmat[i,]=c(invtrans(transth),ll,lastx)
    allweights[i,1:length(unlist(ESSweights))]=unlist(ESSweights)

  }
  
  message("Done with final loop!")
  finaltime<-Sys.time()-begTime
  message(Sys.time()-begTime)
  write.csv(thmat,paste(resultspath,"/FINALthmat.csv",sep=""))
  
  return(thmat)
}


bayessirfinalloop=function(obscholLIST,obsdaysLIST,COVMATLIST,
                  th,trans,invtrans,
                  pmean,psd,
                  betaIIndex,betaWIndex,gammaIndex,kappaIndex,etaIndex,muIndex,alphasIndex,rhoIndex,startmeansIndex,nu1Index,nu2Index,nu3Index,nu4Index,
                  burn,prelim,iters,thin,tune,
                  ll,psigma,
                  deltavalue,critical,PopSize,theMU,
                  numParticles,resultspath,UseGill,UseSIWR,setBetaW,setKappa,setEta,setRatio,setval,setAlpha0,setAlpha0val,setstartmeansval,ssigma,
                  lastxint,usetprior,alphadf,uselaplaceprior,maxWval){

  NumberOfPhases=length(COVMATLIST)
  
  if(UseSIWR){
    numlatentstates=3
  }else{
    numlatentstates=2
  }  
  acc=tot=0
  
  transth<-trans(th)
  p=length(th)+numlatentstates*NumberOfPhases+1
  pthmat=matrix(0,nrow=prelim,ncol=p)
  
  statenames=c("susT","infT","WT")[1:numlatentstates]
  colnames(pthmat)=c(names(th),"ll",rep(statenames,NumberOfPhases))
  
  ptransthmat=matrix(0,nrow=prelim,ncol=p)
  colnames(ptransthmat)=c(names(th),"ll",rep(statenames,NumberOfPhases))
  
  thmat=matrix(0,nrow=iters,ncol=p)
  colnames(thmat)=c(names(th),"ll",rep(statenames,NumberOfPhases))
  
  transthmat=matrix(0,nrow=iters,ncol=p)
  colnames(transthmat)=c(names(th),"ll",rep(statenames,NumberOfPhases))
  
  lastx=lastxint

  maxdays=c()
  for(z in 1:NumberOfPhases){
    maxdays=c(maxdays,length(obsdaysLIST[[z]]))
  }
  allweights=matrix(NA,nrow=iters,ncol=sum(maxdays)-NumberOfPhases)
  
  ##############

  ##############
  
  begTime <- Sys.time()
     
  ##final loop
  for (i in 1:iters) {

#     if(i%%10==0) message(paste(i,""),appendLF=FALSE)
    if(i%%100==0){ 
      message(paste(i,""),appendLF=FALSE)
      write.table(matrix(c(i,difftime(Sys.time(),begTime,unit="hour"),acc/tot),nrow=1),
                  paste(resultspath,"/FINALpmcmctimes.csv",sep=""),
                  append=T,sep=",",col.names=F,row.names=F)
      write.csv(thmat,paste(resultspath,"/FINALthmat.csv",sep=""))
      #write.csv(allweights,paste(resultspath,"/FINALweights.csv",sep=""))
    }
    
    for (j in 1:thin) {
      tot<-tot+1
      
      transthprop=rmvnorm(1,transth,tune*ssigma,method="svd")    
      thprop=invtrans(transthprop) 
      
      if(is.na(nu1Index)){
        nu1prop=1
      }else{
        nu1prop=thprop[nu1Index]
      }
      
      if(is.na(nu2Index)){
        nu2prop=1
      }else{
        nu2prop=thprop[nu2Index]
      }
      
      if(is.na(nu3Index)){
        nu3prop=0
      }else{
        nu3prop=thprop[nu3Index]
      }
      
      if(is.na(nu4Index)){
        nu4prop=0
      }else{
        nu4prop=thprop[nu4Index]
      }
      
      if(any(is.na(startmeansIndex))){
        startmeansprop=setstartmeansval
      }else{
        
        startmeansprop=list()
        for(smI in 1:length(startmeansIndex)) startmeansprop[[smI]]=as.numeric(thprop[startmeansIndex[[smI]]])
        
      }
      
      
      if(is.na(muIndex)){
        MUprop=theMU
      }else{
        MUprop=thprop[muIndex]
      }
      
      
      llprops=rep(NA,NumberOfPhases)
      lastxprop=c()
      ESSweights=list()
      
      for(z in 1:NumberOfPhases){
        
        if(setAlpha0){
          alphaBreakList<-enviforce(COVMATLIST[[z]],obsdaysLIST[[z]],c(setAlpha0val,thprop[alphasIndex]))
        }else{
          alphaBreakList<-enviforce(COVMATLIST[[z]],obsdaysLIST[[z]],thprop[alphasIndex])
        }
        
        if(UseSIWR){
          
          if(setBetaW && setEta && setKappa){
            theBetaW=setval[1]
            theKappa=setval[2]
            theEta=setval[3]
          }else if(setBetaW && setEta){
            theBetaW=setval[1]
            theKappa=thprop[kappaIndex]
            theEta=setval[2]
          }else if(setBetaW && setKappa){
            theBetaW=setval[1]
            theKappa=setval[2]
            theEta=thprop[etaIndex]
          }else if(setKappa && setEta){
            theBetaW=thprop[betaWIndex]
            theKappa=setval[1]
            theEta=setval[2]
          }else if(setBetaW){
            theBetaW=setval
            theKappa=thprop[kappaIndex]
            theEta=thprop[etaIndex]
          }else if(setKappa){
            theBetaW=thprop[betaWIndex]
            theKappa=setval
            theEta=thprop[etaIndex]
          }else if(setEta){
            theBetaW=thprop[betaWIndex]
            theKappa=thprop[kappaIndex]
            theEta=setval
          }else if(setRatio){
            theBetaW=thprop[betaWIndex]
            bdratio=setval
            theKappa=thprop[kappaIndex]
            theEta=theKappa/bdratio
          }else{
            theBetaW=thprop[betaWIndex]
            theKappa=thprop[kappaIndex]
            theEta=thprop[etaIndex]
          }           
          
          
          pf=SMCforSIWR(obscholLIST[[z]],obsdaysLIST[[z]],c(thprop[betaIIndex[z]],theBetaW,thprop[gammaIndex],MUprop,theKappa,theEta),alphaBreakList,
                        startmeansprop[[z]],numParticles,PopSize,thprop[rhoIndex[z]],deltavalue,critical,UseGill,maxWval)    
        }else{
          pf=SMCwithModPossionTL(obscholLIST[[z]],obsdaysLIST[[z]],c(thprop[betaIIndex[z]],thprop[gammaIndex],MUprop,nu1prop,nu2prop,nu3prop,nu4prop),alphaBreakList,
                                 startmeansprop[[z]],numParticles,PopSize,thprop[rhoIndex[z]],deltavalue,critical,UseGill,UseSIWR)    
        }
        
        llprops[z]=pf$ll 
        lastxprop=c(lastxprop,pf$lastX)
        ESSweights[[z]]=pf$weights  
        
      }
      
      llprop=sum(llprops)
      
      if (log(runif(1)) < llprop - ll + prior.fun(transthprop,pmean,psd,usetprior,uselaplaceprior,alphasIndex,alphadf) - prior.fun(transth,pmean,psd,usetprior,uselaplaceprior,alphasIndex,alphadf)) {
        transth=transthprop
        th=thprop
        ll=llprop
        lastx=lastxprop
        acc<-acc+1
      }
    }
    transthmat[i,]=c(transth,ll,lastx)
    thmat[i,]=c(invtrans(transth),ll,lastx)
    allweights[i,1:length(unlist(ESSweights))]=unlist(ESSweights)

  }
  
  message("Done with final loop!")
  finaltime<-Sys.time()-begTime
  message(Sys.time()-begTime)
  write.csv(thmat,paste(resultspath,"/FINALthmat.csv",sep=""))
  
  return(thmat)
}

