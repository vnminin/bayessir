#' @title PMMH algorithm for time-varying SIRS model
#' @param obscholLIST List containing the cholera counts for each phase of data 
#' @param obsdaysLIST List containing the observation days for each phase of data 
#' @param COVMATLIST List containing the matricies of daily covariates for each phase of data 
#' @param th Starting value for parameter vales
#' @param trans function to transform the parameter values 
#' @param invtrans inverse transformation function 
#' @param pmean Prior means for Normal prior distributions
#' @param psd Prior standard deviations for Normal prior distributions
#' @param constParIndex Index of which \code{th} values correspond to \eqn{\beta} and \eqn{\gamma}
#' @param alphasIndex Index of which \code{th} values correspond to \eqn{\alpha} parameters
#' @param rhoIndex Index of which \code{th} values correspond to \eqn{\rho}, the probability of infected individuals seeking treatment
#' @param startmeansIndex Index of which \code{th} values correspond to means of initial distributions for the numbers of susceptible and infected individuals respectively
#' @param burn Number of iterations for burn-in run
#' @param prelim Number of total iterations for preliminary run, preliminary run = burn-in run + secondary run
#' @param iters Number of iterations for final run
#' @param thin Amount to thin the chain; only every \code{thin}th iteration of the PMMH algorithm is saved
#' @param psigma Standard deviation for independent normal proposal distribution in preliminary run
#' @param tune Tuning parameter for the covariance of the multivariate normal proposal distribution in the final run of the PMMH algorithm
#' @param ll Starting value for log-likelihood
#' @param deltavalue Initial value to use for tau in tau-leaping algorithm
#' @param critical Critical size for modified tau-leaping algorithm; if the population of a compartment is lower than this number a single step algorithm is used until the population gets above the critical size.
#' @param PopSize Population size
#' @param theMU The rate at which immunity is lost
#' @param numParticles Number of particles 
#' @param resultspath File path for results
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
#'resultspath<-getwd()
#'
#'deltavalue=1
#'critical=10
#'numParticles=100
#'
#'#Set the population size for inference
#'PopSize=150000
#'#Set the rate immunity is lost
#'theMU=.0009
#' 
#' theseed=13
#' 
#' 
#' # simulate the data
#' SimTimes=seq(0,365*3, by=14)
#' 
#' int<- -7
#' A<- 3.5
#' wave<-pi/(365/2)
#' t<-0:(max(SimTimes))
#' 
#' alpha<-exp(int+A*sin(wave*t))
#' sincov<-sin(wave*t)
#' 
#' pop=150000
#' 
#' phiS=pop*.22
#' phiI=pop*.005
#' 
#' rep=365*4
#' 
#' th1=1.25e-08
#' th2=0.1
#' rho=.0008
#' 
#' set.seed(10)
#' sus0=rpois(1,phiS)
#' inf0=rpois(1,phiI)
#' 
#' shortstart<-as.matrix(c(sus0,inf0))
#' 
#' allcovs<-enviforce(as.matrix(c(sincov)),SimTimes,c(int,A))
#' 
#' simstates<-matrix(NA,nrow=(length(SimTimes)),ncol=2)
#' 
#' simstates[1,]<-shortstart
#' 
#' for (i in 2:length(SimTimes)){
#'   simstates[i,]<-inhomoSIRSGillespie(simstates[i-1,],pop,SimTimes[i-1],SimTimes[i]-SimTimes[i-1],
#'                                     c(th1,th2,theMU),allcovs[[i-1]][,2],allcovs[[i-1]][,1])
#' }
#' 
#' SimData<-c()
#' for(i in 1:dim(simstates)[1]) SimData[i]<-rbinom(1,simstates[i,2],rho)
#' 
#' COVMATLIST=list(as.matrix(c(sincov)))
#' obscholLIST=list(SimData)
#' obsdaysLIST=list(SimTimes)
#' ##
#' 
#' trans=function(p){
#'   c(log(p[1]),  #beta
#'     log(p[2]),  #gamma
#'     p[3],       #alpha0
#'     p[4],       #alpha1
#'     logit(p[5]),#rho
#'     log(p[6]),  #phiS
#'     log(p[7]))  #phiI
#' }
#' 
#' invtrans=function(p){
#'   c(exp(p[1]),  #beta 
#'     exp(p[2]),  #gamma
#'     p[3],       #alpha0
#'     p[4],       #alpha1
#'     expit(p[5]),#rho
#'     exp(p[6]),  #phiS
#'     exp(p[7]))  #phiI
#' }
#' 
#' 
#' #prior means
#' pbeta<-log(.000001) 
#' pgamma<-log(.1)
#' palpha0<-log(.000001) 
#' palpha1<-0
#' prho<-logit(.0001) 
#' pphiS<-log(phiS*1.1)
#' pphiI<-log(phiI*1.1)
#' 
#' pmean=c(pbeta,pgamma,palpha0,palpha1,prho,pphiS,pphiI)
#' 
#' #prior standard deviations
#' psd=c(5,  #beta
#'       .09,#gamma 
#'       5,  #alpha0
#'       5,  #alpha1
#'       2,  #rho
#'       .1, #phiS
#'      .5) #phiI
#' 
#' 
#' constParIndex<-1:2
#' alphasIndex<-3:4
#' rhoIndex<-5
#' startmeansIndex<-6:7
#' 
#' # Iterations set small for example purposes; increase for applications
#' burn = 0
#' prelim = 10
#' iters =10
#' thin =1
#' 
#' tune=1
#' 
#' psigma<-diag(c(0.003, #beta
#'             0.003, #gamma
#'             0.030, #alpha0
#'             0.030, #alpha1
#'             0.030, #rho
#'             0.030, #phiS
#'             0.030)) #phiI
#' 
#' 
#' #start values
#' #Names of th input are used for the column names in the matrix output
#' set.seed(theseed)
#' th=c(
#'   beta=abs(rnorm(1,th1,th1/3)), 
#'   gamma=abs(rnorm(1,th2,th2/10)),
#'   alpha0=rnorm(1,int,1),
#'   alpha1=rnorm(1,A,1),
#'   rho=abs(rnorm(1,rho,rho)),
#'   phiS=.5*PopSize,
#'   phiI=.1*PopSize) 
#' 
#' ll=-50000
#' 
#' 
#' bayessirOUT=bayessir(obscholLIST,obsdaysLIST,COVMATLIST,
#'                      th,trans,invtrans,pmean,psd,constParIndex,alphasIndex,rhoIndex,startmeansIndex,
#'                      burn,prelim,iters,thin,tune,ll,psigma,
#'                      deltavalue,critical,PopSize,theMU,numParticles,resultspath)
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
                  constParIndex,alphasIndex,rhoIndex,startmeansIndex,
                  burn,prelim,iters,thin,tune,
                  ll,psigma,
                  deltavalue,critical,PopSize,theMU,
                  numParticles,resultspath){
  
  
  pacc=ptot=acc=tot=0
  
  transth<-trans(th)
  p=length(th)+3
  pthmat=matrix(0,nrow=prelim,ncol=p)
  colnames(pthmat)=c(names(th),"ll","susT","infT")
  
  ptransthmat=matrix(0,nrow=prelim,ncol=p)
  colnames(ptransthmat)=c(names(th),"ll","susT","infT")
  
  thmat=matrix(0,nrow=iters,ncol=p)
  colnames(thmat)=c(names(th),"ll","susT","infT")
  
  transthmat=matrix(0,nrow=iters,ncol=p)
  colnames(transthmat)=c(names(th),"ll","susT","infT")
  
  
  begTime <- Sys.time()
  
  # prelim pMCMC loop
  for (i in 1:prelim) {
    if(i%%10==0) message(paste(i,""),appendLF=FALSE)
    if(i%%100==0){ 
      write.table(matrix(c(i,difftime(Sys.time(),begTime,unit="hour"),pacc/ptot),nrow=1),
                  paste(resultspath,"/prelimpmcmctimes.csv",sep=""),
                  append=T,sep=",",col.names=F,row.names=F)
      write.csv(pthmat,paste(resultspath,"/prelimthmat.csv",sep=""))
    }
    
    for (j in 1:thin) {    
      ptot<-ptot+1
      
      transthprop=rmvnorm(1,transth,psigma)
      thprop=invtrans(transthprop) 
      
      NumberOfPhases=length(COVMATLIST)
      llprops=rep(NA,NumberOfPhases)
      for(z in 1:NumberOfPhases){
        alphaBreakList<-enviforce(COVMATLIST[[z]],obsdaysLIST[[z]],thprop[alphasIndex])
        pf=SMCwithModPossionTL(obscholLIST[[z]],obsdaysLIST[[z]],c(thprop[constParIndex],theMU),alphaBreakList,
                                thprop[startmeansIndex],numParticles,PopSize,thprop[rhoIndex],deltavalue,critical)    
        llprops[z]=pf$ll       
      }
      
      llprop=sum(llprops)
      
      lastxprop=pf$lastX #only need to save the last state of the last phase
        
      if (log(runif(1)) < llprop - ll + prior.fun(transthprop,pmean,psd) - prior.fun(transth,pmean,psd) +
            dmvnorm(transth,transthprop,psigma,log=T) - dmvnorm(transthprop,transth,psigma,log=T)) {
        transth=transthprop
        th=thprop
        ll=llprop
        lastx=lastxprop
        pacc<-pacc+1
      }
    }
    ptransthmat[i,]=c(transth,ll,lastx)
    pthmat[i,]=c(invtrans(transth),ll,lastx)
  }
  
  message("Done with prelim loop!")
  prelimtime<-Sys.time()-begTime
  message(Sys.time()-begTime)
  write.csv(pthmat,paste(resultspath,"/prelimthmat.csv",sep=""))
  
  ssigma<-cov(ptransthmat[burn:prelim,-c(which(colnames(ptransthmat) %in% c("ll","susT","infT")))])
  
  ##final loop
  for (i in 1:iters) {
    if(i%%10==0) message(paste(i,""),appendLF=FALSE)
    if(i%%100==0){ 
      write.table(matrix(c(i,difftime(Sys.time(),begTime,unit="hour"),acc/tot),nrow=1),
                  paste(resultspath,"/FINALpmcmctimes.csv",sep=""),
                  append=T,sep=",",col.names=F,row.names=F)
      write.csv(thmat,paste(resultspath,"/FINALthmat.csv",sep=""))
    }
    
    for (j in 1:thin) {    
      tot<-tot+1
      
      transthprop=rmvnorm(1,transth,tune*ssigma,method="svd")    
      thprop=invtrans(transthprop) 
      
      NumberOfPhases=length(COVMATLIST)
      llprops=rep(NA,NumberOfPhases)
      for(z in 1:NumberOfPhases){
        alphaBreakList<-enviforce(COVMATLIST[[z]],obsdaysLIST[[z]],thprop[alphasIndex])
        pf=SMCwithModPossionTL(obscholLIST[[z]],obsdaysLIST[[z]],c(thprop[constParIndex],theMU),alphaBreakList,
                               thprop[startmeansIndex],numParticles,PopSize,thprop[rhoIndex],deltavalue,critical)    
        llprops[z]=pf$ll       
      }
      
      llprop=sum(llprops)
      
      lastxprop=pf$lastX #only need to save the last state of the last phase
      
      if (log(runif(1)) < llprop - ll + prior.fun(transthprop,pmean,psd) - prior.fun(transth,pmean,psd)) {
        transth=transthprop
        th=thprop
        ll=llprop
        lastx=lastxprop
        acc<-acc+1
      }
    }
    transthmat[i,]=c(transth,ll,lastx)
    thmat[i,]=c(invtrans(transth),ll,lastx)
  }
  
  message("Done with final loop!")
  finaltime<-Sys.time()-begTime
  message(Sys.time()-begTime)
  write.csv(thmat,paste(resultspath,"/FINALthmat.csv",sep=""))
  
  return(thmat)
}



