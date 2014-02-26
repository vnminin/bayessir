#include <Rcpp.h>
using namespace Rcpp;

//' Sample from initial distributions
//' 
//' Simulates the initial number of susceptible and infected individuals from independent Poisson initial distributions.  
//' @param n Number of particles 
//' @param startState Numeric vector of size two containing the means of initial distributions for the numbers of susceptible and infected individuals respectively
//' @return Matrix of dimension n by 2, where each row is the initial state for a particle trajectory
//' @export
// [[Rcpp::export]]
IntegerMatrix simInit(int n, NumericVector startState){
    using namespace Rcpp;

    IntegerMatrix startMatrix(n,2);

    startMatrix( _, 0) = rpois(n, startState[0]);
    startMatrix( _, 1) = rpois(n, startState[1]);

    return(startMatrix);
}

//' Probability of observed data
//' 
//' Calculates the probability of the observed data given the latent states and the parameters, assuming the number of observed infections has a binomial distribution 
//' @param infected Number of infected individuals
//' @param observed Number of observed individuals
//' @param repRate probability of infected individuals seeking treatment
//' @param givelog boolian; if 1, returns log 
//' @return Probability of observed data 
//' @export
// [[Rcpp::export]]
double dataLik(int infected, int observed, double repRate, bool givelog){ 
    using namespace Rcpp;
    
    IntegerVector RcppObs(1);
    RcppObs[0] = observed;
  
    NumericVector dataloglik = dbinom(RcppObs, infected, repRate, givelog); 
    return(as<double>(dataloglik));
  
}

//' Samples a vector based on weights
//' @param weights input numeric vector
//' @return Index of vector selected
//' @export
// [[Rcpp::export]]
int sampleOnce(NumericVector weights) {
    
    double rUnif = as<double>(runif(1));
    double total = sum(weights);
    double cumProb=0;
 
    int i;
 
    for (i = 0; i < weights.size(); i++) {
  cumProb += weights(i) / total;
	if (rUnif < cumProb) break;
    }
    return(i);
}




//' Modified inhomogeneous Gillespie algorithm 
//' 
//' Uses a modified Gillespie algorithm to simulate the numbers of susceptible and infected individuals forward in time when value of environmental force of infection changes
//' @param startState Integer vector of size two containing the initial number of susceptible and infected individuals respectively
//' @param popSize Population size
//' @param startTime Start time of simulation
//' @param intervalLength Length of simulation interval 
//' @param parameters Numeric vector of size seven containing the current values of the infectious contact rate, the recovery rate, the rate at which immunity is lost, and the powers
//' @param alphas Numeric vector containing the values of alpha
//' @param allbreaks Numeric vector containing the times at which the value of alpha changes
//' @return Vector containing number of susceptible and infected at the end of the simulation interval 
//' @export
// [[Rcpp::export]]
IntegerVector inhomoSIRSGillespie(IntegerVector startState, int popSize, 
  		   double startTime, double intervalLength, 
			   NumericVector parameters,
			   NumericVector alphas, NumericVector allbreaks){

    RNGScope scope;

    double curTime = startTime;
    double endTime = startTime + intervalLength;
    IntegerVector curState(clone(startState));

    int howmanybreaks = allbreaks.size()-1;

    double currentalpha = alphas[0];


    double curNextBreak = allbreaks[1];
    int numChange = 0;

    double h0, h1, h2, h3, u;

    while(curTime < endTime){

      h1 = parameters[0]*pow(curState[1],parameters[3])*pow(curState[0],parameters[4]) + 
            exp(parameters[5]*curState[1])*pow(curState[1],parameters[6])*curState[0]*currentalpha;
	    h2 = parameters[1]*curState[1];
	    h3 = parameters[2]*(popSize - curState[0] - curState[1]);
	    h0 = h1 + h2 + h3;

	    if (h0<1e-10){
	      curTime=1e99;
	    }else{
	      curTime+=as<double>(rexp(1, h0));
	    }

	    if (curTime < endTime){

	      if (curTime <= curNextBreak){ 
		      u = as<double>(runif(1));

		      if (u<h1/h0){
		        curState[0]-=1; 
		        curState[1]+=1;
		      }else{ 
		        if (u<(h1+h2)/h0){ 
			      curState[1]-=1;
		        }else{
			        curState[0]+=1;
		        }
		      }
	      }else{
		      curTime = curNextBreak;		
		      numChange += 1;
		      currentalpha = alphas[numChange];
		      if (numChange < howmanybreaks){
		        curNextBreak = allbreaks[numChange+1];
		      }else{
		        curNextBreak = endTime;
		      }
	      }
	    }
    }
  return(wrap(curState));
}



//' Modified inhomogeneous Poisson tau-leaping algorithm 
//' 
//' Uses a tau-leaping algorithm to simulate the numbers of susceptible and infected individuals forward in time when value of environmental force of infection changes
//' @param startState Integer vector of size two containing the initial number of susceptible and infected individuals respectively
//' @param popSize Population size
//' @param startTime Start time of simulation
//' @param intervalLength Length of simulation interval 
//' @param parameters Numeric vector of size seven containing the current values of the infectious contact rate, the recovery rate, the rate at which immunity is lost, and the powers
//' @param alphas Numeric vector containing the values of alpha
//' @param allbreaks Numeric vector containing the times at which the value of alpha changes
//' @param deltatint Initial value to use for tau in tau-leaping algorithm
//' @param ncrit Critical size for modified tau-leaping algorithm; if the population of a compartment is lower than this number a single step algorithm is used until the population gets above the critical size.
//' @return Vector containing number of susceptible and infected at the end of the simulation interval 
//' @export
// [[Rcpp::export]]
IntegerVector inhomoModPoissonTL(IntegerVector startState, int popSize, 
				    double startTime, double intervalLength, 
				    NumericVector parameters,
				    NumericVector alphas, NumericVector allbreaks,
				    double deltatint,int ncrit){
    
    RNGScope scope;
  
    double curTime = startTime;
    double endTime = startTime + intervalLength;

    IntegerVector curState(clone(startState));
    
    int howmanybreaks = allbreaks.size()-1;

    double currentalpha = alphas[0];

    double curNextBreak = allbreaks[1];
    int numChange = 0;

    IntegerVector transition(3);
    
    IntegerMatrix AllStates(howmanybreaks+2,2);
    
    AllStates(0,0) = curState[0];
    AllStates(0,1) = curState[1];

    double deltat;
   
    double h0, u;	
    NumericVector h(3);		

    while(curTime - endTime < -deltatint/100){

	if(curState[1]<ncrit || curState[0]<ncrit || (popSize-curState[0]-curState[1])<ncrit){
	    
      h[0] = parameters[0]*pow(curState[1],parameters[3])*pow(curState[0],parameters[4]) + 
        exp(parameters[5]*curState[1])*pow(curState[1],parameters[6])*curState[0]*currentalpha;
	    h[1] = parameters[1]*curState[1];
	    h[2] = parameters[2]*(popSize-curState[0]-curState[1]);

	    h0 = h[0] + h[1] + h[2];

	    if( h0 < 1e-10){
		curTime = 1e99;
	    } else {
		curTime += as<double>(rexp(1, h0));
	    }

	    if(curTime <= endTime){

		u = as<double>(runif(1));
    
		if (u < h[0]/h0){
		    curState[0]-=1; 
		    curState[1]+=1;
		} else { 
		    if (u < (h[0]+h[1])/h0){ 
			curState[1] -= 1;
		    } else {
			curState[0] += 1;
		    }
		}
	    }
      
	} else {
	    

	    if(curNextBreak-curTime < deltatint && curNextBreak-curTime > 0){
		deltat=curNextBreak-curTime;
	    }else{
		deltat = deltatint;
	    }

	    transition[0] = as<int>(rpois(1, (parameters[0]*pow(curState[1],parameters[3])*pow(curState[0],parameters[4]) + 
          exp(parameters[5]*curState[1])*pow(curState[1],parameters[6])*curState[0]*currentalpha)*deltat));
	    transition[1] = as<int>(rpois(1, (parameters[1]*curState[1])*deltat));
	    transition[2] = as<int>(rpois(1, (parameters[2]*(popSize-curState[0]-curState[1]))*deltat));

	    curTime += deltat;
      
	    if (
		((curState[0]+transition[2]-transition[0]) >= 0) &&
		((curState[1]+transition[0]-transition[1]) >= 0) &&
		((popSize-curState[0]-curState[1]+transition[1]-transition[2]) >= 0)){

		if(curTime <= curNextBreak){
		    curState[0] += transition[2]-transition[0];
		    curState[1] += transition[0]-transition[1];
		} else {
		    if(curTime <= endTime){
			AllStates(numChange+1,0) = curState[0];
			AllStates(numChange+1,1) = curState[1];

			curTime = curNextBreak;
			numChange += 1;
			currentalpha = alphas[numChange];
			if (numChange < howmanybreaks){
			    curNextBreak = allbreaks[numChange + 1];
			} else {
			    curNextBreak = endTime;
			}

		    }
		}

	    } else {

		//Rcout<<"Warning: Negative population"<<std::endl;
		
		curTime -= deltat;

		while(
		    ((curState[0]+transition[2]-transition[0])<0) ||
		    ((curState[1]+transition[0]-transition[1])<0) ||
		    ((popSize-curState[0]-curState[1]+transition[1]-transition[2])<0)){

		    deltat = deltat/2;
        
		    transition[0] = as<int>(rpois(1, (parameters[0]*pow(curState[1],parameters[3])*pow(curState[0],parameters[4]) + 
          exp(parameters[5]*curState[1])*pow(curState[1],parameters[6])*curState[0]*currentalpha)*deltat));
		    transition[1] = as<int>(rpois(1, (parameters[1]*curState[1])*deltat));
		    transition[2] = as<int>(rpois(1, (parameters[2]*(popSize-curState[0]-curState[1]))*deltat));

		}

		curTime += deltat;
		if (curTime <= curNextBreak){

		    curState[0] += transition[2]-transition[0];
		    curState[1] += transition[0]-transition[1];

		} else {
		    if(curTime <= endTime){
			AllStates(numChange+1,0) = curState[0];
			AllStates(numChange+1,1) = curState[1];

			curTime = curNextBreak;
			numChange += 1;
			currentalpha = alphas[numChange];
			if (numChange < howmanybreaks){
			    curNextBreak = allbreaks[numChange + 1];
			} else {
			    curNextBreak = endTime;
			}


		    }
		}
	    }
	    
	}
    }

    AllStates(howmanybreaks+1,0) = curState[0];
    AllStates(howmanybreaks+1,1) = curState[1];

    return(wrap(curState)); 
}


//' Modified homogeneous Poisson tau-leaping algorithm 
//' 
//' Modified Poisson tau-leaping algorithm. Uses a tau-leaping algorithm to simulate the numbers of susceptible and infected individuals forward in time when value of environmental force of infection remains constant
//' @param startState Integer vector of size two containing the initial number of susceptible and infected individuals respectively
//' @param popSize Population size
//' @param startTime Start time of simulation
//' @param intervalLength Length of simulation interval 
//' @param parameters parameters Numeric vector of size seven containing the current values of the infectious contact rate, the recovery rate, the rate at which immunity is lost, and the powers
//' @param alpha The value of the time-varying environmental force of infection
//' @param deltatint Initial value to use for tau in tau-leaping algorithm
//' @param ncrit Critical number
//' @return Vector containing number of susceptible and infected at the end of the simulation interval 
//' @export
// [[Rcpp::export]]
IntegerVector ModPoissonTL(IntegerVector startState, int popSize, 
			      double startTime, double intervalLength, 
			      NumericVector parameters,
			      double alpha, double deltatint,int ncrit){
    
    RNGScope scope;
  
    double curTime = startTime;
    double endTime = startTime + intervalLength;
    IntegerVector curState(clone(startState));
  
    IntegerVector transition(3);
    
    double deltat;
    double h0, u;	
    NumericVector h(3);		

    while(curTime - endTime < -deltatint/100){
	if(curState[1]<ncrit || curState[0]<ncrit || (popSize-curState[0]-curState[1])<ncrit){

	    h[0] = parameters[0]*pow(curState[1],parameters[3])*pow(curState[0],parameters[4]) + 
        exp(parameters[5]*curState[1])*pow(curState[1],parameters[6])*curState[0]*alpha;
	    h[1] = parameters[1]*curState[1];
	    h[2] = parameters[2]*(popSize-curState[0]-curState[1]);

	    h0 = h[0] + h[1] + h[2];

	    if( h0 < 1e-10){
		curTime = 1e99;
	    } else {
		curTime += as<double>(rexp(1, h0));
	    }

	    if(curTime <= endTime){

		u = as<double>(runif(1));
    
		if (u < h[0]/h0){
		    curState[0]-=1; 
		    curState[1]+=1;
		} else { 
		    if (u < (h[0]+h[1])/h0){ 
			curState[1] -= 1;
		    } else {
			curState[0] += 1;
		    }
		}
	    }


	    
	} else {
	    if(endTime-curTime < deltatint && endTime-curTime > 0){
		deltat=endTime-curTime;
	    }else{
	     	deltat = deltatint;
	    }

	    transition[0] = as<int>(rpois(1, (parameters[0]*pow(curState[1],parameters[3])*pow(curState[0],parameters[4]) + 
        exp(parameters[5]*curState[1])*pow(curState[1],parameters[6])*curState[0]*alpha)*deltat));
	    transition[1] = as<int>(rpois(1, (parameters[1]*curState[1])*deltat));
	    transition[2] = as<int>(rpois(1, (parameters[2]*(popSize-curState[0]-curState[1]))*deltat));
	    
	    curTime+=deltat;

      
	    if (
		((curState[0]+transition[2]-transition[0])>=0) &&
		((curState[1]+transition[0]-transition[1])>=0) &&
		((popSize-curState[0]-curState[1]+transition[1]-transition[2])>=0)){

		if(curTime <= endTime){

		    curState[0]+=transition[2]-transition[0];
		    curState[1]+=transition[0]-transition[1];
		}

	    }else{

		curTime-=deltat;

		while(
		    ((curState[0]+transition[2]-transition[0])<0) ||
		    ((curState[1]+transition[0]-transition[1])<0) ||
		    ((popSize-curState[0]-curState[1]+transition[1]-transition[2])<0)){

		    deltat=deltat/2;

		    transition[0] = as<int>(rpois(1, (parameters[0]*pow(curState[1],parameters[3])*pow(curState[0],parameters[4]) + 
          exp(parameters[5]*curState[1])*pow(curState[1],parameters[6])*curState[0]*alpha)*deltat));
		    transition[1] = as<int>(rpois(1, (parameters[1]*curState[1])*deltat));
		    transition[2] = as<int>(rpois(1, (parameters[2]*(popSize-curState[0]-curState[1]))*deltat));
		}
		curTime+=deltat;
		if (curTime <= endTime){
		    curState[0]+=transition[2]-transition[0];
		    curState[1]+=transition[0]-transition[1];
		}
	    }
	    
	}
    }
    return(wrap(curState)); 
}

//' Sequential Monte Carlo algorithm
//' 
//' 
//' @param observedCounts Integer vector of observed counts
//' @param observedTimes Numeric vector of observation times
//' @param parameters Numeric vector of size seven containing the current values of the infectious contact rate, the recovery rate, the rate at which immunity is lost, and the powers
//' @param alphaBreaks List containing the time-varying environmental force of infection and the times at which the rate changes, broken up by simulation intervals
//' @param startMeans Numeric vector of size two containing the means of initial distributions for the numbers of susceptible and infected individuals respectively
//' @param nparticles Number of particles 
//' @param popSize Population size
//' @param reportingRate Probability of infected individuals seeking treatment
//' @param deltatval Initial value to use for tau in tau-leaping algorithm
//' @param ncrit Critical size for modified tau-leaping algorithm; if the population of a compartment is lower than this number a single step algorithm is used until the population gets above the critical size.
//' @return log-likelihood and estimate of last state
//' @export
// [[Rcpp::export]]
List SMCwithModPossionTL(IntegerVector observedCounts,
			  NumericVector observedTimes,
			  NumericVector parameters,
			  List alphaBreaks,
			  NumericVector startMeans,
			  int nparticles,
			  int popSize,
			  double reportingRate,
			  double deltatval,
			  int ncrit){
  
    RNGScope scope;

    IntegerMatrix xmat(nparticles,2);
    int numIntervals = observedTimes.size()-1;
    NumericVector weights(nparticles);
    double LogLik = 0;

    xmat = simInit(nparticles, startMeans);

    IntegerVector tempStates(2);

    for (int i=0; i<numIntervals; i++){

	NumericMatrix tempAlphaBreaks(as<SEXP>(alphaBreaks[i]));
	int startTime = tempAlphaBreaks(0,0);
	NumericVector alphas = tempAlphaBreaks( _, 1);
	NumericVector tempAlpha = tempAlphaBreaks( _,0);

	int numBreaks = tempAlphaBreaks.nrow();
  
	if(numBreaks > 1){
	    NumericVector allbreaks(tempAlpha.begin(), tempAlpha.end());

	    for (int j=0; j<nparticles; j++){
		tempStates = inhomoModPoissonTL(IntegerVector::create(xmat(j,0),xmat(j,1)), popSize, 
						observedTimes[i], observedTimes[i+1]-observedTimes[i], 
						parameters, alphas,allbreaks,deltatval,ncrit);
		
		xmat(j,0) = tempStates[0];
		xmat(j,1) = tempStates[1];

		weights[j] = dataLik(xmat(j,1), observedCounts[i+1], reportingRate,0);
		
	    }
	}else{ 

	    double thealpha = alphas[0];

	    for (int j=0; j<nparticles; j++){

		tempStates = ModPoissonTL(IntegerVector::create(xmat(j,0),xmat(j,1)), popSize, 
					  observedTimes[i], observedTimes[i+1]-observedTimes[i], 
					  parameters, thealpha,deltatval,ncrit);

		xmat(j,0) = tempStates[0];
		xmat(j,1) = tempStates[1];
		
		weights[j] = dataLik(xmat(j,1), observedCounts[i+1], reportingRate,0);
      
	    }
	}

	bool res = is_true( all( weights == 0 ) );
	if(res){
	    Rcout<<"Warning: all weights zero"<<std::endl;
	    return Rcpp::List::create(Rcpp::Named("xmat") =xmat,
				      Rcpp::Named("ll") =-1e+99);
	}

	LogLik+=log(mean(weights));

	IntegerMatrix tempXmat(nparticles,2);

	for (int j=0; j<nparticles; j++){
	    int tempRowIndex = sampleOnce(weights);
	    tempXmat(j,_)=xmat(tempRowIndex,_);
	}  

	xmat = tempXmat;

    }

    int LastRowIndex = sampleOnce(weights);
    IntegerVector LastX(2);
    LastX(0)=xmat(LastRowIndex,0);
    LastX(1)=xmat(LastRowIndex,1);


    return Rcpp::List::create(Rcpp::Named("lastX") =LastX,
			      Rcpp::Named("ll") =LogLik);
}















