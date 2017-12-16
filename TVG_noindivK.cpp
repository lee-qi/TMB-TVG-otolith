// Can include other random effect relationships using a DATA_FACTOR(Numero); TBD later
#define TMB_LIB_INIT R_init_TVG_noindivK
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_FACTOR(indiv);
  DATA_FACTOR(yrs);
  DATA_FACTOR(sex);

  DATA_VECTOR(widths);
  DATA_VECTOR(age);
  
  DATA_INTEGER(Nfish);											// Total number of individuals
  DATA_INTEGER(Nyears);                     // Total number of years plus 1 because of initial year for initial width
  DATA_INTEGER(Nsex);

  PARAMETER_VECTOR(logmeanK);									// Normally-distributed Kappa
  PARAMETER_VECTOR(logmeanWinf);							// Normally-distributed Winf
  PARAMETER_VECTOR(betaK);
  PARAMETER_VECTOR(betaW);
  PARAMETER_VECTOR(Omega);
//  PARAMETER(Rho);

  PARAMETER_VECTOR(Ke);
  PARAMETER_VECTOR(We);
  PARAMETER(logsigKe);
  PARAMETER(logsigWe);
  PARAMETER(logsigInc);

  int Ndat = widths.size();

  Type nll;														       // Objective function
  Type sigInc = exp(logsigInc);
  Type sigKe = exp(logsigKe);
  Type sigWe = exp(logsigWe);

  vector<Type> meanK(Nsex);
  vector<Type> meanWinf(Nsex);
  vector<Type> predInc(Ndat);
  vector<Type> Kvals(Nfish);
  vector<Type> Winfvals(Nfish);
  vector<Type> Ep(Nyears);
  vector<Type> startInc(Nfish);

  nll = 0.0;

  meanK = exp(logmeanK);
  meanWinf = exp(logmeanWinf);

  for(int fish=0; fish<Nfish; fish++) {
    Kvals(fish) = meanK(sex(fish)) * exp(Ke(fish));
    Winfvals(fish) = meanWinf(sex(fish)) * exp(We(fish));

// Individual sex-specific random effect
//    nll -= dnorm(Ke(fish), Type(0.0), sigKe, true);
    nll -= dnorm(We(fish), Type(0.0), sigWe, true);
  }

  Ep = Omega;

//  for(int yr=1; yr<Nyears; yr++) {
//    Ep(yr) = (Ep(yr-1) * Rho) + (sqrt(1-pow(Rho,2))*Omega(yr));
//  }

  Type tempK = Kvals(indiv(0)) * exp(betaK(sex(indiv(0))) * Ep(yrs(0)));
  Type tempWinf = Winfvals(indiv(0)) * exp(betaW(sex(indiv(0))) * Ep(yrs(0)));

  Type temppreK = Kvals(indiv(0)) * exp(betaK(sex(indiv(0))) * Ep(yrs(0)-1));
  Type temppreWinf = Winfvals(indiv(0)) * exp(betaW(sex(indiv(0))) * Ep(yrs(0)-1));

  startInc(indiv(0)) = temppreWinf * (1.0 - (exp(-temppreK * (age(0)-1))));
  predInc(0) = (exp(-tempK) - 1) * (startInc(indiv(0)) - tempWinf);

  Type tempWid = startInc(indiv(0)) + predInc(0);

  for(int i=1; i<Ndat; i++) {

  	// Time- and individually-varying random effect
  	Type tempK = Kvals(indiv(i)) * exp(betaK(sex(indiv(i))) * Ep(yrs(i)));
  	Type tempWinf = Winfvals(indiv(i)) * exp(betaW(sex(indiv(i))) * Ep(yrs(i)));

  	if(indiv(i) != indiv(i-1)) {
      Type temppreK = Kvals(indiv(i)) * exp(betaK(sex(indiv(i))) * Ep(yrs(i)-1));
      Type temppreWinf = Winfvals(indiv(i)) * exp(betaW(sex(indiv(i))) * Ep(yrs(i)-1)); 
  	  startInc(indiv(i)) = temppreWinf * (1.0 - (exp(-temppreK * (age(i)-1))));
      predInc(i) = (exp(-tempK) - 1) * (startInc(indiv(i)) - tempWinf);

      tempWid = startInc(indiv(i)) + predInc(i);
  	}

  	else {
  	  predInc(i) = (exp(-tempK) - 1.0) * (tempWid - tempWinf);
      tempWid += predInc(i);
  	}
  }

  nll -= sum(dnorm(Omega, Type(0.0), Type(1.0), true));
  nll -= sum(dnorm(widths,predInc,sigInc, true));

  ADREPORT(Kvals);
  ADREPORT(Winfvals);
  ADREPORT(Ep)
  ADREPORT(predInc);

  REPORT(Ke);
  REPORT(We);
  REPORT(Ep);
  REPORT(meanK);
  REPORT(meanWinf);
  REPORT(betaK);
  REPORT(betaW);
  REPORT(Kvals);
  REPORT(Winfvals);

  REPORT(predInc);

  REPORT(sigKe);
  REPORT(sigWe);
  REPORT(sigInc);

  return nll;
}

