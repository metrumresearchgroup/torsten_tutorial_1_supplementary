functions{
  real[] twocptneutmodelode(real t, real[] y, real[] y_pk, real[] theta, real[] rdummy, int[] idummy){
    /* PK variables */
    real V1 = theta[3];

    /* PD variable */
    real mtt      = theta[6];
    real circ0    = theta[7];
    real alpha    = theta[8];
    real gamma    = theta[9];
    real ktr      = 4.0 / mtt;
    real prol     = y[1] + circ0;
    real transit1 = y[2] + circ0;
    real transit2 = y[3] + circ0;
    real transit3 = y[4] + circ0;
    real circ     = fmax(machine_precision(), y[5] + circ0);
    real conc     = y_pk[2] / V1;
    real EDrug    = alpha * conc;

    real dydt[5];

    dydt[1] = ktr * prol * ((1 - EDrug) * ((circ0 / circ)^gamma) - 1);
    dydt[2] = ktr * (prol - transit1);
    dydt[3] = ktr * (transit1 - transit2);
    dydt[4] = ktr * (transit2 - transit3);
    dydt[5] = ktr * (transit3 - circ);

    return dydt;
  }
}

data{
  int<lower = 1> nt;
  int<lower = 1> nObsPK;
  int<lower = 1> nObsPD;
  int<lower = 1> iObsPK[nObsPK];
  int<lower = 1> iObsPD[nObsPD];
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> time[nt];
  vector<lower = 0>[nObsPK] cObs;
  vector<lower = 0>[nObsPD] neutObs;
  real<lower = 0> CLPrior;
  real<lower = 0> QPrior;
  real<lower = 0> V1Prior;
  real<lower = 0> V2Prior;
  real<lower = 0> kaPrior;
  real<lower = 0> CLPriorCV;
  real<lower = 0> QPriorCV;
  real<lower = 0> V1PriorCV;
  real<lower = 0> V2PriorCV;
  real<lower = 0> kaPriorCV;
  real<lower = 0> circ0Prior;
  real<lower = 0> circ0PriorCV;
  real<lower = 0> mttPrior;
  real<lower = 0> mttPriorCV;
  real<lower = 0> gammaPrior;
  real<lower = 0> gammaPriorCV;
  real<lower = 0> alphaPrior;
  real<lower = 0> alphaPriorCV;

  real<lower = 0> rtol;
  real<lower = 0> atol;
}

transformed data{
  real<lower = 0> rate[nt] = rep_array(0.0, nt);
  real<lower = 0> ii[nt] = rep_array(0.0, nt);
  int<lower = 0> addl[nt] = rep_array(0, nt);
  int<lower = 0> ss[nt] = rep_array(0, nt);
  vector[nObsPK] logCObs = log(cObs);
  vector[nObsPD] logNeutObs = log(neutObs);
  int nOde = 5;
  int<lower = 1> nCmt = nOde + 3;
  real biovar[nCmt] = rep_array(1.0, nCmt);
  real tlag[nCmt] = rep_array(0.0, nCmt);
}

parameters{
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  //  real<lower = 0> ka; // ka unconstrained
  real<lower = (CL / V1 + Q / V1 + Q / V2 +
		sqrt((CL / V1 + Q / V1 + Q / V2)^2 -
		     4 * CL / V1 * Q / V2)) / 2> ka; // ka > lambda_1
  real<lower = 0> mtt;
  real<lower = 0> circ0;
  real<lower = 0> alpha;
  real<lower = 0> gamma;
  real<lower = 0> sigma;
  real<lower = 0> sigmaNeut;
}

transformed parameters{
  vector[nt] cHat;
  vector<lower = 0>[nObsPK] cHatObs;
  vector[nt] neutHat;
  vector<lower = 0>[nObsPD] neutHatObs;
  matrix[nCmt, nt] x;
  real<lower = 0> parms[9];

  parms = {CL, Q, V1, V2, ka, mtt, circ0, gamma, alpha};

  x = pmx_solve_twocpt_rk45(twocptneutmodelode, nOde, time, amt, rate, ii, evid, cmt, addl, ss, parms, biovar, tlag, rtol, atol, 1e4);

  cHat = x[2, ]' / V1;
  neutHat = x[8, ]' + circ0;

  cHatObs    = cHat[iObsPK];
  neutHatObs = neutHat[iObsPD];
}

model{
  CL ~ lognormal(log(CLPrior), CLPriorCV);
  Q ~ lognormal(log(QPrior), QPriorCV);
  V1 ~ lognormal(log(V1Prior), V1PriorCV);
  V2 ~ lognormal(log(V2Prior), V2PriorCV);
  ka ~ lognormal(log(kaPrior), kaPriorCV);
  sigma ~ cauchy(0, 1);

  mtt ~ lognormal(log(mttPrior), mttPriorCV);
  circ0 ~ lognormal(log(circ0Prior), circ0PriorCV);
  alpha ~ lognormal(log(alphaPrior), alphaPriorCV);
  gamma ~ lognormal(log(gammaPrior), gammaPriorCV);
  sigmaNeut ~ cauchy(0, 1);

  logCObs ~ normal(log(cHatObs), sigma); // observed data likelihood
  logNeutObs ~ normal(log(neutHatObs), sigmaNeut);
}

generated quantities{
  real cObsPred[nt];
  real neutObsPred[nt];


  for(i in 1:nt){
    if(time[i] == 0){
      cObsPred[i] = 0;
    }else{
      cObsPred[i] = exp(normal_rng(log(fmax(machine_precision(), cHat[i])),
				   sigma));
    }
    neutObsPred[i] = exp(normal_rng(log(fmax(machine_precision(),
					     neutHat[i])), sigmaNeut));
  }
}
