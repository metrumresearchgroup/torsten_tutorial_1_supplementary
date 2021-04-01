functions {
  real[] twoCptNeutModelODE(real t, real[] x, real[] parms, real[] rdummy, int[] idummy){
    real CL = parms[1];
    real Q = parms[2];
    real V1 = parms[3];
    real V2 = parms[4];
    real ka = parms[5];
    real mtt = parms[6];	
    real circ0 = parms[7];
    real gamma = parms[8];
    real alpha = parms[9];
    real k10 = CL / V1;
    real k12 = Q / V1;
    real k21 = Q / V2;
    real ktr = 4 / mtt;
    real dxdt[8];
    real conc = x[2]/V1;
    real EDrug = fmin(1.0, alpha * conc);
    real prol = x[4] + circ0;
    real transit1 = x[5] + circ0;
    real transit2 = x[6] + circ0;
    real transit3 = x[7] + circ0;
    real circ = fmax(machine_precision(), x[8] + circ0);

    dxdt[1] = -ka * x[1];
    dxdt[2] = ka * x[1] - (k10 + k12) * x[2] + k21 * x[3];
    dxdt[3] = k12 * x[2] - k21 * x[3];

    // x[4], x[5], x[6], x[7] and x[8] are differences from circ0.
    dxdt[4] = ktr * prol * ((1 - EDrug) * ((circ0 / circ)^gamma) - 1);
    dxdt[5] = ktr * (prol - transit1);
    dxdt[6] = ktr * (transit1 - transit2);
    dxdt[7] = ktr * (transit2 - transit3);
    dxdt[8] = ktr * (transit3 - circ);

    return dxdt;
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
  int<lower = 0> max_num_step;
}

transformed data{
  real<lower = 0> rate[nt] = rep_array(0.0, nt);
  real<lower = 0> ii[nt] = rep_array(0.0, nt);
  int<lower = 0> addl[nt] = rep_array(0, nt);
  int<lower = 0> ss[nt] = rep_array(0, nt);
  vector[nObsPK] logCObs = log(cObs);
  vector[nObsPD] logNeutObs = log(neutObs);
  int<lower = 1> nCmt = 8;
  real F[nCmt] = rep_array(1.0, nCmt);
  real tLag[nCmt] = rep_array(0.0, nCmt);
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
  vector[nObsPK] cHatObs;
  vector[nt] neutHat;
  vector[nObsPD] neutHatObs;
  matrix[nCmt, nt] x;
  real<lower = 0> parms[9];

  parms = {CL, Q, V1, V2, ka, mtt, circ0, gamma, alpha};

  x = pmx_solve_rk45(twoCptNeutModelODE, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, parms, F, tLag, rtol, atol, max_num_step);

  cHat = x[2, ]' / V1;
  neutHat = x[8, ]' + circ0;

  cHatObs = cHat[iObsPK]; // predictions for observed data records
  neutHatObs = neutHat[iObsPD]; // predictions for observed data records
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
