
data {
  int<lower = 1> nEvent;
  int<lower = 1> nObs;
  int<lower = 1> iObs[nObs];

  // Event schedule
  int<lower = 1> cmt[nEvent];
  int evid[nEvent];
  int addl[nEvent];
  int ss[nEvent];
  real amt[nEvent];
  real time[nEvent];
  real rate[nEvent];
  real ii[nEvent];

  vector<lower = 0>[nObs] cObs;
}

transformed data {
  int nTheta = 5;
  int nCmt = 3;
}

parameters {
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> VC;
  real<lower = 0> VP;
  real<lower = 0> ka;
  real<lower = 0> sigma;
}

transformed parameters {
  real theta[nTheta] = {CL, Q, VC, VP, ka};
  row_vector<lower = 0>[nEvent] concentrationHat;
  matrix<lower = 0>[nCmt, nEvent] mass;

  mass = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, theta);

  concentrationHat = mass[2, ] ./ VC;
}

model {
  // priors
  CL ~ lognormal(log(10), 0.25); 
  Q ~ lognormal(log(15), 0.5);
  VC ~ lognormal(log(35), 0.25);
  VP ~ lognormal(log(105), 0.5);
  ka ~ lognormal(log(2.5), 1);
  sigma ~ normal(0, 1);

  // likelihood
  cObs ~ lognormal(log(concentrationHat[iObs]), sigma);
}

generated quantities {
  real concentrationObsPred[nObs] 
    = lognormal_rng(log(concentrationHat[iObs]), sigma);

  vector[nObs] log_lik;
  for (i in 1:nObs)
    log_lik[i] = lognormal_lpdf(cObs[i] | 
                                log(concentrationHat[iObs[i]]), sigma);
}
