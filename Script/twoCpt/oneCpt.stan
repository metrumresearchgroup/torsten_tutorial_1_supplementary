
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
  vector[nObs] logCObs = log(cObs);
  int nTheta = 3;
  int nCmt = 2;
}

parameters {
  real<lower = 0> CL;
  real<lower = 0> VC;
  real<lower = 0> ka;
  real<lower = 0> sigma;
}

transformed parameters {
  real theta[nTheta] = {CL, VC, ka};
  row_vector<lower = 0>[nEvent] concentration;
  row_vector<lower = 0>[nObs] concentrationObs;
  matrix<lower = 0>[nCmt, nEvent] mass;

  mass = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss, theta);

  concentration = mass[2, ] ./ VC;
  concentrationObs = concentration[iObs];
}

model {
  // priors
  CL ~ lognormal(log(10), 0.25); 
  VC ~ lognormal(log(35), 0.25);
  ka ~ lognormal(log(2.5), 1);
  sigma ~ normal(0, 1);

  logCObs ~ normal(log(concentrationObs), sigma);
}

generated quantities {
  real concentrationObsPred[nObs] 
    = exp(normal_rng(log(concentrationObs), sigma));

  vector[nObs] log_lik;
  for (i in 1:nObs)
    log_lik[i] = normal_lpdf(logCObs[i] | log(concentrationObs[i]), sigma);
}
