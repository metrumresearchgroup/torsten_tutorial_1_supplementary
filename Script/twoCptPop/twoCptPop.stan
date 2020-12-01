
data {
  int<lower = 1> nEvent;
  int<lower = 1> nObs;
  int<lower = 1> iObs[nObs];
  int<lower = 1> nSubjects;
  int<lower = 1> start[nSubjects];
  int<lower = 1> end[nSubjects];

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
  int nTheta = 5;
  int nCmt = 3;
  int nIIV = 5;

  real biovar[nCmt] = {1.0, 1.0, 1.0};
  real tlag[nCmt] = {0.0, 0.0, 0.0};
  
  real prior_sd[nIIV] = {0.25, 0.5, 0.25, 0.5, 1};
}

parameters {
  // Population parameters
  real<lower = 0> CL_pop;
  real<lower = 0> Q_pop;
  real<lower = 0> VC_pop;
  real<lower = 0> VP_pop;
  real<lower = 0> ka_pop;

  // Inter-individual variability
  vector<lower = 0>[nIIV] omega;
  real<lower = 0> theta[nSubjects, nTheta];

  real<lower = 0> sigma;
}

transformed parameters {
  vector<lower = 0>[nTheta] 
    theta_pop = to_vector({CL_pop, Q_pop, VC_pop, VP_pop, ka_pop});
  row_vector<lower = 0>[nEvent] concentration;
  row_vector<lower = 0>[nObs] concentrationObs;
  matrix<lower = 0>[nCmt, nEvent] mass;

  for (j in 1:nSubjects) {
    mass[, start[j]:end[j]] = pmx_solve_twocpt(time[start[j]:end[j]],
                                               amt[start[j]:end[j]],
                                               rate[start[j]:end[j]],
                                               ii[start[j]:end[j]],
                                               evid[start[j]:end[j]],
                                               cmt[start[j]:end[j]],
                                               addl[start[j]:end[j]],
                                               ss[start[j]:end[j]],
                                               theta[j, ], biovar, tlag);

    concentration[start[j]:end[j]] = 
                      mass[2, start[j]:end[j]] / theta[j, 3];
  }

  concentrationObs = concentration[iObs];
}

model {
  // priors
  CL_pop ~ lognormal(log(10), prior_sd[1]); 
  Q_pop ~ lognormal(log(15), prior_sd[2]);
  VC_pop ~ lognormal(log(35), prior_sd[3]);
  VP_pop ~ lognormal(log(105), prior_sd[4]);
  ka_pop ~ lognormal(log(2.5), prior_sd[5]);
  sigma ~ normal(0, 1);
  omega ~ lognormal(1, prior_sd);

  // hierarchical prior
  for (j in 1:nSubjects)
    theta[j, ] ~ lognormal(log(theta_pop), omega);

  // likelihood
  logCObs ~ normal(log(concentrationObs), sigma);
}


generated quantities {
  real concentrationObsPred[nObs] 
    = exp(normal_rng(log(concentrationObs), sigma));
}

