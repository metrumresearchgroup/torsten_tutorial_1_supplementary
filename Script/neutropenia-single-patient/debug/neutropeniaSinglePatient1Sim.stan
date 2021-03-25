functions{

    real[] twoCptNeutModelODE(real t,
			real[] x,
			real[] parms,
			real[] rdummy,
			int[] idummy){
    real k10;
    real k12;
    real k21;
    real CL;
    real Q;
    real V1;
    real V2;
    real ka;
    real mtt;
    real circ0;
    real gamma;
    real alpha;
    real ktr;
    real dxdt[8];
    real conc;
    real EDrug;
    real transit1;
    real transit2;
    real transit3;
    real circ;
    real prol;

    CL = parms[1];
    Q = parms[2];
    V1 = parms[3];
    V2 = parms[4];
    ka = parms[5];
    mtt = parms[6];	
    circ0 = parms[7];
    gamma = parms[8];
    alpha = parms[9];

    k10 = CL / V1;
    k12 = Q / V1;
    k21 = Q / V2;

    ktr = 4 / mtt;
  
    dxdt[1] = -ka * x[1];
    dxdt[2] = ka * x[1] - (k10 + k12) * x[2] + k21 * x[3];
    dxdt[3] = k12 * x[2] - k21 * x[3];
    conc = x[2]/V1;
    EDrug = alpha * conc;
    // x[4], x[5], x[6], x[7] and x[8] are differences from circ0.
    prol = x[4] + circ0;
    transit1 = x[5] + circ0;
    transit2 = x[6] + circ0;
    transit3 = x[7] + circ0;
    circ = fmax(machine_precision(), x[8] + circ0); // Device for implementing a modeled initial condition
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
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> time[nt];
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  real<lower = 0> ka;
  real<lower = 0> mtt;
  real<lower = 0> circ0;
  real<lower = 0> alpha;
  real<lower = 0> gamma;
  real<lower = 0> sigma;
  real<lower = 0> sigmaNeut;
}

transformed data{
  real<lower = 0> rate[nt] = rep_array(0.0, nt);
  real<lower = 0> ii[nt] = rep_array(0.0, nt);
  int<lower = 0> addl[nt] = rep_array(0, nt);
  int<lower = 0> ss[nt] = rep_array(0, nt);
  int<lower = 1> nCmt = 8;
  real F[nCmt] = rep_array(1.0, nCmt);
  real tLag[nCmt] = rep_array(0.0, nCmt);
}

parameters{
}

transformed parameters{
}

model{
}

generated quantities{
  vector[nt] cHat;
  vector[nt] neutHat;
  real cObs[nt];
  real neutObs[nt];
  matrix[nCmt, nt] x;
  real<lower = 0> parms[9];

  parms = {CL, Q, V1, V2, ka, mtt, circ0, gamma, alpha};

  x = pmx_solve_rk45(twoCptNeutModelODE, 8,
                     time, amt, rate, ii, evid, cmt, addl, ss,
                     parms, F, tLag,1e-6, 1e-6, 1e8);

  cHat = x[2, ]' / V1;
  neutHat = x[8, ]' + circ0;

  for(i in 1:nt){
    if(time[i] == 0){
      cObs[i] = 0;
    }else{
      cObs[i] = exp(normal_rng(log(fmax(machine_precision(), cHat[i])),
        sigma));
    }
    neutObs[i] = exp(normal_rng(log(fmax(machine_precision(), neutHat[i])),
      sigmaNeut));
  }

}
