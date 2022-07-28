functions {
  vector system(real time, vector y, real[] theta) {
    real ka = theta[1];
    real V = theta[2];
    real Vm = theta[3];
    real Km = theta[4];
    real C = y[2] / V;
    vector[2] dydt;

    dydt[1] = - ka * y[1];
    dydt[2] = ka * y[1] - Vm * C / (Km + C);
    return dydt;              
  }
}

data {
  int N;
  vector[N] y;
  array[N] real t;
  int stiff_solver;
  vector[2] y0;  // initial drug mass in the gut.
}

transformed data {
  int n_parm = 4;
  int n_cmt = 2;
  real t0 = 0;
}

parameters {
  real<lower = 0> ka;
  real<lower = 0> V;
  real<lower = 0> Vm;
  real<lower = 0> Km;
  real<lower = 0> sigma;
}

transformed parameters {
  array[4] real theta = {ka, V, Vm, Km};
  vector[N] concentration;

  {
    array[N] vector[n_cmt] x;
    if (!stiff_solver) {
      x = ode_rk45(system, y0, t0, t, theta);
    } else {
      x = ode_bdf(system, y0, t0, t, theta);
    }

    concentration = to_vector(x[, 2]) / V;
  }
}

model {
  // priors -- CHECK if these are reasonable values
  ka ~ lognormal(log(2.5), 3);  // 1.85
  V ~ lognormal(log(35), 0.5);
  Vm ~ lognormal(log(10), 0.5);
  Km ~ lognormal(log(2.5), 3);  // 2
  sigma ~ normal(0, 1);

  // likelihood
  y ~ normal(concentration, sigma);
}

generated quantities {
  real y_pred[N] = normal_rng(concentration, sigma);
  vector[N] log_lik;
  for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | concentration[n], sigma);
}
