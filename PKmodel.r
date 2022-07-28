## Clear the environment
rm(list = ls())
gc()

## Set the environment work and import libraries

# Set working path and path to libraries
setwd("C:/Users/stanislas/Documents/R/Columbia/project/pharmacometrics/Internship_Columbia")
.libPaths("C:/Users/stanislas/Documents/R/win-library/4.1")

# necessary package
library(ggplot2)
library(cmdstanr)
library(rjson)
library(posterior)
library(deSolve)
library(bayesplot)
library(tidyverse)
library(gridExtra)
library(comprehenr)
library(loo)
library(docstring)
library(devtools)


## In particular models, use the appropriate cmdstan - Here use Torsten
set_cmdstan_path("C:/Users/stanislas/Documents/R/Torsten/cmdstan")

## Import useful functions to run this Rscript
source("utils/simulations.R")
source("utils/plot_functions.R")

set.seed(12304) ## randomness in the program, set to have the same for all chains


#########################################################################
## Define your model

model_name <- "Michaelis_MentenPK"
mod <- cmdstan_model(paste0("model/",model_name, ".stan"))
saved_name <- "Michaelis_MentenPK_30b"

##########################################################################
## Define the Data and make it stan suitable
## Adapted to Michaelis_Menten model
simulate_data <- FALSE
if (simulate_data) {
  ka <- 3
  V <- 27
  Vm <- 10
  Km <- 14
  sigma <- 0.1
  y0 <- 100
  
  parm <- c(ka, V, Vm, Km)
  state <- c(Y1 = y0, Y2 = 0) # initial condition
  
  dydt <- function(t, state, parm) { # ODE based model, here : Michaelis-Menten model
    with(as.list(c(state, parm)), {
      C = Y2 / V
      dY1 <- - ka * Y1
      dY2 <- ka * Y1 - Vm * C / (Km + C)
      
      list(c(dY1, dY2))
    })
  }
  
  times <- c(seq(0.5, 1.5, by = 0.25), seq(10, 100, by = 10))
  n <- length(times)
  
  mass <- ode(y = state, time = times, func = dydt, parms = parm)
  concentration <- mass[, 3] / V
  y <- rnorm(n, concentration, sigma) #c_obs
  
  stan_data <- list(N = n,  y = y, t = times, y0 = c(y0, 0))
  write_stan_json(stan_data, paste0("data/",model_name,".data.json"))
  
  p <- ggplot(data = data.frame(y = y, times = times),
              aes(y = y, x = times)) +
    geom_point() + theme_bw() + theme(text = element_text(size = 20))
  p
}

# import the data to feed the model
stan_data_rk45 <- fromJSON(file = paste0("data/",model_name,".data.json"))
stan_data_rk45$stiff_solver <- 0
stan_data_bdf <- stan_data_rk45
stan_data_bdf$stiff_solver <- 1

print(stan_data_rk45)

##########################################################################
## Define the parameters and the initialisation

# Draw initial conditions from the prior
pars = c("ka","V","Vm","Km","sigma")


## Draw a markov Chain sample from the model prior distribution
init <- function() {
  list(
    ka = exp(rnorm(1, log(1.5), 3)),
    V = exp(rnorm(1, log(35), 0.5)),
    Vm = exp(rnorm(1, log(10), 0.5)),
    Km = exp(rnorm(1, log(2.5), 3)),
    sigma = abs(rnorm(1, 0, 1))
  )
}


nChains <- 4
parallel_chains <- 4

for (i in 1:nChains) {
  init0 <- init()
  write_stan_json(init0, paste0("init/",saved_name, i, ".json"))
}

init0_files <- paste0("init/",saved_name, 1:nChains, ".json")

##########################################################################
## enter the parameters of the different phases needed to run the MCMC
iter_warmup <- 500
iter_sampling <- 500
stan_seed <- 45678
metric <- "diag_e"
step_size <- NULL
inv_metric <- NULL


###########################################################################
## enter the model condition
run_model <- TRUE
run_posterior <- FALSE
run_parameters <- FALSE

###########################################################################
# run_simulations in utils/simulations.R
## Fit model with rk45 solver.
integrator = "rk45"

fits0 <- run_simulations(mod,saved_name, integrator,run_model,run_posterior,run_parameters, init,stan_data_rk45,
                         iter_warmup, iter_sampling,stan_seed,metric,step_size,inv_metric,
                         nChains,parallel_chains) 

## Fit model with bdf solver. 

bdf = "bdf"
fits1 <- run_simulations(mod,saved_name,  integrator = bdf,run_model,run_posterior,run_parameters, init,stan_data_bdf,
                         iter_warmup, iter_sampling,stan_seed,metric,step_size,inv_metric,
                         nChains,parallel_chains)

###########################################################################
## parameters to analyse
parms <- c("ka","V","Vm","Km","sigma", "lp__")
lp_index <- match("lp__", parms)
integrators <- c("rk45","bdf")
n_drop <- 10
min_chains <- 12

# from utils/simulations.R
fit_objects <- list(fit_object(fits0,"rk45"),fit_object(fits1,"bdf"))
fit_objects0 = fits_all_ini(saved_name,"rk45")
fit_objects1 = fits_all_ini(saved_name,"bdf")

for (fitted in fit_objects){
  print(fitted$weights)
}

###################################################################################
## Posterior analysis and influence of initialisation on performances
######################################################################################



## plot analyses - utils/plot_functions.R
?plot_posterior
plot_posterior(saved_name, integrators, run_posterior, run_parameters,fit_objects,parms, iter_sampling)



## plot the markov chains samples - from Bayesplot package
color_scheme_set("mix-blue-red")
mcmc_trace(fits0$fit$draws(), pars = parms)
mcmc_trace(fits1$fit$draws(), pars = parms)


#####################################################################
# Examine global efficiency of parameters - utils/plot_functions.R

plot_ess(saved_name, run_posterior, run_parameters, fit_objects,parms)


# Compare global efficiency of parameters for all type of initialisations - utils/plot_function.R
plot_ess_ini(fit_objects0, "rk45",parms)
plot_ess_ini(fit_objects1, "bdf",parms)

## RK45 and BDF full review - utils/plot_function.R
plot_phase_time(saved_name,fit_objects, integrators, run_parameters, run_posterior)


## Dropping chains - utils/plot_function.R
plot_ess_chains_para(saved_name, run_posterior, run_parameters, fit_objects,integrators,parms,n_drop)



#### ESS/s en fonction des chaines extraites - utils/plot_function.R
plot_integrator_lp_summary(saved_name,integrators, parameter= "lp__")


## plot relaxation time - utils/plot_function.R
plot_lp_tau(saved_name,integrators,parameter = "lp__")


## plot efficacity and ess_bulk of chains kept - utils/plot_function.R
plot_drop_chain_distribution(saved_name, fit_objects,run_posterior,run_parameters,min_chains)
plot_ess_chain_distribution(saved_name, fit_objects,run_posterior,run_parameters,min_chains)


## plot Rhat of chains kept - utils/plot_function.R
plot_Rhat(saved_name, fit_objects,run_posterior,run_parameters,min_chains)


## plot the weighted markov chain samples - from Bayesplot package
mcmc_trace(fit_object(fits0,"rk45")$weight_draws, pars =parms)



#############################################################################
## pseudo-algorithm 
#############################################################################


fit_object0 <- fit_object(fits0,"rk45")
fastest_chains <- fit_object0$fastest_chain

# extract divergent chains
diagnostic = fits0$fit$sampler_diagnostics(format = "df")
non_mixing_chains = unique(diagnostic[(diagnostic$"divergent__" == 1),]$.chain)
non_mixing_chains


n<-15 # number of observations
n_min <- 5 #number min of chains we want in final stacking
ESS_target = 1000 # ESS target

new_ess <- NULL
log_liks = paste0("log_lik[",c(1:n),"]") # to extract the pointwise log likelihood
ess_weight_chains <- c()
all_ess = rep(0,30-n_min)


for (i in n_min:30){
  study_chains <- as_draws_array(fits0$fit$draws("log_lik")[,fastest_chains[1:i],]) # extract pointwise log likelihood
  new_inv_ess <- 1/dplyr::pull(summarize_draws(fits0$fit$draws("lp__")[,fastest_chains[i],], measure = c("ess_bulk"))) #extract new ESS
  
  if (is.na(new_inv_ess)) new_inv_ess = 0
  ess_weight_chains <- c(ess_weight_chains, new_inv_ess)
  weights <- chain_stacking(fits0,"log_lik",fastest_chains[1:i]) ## realise chain stacking - utils/simulations.R
  weight_chains <- mixture_draws(study_chains,weight = weights) ## extract mixture draw - utils/simulations.R
  ESS_measure <- 1/sum( weights**2 * ess_weight_chains) ## calculate new ESS (see stacking paper for formula)
  pareto_k_measure <- loo(weight_chains)$diagnostics$pareto_k ## realise pareto-k diagnostic
  all_ess[i-n_min] = ESS_measure
  if (ESS_measure> ESS_target && sum((pareto_k_measure>0.7))==0 ) break
}


print("###########################")
print(paste0("we stop at ",i," chains"))
print(paste0("the result ESS is ", ESS_measure))
print(pareto_k_measure )
total_running_time =  fits0$fit$time()$chains[fastest_chains[i],4]+
  fits0$fit_w1$time()$chains[fastest_chains[i],4]+
  fits0$fit_w2$time()$chains[fastest_chains[i],4]+
  fits0$fit_w3$time()$chains[fastest_chains[i],4]
print(paste0("total running time : ",total_running_time))


####### Uncomment and Run all the pseudo algorithm without stopping condition to plot the evolution of ESS through it
# ggplot(data.frame(chain_finished = c(n_min:29), ESS_stacked = all_ess), aes(x= chain_finished, y = ESS_stacked))+
#   geom_point() + theme_bw() + geom_hline(yintercept = ESS_measure)



##########################################################################
################### stacking on all chains ###############################

weights_pseudo = chain_stacking(fits0, study_chains = fastest_chains, method = "pseudobma")
weights_stack = chain_stacking(fits0, study_chains = fastest_chains, method = "stacking")

## plot of weights from pseudo-bma
ggplot(data = data.frame(fastest_chains = c(1:30), weights = c(weights_pseudo)), aes(x = fastest_chains,y=weights))+
  geom_bar(stat = "identity") + coord_flip() + 
  theme(text = element_text(size = 40)) + theme_bw()
#mcmc_trace associated

mixture_draws((fits0$fit$draws(parms)),weights_pseudo)


## plot of weights from stacking
ggplot(data = data.frame(fastest_chains = c(1:30), weights = c(weights_stack)), aes(x = fastest_chains,y=weights))+
  geom_bar(stat = "identity") + coord_flip()+
  theme(text = element_text(size = 40)) + theme_bw()



############ distribution of weighted chains ################
samples <- weight_draws(fits0, study_chains[1,,1])$weight_draws[1:length(parms)]
samples_list = c()
for (par in names(samples)){
  samples_list = c(samples_list, samples[[par]])
}
parm_name = rep(parms, each = dim(samples)[1])

plot_data = data.frame(sample = samples_list, parm = parm_name)

ggplot(plot_data, aes(x=sample))+
  geom_histogram(alpha = 0.25, position = "identity", bins = 30) +
  theme_bw() + facet_wrap(~ parm, scale = "free")



############## Analyse of the weighted chains using Bayesplot package ########################
mcmc_hist(weight_draws(fits0, study_chains[1,,1])$weight_draws[parms])
mcmc_areas(weight_draws(fits0, study_chains[1,,1])$weight_draws["lp__"], prob = 0.8)


#################### Predictive Posterior Checks   ################################"
## extract y_pred
y_preds = paste0("y_pred[",c(1:15),"]")
y_pred_rk45 = as.matrix(weight_draws(fits0, study_chains[1,,1])$weight_draws[y_preds])
ppc_hist(stan_data_rk45$y,y_pred_rk45[1:8], binwidth = 1)
ppc_dens_overlay(stan_data_rk45$y, y_pred_rk45[1:50, ])


## ppc plot
yrep <- as.matrix(posterior::as_draws_df(weight_draws(fits0, study_chains[1,,1])$weight_draws[y_preds]))[, -(16:18)]
yobs <- stan_data_rk45$y
t <- stan_data_rk45$t
p <- bayesplot::ppc_ribbon(y = yobs, yrep = yrep, x = t) + ggplot2::xlab("time (h)") + ggplot2::ylab("plasma concentration (mg/L)")
p



##########################################################################
## without stacking on the 14 fastest chains -- Posterior distribution
samples <- as_draws_df(fits0$fit$draws(parms)[,fastest_chains[1:14],])[1:length(parms)]
samples_list = c()
for (par in names(samples)){
  samples_list = c(samples_list, samples[[par]])
}
parm_name = rep(parms, each = dim(samples)[1])
plot_data = data.frame(sample = samples_list, parm = parm_name)

ggplot(plot_data, aes(x=sample))+
  geom_histogram(alpha = 0.25, position = "identity", bins = 30) +
  theme_bw() + facet_wrap(~ parm, scale = "free")


################################################################################
## Comparing results with BDF for the same computational time 

## extract chains finished in computational time
bdf_fit_model <- fit_object(fits1,"bdf")
chains_below_rk45_runtime = bdf_fit_model$chain_time[bdf_fit_model$chain_time< total_running_time]
study_chains_bdf <- bdf_fit_model$fastest_chain[1:length(chains_below_rk45_runtime)]

## posterior draw
samples <- as_draws_df(fits1$fit$draws(parms)[,study_chains_bdf,])[1:length(parms)]
samples_list = c()
for (par in names(samples)){
  samples_list = c(samples_list, samples[[par]])
}
parm_name = rep(parms, each = dim(samples)[1])
plot_data = data.frame(sample = samples_list, parm = parm_name)

ggplot(plot_data, aes(x=sample))+
  geom_histogram(alpha = 0.25, position = "identity", bins = 30) +
  theme_bw() + facet_wrap(~ parm, scale = "free")


## ESS
dplyr::pull(summarize_draws(fits1$fit$draws("lp__")[,study_chains_bdf,], measure = c("ess_bulk")))

## pareto-k
loo(as_draws_array(fits1$fit$draws("log_lik")[,study_chains_bdf,]))

## running time 
fit_object(fits1,"bdf", study_chains_bdf)$total_time

## trace
mcmc_trace(fits1$fit$draws()[,study_chains_bdf,], pars = parms)


############# plot the evolution of the pareto-k during stacking ############################"
n<-0
n_min <- 5
ESS_min = 1000

log_liks = paste0("log_lik[",c(1:n),"]")
pareto_k <- matrix(0,30,15)

for (i in 2:30){
  study_chains <- as_draws_array(fits0$fit$draws("log_lik")[,fastest_chains[1:i],])
  weights <- chain_stacking(fits0,"log_lik",fastest_chains[1:i])
  weight_chains <- mixture_draws(study_chains,weight = weights)
  print(loo(weight_chains))
  
}

data_frame_pareto <- as.data.frame(pareto_k)
data_frame_pareto$id <- 1:nrow(data_frame_pareto)
library(reshape2)
pareto_long = melt(data_frame_pareto, id.vars = "id")

ggplot(pareto_long, aes(x = value, y = id)) + 
  ylab("n fastest chains") + xlab("Pareto k diagnostic value") +
  geom_point() + geom_vline(xintercept=0.7, linetype="dashed", color = "red")+
  facet_wrap(~ variable) + annotate(geom = "text",label = "k = 0.7", x = 0.72,y =15,angle = 90, vjust = 1)


##########################################################
####################  convergence test - posterior distribution   ######################
eps <- 0.05
measure <- 0 
new_measure <- 0
step = 5
n_min = 8
for (i in seq(n_min, nChains, step)){
  new_draw = weight_draws(fits0,fastest_chains[1:i],"log_lik")$weight_draws["lp__"][,1]
  new_measure <- quantile(new_draw,probs = seq(0.2,0.8,0.1))
  print(paste0("###### iteration",i," ######"))
  print(paste0("max relative error ",max(abs(measure - new_measure)/measure)))
  
  print(paste0("kolmogorov-smirnov test ",ks.test(measure, new_measure)))
  measure = new_measure
}


##########################################################

#################  Rhat wagon evolution inde chains  ######################
## idea : because the rhat is small, it might be because there chains are too dependent
## we delete this dependent and realise stacking on different chains

n_min <- 4
n_wagon <- 4
rhat_evolution <- rep(0,nChains)
k=1
wagon <- as_draws_array(fits0$fit$draws("lp__")[,c(1:n_wagon),])
for (i in seq(2*n_min,30,n_wagon)){
  for (j in 1:n_wagon){
    wagon[,j,]= weight_draws(fits0,fastest_chains[seq(j,i,n_wagon)],"log_lik")$weight_draws["lp__"][,1]
  }
  if (i >= n_wagon + n_min){
    rhat_evolution[k] <- dplyr::pull(summarize_draws(wagon,measure = c("rhat")),2)
    k=k+1
  }
}

rhat_evolution = rhat_evolution[(rhat_evolution!=0)]
plot_data = data.frame(evaluation = c(1:length(rhat_evolution)), rhat = rhat_evolution)
ggplot(plot_data, aes(x= rhat, y = evaluation)) + geom_point() + theme_bw()

mcmc_trace(wagon)


###################################################################
## convergence test of the stability of the log predictive density
lpd_loo  <- rep(0,29)
for (i in 2:30){
  stacked = weight_draws(fits0,fastest_chains[1:i],"log_lik")$weight_draws[log_liks]
  lpd_loo[i-1] = (loo(data.matrix(stacked))$estimate[1,1])
}

ggplot(data = data.frame(chain_finished = c(2:30), lpd_loo = lpd_loo), aes(x = chain_finished, y=lpd_loo))+
  geom_point() + theme_bw()

#############################################################

## Rhat of the slowest chains
slow_chains <- as_draws_array(fits0$fit$draws(c("ka","lp__"))[,,])
dplyr::pull(summarize_draws(slow_chains,measure = c("ess_bulk")),2)
color_scheme_set("mix-blue-red")
mcmc_trace(slow_chains)


############### step - size plotting with red for divergent chains ##############
step_sizes = fits0$fit$metadata()$step_size_adaptation
ggplot()+
  geom_point(data = data.frame(chain_id = fastest_chains, step_size = step_sizes[fastest_chains]),aes(x = step_size, y = chain_id, colour = "Mixing chains")) + 
  geom_point(data = data.frame(chain_id = non_mixing_chains, step_size =step_sizes[non_mixing_chains]), aes(x = step_size, y = chain_id, colour = "Biased chains"))+
  scale_colour_manual(values = c("red","black"))+ theme_bw()



############### mass metrix with red for divergent chains ########################
fits0$fit$inv_metric(F)
