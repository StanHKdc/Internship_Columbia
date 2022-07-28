
library(dplyr)
library(wesanderson)
library(docstring)

# TODO: Create a saved_last_iter function which works for every model.
# A function to extract the final sample of warmup (first sample of
# "sampling" phases) and save it.
# Designed for the Michaelis-Menten PK model for one individual.


save_last_iter <- function(samples, nChains, simulation_name, pars) {
  #' @title intermediate function
  #' @description Save the last iteration 
  
  
  saved_file_name <- paste0("init/", simulation_name, 1:nChains, ".json")
  for (i in 1:nChains) {
    init_saved <- c(as_draws_df(samples[1,i,]))[1:length(pars)+1]
    write_stan_json(init_saved, saved_file_name[i])
  }
  saved_file_name
}


fits <- function(model_name,integrator, run_posterior,run_parameters){
  #' @description Returns a list containing all the fits for a simulation previously done for a particular initialisation
  #' @param model_name String - model we study
  #' @param integrator String - "rk45" or "bdf"
  #' @param run_posterior Boolean - if FALSE : the Markov Chain initializations are drawn from the prior; 
  #' if TRUE : the Markov Chain initializations are drawn from a posterior distribution of the model ran (necessity to have a model already ran)
  #' @param run_parameters Boolean - if FALSE : HMC parameters are randomly initialized
  #' if TRUE : HMC parameters are initialized with the tuned parameters of a previously ran model
  
  saved_fit_file <- paste0("output/",model_name,integrator,"_post",run_posterior,"_tuned", run_parameters)
  fit_w1 <- readRDS(paste0(saved_fit_file,".phase1.fit.RDS"))
  fit_w2 <- readRDS(paste0(saved_fit_file,".phase2.fit.RDS"))
  fit_w3 <- readRDS(paste0(saved_fit_file,".phase3.fit.RDS"))
  fit <- readRDS(paste0(saved_fit_file,".RDS"))
  
  list(fit = fit,
       fit_w1 = fit_w1,
       fit_w2 = fit_w2,
       fit_w3= fit_w3)
}

fit_object <- function(fitted, integrator, study_chains = NULL){
  #' @description Returns all the basic information of a fitted model for an integrator and chains wanted
  #' @param fitted list of fitted model - output of fits function just above
  #' @param integrator String - "rk45" or "bdf"
  #' @param study_chains vector of integer - chains we want to extract the weights
  #' @return list with the following information : 
  #' - draw = sample draw
  #' - total_time = method total computational time
  #' - chain_time = computational time of each chain
  #' - phase_chain_time = computational time of each chain for each phase
  #' - integrator = "rk45" or "bdf"
  #' - num_chains = number of chains to be studied
  #' - weight = weights on the chains we want to study
  #' - weights_draw = mixture draw on the study chains according to the weights
  #' - fastest_chain = chains ordered by their computational time
  #' 
  if (is.null(study_chains)) study_chains = c(1:fitted$fit$num_chains())
  weight_draws = weight_draws(fitted, study_chains)
  chain_time = (fitted$fit$time()$chains[, 4] +
                  fitted$fit_w1$time()$chains[, 4] +
                  fitted$fit_w2$time()$chains[, 4] +
                  fitted$fit_w3$time()$chains[, 4])
  list(draw = fitted$fit$draws(), total_time = (fitted$fit$time()$total +
                                                  fitted$fit_w1$time()$total +
                                                  fitted$fit_w2$time()$total +
                                                  fitted$fit_w3$time()$total),
              chain_time = chain_time, 
       phase_chain_time = c(fitted$fit_w1$time()$chains[, 4] ,
                            fitted$fit_w2$time()$chains[, 4] ,
                            fitted$fit_w3$time()$chains[, 4],
                            fitted$fit$time()$chains[, 4]),
       integrator = integrator,num_chains = fitted$fit$num_chains(),
       weight_draws = weight_draws$weight_draws ,
       weights = weight_draws$weights,
       fastest_chain = order(chain_time)
)
  
}

fits_all_ini <- function(saved_name,integrator){
  #' @description extract the output of fit_object for all different initialisation done
  #' @param saved_name string - model we study
  #' @param integrator string - "rk45" or "bdf"
  #' @return list of the output of fit_object in the following order ("FF","TF","FT","TT") @seealso utils/plot_function/plot_lp_tau for more information on these acronym
  
  fitsFF <- fit_object(fits(saved_name,integrator,FALSE,FALSE),integrator)
  fitsVF <- fit_object(fits(saved_name,integrator,TRUE,FALSE),integrator)
  fitsFV <- fit_object(fits(saved_name,integrator,FALSE,TRUE),integrator)
  fitsVV <- fit_object(fits(saved_name,integrator,TRUE,TRUE),integrator)
  
  list(fitsFF,fitsVF,fitsFV,fitsVV)
}





# A function to run the warmup and save the fit for each phase.
# save_function should either be save_last_iter or save_last_iter_pop.
fit_warmup_phases <- function(mod, simulation_name, data, 
                              init = NULL,
                              iter_warmup = 500, 
                              init_buffer = 75,
                              term_buffer = 50, nChains = 4, parallel_chains = 4,
                              seed = 123,
                              metric = "diag_e",
                              adapt_delta = 0.8,
                              inv_metric = NULL,
                              step_size = NULL,
                              ...) {
  
  #' @title fitting of the warm up phases
  #' @description This function highlights the 3 subphases of the warm-up phase described in the stan documentation
  #' We save the results at the end of each subphase so we can study them later and understand better how the warm-up phase goes on
  #' 
  #' @param mod - stan model already computated
  #' @param simulation_name - name of the simulation (model_saved + info on initialisation)
  #' @param data json - stan data   #' @param init - markov chain initialisation
  #' @return a list of the sampling of the 3 subphases of the warm-up phase

  
  pars <- names(init())
  
  ## sample for the 1st subphase of warmup
  fit_w1 <- mod$sample(data = data, init = init,
                       chains = nChains, parallel_chains = parallel_chains,
                       iter_warmup = init_buffer,
                       iter_sampling = 1,
                       init_buffer = init_buffer,
                       term_buffer = 0, window = 0,
                       metric = metric,
                       seed = seed, refresh = 10, adapt_delta = adapt_delta,
                       inv_metric = inv_metric,step_size = step_size)
  
  samples_w1 <- fit_w1$draws(variables = pars)
  nChains <- length(samples_w1[1,,1])
  
  ## save the results
  init_files <- save_last_iter(samples_w1,nChains ,
                              paste0(simulation_name, "_w1."),pars)
  
  if (metric == "diag_e") {
    metric_is_matrix <- FALSE
  } else {
    metric_is_matrix <- TRUE
  }
  
  ## sample the 2nd subphase of warmup, initialise at the end of 1st phase
  window_size <- iter_warmup - (init_buffer + term_buffer)
  fit_w2 <- mod$sample(data = data, init = init_files,
                       chains = nChains, parallel_chains = parallel_chains,
                       iter_sampling = 1,
                       iter_warmup = window_size,
                       init_buffer = 0, term_buffer = 0,
                       metric = metric,
                       seed = seed,
                       step_size = fit_w1$metadata()$step_size_adaptation,
                       inv_metric = fit_w1$inv_metric(matrix = metric_is_matrix),
                       refresh = 10, adapt_delta = adapt_delta)
  ## save the results
  
  samples_w2 <- fit_w2$draws(variables = pars)
  nChains <- length(samples_w2[1,,1])
  init_files2 <- save_last_iter(samples_w2, nChains, paste0(simulation_name, "_w2."),
                               pars)
  
  ## sample the 3rd subphase of warmup, initialise at the end of 2nd phase
  fit_w3 <- mod$sample(data = data, init = init_files2,
                       chains = nChains, parallel_chains = parallel_chains,
                       iter_sampling = 1,
                       iter_warmup = term_buffer, window = 0,
                       init_buffer = 0, term_buffer = term_buffer,
                       metric = metric,
                       seed = seed,
                       step_size = fit_w2$metadata()$step_size_adaptation,
                       inv_metric = fit_w2$inv_metric(matrix = metric_is_matrix),
                       refresh = 10, adapt_delta = adapt_delta)
  
  # Return the results
  list(fit_w1 = fit_w1, fit_w2 = fit_w2, fit_w3 = fit_w3)
}




run_simulations <- function (model,model_name, integrator,run_model = FALSE,run_posterior = FALSE,run_parameters = FALSE, 
                             inits= NULL, data,iter_warmup = 500,iter_sampling = 500,seed = 2022,metric = "diag_e",step_size = NULL,inv_metric = NULL, 
                             nChains = 8, parallel_chains = 8){
  #' @title Run the HMC algorithm on a model
  #' @param model stan model we want to run
  #' @param model_name string - name of the model
  #' @param run_model Boolean - FALSE: load a previous simulation // TRUE: run simulation
  #' @param run_posterior Boolean - if FALSE : the Markov Chain initializations are drawn from the prior; 
  #' if TRUE : the Markov Chain initializations are drawn from a posterior distribution of the model ran (necessity to have a model already ran)
  #' @param run_parameters Boolean - if FALSE : HMC parameters are randomly initialized
  #' if TRUE : HMC parameters are initialized with the tuned parameters of a previously ran model
  #' @param inits JSON - Markov Chain initialisations
  #' @param data JSON - model observed parameters
  #' 
  #' @description This function runs the HMC algorithm with the settings entered in the parameters. 
  #' The initialisation can be modified and taken from a previous simulation thanks to the arguments run_posterior and run_parameters.
  #' 
  #' 
  #' @return a list containing : 
  #' - fit = HMC sampling
  #' - fit_w1 = 1st warm-up phase
  #' - fit_w2 = 2nd warm-up phase
  #' - fit_w3 = 3rd warp-up phase
  
  pars = names(inits())
  saved_fit_file <- paste0("output/",model_name,integrator,"_post",run_posterior,"_tuned",run_parameters)
  simulation_name <- paste0(model_name,integrator,"_post",run_posterior,"_tuned",run_parameters)
  
  if (run_model){
    init <- inits
    if (run_posterior){ # dram MC samples from the posterior distribution of a previous simulation
      samples <- fits(model_name,integrator, FALSE,FALSE)$fit$draws(variables= pars)
      nChains <- length(samples[1,,1])
      init <- save_last_iter(samples, nChains, simulation_name, pars) 
    }
    
    if (run_parameters) { # HMC parameters taken from tuned parameters of a previous simulation
      fit_object <- fits(model_name,integrator, FALSE,FALSE)
      inv_metric <- fit_object$fit$inv_metric(matrix = F)
      step_size <- fit_object$fit$metadata()$step_size_adaptation
      nChains <- fit_object$fit$num_chains()
    }
    
    # realise the warmup phase
    fit_warmups <- fit_warmup_phases(mod=model,
                                     simulation_name = simulation_name,
                                     data = data,
                                     init = init,
                                     iter_warmup = iter_warmup,
                                     init_buffer = 75,
                                     term_buffer = 50,
                                     nChains = nChains,
                                     seed = seed,
                                     metric = metric,                                  
                                     inv_metric = inv_metric,
                                     step_size = step_size)
    
    fit_warmups$fit_w1$save_object(paste0(saved_fit_file, ".phase1.fit.RDS"))
    fit_warmups$fit_w2$save_object(paste0(saved_fit_file, ".phase2.fit.RDS"))
    fit_warmups$fit_w3$save_object(paste0(saved_fit_file,".phase3.fit.RDS"))
    
    ## sample the MC, initialise at the end of warm up phase
    samples_w3 <- fit_warmups$fit_w3$draws(variables = pars)
    nChains <- length(samples_w3[1,,1])
    init_files <- save_last_iter(samples_w3, nChains, paste0(simulation_name, "_w3."),pars)
    
    fit <- model$sample(
      data = data, chains = nChains,
      parallel_chains = parallel_chains,
      iter_warmup = 0, iter_sampling = iter_sampling,
      seed = seed, adapt_delta = 0.8,
      init = init_files,
      metric = metric,
      step_size = fit_warmups$fit_w3$metadata()$step_size_adaptation,
      inv_metric = fit_warmups$fit_w3$inv_metric(matrix = F),
      adapt_engaged = FALSE)
    
    #save the results and return the fits
    fit$save_object(paste0(saved_fit_file,".RDS"))
    fit_object <- list(fit = fit,
                       fit_w1 = fit_warmups$fit_w1,
                       fit_w2 = fit_warmups$fit_w2,
                       fit_w3 = fit_warmups$fit_w3)
    fit_object
  }
  
  else{
    fits(model_name,integrator, run_posterior,run_parameters)
  }
}


chain_stacking <- function(fits,log_lik ="log_lik",study_chains, method = "pseudobma"){
  #' @title Chain stacking
  #' @param fits - a model sampled
  #' @param study_chains vector of integers - id of chains to be stacked
  #' @param method string - way of sampling
  #' 
  #' @description This function uses the loo package to realise a chain stacking as done in Andrew's paper
  #' @return vector of weights for the study_chains entered
  
  y_pred_mat <- as.array(fits$fit$draws(log_lik)[,study_chains,])
  n = dim(y_pred_mat)[3]
  K = dim(y_pred_mat)[2]
  S = dim(y_pred_mat)[1]
  log_lik_list = list()
  for (i in 1:K){
    log_lik_list[[i]] = y_pred_mat[,i,]
  }
  r_eff_list <- lapply(log_lik_list,function(x){
    relative_eff(exp(x))
  })
  
  loo_list <- lapply(1:length(log_lik_list), function(j) {
    loo(log_lik_list[[j]], r_eff = r_eff_list[[j]])
  })
  
  loo_model_weights(loo_list,method = method,optim_control = list(reltol=1e-10))
}

mixture_draws= function (individual_draws,  weight, random_seed=1, S=NULL, permutation=TRUE){
  #' @title extract a draw from the chains
  #' @param individual_draws draw of chains - draw that we will weight
  #' @param weights vector of float - same size of the number of chains in the draw
  #' 
  #' @description this function is returning a mixture draw containing samples from each chain in individual draws
  #' according to the weight they have
  #' 
  #' @return a matrix of the size of n_samples x parameters in the draw
  #' 
  #' 
  set.seed(random_seed)
  S_sample=dim(individual_draws)[1]
  K=dim(individual_draws)[2]
  
  if(is.null(S))
    S=S_sample
  
  if(permutation==TRUE)
    individual_draws=individual_draws[sample(1:S_sample),, ]	 # random permutation of draws
  
  integer_part=floor(S*weight)
  existing_draws=sum(integer_part)
  
  if(existing_draws<S){
    remaining_draws=S-existing_draws
    update_w=(weight- integer_part/S)*  S / remaining_draws
    remaining_assignment=sample(1:K, remaining_draws, prob =update_w , replace = F)
    integer_part[remaining_assignment] =integer_part[remaining_assignment]+1
  }
  integer_part_index=c(0,cumsum(integer_part))

  mixture_vector= matrix(NA,S,dim(individual_draws)[3])
  for(k in 1:K){
    if((1+integer_part_index[k])<=integer_part_index[k+1]){
      mixture_vector[(1+integer_part_index[k]):integer_part_index[k+1],]=individual_draws[1:integer_part[k],k,]
    }
    
  }
  return(mixture_vector)
}


weight_draws <- function (fits, study_chains, log_lik = "log_lik"){
  #' @title recap of the whole stacking method in one function
  #' @param fits model sampled
  #' @param study_chains vector of integer - chains we want to stack
  #' 
  #' @description This function realise the stacking
  #' 
  #' @return Returns a list :
  #' - weight_draws = mixture draw same form as the result of model$sample()
  #' - weights = weights of the stacking of the study_chains
  #' 
  individual_draws <- as_draws_array(fits$fit$draws()[,study_chains,])
  weights <- chain_stacking(fits,log_lik,study_chains)
  new_draws <- as.data.frame(mixture_draws(individual_draws,weight = weights))
  names(new_draws) = names(as_draws_df(individual_draws))[1:length(names(new_draws))]
  list(weight_draws = new_draws, weights = weights)
}
