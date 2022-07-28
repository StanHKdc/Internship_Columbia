library(wesanderson)
library(docstring)
source("utils/simulations.R")


fits <- function(saved_name,integrator, run_posterior,run_parameters){
  #' @title Retrieve a model saved
  #' @description After simulating a model and saving the phases of the warm up 
  #' and its sampling phase, we retrieve it to do analysis on them. 
  #' As convention, the models are saved in a folder "output/" and have the following name : 
  #' "{saved_name}{integrator}_post{boolean}_tuned{boolean}"
  #' @param saved_name string - name of the model ran
  #' @param integrator string : "rk45" or "bdf" - Runge Kutta 4/5 (rk45) or Backward Differentiate Formula (bdf)
  #' @param run_posterior Boolean - if FALSE : the Markov Chain initializations are drawn from the prior; 
  #' if TRUE : the Markov Chain initializations are drawn from a posterior distribution of the model ran (necessity to have a model already ran)
  #' @param run_parameters Boolean - if FALSE : HMC parameters are randomly initialized
  #' if TRUE : HMC parameters are initialized with the tuned parameters of a previously ran model
  #' @return This function returns a list(fit,fit_w1,fit_w2,fit_w3) where fit corresponds to the sampling phase and 
  #' (fit_w1,fit_w2,fit_w3) the 3 subphases of the sampling phase
  
  saved_fit_file <- paste0("output/",saved_name,integrator,"_post",run_posterior,"_tuned", run_parameters)
  fit_w1 <- readRDS(paste0(saved_fit_file,".phase1.fit.RDS"))
  fit_w2 <- readRDS(paste0(saved_fit_file,".phase2.fit.RDS"))
  fit_w3 <- readRDS(paste0(saved_fit_file,".phase3.fit.RDS"))
  fit <- readRDS(paste0(saved_fit_file,".RDS"))
  
  list(fit = fit,
       fit_w1 = fit_w1,
       fit_w2 = fit_w2,
       fit_w3= fit_w3)
}

ess_summary <- function(fit_object, parms) {
  #' @title Ess and Eff summary
  #' @return Returns a list containing the ESS and the Efficiency (ESS/running time) for each chain and the parameters wanted for a
  #' particular fit_model object
  #' @param fit_object object(utils/simulations.R) - informations of the simulations we want to extract the data
  #' @param parms string or vector of strings - parameters we want to study
  #' @description We extract for each chain its ess_bulk (posterior package) and then we stock it and calculate the efficiency of the chain 
  time <-fit_object$chain_time
  n_parm <- length(parms)
  draws <- fit_object$draw[,,parms]
  nChains <- dim(fit_object$draw)[2]
  chain_ess <- matrix(NA, n_parm, nChains)
  chain_eff <- matrix(NA, n_parm, nChains)
  for (i in 1:nChains) {
    chain_draw <- as_draws_df(draws[, i, ])
    chain_ess[, i] <- dplyr::pull(summarize_draws(chain_draw, 
                                                  measure = c("ess_bulk")), 2) 
    chain_eff[, i] <- chain_ess[, i] / time[i]
  }
  
  list(chain_ess = chain_ess, chain_eff = chain_eff)
}



plot_phase_time <- function(saved_name,fit_objects, integrators,run_posterior,run_parameters){
  #' @title Plot the running time of chains
  #' @description This function draw the running time of each chain for both integrators and a particular initialisation,
  #' highligthing the phases
  #' 
  #' @param saved_name string - name of the model ran
  #' @param fit_object list() - list of fits returned by fits() /// see also à rajouter
  #' @param integrators vector - integrators to be studied
  #' @param run_posterior Boolean - Initialisation of Markov chains from the posterior distribution of a previous model
  #' @param run_parameters Boolean - Initialisation of HMC parameters from a tuned parameters of a previous model

  
  time_all_data <- c()
  phase <- c()
  chain_id <- c()
  integrator <- c()
  for (fit_object in fit_objects) {
    time <- fit_object$phase_chain_time
    time_all_data <- c(time_all_data,time)
    phase <- c(phase,rep(c("phase 1", "phase 2", "phase 3", "sampling"), each = fit_object$num_chains))
    chain_id <- c(chain_id,rep(1:fit_object$num_chains,4))
    integrator <- c(integrator,rep(fit_object$integrator,length(time)))
  }
  
  time_all_data[time_all_data>1000] = 1000
  phase <- factor(phase, levels = c("sampling", "phase 3", "phase 2", "phase 1"))
  plot_data <- data.frame(run_times = time_all_data,
                          chain_id = chain_id,
                          phase = phase,
                          integrator = integrator)
  
  plot <- ggplot(data = plot_data, 
                 aes(x = chain_id, y = run_times, fill = phase))+
    geom_bar(stat = "identity", position = "stack", width = 0.3) +
    coord_flip() + theme_bw() +
    ylab("Run time (s)") + xlab("Chain ID") +
    theme(text = element_text(size = 12)) +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "coral")) + 
    facet_wrap(~integrator)
  print(plot)
  ggsave(plot, filename =  paste0("figures/",saved_name,"phase time ",run_posterior,"_tuned",run_parameters,".png"))
  
}

plot_lp_tau <- function(saved_name,integrators,parameter){
  #' @title Plot the relaxation time of 'lp__' parameter for each chain for both integrators and types of initialisation
  #' @description This function is drawing the relaxation time (total running time / Effective Sample Size) for each chain.
  #' To be run, this function requires that the simulations for all type of initialisation have previously been done. The different types of initialisation are : 
  #' - MC sample from the Prior Distribution + Default HMC parameters (FF)
  #' - MC sample from a Posterior Distribution + Default HMC parameters (TF)
  #' - MC sample from the Prior Distribution + Tuned HMC parameters (FT)
  #' - MC sample from a Posterior Distribution + Tuned HMC parameters (TT)
  #' 
  #' For the last 3 initialisation configuration, we need to first run the first simulation to realise them. 
  #' 
  #' @param saved_name string - name of the model ran
  #' @param integrators vector - integrators to be studied, usually BDF and RK45
  #' @param parameter string - Parameter to be studied
  #' @seealso plot_integrator_lp_summary
  #' 
  #' @return Plots of the relaxation time and log relaxation time (great if outliers make the graph difficult to read)
  #' 
  recorded_tau <- c()
  method_names <- c("VV","VF","FV","FF")
  method <- c()
  weight_tau <- c()
  weight_method <- c()
  recorded_tau <- c()
  
  for (integrator in integrators){
    fitsFF <- fit_object(fits(saved_name,integrator,FALSE,FALSE),integrator)
    fitsVF <- fit_object(fits(saved_name,integrator,TRUE,FALSE),integrator)
    fitsFV <- fit_object(fits(saved_name,integrator,FALSE,TRUE),integrator)
    fitsVV <- fit_object(fits(saved_name,integrator,TRUE,TRUE),integrator)
    
    essFF <- ess_summary(fitsFF,parameter)
    essFV <- ess_summary(fitsFV,parameter)
    essVF <- ess_summary(fitsVF,parameter)
    essVV <- ess_summary(fitsVV,parameter)
    recorded_tau <- c(recorded_tau, 1/essFF$chain_eff[1,],1/essVF$chain_eff[1,],1/essFV$chain_eff[1,],1/essVV$chain_eff[1,])
    
    w_effFF <- w_eff(fitsFF,essFF$chain_ess[1,])
    w_effFV <-  w_eff(fitsFF,essFV$chain_ess[1,])
    w_effVF <-  w_eff(fitsFF,essVF$chain_ess[1,])
    w_effVV <-  w_eff(fitsFF,essVV$chain_ess[1,])
    
    weight_tau <- c(weight_tau, 1/w_effFF,1/w_effVF,1/w_effFV, 1/w_effVV)
    
    method <- c(method,rep(paste0(integrator,"-", method_names),each =fitsFF$num_chains))
    weight_method <- c(weight_method, paste0(integrator,"-", method_names))    
    
  }
  median_tau <- median(recorded_tau[!is.na(recorded_tau)])

  plot_data <- data.frame(method = method, tau = recorded_tau)
  
  plot_tau <- ggplot(NULL,aes(x = method, y = tau)) + theme_bw() +
    geom_point(data = subset(plot_data, tau < 50*median_tau))+
    geom_point(data = data.frame(method = weight_method, tau = weight_tau), color = "red")+ ylab("Relaxation time (s)") + coord_flip() + 
    ggtitle(paste0("relaxation time for lp per chain \n model ",saved_name))+
    theme(text = element_text(size = 9),plot.title = element_text(hjust = 0.5) )
  
  plot_log_tau <- ggplot(NULL,aes(x = method, y = tau)) + theme_bw() +
    geom_point(data = subset(plot_data, tau < 50*median_tau))+
    geom_point(data = data.frame(method = weight_method, tau = weight_tau),color = "red")+ scale_y_continuous(trans = 'log10') +
    ylab("Relaxation time (s)") + xlab("") + coord_flip() +
    ggtitle(paste0("log relaxation time for lp per chain \n model ",saved_name)) +
    theme(text = element_text(size = 9),plot.title = element_text(hjust = 0.5))
  
  print(grid.arrange(plot_tau,plot_log_tau,ncol = 2))
  ggsave(plot =grid.arrange(plot_tau,plot_log_tau,ncol = 2) , filename = paste0("figures/",saved_name,"lp_summary_log_tau.png") )
}

w_eff <- function(fit_object,ess){
  #' @title ESS of weighted object
  #' @description We took the formula from Andrew Gelman stacking paper to calculate the ESS of the weighted chains
  #' @param fit_object object(utils/simulation.R) - object with all the basic informations of the object we want
  #' @param ess vector - containing the ess of the chains studied but non weighted
  inv_ess = rep(0,length(ess))
  for (i in 1:length(ess)){
    if (!is.na(ess[i])) inv_ess[i] = 1/ess[i]
  }
  ESS_measure <- 1/sum( fit_object$weights**2 * inv_ess)
  ESS_measure/fit_object$chain_time[tail(order(fit_object$chain_time[which(fit_object$weights> 0.01)]),1)]
}


plot_integrator_lp_summary <- function(saved_name,integrators, parameter = "lp__"){
  #' @title Plot the efficiency of the entered parameter for each chain for both integrators and types of initialisation
  #' @description This function is drawing the efficiency (Effective Sample Size / total running time) for each chain.
  #' To be run, this function requires that the simulations for all type of initialisation have previously been done. The different types of initialisation are : 
  #' - MC sample from the Prior Distribution + Default HMC parameters (FF)
  #' - MC sample from a Posterior Distribution + Default HMC parameters (TF)
  #' - MC sample from the Prior Distribution + Tuned HMC parameters (FT)
  #' - MC sample from a Posterior Distribution + Tuned HMC parameters (TT)
  #' 
  #' For the last 3 initialisation configuration, we need to first run the first simulation to realise them. 
  #' 
  #' @param saved_name string - name of the model ran
  #' @param integrators vector - integrators to be studied, usually BDF and RK45
  #' @param parameter string - Parameter to be studied
  #' @seealso plot_lp_tau
  #' @return Plots of the relaxation time and log relaxation time (great if outliers make the graph difficult to read)
  #'  
  recorded_eff <- NULL
  method_names <- c("FF","TF","FT","TT")
  method <- c()
  weight_eff <- c()
  weight_method <- c()
  for (integrator in integrators){
    fitsFF <- fit_object(fits(saved_name,integrator,FALSE,FALSE),integrator)
    fitsVF <- fit_object(fits(saved_name,integrator,TRUE,FALSE),integrator)
    fitsFV <- fit_object(fits(saved_name,integrator,FALSE,TRUE),integrator)
    fitsVV <- fit_object(fits(saved_name,integrator,TRUE,TRUE),integrator)
    
    essFF <- ess_summary(fitsFF,parameter)
    essFV <- ess_summary(fitsFV,parameter)
    essVF <- ess_summary(fitsVF,parameter)
    essVV <- ess_summary(fitsVV,parameter)
    recorded_eff <- c(recorded_eff, essFF$chain_eff[1,],essVF$chain_eff[1,],essFV$chain_eff[1,],essVV$chain_eff[1,])

    w_effFF <- w_eff(fitsFF,essFF$chain_ess[1,])
    w_effFV <-  w_eff(fitsFF,essFV$chain_ess[1,])
    w_effVF <-  w_eff(fitsFF,essVF$chain_ess[1,])
    w_effVV <-  w_eff(fitsFF,essVV$chain_ess[1,])
    
    weight_eff <- c(weight_eff, w_effFF,w_effVF,w_effFV,w_effVV)
    
    method <- c(method,rep(paste0(integrator,"-", method_names),each =fitsFF$num_chains))
    weight_method <- c(weight_method, paste0(integrator,"-", method_names))
  }
  
  plot <- ggplot() + theme_bw() + geom_point(aes(x = method, y = eff, colour = "Parallel chains"),data = data.frame(method = method, eff = recorded_eff)) +
    geom_point(aes(x = method, y = eff, colour = "Weighted chain"),data = data.frame(method = weight_method, eff= weight_eff ))+
    xlab("Integrator - Initialisation configuration") + ylab("ESS/s of log posterior density")+
    scale_colour_manual(values = c("black","red"))
  print(plot)
  
  ggsave(plot, filename = paste0("figures/",saved_name,"lp_summary_ess.png"))
}



plot_posterior <- function(saved_name,integrators, run_posterior, run_parameters,fit_objects, parms, iter_sampling){
  #' @title Plot Posterior Distribution of parms for both integrator in 
  #' @description In this plot, we calculate the Efficiency (ESS/computational time) for each chain and each parameter for both integrators
  #' We highlight in red the n_drop slowest chains.
  #' @param saved_name string - name of the model ran
  #' @param integrators String or vector of strings - integrators to be studied
  #' @param run_posterior Boolean - if FALSE : the Markov Chain initialisations are drawn from the prior; 
  #' if TRUE : the Markov Chain initialisations are drawn from a posterior distribution of the model ran (necessity to have a model already ran)
  #' @param run_parameters Boolean - if FALSE : HMC parameters are randomly initialized
  #' if TRUE : HMC parameters are initialized with the tuned parameters of a previously ran model
  #' @param fit_objects vector of fit_object - fit_object are an object created with a function in utils/simulation.r
  #' @param parms vector of strings - parameters we want the posterior distribution to be plotted
  #' @param iter_sampling Integer - Number of samples in a fit
  #' @return Return a graph for each parameter drawing the posterior distribution for both integrator in a particular case of initialisation
  
  samples_all <- c()
  samples_number <- c()
  parm_name <- c()
  for (fit_object  in fit_objects){
    samples <- as_draws_df(fit_object$draw[,,parms])[1:length(parms)]
    samples_list = c()
    for (par in names(samples)){
      samples_list = c(samples_list, samples[[par]])
    }
    samples_all <- c(samples_all,samples_list)
    samples_number <- c(samples_number,length(samples_list))
    parm_name <- c(parm_name,rep(parms, each = iter_sampling * dim(fit_object$draw)[2]))
  }
  
  
  method <- rep(integrators, samples_number)

  plot_data <- data.frame(samples = samples_all,
                          parm = parm_name,
                          method = method)
  plot <- ggplot(data = plot_data,
                 aes(x = samples, color = method, fill = method)) +
    geom_histogram(alpha = 0.25, position = "identity", bins = 30) +
    theme_bw() + facet_wrap(~ parm, scale = "free")
    
  print(plot)
  ggsave(plot = plot, filename =paste0("figures/",saved_name,"_post",run_posterior,"_tuned",run_parameters,"plot_posterior.png" ))
  
}


plot_ess <- function (saved_name, run_posterior, run_parameters, fit_objects, parms){
  #' @title Plot the ESS/s for each chain, each parameter and each integrator for a specific initialisation
  #' @description In this plot, we calculate the Efficiency (ESS/computational time) for each chain and each parameter for both integrators
  #' We highlight in red the n_drop slowest chains.
  #' @param saved_name string - name of the model ran
  #' @param fit_objects vector of fit_object - fit_object are an object created with a function in utils/simulation.r
  #' @param run_posterior Boolean - if FALSE : the Markov Chain initialisations are drawn from the prior; 
  #' if TRUE : the Markov Chain initialisations are drawn from a posterior distribution of the model ran (necessity to have a model already ran)
  #' @param run_parameters Boolean - if FALSE : HMC parameters are randomly initialized
  #' if TRUE : HMC parameters are initialized with the tuned parameters of a previously ran model
  #' @param n_drop Integer - Number of slowest chains to be highlighted in red
  #' @param integrators String or vector of strings - integrators to be studied
  #' @return Return a graph for each fit_object describing the efficiency for each chain sampled for a particular initialisation
  eff <- NULL
  for (fit_object in fit_objects){
    ess_bulk <-  summarize_draws(fit_object$draw)$ess_bulk[1:length(parms)]
    time_total <- fit_object$total_time
    eff <- c(eff, ess_bulk / time_total)
  }
  parm <- rep(parms, 2)
  method <- rep(c("rk45", "bdf"), each = length(parms))
  plot_data <- data.frame(eff = eff, parm = parm, method = method)
  
  plot <- ggplot(data = plot_data,
                 aes(x = parm, y = eff, fill = method)) +
    geom_bar(stat = "identity", width = 0.3, alpha = 0.8,
             position = "dodge") + 
    theme_bw() + theme(text = element_text(size = 10)) + coord_flip() +
    ylab("ESS / s") + xlab("model parameters")
  print(plot)
  ggsave(plot = plot, filename = paste0("figures/",saved_name,"post_",run_posterior,"tuned_",run_parameters,"global_efficiency.png"))
} 

plot_ess_chains_para <- function(saved_name, run_posterior, run_parameters, fit_objects,integrators,parms,n_drop){
  #' @title Plot the ESS/s for each chain, each parameter and each integrator for a specific initialisation
  #' @description In this plot, we calculate the Efficiency (ESS/computational time) for each chain and each parameter for both integrators
  #' We highlight in red the n_drop slowest chains.
  #' @param saved_name string - name of the model ran
  #' @param fit_objects vector of fit_object - fit_object are an object created with a function in utils/simulation.r
  #' @param run_posterior Boolean - if FALSE : the Markov Chain initialisations are drawn from the prior; 
  #' if TRUE : the Markov Chain initialisations are drawn from a posterior distribution of the model ran (necessity to have a model already ran)
  #' @param run_parameters Boolean - if FALSE : HMC parameters are randomly initialized
  #' if TRUE : HMC parameters are initialized with the tuned parameters of a previously ran model
  #' @param n_drop Integer - Number of slowest chains to be highlighted in red
  #' @param integrators String or vector of strings - integrators to be studied
  #' @return Return a graph for each fit_object describing the efficiency for each chain sampled for a particular initialisation
  
  method_fast <- c()
  method_slow <- c()
  eff_fast <- c()
  eff_slow <- c()
  chain_number <- c()
  for (fit_object in fit_objects){
    time <- fit_object$chain_time
    
    fastest_chains <- order(time)
    eff_chains <- ess_summary(fit_object,parms)$chain_eff
    num_chains <- dim(fit_object$draw)[2]
    
    chain_number <- c(chain_number,num_chains - n_drop)
    
    method_fast <- c(method_fast, rep(parms, num_chains-n_drop))
    method_slow <- c(method_slow, each = rep(parms,n_drop))
    
    eff_fast <- c(eff_fast,c(eff_chains[,head(fastest_chains,n=num_chains-n_drop)]))
    eff_slow <- c(eff_slow, c(eff_chains[,tail(fastest_chains,n=n_drop)]))
  }
  
  
  method_fast <- paste0(rep(integrators, length(parms)*(chain_number))," ", method_fast)
  method_slow <- paste0(rep(integrators, each = length(parms)*n_drop)," ", method_slow)
  
  fast_data <- data.frame(method =method_fast , eff =eff_fast)
  slow_data <- data.frame(method = method_slow, eff=  eff_slow)
  
  plot <- ggplot()+
    geom_point(aes(x=method, y= eff, colour = "Fastest 20 chains"),data = fast_data)+
    geom_point(aes(x=method, y= eff, colour = "Slowest 10 chains"),data = slow_data)+
    theme_bw()+scale_y_continuous(trans = 'log10') +
    coord_flip() + ylab("Integrator - parameter") + xlab("Efficiencey (ESS/s)")+
    scale_colour_manual(values=c("black","red"))
  print(plot)
  ggsave(plot,filename = paste0("figures/",saved_name,"post_",run_posterior,"tuned_",run_parameters,"drop_chains.png"))
}


plot_drop_chain_distribution <- function(saved_name,fit_objects,run_posterior,run_parameters,parameter = "lp__",min_chains){
  #' @title Plot the ESS/s for the n fastest chains
  #' @description To plot this graph, we sort the chains by their computational time, and note them (1,...,30) and draw the Efficiency of the 
  #' n fastest chains, for each object fit_object in fit_objects. 
  #' The Efficiency is defined by the Effective Sample Size / running time of the whole method
  #' 
  #' @param saved_name string - name of the model ran
  #' @param fit_objects vector of fit_object - fit_object are an object created with a function in utils/simulation.r
  #' @param run_posterior Boolean - if FALSE : the Markov Chain initializations are drawn from the prior; 
  #' if TRUE : the Markov Chain initializations are drawn from a posterior distribution of the model ran (necessity to have a model already ran)
  #' @param run_parameters Boolean - if FALSE : HMC parameters are randomly initialized
  #' if TRUE : HMC parameters are initialized with the tuned parameters of a previously ran model
  #' @param min_chains Integer - Number minimum of chains to be taken in account
  #' @param parameter String or vector of strings - parameters to be studied
  #' @return Return a graph for each fit_object describing the evolution of the efficiency
  
   inte_plots = list()
  
  for (fit_object in fit_objects){
    time <- fit_object$chain_time
    fastest_chains <- order(time)
    draws <- fit_object$draw[,, parameter]
    
    nChains <- fit_object$num_chains
    eff_group <- c()
    ess_group <- c()
    for (i in (nChains):min_chains){ 
      study_chains <- head(fastest_chains, n = i)
      new_ess <- dplyr::pull(summarize_draws(as_draws_df(draws[,study_chains,]),measure = c("ess_bulk")), 2)
      
      ess_group <- c(ess_group, new_ess)
      eff_group <- c(eff_group,new_ess/ time[fastest_chains[i]])
    }
    

    inte_plots[[fit_object$integrator]] =  ggplot(data = data.frame(eff = c(eff_group), chains = c(nChains:min_chains)),
                  aes(x = chains, y = eff)) + geom_point() + coord_flip() + theme_bw()+
      xlab("n fastest chains") + ylab("Efficiency (ESS associated/ run time of the last n chains)")
    print(inte_plots[[fit_object$integrator]])
    ggsave(inte_plots[[fit_object$integrator]], file =  paste0("figures/",saved_name, fit_object$integrator,"_post",run_posterior,"_tuned",run_parameters,"lp_eff.png"))

  }
  
}

plot_ess_chain_distribution <- function(saved_name, fit_objects,run_posterior,run_parameters,parameter = "lp__",min_chains){
  #' @title Plot the ESS for each chain sorted in their order of arrival
  #' @description To plot this graph, we sort the chains by their computational time, and note them (1,...,30) and draw the ESS of the n^{th} fastest chain, for each object fit_object in fit_objects 
  #' 
  #' @param saved_name string - name of the model ran
  #' @param fit_objects vector of fit_object - fit_object are an object created with a function in utils/simulation.r
  #' @param run_posterior Boolean - if FALSE : the Markov Chain initializations are drawn from the prior; 
  #' if TRUE : the Markov Chain initializations are drawn from a posterior distribution of the model ran (necessity to have a model already ran)
  #' @param run_parameters Boolean - if FALSE : HMC parameters are randomly initialized
  #' if TRUE : HMC parameters are initialized with the tuned parameters of a previously ran model
  #' @param min_chains Integer - Number minimum of chains to be taken in account
  #' @param parameter String or vector of strings - parameters to be studied
  #' @return Return a graph for each fit_object describing the ESS brought by each chain
  
  inte_plots = list()
  
  for (fit_object in fit_objects){
    time <- fit_object$chain_time
    fastest_chains <- order(time)
    draws <- fit_object$draw[,,parameter]
    
    nChains <- fit_object$num_chains
    ess_group <- c()
    for (i in (nChains):1){
      new_ess <- dplyr::pull(summarize_draws(as_draws_df(draws[,fastest_chains[i],]),measure = c("ess_bulk")), 2)
      ess_group <- c(ess_group, new_ess)
    }
    

    inte_plots[[fit_object$integrator]] =  ggplot(data = data.frame(ess = c(ess_group), chains = c(nChains:1)),
                                       aes(x = chains, y = ess)) + geom_point() + coord_flip()+ theme_bw() +
      xlab("Chains (from fastest to slowest)") + ylab("ESS")
    print(inte_plots[[fit_object$integrator]])
    ggsave(inte_plots[[fit_object$integrator]], file =  paste0("figures/",saved_name, fit_object$integrator,"_post",run_posterior,"_tuned",run_parameters,"lp_ess.png"))
    
  }
  
}



plot_Rhat <- function(saved_name, fit_objects,run_posterior,run_parameters,min_chains){
  #' @title Plot of Rhat of the fastest chains to finish
  #' @description To plot this graph, we sort the chains by their computational time, and note them (1,...,30) and draw the Rhat(1,...,n), the Rhat of the 
  #' n fastest chains, for each object fit_object in fit_objects 
  #' 
  #' @param saved_name string - name of the model ran
  #' @param fit_objects vector of fit_object - fit_object are an object created with a function in utils/simulation.r
  #' @param run_posterior Boolean - if FALSE : the Markov Chain initializations are drawn from the prior; 
  #' if TRUE : the Markov Chain initializations are drawn from a posterior distribution of the model ran (necessity to have a model already ran)
  #' @param run_parameters Boolean - if FALSE : HMC parameters are randomly initialized
  #' if TRUE : HMC parameters are initialized with the tuned parameters of a previously ran model
  #' @param min_chains Integer - Number minimum of chains to be taken in account
  #' @return Return a graph for each fit_object describing the evolution of Rhat
  #'  
  inte_plots = list()
  for (fit_object in fit_objects){
    rhat_group <- c()
    time <- fit_object$chain_time
    fastest_chains <- order(time)
    draws <- fit_object$draw[,,"lp__"]
    nChains <- fit_object$num_chains
    for (i in (nChains):min_chains){
      study_chains <- head(fastest_chains, n = i)
      new_rhat <- dplyr::pull(summarize_draws(as_draws_df(draws[,study_chains,]),measure = c("rhat")), 2)
      rhat_group <- c(rhat_group, new_rhat)
    }
    #rhat_weight =dplyr::pull(summarize_draws(fit_object$weight_draws["lp__"],measure = c("rhat")), 2)
    inte_plots[[fit_object$integrator]] = ggplot(data = data.frame(rhat = c(rhat_group), chains = c(nChains:min_chains)),
                                      aes(x = chains, y = rhat)) + geom_point()+
      geom_hline(yintercept = 1.1, linetype = "dashed", color = "red") +
      geom_hline(yintercept = 1, linetype = "blank")+
      coord_flip()+ theme_bw() +
      xlab("n fastest chains") + ylab("Rhat associated") +  theme(text = element_text(size = 10))

    print(inte_plots[[fit_object$integrator]])
    ggsave(inte_plots[[fit_object$integrator]],file = paste0("figures/",saved_name, fit_object$integrator,"_post",run_posterior,"_tuned",run_parameters,"lp_rhats.png"))
  }
}

plot_ess_ini <- function(fit_objects, integrator, parms){
  #' @title Efficiency for each parameter and each initialisation
  #' @description We extract the efficiency for each parameter for the different simulations with different initialisation
  #' @param fit_objects list of fit_object - contains all the informations of the simulations with their own initialisation
  #' in the following order ("FF","TF","FT","TT") @seealso plot_integrator_lp_summary for details on the meaning of ("FF","TF","FT","TT")
  #' @param integrator string - integrator studied
  #' @param parms vector of strings - parameters to be studied
  #' @return the evolution of the efficiency of parameters throughout different initialisation tested
  #' 
  eff <- NULL
  method_names <- c("FF","TF","FT","TT")
  for (fitted in fit_objects){
    ess_bulk <-  summarize_draws(fitted$draw)$ess_bulk[1:length(parms)]
    time_total <- fitted$total_time
    eff <- c(eff, ess_bulk / time_total)
  }
  parm <- rep(parms, 4)
  method <- rep(paste0(integrator," ",method_names), each = length(parms))
  plot_data <- data.frame(eff = eff, parm = parm, method = method)
  
  plot <- ggplot(data = plot_data,
                 aes(x = parm, y = eff, fill = method)) +
    geom_bar(stat = "identity", width = 0.5, alpha = 0.8,
             position = "dodge") + 
    theme_bw() + theme(text = element_text(size = 10)) + coord_flip() +
    ylab("ESS / s") + xlab("model parameters")
  print(plot)
}

