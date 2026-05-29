# R script to infer the parameters of the logistic growth model for the 
# simulated data sets.

################################################################################
################################ Preparations ##################################
################################################################################

##### Load R packages #####
library("dMod")
library("dplyr")
library("foreach")
library("doParallel")
library("parallel")

##### Folder path of this file to load data and to save results #####
# folder.path <- dirname(rstudioapi::getActiveDocumentContext()$path)

###### Load functions #####
# functions.path <- ".../utils.R"
# source(functions.path)

##### Generate logistic growth model for dMod #####
# ODE model
ode_equation <- eqnvec(
  x1 = "lambda1*x1*(1-(x1+x2)/K)", x2 = "lambda2*x2*(1-(x1+x2)/K)")
model <- odemodel(ode_equation, modelname = "logistic_model") 
x_dMod <- Xs(model)

# Observable model
observable <- eqnvec(y = "x1/(x1+x2)")
g_dMod <- Y(observable, x_dMod, compile = T, attach.input = F)

# Error model
error_pars <- eqnvec(y = "sigma")
error_dMod <- Y(error_pars, f = c(observable, ode_equation), attach.input = T, compile = T)

# Global parameters
p_global <- P(
  eqnvec(
    lambda1 = "lambda1", lambda2 = "exp(lambda2)", K = "exp(K)", x1 = "exp(x1)", 
    x2 = "exp(x2)", sigma = "exp(sigma)"), condition = "global")


##### Simulated data ##### 
# Specify data sets by file names
filenames_vector <- c(
  "simulated_data_size=8.rds",
  "simulated_data_size=16.rds",
  "simulated_data_size=32.rds")


################################################################################
##################### Parameter estimation #####################################
################################################################################

# Function which is parallelised by foreach
estimation_log_model <- function(data, pars_sim){
  
  # Prepare mouse-specific data sets and parameter transformations for dMod 
  unique_mice <- unique(data$mouse)
  data_dMod <- NULL
  p_dMod <- NULL
  
  for (j in 1:length(unique_mice)){
    data_m <- data[data$mouse == unique_mice[j], ]
    data_m <- data_m[, c("name", "time", "value")]
    data_m <- data_m[order(data_m$time), ]
    data_m$sigma <- NA
    data_m <- datalist(data_m)
    names(data_m) <- unique_mice[j]
    data_dMod <- data_dMod + data_m
    
    trafo <- getEquations(p_global, conditions = "global")
    trafo["x1"] <- paste0("exp(x1_", unique_mice[j], ")")
    trafo["x2"] <- paste0("exp(x2_", unique_mice[j], ")")
    p_dMod <- p_dMod + P(trafo, condition = unique_mice[j])
  }
  
  # Center of starting values for parameter estimation
  outerpars <- getParameters(p_dMod)
  outerpars_x1 <- outerpars[grepl("x1_", outerpars)]
  outerpars_x2 <- outerpars[grepl("x2_", outerpars)]
  fixed_par <- c("K", outerpars_x1)
  pouter <- structure(
    c(as.numeric(pars_sim[1, "b1"] - pars_sim[1, "d1"]), 
      log(as.numeric(pars_sim[1, "b2"] - pars_sim[1, "d2"])), 
      log(as.numeric(pars_sim[1, "K"])), 
      rep(log(pars_sim[1, "x1"]), length(outerpars_x1)), 
      rep(log(pars_sim[1, "x2"]), length(outerpars_x2)), 
      log(pars_sim[1, "sigma"])), 
    names = c("lambda1", "lambda2", "K", outerpars_x1, outerpars_x2, "sigma"))
  
  # Objective function for parameter estimation with L2-constraint
  obj <- normL2(data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(
    pouter, sigma = 3)
  
  # Multiple parameter estimation with randomly chosen starting values
  fits <- mstrust(obj, pouter[!names(pouter) %in% fixed_par], fixed = pouter[fixed_par], 
                  studyname = "multiple_fits", rinit = 0.1, rmax = 10, fits = 30, 
                  samplefun = "runif", iterlim = 500)
  
  # Extract estimated parameters of best estimation run
  bestfit <- as.parvec(as.parframe(fits))
  
  # Objective value for best estimation run without L2-contribution
  obj_value <- normL2(data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod)(
    c(bestfit, pouter[fixed_par]))$value
  
  # Profile likelihoods
  profiles <- profile(
    obj = obj, pars = bestfit, fixed = pouter[fixed_par], 
    whichPar = names(pouter)[!names(pouter) %in% c(fixed_par)], 
    alpha = 0.001, method = "integrate")
  
  # Profile likelihood-based confidence intervals (for 0.95 quantile of 
  # chi-squared distribution and for adapted threshold based on Cantelli's 
  # inequality)
  ci_0.99255_confint <- confint(profiles, val.column = "data", level = 0.99255)
  ci_0.95_confint <- confint(profiles, val.column = "data", level = 0.95)
  
  # If a confidence interval crosses confidence threshold more than twice,
  # adapt confidence interval; otherwise, just transform CIs to list
  ci_0.99255_list <- ci_extraction(ci_0.99255_confint, profiles, obj_value, 0.99255)
  ci_0.95_list <- ci_extraction(ci_0.95_confint, profiles, obj_value, 0.95)
  
  # Prepare output
  result_list <- list(
    mllik = obj_value, bestfit = bestfit, profiles = profiles, 
    ci_0.95_list = ci_0.95_list, ci_0.99255_list = ci_0.99255_list)
  
  return(result_list)
}


##### Run parameter estimations in parallel #####
set.seed(210)

# List for results
results_list <- list()

# Loop through sample sizes
for (n in 1:length(filenames_vector)) {
  # Load data sets for sample size
  read_list <- readRDS(paste0(folder.path, "/RDS/", filenames_vector[n]))
  data_list <- read_list$data_list
  experiment_pars <- read_list$experiment_pars
  result_list <- list()
  
  # Loop through modification scenarios
  for (s in 1:length(data_list)){ 
    # Extract data sets for scenario
    data_list_exp <- data_list[[s]]
    pars_sim <- experiment_pars[s, 2:ncol(experiment_pars)]
    num_sim <- length(data_list_exp)

    # Detect the number of available cores and activate cluster
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)
    
    # Parallelise esimation for the simulated data sets
    result <- foreach(
      i = 1:num_sim, 
      .packages=c("dMod")) %dopar% {
        estimation_log_model(
          data = data_list_exp[[i]], 
          pars_sim = pars_sim)
      }
    # Stop cluster
    stopCluster(cl)
    
    result_list <- c(result_list, list(result))
    names(result_list)[length(result_list)] <- paste0("experiment_", s)
  }
  
  results_list <- c(results_list, list(result_list))
  names(results_list)[length(results_list)] <- paste0("sample_size=",nrow(data_list_exp[[1]]))
}

# Save results
# filename <- "results_log_model"
# file.path <- paste0(folder.path, "/RDS/", filename, ".rds")
# saveRDS(results_list, file = file.path)
# Due to storage reasons, profile likelihoods are removed from the rds file


################################################################################
#################### Bootstrap for lambda2 - lambda1 ###########################
################################################################################

##### Function which is parallelised by foreach #####
bootstrap_log_model <- function(data, pars_sim){
  
  # Set design for bootstrap samples
  groups <- as.character(data[, "time"])
  sample_sizes <- sapply(unique(groups), function(x) 
    return(sum(data[, "time"] == as.numeric(x))))
  
  # Sample indices for bootstrap sample
  data_b_index <- unlist(lapply(unique(groups), function(gr) {
    group_indices <- which(groups == gr)
    n_samples <- sample_sizes[gr]
    sample(group_indices, n_samples, replace = TRUE)
  }))
  data_b <- data[data_b_index,]
  
  # Label measurements of bootstrap sample with fictitious mice; 
  # numbering does not necessarily coincide with numbering of original data set
  num_mice <- sum(data_b$time == 0)
  mouse_numbers <- rep(1:num_mice, 2)
  mouse_numbers <- paste0("m", mouse_numbers)
  data_b <- cbind(mouse_numbers, data_b)
  colnames(data_b)[1] <- "mouse"
  rownames(data_b) <- NULL
  
  # Prepare mouse-specific data sets and parameter transformations for dMod 
  unique_mice <- unique(data_b$mouse)
  data_dMod <- NULL
  p_dMod <- NULL
  
  for (j in 1:length(unique_mice)){
    data_m <- data_b[data_b$mouse == unique_mice[j], ]
    data_m <- data_m[, c("name", "time", "value", "sigma")]
    data_m <- data_m[order(data_m$time), ]
    data_m <- datalist(data_m)
    names(data_m) <- unique_mice[j]
    data_dMod <- data_dMod + data_m
    
    # add initial parameters for each mouse individually to the parameter transformations
    trafo <- getEquations(p_global, conditions = "global")
    trafo["x1"] <- paste0("exp(x1_", unique_mice[j], ")")
    trafo["x2"] <- paste0("exp(x2_", unique_mice[j], ")")
    p_dMod <- p_dMod + P(trafo, condition = unique_mice[j])
  }
  
  # Center of starting values for parameter estimation
  outerpars <- getParameters(p_dMod)
  outerpars_x1 <- outerpars[grepl("x1_", outerpars)]
  outerpars_x2 <- outerpars[grepl("x2_", outerpars)]
  fixed_par <- c("K", outerpars_x1)
  pouter <- structure(
    c(as.numeric(pars_sim[1, "b1"] - pars_sim[1, "d1"]), 
      log(as.numeric(pars_sim[1, "b2"] - pars_sim[1, "d2"])), 
      log(as.numeric(pars_sim[1, "K"])), 
      rep(log(pars_sim[1, "x1"]), length(outerpars_x1)), 
      rep(log(pars_sim[1, "x2"]), length(outerpars_x2)), 
      log(pars_sim[1, "sigma"])), 
    names = c("lambda1", "lambda2", "K", outerpars_x1, outerpars_x2, "sigma"))
  
  # Objective function for parameter estimation with L2-constraint
  obj <- normL2(data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(
    pouter, sigma = 3)
  
  # Multiple parameter estimation with randomly chosen starting values
  fits <- mstrust(
    obj, pouter[!names(pouter) %in% fixed_par], fixed = pouter[fixed_par], 
    studyname = "multiple_fits", rinit = 0.1, rmax = 10, fits = 30, 
    samplefun = "runif", iterlim = 500)
  
  # If at least one estimation converged, calculate difference of lambda_2 and 
  # lambda_1 for estimation with highest likelihood 
  if (!inherits(try(as.parframe(fits)), "try-error")){
    bestfit_b <- as.parvec(as.parframe(fits))
    diff_lambda <- exp(as.numeric(bestfit_b["lambda2"])) - as.numeric(bestfit_b["lambda1"])
  } else {
    diff_lambda <- NA
  }
  
  return(diff_lambda)
}


##### Run bootstraps in parallel #####
set.seed(211)

# Number of bootstrap samples 
B <- 499

# List for results
boot_results <- list()

# Loop through sample sizes
for (n in 1:length(filenames_vector)) { 
  # Load data sets for sample size
  read_list <- readRDS(paste0(folder.path, "/RDS/", filenames_vector[n]))
  data_list <- read_list$data_list
  experiment_pars <- read_list$experiment_pars
  boot_results_sample_size <- list()
  
  # Loop through modification scenarios
  for (s in 1:length(data_list)){ 
    # Extract data sets for scenario
    data_list_exp <- data_list[[s]]
    pars_sim <- experiment_pars[s, 2:ncol(experiment_pars)]
    num_sim <- length(data_list_exp)
    boot_results_exp <- list()
    
    # Loop through simulated data sets for experimental scenario 
    for (j in 1:length(data_list_exp)){ 
      # Extract data set and prepare it for dMod
      data <- data_list_exp[[j]]
      data <- data[, c("name", "time", "value")]
      data <- data[order(data$time),]
      data$sigma <- NA
      
      # Detect the number of available cores and activate cluster
      cl <- makeCluster(detectCores() - 1)
      registerDoParallel(cl) 
      
      # Loop through all data sets of a scenario via parallelisation
      boot_results_sim <- foreach(
        i = 1:B,
        .packages = c("dMod")) %dopar% {
          bootstrap_log_model(
            data = data,
            pars_sim = pars_sim)
        }
      
      # Stop cluster
      stopCluster(cl)
      
      # Store results for this simulation of the experimental scenario
      boot_results_unlist <- unlist(boot_results_sim)
      print(paste("Number of diverged bootstrap estimations:", sum(is.na(boot_results_unlist))))
      boot_results_unlist <- boot_results_unlist[!is.na(boot_results_unlist)]
      boot_results_exp <- c(boot_results_exp, list(boot_results_unlist))
    }
    
    # Store results for experimental scenario
    boot_results_sample_size <- c(boot_results_sample_size, list(boot_results_exp))
    names(boot_results_sample_size)[length(boot_results_sample_size)] <- names(data_list)[s]
  }
  
  # Store results for sample size
  boot_results <- c(boot_results, list(boot_results_sample_size))
  names(boot_results)[length(boot_results)] <- paste0(
    "sample_size=", nrow(data_list[[1]][[1]]))
  
}

# Save results
# filename <- "results_log_model_bootstrap"
# file.path <- paste0(folder.path, "/RDS/", filename, ".rds")
# saveRDS(boot_results, file = file.path)
