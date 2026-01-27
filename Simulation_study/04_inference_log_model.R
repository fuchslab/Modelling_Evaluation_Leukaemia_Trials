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
ode_equation <- eqnvec(x1 = "lambda1*x1*(1-(x1+x2)/K)",
                       x2 = "lambda2*x2*(1-(x1+x2)/K)")
model_dMod <- odemodel(ode_equation, modelname = "logistic_model") 
x_dMod <- Xs(model_dMod) 

# Observable model
observable <- eqnvec(y = "log(x1/(x1+x2))-(sigma^2)/2")
g_dMod <- Y(observable, x_dMod, compile = T, attach.input = F)

# Error model
error_pars <- eqnvec(y = "sigma")
error_dMod <- Y(error_pars, f = c(observable, ode_equation), attach.input = TRUE, 
                compile = TRUE)

# Parameter transformations
p_dMod <- P(trafo = eqnvec(
  lambda1 = "lambda1", lambda2 = "lambda2", K = "exp(K)", x1 = "exp(x1)", 
  x2 = "exp(x2)", sigma = "exp(sigma)"), condition = "C1")
outerpars <- getParameters(p_dMod)


##### Simulated data ##### 
# Specify data sets by file names
filenames_vector <- c(
  "simulated_data_size=8_seed=8.rds",
  "simulated_data_size=16_seed=16.rds",
  "simulated_data_size=32_seed=32.rds",
  "simulated_data_size=64_seed=64.rds")


################################################################################
##################### Parameter estimation #####################################
################################################################################

# Function which is parallelised by foreach
estimation_log_model <- function(data, pouter, fixed_par){
  
  # Extract time and value column of data set 
  data <- data[, c("time", "value")]
  
  # Prepare data for dMod and log-transform observations
  data <- data[order(data$time),]
  data$value <- log(data$value)
  data$name <- "y"
  data$sigma <- NA
  data <- data[, c(3, 1, 2, 4)]
  data <- datalist(C1 = data)
  
  # Objective function for parameter estimation
  obj <- normL2(
    data, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(
      pouter, sigma = c(3, 3, 5, 5, 5, 3))
  
  # Multiple parameter estimation with randomly chosen starting values
  fits <- mstrust(obj, pouter[!names(pouter) %in% fixed_par], fixed = pouter[fixed_par], 
                  studyname = "multiple_fits", rinit = 0.1, rmax = 10, fits = 30, 
                  samplefun = "runif", iterlim = 500)
  
  # Extract estimated parameters of best estimation run
  bestfit <- as.parvec(as.parframe(fits))
  
  # Objective value for best estimation run without L2-contribution
  obj_value <- normL2(data, g_dMod*x_dMod*p_dMod, errmodel = error_dMod)(
    c(bestfit, pouter[fixed_par]))$value
  
  # Profile likelihoods
  profiles <- profile(obj = obj, pars = bestfit, fixed = pouter[fixed_par],
                      method = "integrate", whichPar = names(bestfit), alpha = 0.001)
  
  # Profile likelihood-based confidence intervals with 0.95 quantile of 
  # chi-squared distribution
  ci_confint <- confint(profiles, val.column = "data", level = 0.95) 
  
  # If a confidence interval crosses confidence threshold more than twice,
  # adapt confidence interval; otherwise, just transform CIs to list
  ci_0.95_list <- ci_extraction(ci_confint, profiles, obj_value, 0.95)
  
  # Prepare output
  bestfit_df <- data.frame(
    mllik = obj_value, lambda1 = bestfit["lambda1"], lambda2 = bestfit["lambda2"], 
    K = exp(pouter["K"]), x1 = exp(pouter["x1"]), x2 = exp(bestfit["x2"]), 
    sigma = exp(bestfit["sigma"]))
  
  # Both bestfit_df and bestfit contain estimated (reparametrised) parameter values. 
  # bestfit_df is a data frame with one row for further processing and bestfit
  # is a parvec (a dMod class).
  output_list <- list(bestfit_df = bestfit_df, bestfit = bestfit, 
                      ci_0.95_list = ci_0.95_list)
  
  return(output_list)
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
    num_sim <- length(data_list_exp)
    
    # Center of starting values for parameter estimation and fix K and x1
    pouter <- structure(c(
      experiment_pars[s, "b1"] - experiment_pars[s, "d1"],
      experiment_pars[s, "b2"] - experiment_pars[s, "d2"], 
      log(1e9),
      log(experiment_pars[s, "x1"]), 
      log(experiment_pars[s, "x2"]),
      log(experiment_pars[s, "sigma"])), 
      names = outerpars)
    fixed_par <- c("K", "x1")
    
    # Detect the number of available cores and activate cluster
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)
    
    # Parallelise esimation for the simulated data sets
    result <- foreach(
      i = 1:num_sim, 
      .packages=c("dMod")) %dopar% {
        estimation_log_model(
          data = data_list_exp[[i]], 
          pouter = pouter,
          fixed_par = fixed_par)
      }
    # Stop cluster
    stopCluster(cl)
    
    result_list <- c(result_list, list(result))
    names(result_list)[length(result_list)] <- paste0("experiment_", s)
  }
  
  results_list <- c(results_list, list(result_list))
  names(results_list)[length(results_list)] <- paste0(
    "sample_size=", nrow(data_list[[1]][[1]]))
}

# Save results
# filename <- "results_log_model_server_fixed=(K,x1)_L2=(3,3,5,5,5,3)_seed=210"
# file.path <- paste0(folder.path, "/RDS/", filename, ".rds")
# saveRDS(results_list, file = file.path)


################################################################################
#################### Bootstrap for lambda2 - lambda1 ###########################
################################################################################

# Function which is parallelised by foreach
bootstrap_log_model <- function(
    data, groups, sample_sizes, pouter, fixed_par){
  
  # Sample indices for bootstrap sample
  sampled_data <- unlist(lapply(unique(groups), function(gr) {
    group_indices <- which(groups == gr)
    n_samples <- sample_sizes[gr]
    sample(group_indices, n_samples, replace = TRUE)
  }))
  
  # Create datalist for dMod of bootstrap sample
  data_b <- data[sampled_data,]
  rownames(data_b) <- NULL
  data_b <- datalist(C1 = data_b[order(data$time),])
  
  # Objective function for parameter estimation
  obj <- normL2(data_b, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(
    pouter, sigma = c(3, 3, 5, 5, 5, 3))
  
  # Multiple parameter estimation with randomly chosen starting values
  fits <- mstrust(obj, pouter[!names(pouter) %in% fixed_par], fixed = pouter[fixed_par], 
                  rinit = 0.1, rmax = 10, fits = 20, samplefun = "runif", iterlim = 500,
                  studyname = "multiple_fits")
  
  # Run multiple estimations until at least one run converges
  while (inherits(try(as.parframe(fits)), "try-error")) {
    fits <- mstrust(obj, pouter[!names(pouter) %in% fixed_par], fixed = pouter[fixed_par], 
                    rinit = 0.1, rmax = 10, fits = 10, samplefun = "runif", iterlim = 1000,
                    studyname = "multiple_fits")
  }
  
  # Calculate difference of lambda_2 and lambda_1 for best estimation run
  bestfit_b <- as.parvec(as.parframe(fits))
  diff_lambda <- as.numeric(bestfit_b["lambda2"] - bestfit_b["lambda1"])
  
  return(diff_lambda)
}


##### Run bootstraps in parallel #####
set.seed(211)

# Number of bootstrap samples 
B <- 999

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
    num_sim <- length(data_list_exp)
    boot_results_exp <- list()
    
    # Center of starting values for parameter estimation and fix K and x1
    pouter <- structure(c(
      experiment_pars[s, "b1"] - experiment_pars[s, "d1"],
      experiment_pars[s, "b2"] - experiment_pars[s, "d2"], 
      log(1e9),
      log(experiment_pars[s, "x1"]), 
      log(experiment_pars[s, "x2"]),
      log(experiment_pars[s, "sigma"])), 
      names = outerpars)
    fixed_par <- c("K", "x1")
    
    # Set design for bootstrap samples
    groups <- as.character(data_list_exp[[1]][, "time"])
    sample_sizes <- sapply(unique(groups), function(x) 
      return(sum(data_list_exp[[1]][, "time"] == as.numeric(x))))
    
    # Loop through simulated data sets for experimental scenario 
    for (j in 1:length(data_list_exp)){ 
      # Extract data set and prepare it for dMod
      data <- data_list_exp[[j]]
      data <- data[, c("time", "value")]
      data$value <- log(data$value)
      data$name <- "y"
      data$sigma <- NA
      data <- data[, c(3, 1, 2, 4)]
      
      # Detect the number of available cores and activate cluster
      cl <- makeCluster(detectCores() - 1)
      registerDoParallel(cl) 
      
      # Loop through all data sets of a scenario via parallelisation
      boot_results_sim <- foreach(
        i = 1:B,
        .packages = c("dMod")) %dopar% {
          bootstrap_log_model(
            data = data,
            groups = groups,
            sample_sizes = sample_sizes,
            pouter = pouter,
            fixed_par = fixed_par
          )
        }
      
      # Stop cluster
      stopCluster(cl)
      
      # Store results for this simulation of the experimental scenario
      boot_results_exp <- c(boot_results_exp, list(unlist(boot_results_sim)))
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
# filename <- "results_log_model_server_bootstrap_seed=211"
# file.path <- paste0(folder.path, "/RDS/", filename, ".rds")
# saveRDS(boot_results, file = file.path)