# R script to infer the parameters of the exponential growth model for the 
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

###### Folder path of this file to load data and to save results #####
# folder.path <- dirname(rstudioapi::getActiveDocumentContext()$path)

###### Load functions #####
# functions.path <- ".../utils.R"
# source(functions.path)

##### Generate exponential growth model for dMod #####
# ODE model
ode_equation <- eqnvec(x1 = "beta1*x1", x2 = "beta2*x2")
model_dMod <- odemodel(ode_equation, modelname = "exponential_model")
x_dMod <- Xs(model_dMod)

# Observable model
observable <- eqnvec(y = "log(x1/(x1+x2))-(sigma^2)/2")
g_dMod <- Y(observable, x_dMod, compile = T, attach.input = F)

# Error model
error_pars <- eqnvec(y = "sigma")
error_dMod <- Y(
  error_pars, f = c(observable, ode_equation), attach.input = TRUE, compile = TRUE)

# Parameter transformations for original full model 
p_full <- P(trafo=eqnvec(
  x1 = "exp(x1)", x2 = "exp(x2)", beta1 = "beta1", beta2 = "beta2", 
  sigma = "exp(sigma)"), condition = "C1")
outerpars_full <- getParameters(p_full)

# Reparametrisations to make model structurally identifiable
equations <- getEquations(p_full)
equations[["C1"]]["beta2"] <- "theta1+beta1"
equations[["C1"]]["x2"] <- "exp(theta2)*exp(x1)"
p_repar <- P(equations[["C1"]], condition = "C1")
outerpars <- getParameters(p_repar)
theta <- c("theta1", "theta2", "sigma")


##### Simulated data ##### 
# Specify data sets by file names
filenames_vector <- c(
  "simulated_data_size=8_seed=8.rds",
  "simulated_data_size=16_seed=16.rds",
  "simulated_data_size=32_seed=32.rds",
  "simulated_data_size=64_seed=64.rds")


################################################################################
########### Function to infer parameters of exponential growth model ###########
################################################################################
# Function which is parallelised by foreach
estimation_exp_model <- function(data, pouter){
  
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
  obj <- normL2(data, g_dMod*x_dMod*p_repar, errmodel = error_dMod) + constraintL2(pouter, sigma = 5)
  # Multiple parameter estimation with randomly chosen starting values
  fits <- mstrust(obj, pouter, studyname = "multiple_fits", rinit = 0.1, rmax = 10, 
                  fits = 30, samplefun = "runif", iterlim = 500)
  
  # Extract estimated parameters of best estimation run
  bestfit <- as.parvec(as.parframe(fits))
  
  # Objective value for best estimation run without L2-contribution
  obj_value <- normL2(data, g_dMod*x_dMod*p_repar, errmodel = error_dMod)(bestfit)$value
  
  # Profile likelihoods
  profiles <- profile(obj = obj, pars = bestfit, whichPar = theta, alpha = 0.0001,
                      method = "integrate")
  
  # Profile likelihood-based confidence intervals (for 0.95 quantile of 
  # chi-squared distribution and for adapted threshold based on Cantelli's 
  # inequality)
  ci_0.99255_confint <- confint(profiles, val.column = "data", level = 0.99255) 
  ci_0.95_confint <- confint(profiles, val.column = "data", level = 0.95)
  
  # If a confidence interval crosses confidence threshold more than twice,
  # adapt confidence interval; otherwise, just transform CIs to list
  ci_0.99255_list <- ci_extraction(ci_0.99255_confint, profiles, obj_value, 0.99255)
  ci_0.95_list <- ci_extraction(ci_0.95_confint, profiles, obj_value, 0.95)
  
  # Create output list with parameter estimates and confidence intervals
  output_list <- list(
    bestfit = bestfit, ci_0.95_list = ci_0.95_list, ci_0.99255_list = ci_0.99255_list)
  
  return(output_list)
}


################################################################################
######################### Estimations ##########################################
################################################################################
set.seed(110)

# list for results
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
    
    # Center of starting values for parameter estimation
    pouter <- structure(
      c(log(experiment_pars[s, "x1"]),
        log(experiment_pars[s, "x2"]/experiment_pars[s, "x1"]),
        experiment_pars[s, "b1"] - experiment_pars[s, "d1"],
        experiment_pars[s, "b2"] - experiment_pars[s, "d2"] - (experiment_pars[s, "b1"]-experiment_pars[s, "d1"]),
        log(experiment_pars[s, "sigma"])), names = outerpars)
    
    # Detect the number of available cores and activate cluster
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)
    
    # Parallelise esimations for the simulated data sets
    result <- foreach(
      i = 1:num_sim,
      .packages = c("dMod")) %dopar% {
        estimation_exp_model(
          data = data_list_exp[[i]], 
          pouter = pouter)
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
# filename <- "results_exp_model_L2=5_seed=110"
# file.path <- paste0(folder.path, "/RDS/", filename, ".rds")
# saveRDS(results_list, file = file.path)