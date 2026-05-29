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
model <- odemodel(ode_equation, modelname = "exponential_model") 
x_dMod <- Xs(model)

# Observable model
observable <- eqnvec(y = "x1/(x1+x2)")
g_dMod <- Y(observable, x_dMod, compile = T, attach.input = F)

# Error model
error_pars <- eqnvec(y = "sigma")
error_dMod <- Y(error_pars, f = c(observable, ode_equation), attach.input = T, compile = T)

# Parameter transformations for original full model 
# Global parameters
p_full <- P(
  eqnvec(
    beta1 = "beta1", beta2 = "beta2", x1 = "exp(x1)", x2 = "exp(x2)", 
    sigma = "exp(sigma)"), condition = "full") 

# Reparametrisations to make model structurally identifiable
equations <- getEquations(p_full)
equations[["full"]]["beta2"] <- "theta1+beta1"
equations[["full"]]["x2"] <- "exp(theta2)*exp(x1)"
p_global <- P(equations[["full"]], condition = "global")


##### Simulated data ##### 
# Specify data sets by file names
filenames_vector <- c(
  "simulated_data_size=8.rds",
  "simulated_data_size=16.rds",
  "simulated_data_size=32.rds")


################################################################################
########### Function to infer parameters of exponential growth model ###########
################################################################################
# Function which is parallelised by foreach
estimation_exp_model <- function(data, pars_sim){
  
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
    trafo["x2"] <- paste0("exp(x1)/exp(theta2_", unique_mice[j], ")")
    p_dMod <- p_dMod + P(trafo, condition = unique_mice[j])
  }
  
  # Center of starting values for parameter estimation
  outerpars <- getParameters(p_dMod)
  outerpars_theta2 <- outerpars[grepl("theta2_", outerpars)]
  pouter <- structure(c(
    pars_sim[1, "b2"] - pars_sim[1, "d2"] - (pars_sim[1, "b1"] - pars_sim[1, "d1"]),
    pars_sim[1, "b1"] - pars_sim[1, "d1"], 
    log(pars_sim[1, "x1"]),
    log(pars_sim[1, "sigma"]),
    rep(log(1), length(outerpars_theta2))), 
    names = c("theta1", "beta1", "x1", "sigma", outerpars_theta2))
  
  # Objective function for parameter estimation with L2-constraint
  obj <- normL2(data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(
    pouter, sigma = 3)
  
  # Multiple parameter estimation with randomly chosen starting values
  fits <- mstrust(obj, pouter, studyname = "multiple_fits", rinit = 0.1, 
                  rmax = 10, fits = 30, samplefun = "runif", iterlim = 500)
  
  # Extract estimated parameters of best estimation run
  bestfit <- as.parvec(as.parframe(fits))
  
  # Objective value for best estimation run without L2-contribution
  obj_value <- normL2(data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod)(
    bestfit)$value
  
  # Profile likelihoods
  profiles <- profile(
    obj = obj, pars = bestfit, whichPar = c("theta1", outerpars_theta2, "sigma"),
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


################################################################################
######################### Estimations ##########################################
################################################################################
set.seed(110)

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
    pars_sim <- experiment_pars[s, 2:ncol(experiment_pars)]
    
    # Detect the number of available cores and activate cluster
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)
    
    # Parallelise esimations for the simulated data sets
    result <- foreach(
      i = 1:num_sim, 
      .packages = c("dMod")) %dopar% {
        estimation_exp_model(
          data = data_list_exp[[i]], 
          pars_sim = pars_sim)
      }
    # Stop cluster
    stopCluster(cl)
    
    result_list <- c(result_list, list(result))
    names(result_list)[length(result_list)] <- paste0("experiment_",s)
  }
  
  results_list <- c(results_list, list(result_list))
  names(results_list)[length(results_list)] <- paste0("sample_size=", nrow(data_list_exp[[1]]))
}

# Save results
# filename <- "results_exp_model"
# file.path <- paste0(folder.path, "/RDS/", filename, ".rds")
# saveRDS(results_list, file = file.path)
# Due to storage reasons, profile likelihoods are removed from the rds file.
