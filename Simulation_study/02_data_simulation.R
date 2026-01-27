# R script to simulate the data sets for the simulation study.

################################################################################
################################# Preparation ##################################
################################################################################

# Load R packages 
library("GillespieSSA")
library("foreach")
library("doParallel")
library("parallel")
library("dMod")

###### Folder path of this file to save data #####
#folder.path <- dirname(rstudioapi::getActiveDocumentContext()$path) 

###### Load functions stored in utils.R #####
#functions.path <- ".../utils.R"
#source(functions.path)


################################################################################
############# Simulation of different sample size settings #####################
################################################################################

##### Simulation framework #####
### Framework for Gillespie's algorithm
# Stoichiometry matrix
smat <- matrix(0, nrow = 2, ncol = 4)
smat[1,1] <- 1
smat[1,2] <- -1
smat[2,3] <- 1
smat[2,4] <- -1
rownames(smat) <- c("X1", "X2")

# Propensity vector
props <- c("b1*X1", "d1*X1", "b2*X2", "d2*X2")

### General parameters
# Parameters for the four modification scenarios
p1 <- c(experiment = 1, b1 = 0.2, d1 = 0.11, b2 = 0.2, d2 = 0.11, x1 = 50000, 
        x2 = 50000, sigma = 0.2)
p2 <- c(experiment = 2, b1 = 0.2, d1 = 0.13, b2 = 0.2, d2 = 0.11, x1 = 50000, 
        x2 = 50000, sigma = 0.2)
p3 <- c(experiment = 3, b1 = 0.2, d1 = 0.15, b2 = 0.2, d2 = 0.11, x1 = 50000, 
        x2 = 50000, sigma = 0.2)
p4 <- c(experiment = 4, b1 = 0.2, d1 = 0.21, b2 = 0.2, d2 = 0.11, x1 = 50000, 
        x2 = 50000, sigma = 0.2)
experiment_pars <- as.data.frame(rbind(p1,p2,p3,p4))

# Number of data sets to be simulated per scenario
num_sim <- 100


##### Simulation for sample size = 8 #####
set.seed(8)

# Allocation of measurements
design <- data.frame(input_day = c(0, 0), output_day = c(14, 40), num_pdx = c(2, 2)) 

# List to save simulation
data_list <- list()

# Loop through modification scenarios
for (s in 1:nrow(experiment_pars)) {
  
  # Detect the number of available cores and activate cluster
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  
  # Simulate data sets in parallel 
  data <- foreach(
    i = 1:num_sim,
    .packages = c("GillespieSSA")) %dopar% {
      gillespie_sim_ko_exp(
        theta = experiment_pars[s, c("b1", "d1", "b2", "d2")], 
        x0 = c(X1 = experiment_pars[s, "x1"], X2 = experiment_pars[s, "x2"]),
        sigma = experiment_pars[s, "sigma"], 
        smat = smat, 
        props = props,
        design = design)}
  
  # Stop cluster
  stopCluster(cl)
  
  data_list <- c(data_list, list(data))
  names(data_list)[length(data_list)] <- paste0("experiment_", experiment_pars[s, 1])
}

# Save simulation
# save_list <- list(data_list=data_list, design=design, experiment_pars=experiment_pars)
# saveRDS(save_list, file=paste0(folder.path, "/RDS/simulated_data_size=8_seed=8.rds"))


##### Simulation for sample size = 16 #####
set.seed(16)

# Allocation of measurements
design <- data.frame(input_day = c(0, 0), output_day = c(14, 40), num_pdx = c(4, 4)) 

# List to save simulation
data_list <- list()

# Loop through modification scenarios
for (s in 1:nrow(experiment_pars)) {
  
  # Detect the number of available cores and activate cluster
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  
  # Simulate data sets in parallel 
  data <- foreach(
    i = 1:num_sim,
    .packages = c("GillespieSSA")) %dopar% {
      gillespie_sim_ko_exp(
        theta = experiment_pars[s, c("b1", "d1", "b2", "d2")], 
        x0 = c(X1 = experiment_pars[s, "x1"], X2 = experiment_pars[s, "x2"]),
        sigma = experiment_pars[s, "sigma"], 
        smat = smat, 
        props = props,
        design = design)}
  
  # Stop cluster
  stopCluster(cl)
  
  data_list <- c(data_list, list(data))
  names(data_list)[length(data_list)] <- paste0("experiment_", experiment_pars[s,1])
}

# Save simulation
# save_list <- list(data_list=data_list, design=design, experiment_pars=experiment_pars)
# saveRDS(save_list, file=paste0(folder.path, "/RDS/simulated_data_size=16_seed=16.rds"))


##### Simulation for sample size = 32 #####
set.seed(32)

# Allocation of measurements
design <- data.frame(input_day = c(0, 0), output_day = c(14, 40), num_pdx = c(8, 8)) 

# List to save simulation
data_list <- list()

# Loop through modification scenarios
for (s in 1:nrow(experiment_pars)) {
  
  # Detect the number of available cores and activate cluster
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  
  # Simulate data sets in parallel 
  data <- foreach(
    i = 1:num_sim,
    .packages = c("GillespieSSA")) %dopar% {
      gillespie_sim_ko_exp(
        theta = experiment_pars[s, c("b1", "d1", "b2", "d2")], 
        x0 = c(X1 = experiment_pars[s, "x1"], X2 = experiment_pars[s, "x2"]),
        sigma = experiment_pars[s, "sigma"], 
        smat = smat, 
        props = props,
        design = design)}
  
  # Stop cluster
  stopCluster(cl)
  
  data_list <- c(data_list, list(data))
  names(data_list)[length(data_list)] <- paste0("experiment_", experiment_pars[s,1])
}

# Save simulation
# save_list <- list(data_list=data_list, design=design, experiment_pars=experiment_pars)
# saveRDS(save_list, file=paste0(folder.path, "/RDS/simulated_data_size=32_seed=32.rds"))


##### Simulation for sample size = 64 #####
set.seed(64)

# Allocation of measurements
design <- data.frame(input_day = c(0, 0), output_day = c(14, 40), num_pdx = c(16, 16)) 

# List to save simulation
data_list <- list()

# Loop through modification scenarios
for (s in 1:nrow(experiment_pars)) {
  
  # Detect the number of available cores and activate cluster
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  
  # Simulate data sets in parallel 
  data <- foreach(
    i = 1:num_sim,
    .packages = c("GillespieSSA")) %dopar% {
      gillespie_sim_ko_exp(
        theta = experiment_pars[s, c("b1", "d1", "b2", "d2")], 
        x0 = c(X1 = experiment_pars[s, "x1"], X2 = experiment_pars[s, "x2"]),
        sigma = experiment_pars[s, "sigma"], 
        smat = smat, 
        props = props,
        design = design)}
  
  # Stop cluster
  stopCluster(cl)
  
  data_list <- c(data_list, list(data))
  names(data_list)[length(data_list)] <- paste0("experiment_", experiment_pars[s,1])
}

# Save simulation
# save_list <- list(data_list=data_list, design=design, experiment_pars=experiment_pars)
# saveRDS(save_list, file=paste0(folder.path, "/RDS/simulated_data_size=64_seed=64.rds"))
