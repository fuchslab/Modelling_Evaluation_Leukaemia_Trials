# R script to investigate structural identifiability of the exponential and the
# logistic growth model. 

################################################################################
################################# Preparation ##################################
################################################################################

##### Load R packages #####
library("dMod")
library("ggplot2")
library("latex2exp")


################################################################################
###################### Exponential growth model ################################
################################################################################

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

# Parameter transformations
p_full <- P(eqnvec(
  beta1 = "beta1", beta2 = "beta2", x1 = "exp(x1)", x2 = "exp(x2)", 
  sigma = "exp(sigma)"), condition = "standard") 

# Reparametrisations to make model structurally identifiable
equations <- getEquations(p_full)
equations[["standard"]]["beta2"] <- "theta1+beta1"
equations[["standard"]]["x2"] <- "exp(theta2)*exp(x1)"
p_dMod <- P(equations[["standard"]], condition = "standard")


##### Simulate data for structural identifiability analysis #####
set.seed(300)

# Parameter values
pars <- c(x1 = 50, x2 = 50, beta1 = 0.3, beta2 = 0.5)
sigma <- 0.05

# Time range
times <- seq(0, 50, 1)

# Simulation
out <- (g_dMod*x_dMod)(times, pars)

data_sim <- do.call(rbind, replicate(20, out[[1]], simplify = FALSE))
data_sim <- data_sim[order(data_sim[, 1]),]

# Add normally distributed error to simulated observable
data_sim[, 2] <- data_sim[,2] + rnorm(nrow(data_sim), mean = 0, sd = sigma)
data <- data.frame(
  name = rep("y", nrow(data_sim)), time = data_sim[, 1], value = data_sim[,2], 
  sigma = rep(NA, nrow(data_sim)))
data_dMod <- datalist(standard = data)


##### Identifiability of full model #####
### Parameter estimation
set.seed(33)

# Center of starting values for parameter estimation
pouter <- structure(c(
  pars["beta1"], pars["beta2"], log(pars["x1"]), log(pars["x2"]), log(sigma)),
  names = c("beta1", "beta2", "x1", "x2", "sigma"))

# Objective function
obj <- normL2(data_dMod, g_dMod*x_dMod*p_full, errmodel = error_dMod) + constraintL2(
  pouter, sigma = 3)

# Multiple parameter estimations with randomly chosen starting values
fits <- mstrust(
  obj, pouter, studyname = "multiple_fits", rinit = 0.1, rmax = 10, cores = 4, 
  fits = 50, samplefun = "rnorm", iterlim = 500)
summary(fits)

# Extract estimated parameters of best estimation run
bestfit <- as.parvec(as.parframe(fits))

### Profile likelihoods
profiles <- profile(
  obj = obj, pars = bestfit, whichPar = names(bestfit), alpha = 0.05, cores = 5)
plotProfile(profiles, mode == "data")


##### Identifiability of reparametrised model #####
### Parameter estimation
set.seed(98)

# Center of starting values for parameter estimation
pouter <- structure(
  c(pars["beta2"]-pars["beta1"], pars["beta1"], log(pars["x2"]/pars["x1"]), 
    log(pars["x1"]), log(sigma)),
  names = c("theta1", "beta1", "theta2", "x1", "sigma"))

# Objective function
obj <- normL2(data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(
  pouter, sigma = 3)

# Multiple parameter estimations with randomly chosen starting values
fits <- mstrust(
  obj, pouter, studyname = "multiple_fits", rinit = 0.1, rmax = 10, cores = 5,
  fits = 50, samplefun = "runif", iterlim = 500)
summary(fits)

# Extract estimated parameters of best estimation run
bestfit <- as.parvec(as.parframe(fits))

### Profile likelihoods
profiles <- profile(
  obj = obj, pars = bestfit, whichPar = c("theta1", "theta2", "sigma"), 
  alpha = 0.05, cores = 3)
plotProfile(profiles, mode == "data")


################################################################################
######################### Logistic growth model ################################
################################################################################

##### Generate logistic growth model for dMod #####
# ODE model
ode_equation <- eqnvec(
  x1 = "lambda1*x1*(1-(x1+x2)/K)", x2 = "lambda2*x2*(1-(x1+x2)/K)")
model <- odemodel(ode_equation, modelname = "logistic_model") 
x_dMod <- Xs(model) 

# Observable
observable <- eqnvec(y = "x1/(x1+x2)")
g_dMod <- Y(observable, x_dMod, compile = T, attach.input = F)

# Error model
error_pars <- eqnvec(y = "sigma")
error_dMod <- Y(error_pars, f = c(observable, ode_equation), attach.input = T, compile = T)

# Parameter transformations
p_dMod <- P(
  eqnvec(
    lambda1 = "lambda1", lambda2 = "lambda2", K = "exp(K)", x1 = "exp(x1)", 
    x2 = "exp(x2)", sigma = "exp(sigma)"), condition = "standard")


##### Simulate data for structural identifiability analysis #####
set.seed(59)

# Parameter values
pars <- c(x1 = 50, x2 = 50, lambda1 = 0.1, lambda2 = 0.3, K = 1e3)
sigma <- 0.05

# Time range
times <- seq(0, 50, 1)

# Simulation with dMod 
out <- (g_dMod*x_dMod)(times, pars)
data_sim <- do.call(rbind, replicate(20, out[[1]], simplify = F))
data_sim <- data_sim[order(data_sim[,1]),]

# Add normally distributed error to simulated observable
data_sim[,2] <- data_sim[,2] + rnorm(nrow(data_sim), mean = 0, sigma)
data <- data.frame(name = rep("y",nrow(data_sim)), time = data_sim[,1], 
                   value = data_sim[,2], sigma = rep(NA,nrow(data_sim)))
data_dMod <- datalist(standard = data)


##### Identifiability of full model #####
### Parameter estimation
set.seed(100)


 #sigma=10

#####

# Center of starting values for parameter estimation
pouter <- structure(
  c(pars["lambda1"], pars["lambda2"], log(pars["K"]), log(pars["x1"]), 
    log(pars["x2"]), log(sigma)), 
  names = c("lambda1", "lambda2", "K", "x1", "x2", "sigma"))

# Objective function
obj <- normL2(data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(
  pouter, sigma=10)

# Multiple parameter estimations with randomly chosen starting values
fits <- mstrust(obj, pouter, studyname = "multiple_fits", rinit = 0.1, rmax = 10, 
                cores = 5, fits = 50, samplefun = "rnorm", iterlim = 500)
summary(fits)

# Extract estimated parameters of best estimation run 
bestfit <- as.parvec(as.parframe(fits))

### Profile likelihoods
profiles <- profile(obj = obj, pars = bestfit, whichPar = names(bestfit), 
                    alpha = 0.05, cores = 6)
plotProfile(profiles, mode == "data")


##### Identifiability with fixed x1(0) #####
### Parameter estimation
set.seed(100)


 #sigma=10

#####

# Specify fixed parameter x1
fixed_par <- c("x1") 

# Center of starting values for parameter estimation
pouter <- structure(
  c(pars["lambda1"], pars["lambda2"], log(pars["K"]), log(pars["x1"]), 
    log(pars["x2"]), log(sigma)), 
  names = c("lambda1", "lambda2", "K", "x1", "x2", "sigma"))

# Objective function
obj <- normL2(data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(
  pouter, sigma=10) 

# Multiple parameter estimations with randomly chosen starting values
fits <- mstrust(
  obj, pouter[!names(pouter) %in% fixed_par], fixed = pouter[fixed_par],  
  studyname = "multiple_fits", rinit = 0.1, rmax = 10, cores = 5, fits = 50, 
  samplefun = "rnorm", iterlim = 500)
summary(fits)

# Extract estimated parameters of best estimation run 
bestfit <- as.parvec(as.parframe(fits))

### Profile likelihoods
profiles <- profile(obj = obj, pars = bestfit, fixed = pouter[fixed_par],
                    whichPar = names(bestfit), alpha = 0.05, cores = 5)
plotProfile(profiles, mode == "data")

