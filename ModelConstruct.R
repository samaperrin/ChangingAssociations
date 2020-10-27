#-------------------------#
### 1. MODEL CONSTRUCT ####
#-------------------------#

# Set parameters
reduced_species <- TRUE                # All species in model or only those in >10% of lakes?
spatial_component <- FALSE             # Introduce spatial component?
association <- TRUE                    # Estimate associations over temperature (FALSE implies 
                                       # just measured at mean temperature)?

#---------------------------#
### A. Data Construction ####
#---------------------------#

# Data import
library(greta)
library(dplyr)
load("./CodesAndData/NOFA_Data_WC.rda")
source("JSAM_Functions.R")

# Colinnearity check
source("http://www.sthda.com/upload/rquery_cormat.r")
colin <- rquery.cormat(NOFA_Data$X, type = "flatten", graph = FALSE)
colin$r %>% filter(abs(cor) > 0.4)

# Check to see data distribution of anything suspicious on collinearity
with(NOFA_Data$X,plot(HFP_sc,temperature))

# Remove elevation data
NOFA_Data$X <- NOFA_Data$X[,-1]

# The first task here is to create a matrix which shows numbers of cooccurrences between different species.
# The idea being that if they only cooccur below a certain limit, there is insufficient data to
# gauge associations (Tikhonov et al)
n_species <- ncol(NOFA_Data$Y)

# Isolate environmental data, get rid of elevation (unnecessary)
X_base <- NOFA_Data$X %>%
  mutate(temp_sq = temperature^2)

# Take species data
if (reduced_species == TRUE) {
  valid_species <- which(apply(NOFA_Data$Y, 2, mean) > 0.1)
  Y_base <- NOFA_Data$Y[,valid_species]
  n_species <- length(valid_species)
}

# Load spatial component if necessary  
if (spatial_component == TRUE) {
  load("./CodesAndData/sites.RDA")
  sites <- sites[,c("waterBodyID","inout")]
  site_no <- as_data(sites["inout"])
  n_valleys <- length(unique(site_no))
}

# Convert main tables to greta
X <- as_data(as.matrix(X_base))
Y <- as_data(as.matrix(Y_base))

# Set some parameters
nofa_env_names <- colnames(X_base)
nofa_species_names <- colnames(Y_base)

# Set species association vector up
X_assoc <- as_data(X_base$temperature)

if (association == TRUE) {
  minXassoc <- floor(min(X_base$temperature)*10)/10
  maxXassoc <- ceiling(max(X_base$temperature)*10)/10
  X_assoc_pred <- seq(minXassoc,maxXassoc,length.out=11)
}

# Define dimensions
n_sites<- nrow(X)
n_env <- ncol(X)
n_latent <- 3
n_species <- ncol(Y)

# Finally, group parameters for later
parameters <- list(n_sites = n_sites, n_env = n_env, n_latent=n_latent, n_species=n_species, 
                   species_names = nofa_species_names, env_names = nofa_env_names,
                   valid_species = valid_species)
if (association == TRUE) {parameters[["X_assoc_pred"]] <-  X_assoc_pred}
if (spatial_component == TRUE) {parameters[["site_no"]] <- site_no}

#-----------------------#
### B. Model fitting ####
#-----------------------#

alpha <- normal(-2, 1, dim = n_species)
beta <- normal(0, 10, dim = c(n_env, n_species))

z <- normal(0, 1, dim = c(n_latent, n_sites))
  
# Define our priors for lambda. Lambda_int will create the correlation matrix later on.
lambda_int <- create_lambda_wide(n_species,n_latent, 10)
lambda_coef <- create_lambda_wide(n_species,n_latent, 10)

# Hierarchical component
if (spatial_component == TRUE) {
  valley_sd <- normal(0, 10, truncation = c(0,Inf))
  valley_offset_raw <- normal(0, 1, dim = c(n_valleys,n_species))
  valley_offset <- valley_offset_raw * valley_sd
  valley_effect <- valley_offset[as.factor(site_no),]
}

# Now produce our prediction function for our later correlation matrix
if (association == TRUE) {
  R_pred_ident <- pred_R_ident(X_assoc_pred,lambda_int=lambda_int, lambda_coef = lambda_coef)
}

# Now just put the model together. THis can get tricky.
if (association == TRUE & spatial_component == TRUE) {
  eta <- pred_eta(X_newdata = X, 
                  X_assoc_newdata = X_assoc,alpha = alpha,beta=beta,
                  lambda_int = lambda_int, lambda_coef = lambda_coef,z=z,
                  catchment_effect = valley_effect, site_id_newdata = site_no)
} else if (association == TRUE & spatial_component == FALSE) {
  eta <- pred_eta(X_newdata = X, 
                  X_assoc_newdata = X_assoc,alpha = alpha,beta=beta,
                  lambda_int = lambda_int, lambda_coef = lambda_coef,z=z)
} else if (association == FALSE & spatial_component == TRUE) {
  eta <- pred_eta_null(X_newdata = X,
                       X_assoc_newdata = X_assoc, alpha = alpha,beta=beta,
                       lambda_int = lambda_int,z=z,
                       catchment_effect = valley_effect, site_id_newdata = site_no)
} else {
  eta <- pred_eta_null(X_newdata = X,
                       X_assoc_newdata = X_assoc, alpha = alpha,beta=beta,
                       lambda_int = lambda_int,z=z)
}


# Define the distribution
p <- iprobit(eta)
distribution(Y) <- bernoulli(p)

# Define our correlation matrix
R_ident <- greta::cov2cor(lambda_int %*% t(lambda_int) + diag(n_species))
lower_idx_R <- lower.tri(R_ident)
R_lower <- R_ident[lower_idx_R]

# Define and plot  model
nofa_gmodel <- model(beta)
plot(nofa_gmodel)

# Aaaaand perform MCMC
Lmin <- 40
Lmax <- round(Lmin*1.5)
nofa_Gdraws <- mcmc(nofa_gmodel,chains = 1, n_samples = 3000, warmup = 2000 
                    ,sampler = hmc(Lmin = Lmin, Lmax = Lmax))

#---------------------------------#
### C. Model convergence check ####
#---------------------------------#

calc_R_lower <- calculate(R_lower,values=nofa_Gdraws)
gelman_sims_R_lower <- coda::gelman.diag(calc_R_lower,multivariate = FALSE)
nrow(gelman_sims_R_lower$psrf[gelman_sims_R_lower$psrf[,1] < 1.1,])/nrow(gelman_sims_R_lower$psrf)

calc_beta <- calculate(beta,values=nofa_Gdraws)
gelman_sims_beta <- coda::gelman.diag(calc_beta,multivariate = FALSE)
nrow(gelman_sims_beta$psrf[gelman_sims_beta$psrf[,1] < 1.1,])/nrow(gelman_sims_beta$psrf)

calc_lambda <- calculate(lambda_coef,values=nofa_Gdraws)
gelman_sims_lambda <- coda::gelman.diag(calc_lambda,multivariate = FALSE)
nrow(gelman_sims_lambda$psrf[gelman_sims_lambda$psrf[,1] < 1.1,])/nrow(gelman_sims_lambda$psrf)



#------------------------------#
### D. Consolidate and save ####
#------------------------------#


NOFA_JSAM_model_subset <- list(draws = nofa_Gdraws, beta = beta, alpha = alpha, 
                               R_ident = R_ident, R_lower = R_lower, p = p, 
                               eta = eta, z = z, lambda_int = lambda_int
                               )
if (association == TRUE) {
  NOFA_JSAM_model_subset[["R_pred_ident"]] <- R_pred_ident
  NOFA_JSAM_model_subset[["lambda_coef"]] <- lambda_coef
}

if (spatial_component == TRUE) {
  NOFA_JSAM_model_subset[["valley_offset"]] <- valley_offset
}

NOFA_JSAM_subset <- list(model=NOFA_JSAM_model_subset, data=NOFA_Data, parameters=parameters)

file_name <- getPass("Name the file")
saveRDS(NOFA_JSAM_subset,paste0("./SelectYourRepository/",file_name,".RDS"))
