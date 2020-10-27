#----------------------#
### 1. MODEL IMPORT ####
#----------------------#

# Set parameters
reduced_species <- TRUE                # All species in model or only those in >10% of lakes?
spatial_component <- FALSE             # Introduce spatial component?
association <- TRUE                    # Estimate associations over temperature?

# Load libraries
library(greta)
library(tensorflow)
library(dplyr)
library(corrplot)
library(ggplot2)
library(magrittr)
library(gridExtra)


# Import model and scripts
NOFA_JSAM_results <- readRDS("ModelFile.RDS")
source("./CodesAndData/greta_pred_fns.R")

# Define data and parameters
NOFA_Data <- NOFA_JSAM_results$data

model <- NOFA_JSAM_results$model
parameters <- NOFA_JSAM_results$parameters

Y <- NOFA_Data$Y[,parameters$species_names]

# Define model components
draws <- model$draws
R_lower <- model$R_lower
R <- model$R_ident
beta <- model$beta
alpha <- model$alpha
p <- model$p
lambda_int <- model$lambda_int
R_pred_ident <- model$R_pred_ident
valley_offset <- model$valley_offset
  
  
n_species <- parameters$n_species
valid_species <- parameters$valid_species

temp_sd <- 10/sd(NOFA_Data$full_data$eurolst_bio10)
X_assoc_pred <- parameters$X_assoc_pred
X_assoc_pred_nonSD <- round(mean(NOFA_Data$full_data$eurolst_bio10) + 
                              X_assoc_pred * sd(NOFA_Data$full_data$eurolst_bio10),1)/10
#X_assoc_pred_test <- seq(-3.5,2.5,length.out=11)

#---------------------#
### 2. MODEL CHECK ####
#---------------------#

# Check convergence

calc_R <- calculate(R_lower,values=draws)
gelman_sims_R <- coda::gelman.diag(calc_R,multivariate = FALSE)
nrow(gelman_sims_R$psrf[gelman_sims_R$psrf[,1] < 1.1,])/nrow(gelman_sims_R$psrf)

calc_beta <- calculate(beta,values=draws)
gelman_sims_beta <- coda::gelman.diag(calc_beta,multivariate = FALSE)
nrow(gelman_sims_beta$psrf[gelman_sims_beta$psrf[,1] < 1.1,])/nrow(gelman_sims_beta$psrf)

if (spatial_component == TRUE) {
  calc_catchment <- calculate(valley,values=draws)
  gelman_sims_catchment <- coda::gelman.diag(calc_catchment,multivariate = FALSE)
  nrow(gelman_sims_catchment$psrf[gelman_sims_catchment$psrf[,1] < 1.1,])/nrow(gelman_sims_catchment$psrf)
}


#----------------------------#
### 3. Parameter Analysis ####
#----------------------------#

nofa_beta_list <- get_beta_list(draws,beta_shared = beta,species_names=parameters$species_names,
                                env_names=parameters$env_names)

nofa_beta_list$significance <- as.factor(ifelse(nofa_beta_list$lower > 0 | nofa_beta_list$upper < 0, 0,1))
nofa_beta_sig <- filter(nofa_beta_list,
                        significance == 0)
filter(nofa_beta_sig,variable=="temperature")


ggplot(nofa_beta_list, aes(x=variable, y=mean,col=variable, linetype = significance)) + 
  ylim(-4,4) +
  geom_linerange(aes(ymin=lower, ymax=upper)) +
  facet_wrap(~ species) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  geom_hline(yintercept=0)

ggplot(nofa_beta_list, aes(x=species, y=mean,col=species, linetype = significance)) + 
  ylim(-2,2) +
  geom_linerange(aes(ymin=lower, ymax=upper)) +
  facet_wrap(~ variable) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank()) + 
  geom_hline(yintercept=0)

# # Check if and valley values are significant
# valley_mean <- matrix(apply(as.matrix(calc_catchment),2, mean), 
#                       nrow = nrow(valley), ncol = ncol(valley))
# valley_lower <- matrix(apply(as.matrix(calc_catchment),2, quantile, probs = 0.025), 
#                        nrow = nrow(valley), ncol = ncol(valley))
# valley_upper <- matrix(apply(as.matrix(calc_catchment),2, quantile, probs = 0.975), 
#                        nrow = nrow(valley), ncol = ncol(valley))
# valley_values <- valley_mean
# 
# valley_values[valley_lower < 0 & valley_upper > 0] <- 0

#------------------------------#
### C. Corellation Analysis ####
#-------------------..---------#

# Get correlaton matrix at mean temperature first
nofa_cooccurrence_R <- get_cooc_ints(draws,parameters$species_names, R=R)
diag(nofa_cooccurrence_R) <- 0
nofa_cooccurrence_R_lower <- get_cooc_ints(draws,parameters$species_names,
                                           R, quantile=0.025)
nofa_cooccurrence_R_upper <- get_cooc_ints(draws,parameters$species_names,
                                           R, quantile=0.975)
# Stipulate species order early on
species_order <- row.names(nofa_cooccurrence_R)
common_names <- c("Burbot", "Roach", "Zander", "Bleak", "Crucian carp", "Perch","Vendace", "Bream",
                  "Pike", "Whitefish", "Rudd", "Tench", "Ruffe", "Brown trout", "Artic charr")

# Now let's show species with significant associations
nofa_cooccurrence_sig <- nofa_cooccurrence_R
nofa_cooccurrence_sig[nofa_cooccurrence_R_lower < 0 & nofa_cooccurrence_R_upper > 0] <- 0
corrplot(nofa_cooccurrence_R, method = "color",type="lower",order="original", diag=FALSE,
         cl.lim=c(-1,1))
corrplot(nofa_cooccurrence_sig, method = "color",type="lower", order="original",diag=FALSE,
         cl.lim=c(-1,1))



# Now map correlations in species over the temperature gradients
nofa_correlation_temp <- pred_correlation(draws,X_assoc_pred,parameters$species_names,R_pred_ident)
nofa_correlation_temp_lower <- pred_correlation(draws,X_assoc_pred,
                                                parameters$species_names,R_pred_ident,
                                                quantile_n=0.025)
nofa_correlation_temp_upper <- pred_correlation(draws,X_assoc_pred,
                                                parameters$species_names,R_pred_ident,
                                                quantile_n=0.975)

# Save the above 3 lists to speed up future analysis
correlation_CIs <- list(lower = nofa_correlation_temp_lower, mean = nofa_correlation_temp,
                        upper = nofa_correlation_temp_upper)

# Isolate correlation matrices for associations at specific temperatures
# for easier visualisation
NCT_mean <- nofa_correlation_temp[[7]][species_order,species_order] %>%
  set_colnames(common_names) %>%
  set_rownames(common_names)
NCT_low <- nofa_correlation_temp[[1]][species_order,species_order] %>%
  set_colnames(common_names) %>%
  set_rownames(common_names)
NCT_high <- nofa_correlation_temp[[11]][species_order,species_order] %>%
  set_colnames(common_names) %>%
  set_rownames(common_names)

# Plto associations at low, high and mean values
corrplot(NCT_mean, method = "color",type="upper",order="original",diag=FALSE, tl.col="black",
         cl.lim=c(-1,1),is.corr=FALSE,mar=c(0,0,1,1))
corrplot(NCT_low, method = "color",type="upper",order="original",diag=FALSE, tl.col="black",
         cl.lim=c(-1,1),is.corr=FALSE,mar=c(0,0,1,1))
corrplot(NCT_high, method = "color",type="upper",order="original",diag=FALSE, tl.col="black",
         cl.lim=c(-1,1),is.corr=FALSE,mar=c(0,0,1,1))



# The following function shows associations between 2 species
view_interaction <- function(interacting_species, correlation_CIs, X_assoc_pred, ylim_p = c(-1,1), title.int=NA,
                             xlab.int = NA, ylab.int=NA) {
  
  all_assocs <- create_all_association_gradient(correlation_CIs$mean)
  all_assocs_lower <- create_all_association_gradient(correlation_CIs$lower)
  all_assocs_higher <- create_all_association_gradient(correlation_CIs$upper)
  
  interacting_species_frame <- as.data.frame(matrix(NA, nrow = length(X_assoc_pred), ncol = 4))
  interacting_species_frame[,1] <- X_assoc_pred_nonSD
  interacting_species_frame[,2] <- all_assocs_lower %>%
    filter(species1 == interacting_species[1] & species2 == interacting_species[2]) %>%
    select(1:length(X_assoc_pred)) %>% t()
  interacting_species_frame[,3] <- all_assocs %>%
    filter(species1 == interacting_species[1] & species2 == interacting_species[2]) %>%
    select(1:length(X_assoc_pred)) %>% t()
  interacting_species_frame[,4] <- all_assocs_higher %>%
    filter(species1 == interacting_species[1] & species2 == interacting_species[2]) %>%
    select(1:length(X_assoc_pred)) %>% t()
  colnames(interacting_species_frame) <- c("temp","low","mean","high")
  
  D0 <- ggplot(interacting_species_frame, aes(temp, mean)) +
    geom_line(aes(colour = "red"), # Line type depends on cond
              size = 1.5) +
    ylim(ylim_p) +
    geom_ribbon(aes(ymin = low, ymax = high, fill = "red"), alpha = .35) +
    theme_bw() +
    theme(legend.position = "none") + 
    scale_x_continuous(breaks=seq(7.5,16.5,3))
  if (is.na(title.int)) {D1 <- D0 + 
    labs(title = paste0("Relationship between ",interacting_species[1],
                        " and ",interacting_species[2]))} else {
                          D1 <- D0 + labs(title = title.int)
                        }
  
  if (is.na(xlab.int)) {
    D2 <- D1 +
      xlab("Standardised temperature")
  } else {
    D2 <- D1 +
      xlab(xlab.int) 
  }
  
  if (is.na(ylab.int)) {
    D3 <- D2 +
      ylab("Association")
  } else {
    D3 <- D2 + ylab(ylab.int)
  }
  D3
}

# Example
view_interaction(c("Rutilus_rutilus","Lota_lota"), correlation_CIs, X_assoc_pred, title = "")



# Calculate correlation in environmental responses
beta_mat <- as.matrix(calc_beta)
all_enviro_cor_mat <- array(0, dim = c(nrow(beta_mat), n_species, n_species))

env_matrix <- as.matrix(cbind(NOFA_Data$X[,-1], NOFA_Data$X$temperature^2)) # all env

# This uses the same method of calculation as Francis Hui's boral package.
for (i in 1:nrow(draws$`11`)) {
  cw_X_coefs <- matrix(beta_mat[i,],nrow=n_species)
  enviro.linpreds <- tcrossprod(env_matrix, as.matrix(cw_X_coefs[,c(1:7,8)]))
  #  enviro.linpreds <- tcrossprod(env_matrix, as.matrix(cw_X_coefs[,c(2,9)]))
  all_enviro_cor_mat[i,,] <- cor(enviro.linpreds)
}

enviro_cor_mat <- matrix(0, n_species, n_species)

for (j in 1:n_species) {
  for (j2 in 1:n_species) {
    enviro_cor_mat[j, j2] <- mean(all_enviro_cor_mat[, j, j2])
  }
}
colnames(enviro_cor_mat) <- colnames(Y)
rownames(enviro_cor_mat) <- colnames(Y)

enviro_cor_mat_ordered <- enviro_cor_mat[species_order, species_order] %>%
  set_colnames(common_names) %>%
  set_rownames(common_names)

enviro_cor_plot <- corrplot(enviro_cor_mat_ordered, method = "color",type="upper",order="original",diag=FALSE,
                            tl.col='black',mar=c(0,0,1,1))


#-------------------#
### D. Model fit ####
#-------------------#



### Calculating deviance

# Now we look at deviance
Y_dev <- NOFA_Data$full_data[,parameters$species_names]

# intercept-only NULL deviance
p_null <- colMeans(Y_dev)
p_null <- kronecker(rep(1, parameters$n_sites), t(p_null))
deviance_yp(as.matrix(Y_dev),p_null)

# environment-only deviance
dev_env <- 0
for (i in 1:n_species) {
  dat <- data.frame(y = Y_dev[, i],
                    NOFA_Data$full_data[,parameters$env_names])
  m <- glm(y ~ ., data = dat,
           family = stats::binomial("probit"))
  pred <- predict(m, type = "response")
  dev_env <- dev_env + deviance_yp(dat$y, pred)
}
dev_env



# Model deviance
# Note: THis is a fairly time-consuming process, and it takes a lot of memory.
# I recommend deleting any large files you no longer need from your R workspace before
# starting it.
p_mat <- matrix(ncol=ncol(p),nrow=nrow(p))

time1 <- Sys.time()
for (i in 1:nrow(p_mat)) {
  p_temp <- calculate(p[,i], values = draws)
  p_out <- as.matrix(p_temp)
  rm(p_temp)
  p_out_medians <- apply(p_out,2,median)
  rm(p_out)
  p_mat[,i] <-p_out_medians
  
  # Quick function to let us know how the loop is progressing
  time2 <- Sys.time()
  difftime_1 <- round(as.numeric(difftime(time2, time1,
                                          units = "mins")),4)
  if (i %% 3 == 0) {print(paste0("Run ", i, " of ",  nrow(p_reduced), " finished in ",difftime_1, " minutes. Estimated ", round(difftime_1*nrow(p_reduced)/i-difftime_1,3), " minutes left.") )}
}

deviance_yp(as.matrix(Y_dev),p_mat)

saveRDS(p_mat, file = "./CodesAndData/reducedmodel_devianceValues.RDS")
