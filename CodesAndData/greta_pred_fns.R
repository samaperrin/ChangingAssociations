# 
# library(greta)
# X_newdata <- X_full
#  X_assoc_newdata=X_assoc
# site_id_newdata <- as_data(Spatial)

# X_newdata = X_full,X_assoc_newdata = X_assoc,alpha = alpha,beta=beta,lambda_int = lambda_int, lambda_coef = lambda_coef,z=z)

pred_eta <- function(X_newdata, site_id_newdata = NULL, X_assoc_newdata = NULL,alpha=alpha,beta=beta,
                     lambda_int=lambda_int, lambda_coef = lambda_coef,z=z, catchment_effect = NULL) {
  n_sites <- nrow(X_newdata)
  n_species <- ncol(beta)
  
  
  eta <- pred_env(X_newdata,alpha,beta)
  
  if (!is.null(site_id_newdata)) {
    eta <- eta + pred_site_eachspecies(site_id_newdata, catchment_effect, n_species)
  }
  
  if (!is.null(X_assoc_newdata)) {
    eta <- eta + pred_assoc(X_assoc_newdata,lambda_int,lambda_coef,z)
  }
  
  return(eta)
}


pred_eta_null <- function(X_newdata, site_id_newdata = NULL, X_assoc_newdata = NULL,alpha=alpha,beta=beta,
                          lambda_int=NULL,z=NULL, catchment_effect = NULL) {
  n_sites <- nrow(X_newdata)
  n_species <- nrow(lambda_int)
  
  
  eta <- pred_env(X_newdata,alpha,beta)
  
  if (!is.null(site_id_newdata)) {
    eta <- eta + pred_site_eachspecies(site_id_newdata, catchment_effect, n_species)
  }
  
  if (!is.null(X_assoc_newdata)) {
    eta <- eta + pred_assoc_null(X_assoc_newdata,lambda_int,z)
  }
  
  return(eta)
}




pred_env <- function(X_newdata,alpha,beta) {
  XB <- X_newdata %*% beta
  eta <- sweep(XB, 2, alpha, "+")
  return(eta)
}


pred_site <- function(site_id_newdata, catchment_effect, n_species) {
  n_sites <- length(site_id_newdata)
  catchmentB <- sweep(ones(nrow(catchment_effect),n_species), 1, catchment_effect, "+")
  return(catchmentB)
}

pred_site_eachspecies <- function(site_id_newdata, catchment_effect, n_species) {
  n_sites <- length(site_id_newdata)
  catchmentB <- catchment_effect
  return(catchmentB)
}


# X_assoc_sim <- X_assoc

pred_assoc <- function(X_assoc,lambda_int,lambda_coef,z=z) {
  n_sites <- length(X_assoc)
  n_species <- nrow(lambda_coef)
  lambda_effect_rep <- kronecker(as_data(X_assoc), lambda_coef) # (species-lambdas stacked on top of one another)
  lambda_int_rep <- kronecker(ones(n_sites, 1), lambda_int)
  lambda_rep <- lambda_int_rep + lambda_effect_rep
  
  #Do the same thing for our zs so we can perform multiplication
  z_rep <- kronecker(t(z),ones(n_species, 1))
  
  t_e <- rowSums(lambda_rep * z_rep)
  dim(t_e) <- c(n_species, n_sites)
  e <- t(t_e)
  return(e)
}

# a <- matrix(c(1,2,3,4,1,2,3,4),4,2)
# b <- matrix(c(1,2),2,1)
# kronecker(a,b)
# dim(lambda_rep)
# dim(z_rep)
# 
# ex_list <- 1:15
# # 3 species, 5 sites
# dim(ex_list) <- c(3,5)  



pred_assoc_null <- function(X_assoc,lambda_int,z=z) {
  n_sites <- length(X_assoc)
  n_species <- nrow(lambda_int)

  t_e <- lambda_int %*% z
  dim(t_e) <- c(n_species, n_sites)
  e <- t(t_e)
  return(e)
}

# This SHOULD give us our new R matrix, but it needs work.
#X_assoc_pred <- X_assoc_pred[1,]
# lambda_int <- create_lambda(n_species,n_latent)

pred_one_R <- function(X_assoc_pred = NULL, lower_only = FALSE, lambda_int, lambda_coef){
  
  lambda <- lambda_int
  if (!is.null(X_assoc_pred)) {
    lambda <- lambda + lambda_coef * X_assoc_pred
  }
  
  R <- greta::cov2cor(lambda %*% t(lambda))
  
  if (lower_only) {
    lower_idx_R <- lower.tri(R)
    R <- R[lower_idx_R]
  }
  
  return(R)
}


# This SHOULD give us our new R matrix, but it needs work.
#X_assoc_pred <- X_assoc_pred[1,]
# lambda_int <- create_lambda(n_species,n_latent)

pred_one_R_ident <- function(X_assoc_pred = NULL, lower_only = FALSE, lambda_int, lambda_coef){
  
  lambda <- lambda_int
  if (!is.null(X_assoc_pred)) {
    lambda <- lambda + lambda_coef * X_assoc_pred
  }
  
  R <- greta::cov2cor(lambda %*% t(lambda) + diag(n_species))
  
  if (lower_only) {
    lower_idx_R <- lower.tri(R)
    R <- R[lower_idx_R]
  }
  
  return(R)
}


# This SHOULD give us our new R matrix, but it needs work.
#X_assoc_pred <- X_assoc_pred[1,]
# lambda_int <- create_lambda(n_species,n_latent)

pred_one_R_ident_v2 <- function(X_assoc_pred = NULL, lower_only = FALSE, lambda_coef){
  
  if (!is.null(X_assoc_pred)) {
    lambda <- lambda_coef * X_assoc_pred
  }
  
  R <- greta::cov2cor(lambda %*% t(lambda) + diag(n_species))
  
  if (lower_only) {
    lower_idx_R <- lower.tri(R)
    R <- R[lower_idx_R]
  }
  
  return(R)
}



# This SHOULD give us our new R matrix, but it needs work. new versioN

#X_assoc_pred <- X_assoc_pred_sim
#lambda_int=lambda_true_int


pred_R <- function(X_assoc_pred, lower_only = FALSE, lambda_int, lambda_coef) {
  
  if (!inherits(X_assoc_pred, "greta_array")) {
    X_assoc_pred <- as.matrix(X_assoc_pred)
  }
  
  n <- nrow(X_assoc_pred)
  
  result <- list()
  
  for (i in seq_len(n)) {
    result[[i]] <- pred_one_R(X_assoc_pred[i, ], lower_only = lower_only, lambda_int, lambda_coef)
  }
  
  return(result)
  
}



# This SHOULD give us our new R matrix, but it needs work. new versioN

#X_assoc_pred <- X_assoc_pred_sim
#lambda_int=lambda_true_int


pred_R_ident <- function(X_assoc_pred, lower_only = FALSE, lambda_int, lambda_coef) {
  
  if (!inherits(X_assoc_pred, "greta_array")) {
    X_assoc_pred <- as.matrix(X_assoc_pred)
  }
  
  n <- nrow(X_assoc_pred)
  
  result <- list()
  
  for (i in seq_len(n)) {
    result[[i]] <- pred_one_R_ident(X_assoc_pred[i, ], lower_only = lower_only, lambda_int, lambda_coef)
  }
  
  return(result)
  
}


# Bert's curiosity
pred_R_ident_v2 <- function(X_assoc_pred, lower_only = FALSE, lambda_coef) {
  
  if (!inherits(X_assoc_pred, "greta_array")) {
    X_assoc_pred <- as.matrix(X_assoc_pred)
  }
  
  n <- nrow(X_assoc_pred)
  
  result <- list()
  
  for (i in seq_len(n)) {
    result[[i]] <- pred_one_R_ident_v2(X_assoc_pred[i, ], lower_only = lower_only, lambda_coef)
  }
  
  return(result)
  
}


# This SHOULD give us our new R matrix, but it needs work.
#X_assoc_pred <- X_assoc_pred[1,]
# lambda_int <- create_lambda(n_species,n_latent)

pred_one_R_covariance <- function(X_assoc_pred = NULL, lower_only = FALSE, lambda_int, lambda_coef){
  
  lambda <- lambda_int
  if (!is.null(X_assoc_pred)) {
    lambda <- lambda + lambda_coef * X_assoc_pred
  }
  
  R <- lambda %*% t(lambda)
  
  if (lower_only) {
    lower_idx_R <- lower.tri(R)
    R <- R[lower_idx_R]
  }
  
  return(R)
}

pred_R_covariance <- function(X_assoc_pred, lower_only = FALSE, lambda_int, lambda_coef) {
  
  if (!inherits(X_assoc_pred, "greta_array")) {
    X_assoc_pred <- as.matrix(X_assoc_pred)
  }
  
  n <- nrow(X_assoc_pred)
  
  result <- list()
  
  for (i in seq_len(n)) {
    result[[i]] <- pred_one_R_covariance(X_assoc_pred[i, ], lower_only = lower_only, lambda_int, lambda_coef)
  }
  
  return(result)
  
}




## assoc ----
assoc <- function(temp) {
  lambda <- lambda_int + temp * lambda_coef
  cov2cor(lambda %*% t(lambda))
}

## Produces beta estimates, means and quantiles
get_param_ints <- function(parameter,draws,species_names,env_names,ssv=FALSE,quantile_n="mean") {
  n_species <- length(species_names)
  n_env_shared <- length(env_names)
  initial_matrix <- as.matrix(calculate(parameter,values=draws))
  if (quantile_n=="mean") { result_matrix <- apply(initial_matrix, 2, mean)
  } else {
    result_matrix <- apply(initial_matrix, 2, quantile,
                           probs = quantile_n)
  }
  if (ssv==TRUE) {
    ssv_matrix <- matrix(result_matrix,nrow = n_species, ncol=n_species)
    betas <- matrix(diag(ssv_matrix),nrow = n_species,ncol=1)
    colnames(betas) <- "ssv"
  } else {
    betas <- matrix(result_matrix,nrow = nrow(parameter), ncol=ncol(parameter))
    rownames(betas) <- env_names
    colnames(betas) <- species_names
  }
  return(betas)
}


### Produces cooccurrence intervals
get_cooc_ints <- function(draws,species_names,R , quantile_n="mean") {
  n_species <- length(species_names)
  initial_matrix <- as.matrix(calculate(R,values=draws))
  if (quantile_n=="mean") { result_matrix <- apply(initial_matrix, 2, mean)
  } else {
    result_matrix <- apply(initial_matrix, 2, quantile,
                           probs = quantile_n)
  }
  cooc <- matrix(result_matrix,nrow = n_species, ncol=n_species)
  rownames(cooc) <- species_names
  colnames(cooc) <- species_names
  return(cooc)
}

# Creates lambdas with contrained priors
create_lambda <- function(n_species,n_latent) {
  lambda <- zeros(n_species, n_latent)
  diag(lambda) <- normal(0, 1, dim = n_latent, truncation = c(0,Inf))
  lower_idx <- lower.tri(lambda)
  lambda[lower_idx] <- normal(0, 1, dim = sum(lower_idx))
  return(lambda)
}

# Creates lambdas with contrained priors
create_lambda_wide <- function(n_species,n_latent, sd) {
  lambda <- zeros(n_species, n_latent)
  diag_raw <- normal(0, 1, dim = n_latent, truncation = c(0,Inf))
  diag(lambda) <- diag_raw * sd
  lower_idx <- lower.tri(lambda)
  lower_raw <- normal(0, 1, dim = sum(lower_idx))
  lambda[lower_idx] <- lower_raw * sd
  return(lambda)
}

# Produce series of associations matrices based on incremental average temperature increases
pred_correlation <- function(draws,X_assoc_pred,species_names, R_pred, quantile_n="mean") {
  n_species <- length(species_names)
  returned <- list()
  for (i in 1:length(X_assoc_pred)) {
    if(quantile_n == 'mean') {
      temp_R <- colMeans(as.matrix(calculate(R_pred[[i]],values=draws)))
    } else {
      temp_R <- apply(as.matrix(calculate(R_pred[[i]],values=draws)), 2, quantile,
                      probs = quantile_n)
    }
    temp_R_matrix <- matrix(temp_R,n_species,n_species)
    colnames(temp_R_matrix) <- species_names
    rownames(temp_R_matrix) <- species_names
    diag(temp_R_matrix) <- 0
    returned[[i]] <- temp_R_matrix
  }
  return(returned)
}

### Deviance function to calculate deviance of given model
deviance_yp <- function (y, p, log = FALSE) {
  -2 * sum(dbinom(y, 1, p, log))
}


# Produces list of beta CIs and mean with descriptors
get_beta_list <- function(draws,beta_shared,beta_ssv,species_names,env_names,ssv=FALSE){
  demo_betas_upper <- t(as.data.frame(get_param_ints(beta_shared,draws,species_names,env_names,quantile_n = 0.975)))
  demo_betas_lower <- t(as.data.frame(get_param_ints(beta_shared,draws,species_names,env_names,quantile_n = 0.025)))
  demo_betas_mean <- t(as.data.frame(get_param_ints(beta_shared,draws,species_names,env_names)))
  if (ssv==TRUE) {
    demo_betas_ssv_upper <- get_param_ints(beta_ssv,draws,species_names,env_names,ssv=TRUE,quantile_n = 0.975)
    demo_betas_upper <- as.data.frame(cbind(demo_betas_upper,demo_betas_ssv_upper))
    demo_betas_ssv_lower <- get_param_ints(beta_ssv,draws,species_names,env_names,ssv=TRUE,quantile_n = 0.025)
    demo_betas_lower <- as.data.frame(cbind(demo_betas_lower,demo_betas_ssv_lower))
    demo_betas_ssv_mean <- get_param_ints(beta_ssv,draws,species_names,env_names,ssv=TRUE)
    demo_betas_mean <- as.data.frame(cbind(demo_betas_mean,demo_betas_ssv_mean))
  }
  library(reshape2)
  upper_betas_reshaped <- melt(demo_betas_upper)
  lower_betas_reshaped <- melt(demo_betas_lower)
  mean_betas_reshaped <- melt(demo_betas_mean)
  full_betas1 <- merge(lower_betas_reshaped,mean_betas_reshaped,by=c("Var1","Var2"))
  full_betas <- merge(full_betas1,upper_betas_reshaped,by=c("Var1","Var2"))
  names(full_betas) <- c("species","variable","lower","mean","upper")
  return(full_betas)
}


# draws <- demo_draws_extra
# X_newdata <- X_full_new

# This predicts the likelihood of presence/absense of one species at one site, given the environmental values and the presence of other species.


predict_species_perc <- function(draws,X_newdata, site_id_newdata = NULL, X_assoc_newdata = NULL,species_list) {
  
  eta_new <- pred_eta(X_newdata, site_id_newdata,X_assoc_newdata)
  p_new <- ilogit(eta_new)
  
  fill_preds <- matrix(NA,nrow=nrow(p_new),ncol=ncol(p_new))
  
  for (i in 1:length(species_list)) {
    p_new1_draws1 <- calculate(p_new[, i], values=draws)
    p_new1_mn1 <- c(colMeans(as.matrix(p_new1_draws1)))
    fill_preds[,i] <- p_new1_mn1
  }
  colnames(fill_preds) <- species_list
  return(fill_preds)
}


# 
# 
# 
# eta_new <-

# betas_shared <- kauto_mean_betas
# site_env=t(kauto_env[1,])
# occupancy=kauto_occ[1,]
# cooccurrence=kauto_cooccurrence
# focal_species=4

# This predicts the likelihood of presence/absense of one species at one site, given the environmental values and the presence of other species.
#pred_cond(betas_shared=nofa_mean_betas,site_env=orp_env,occupancy=orp_species,cooccurrence=one_row_correlation_temp[[1]],focal_species=2)

pred_cond <- function(betas_shared,betas_ssv=NULL,site_env,occupancy,cooccurrence,focal_species,ssv=FALSE){
  library(tmvtnorm)
  if (ssv==TRUE) {
    betas_ssv_matrix <- matrix(0,length(betas_ssv),length(betas_ssv))
    diag(betas_ssv_matrix) <- betas_ssv
    betas <- cbind(betas_shared,betas_ssv_matrix)
  } else {betas <- betas_shared}
  mean1 <- c(betas %*% site_env)
  sigma1 <- cooccurrence # Need to edit this so it generates cooccurence matrix based on our temeprature
  diag(sigma1) <- 1
  
  lower1 <- c(ifelse(occupancy==0,-Inf,0))
  lower1[focal_species] <- -Inf
  upper1 <- c(ifelse(occupancy==0,0,Inf))
  
  lowerx1 <- c(lower1)
  lowerx1[focal_species] <- 0
  upperx1 <- c(upper1)
  upperx1[focal_species] <- Inf
  
  result <- ptmvnorm(mean=mean1,sigma=sigma1,lower=lower1,upper=upper1,lowerx=lowerx1,upperx=upperx1)
  return(result[1])
}


# This predicts the likelihood of presence/absense of one species at a series of sites, given the environmental values and the presence of other species.
pred_cond_entire <- function(draws,betas_shared,betas_ssv=NULL,site_matrix,occupancy_matrix,temp_vector,focal_species,ssv=FALSE,R_pred) {
  presences <- matrix(NA,ncol=1,nrow=nrow(site_matrix))
  for (i in 1:nrow(site_matrix)) {
    site_env <- matrix(site_matrix[i,],nrow=ncol(site_matrix),ncol=1)
    occupancy <- as.vector(occupancy_matrix[i,])
    cooccurrence <- pred_correlation(draws,temp_vector[i],colnames(occupancy_matrix),R_pred)
    diag(cooccurrence[[1]]) <- 1
    
    presences[i] <- pred_cond(betas_shared,betas_ssv,site_env,occupancy,cooccurrence[[1]],focal_species,ssv)
    if (i %% 10 == 0) {print(paste0(i," runs are complete.")) }
  }
  return(presences)
}
# conditional_entire_value <- pred_cond_entire(nofa_draws,nofa_mean_betas,site_matrix=as.matrix(NOFA_Data$X_val),occupancy_matrix=as.matrix(NOFA_Data$Y_val),temp_vector = as.matrix(NOFA_Data$X_val)[,2],focal_species=4)

# This takes a species correlation matrix and accentuates the values so it's easier to see
accentuate <- function(correlation_matrix_list) {
  unlisted <- unlist(correlation_matrix_list)
  
  times_factor <- 1/max(abs(unlisted))
  accentuated <- lapply(correlation_matrix_list, "*",times_factor)
  return(list(accentuated=accentuated,factor=times_factor))
}

# This narrows down our matrices to the species we want
narrowing <- function(correlation_matrix_list, focal_species_list) {
  narrowed_correlation_matrix_list <- vector("list", length(correlation_matrix_list))
  for(i in 1:length(correlation_matrix_list)) {
    narrowed_correlation_matrix_list[[i]] <- correlation_matrix_list[[i]][focal_species_list,focal_species_list]
  }
  return(narrowed_correlation_matrix_list)
}

# This creates a matrix with the mean associations given different temepratures gradients. Plans are to expand this to include credible intervals.

#association_matrix_list <-  demo_acc_narr
#focal_species <- "Golden_perch"

create_association_gradient <- function(association_matrix_list,focal_species) {
  ascmat_length <- length(association_matrix_list)
  n_species <- nrow(association_matrix_list[[1]])
  asc_map <- matrix(NA,n_species,ascmat_length)
  for(i in 1:ascmat_length) {
    asc_ind <- association_matrix_list[[i]][focal_species,]
    asc_map[,i] <- asc_ind
  }
  rownames(asc_map) <- colnames(association_matrix_list[[1]])
  return(asc_map)
}

# function to plot all associations
create_all_association_gradient <- function(association_matrix_list) {
  ascmat_length <- length(association_matrix_list)
  n_species <- nrow(association_matrix_list[[1]])
  asc_map <- as.data.frame(matrix(NA,n_species*n_species,ascmat_length+2))
  for(i in 1:ascmat_length) {
    asc_ind <- c(association_matrix_list[[i]])
    asc_map[,i] <- asc_ind
  }
  colnames(asc_map) <- c(paste0("temp",seq(1,ascmat_length)),"species1","species2")
  if (is.null(colnames(association_matrix_list[[1]]))) {
    asc_map$species1 <- rep(1:n_species, each=n_species)
    asc_map$species2 <- rep.int(1:n_species, times=n_species)
  } else { 
    asc_map$species1 <- rep(colnames(association_matrix_list[[1]]), each=n_species)
    asc_map$species2 <- rep.int(colnames(association_matrix_list[[1]]), times=n_species)
  }
  asc_map <- asc_map %>% filter(species1 != species2)
  return(asc_map)
}


# Easy function to plot associations given by "create_association_gradient".
#association_gradient_change_matrix <- asc_map
plot_change_asc <- function(association_gradient_change_matrix,legend_position=c(2,1)) {
  n_temps <- ncol(association_gradient_change_matrix)
  n_species <- nrow(association_gradient_change_matrix)
  plot(seq(-1,1,length.out=n_temps) ~ c(1:n_temps),type='n',xlab="Temperature",ylab="Association")
  cl <- brewer.pal(n_species, "Set1")
  for(i in 1:n_temps) {
    lines(seq(1,n_temps,1) ,association_gradient_change_matrix[i,],col=cl[i])
  }
  legend(legend_position[1], legend_position[2], legend=rownames(association_gradient_change_matrix), fill = cl)
}


###############################################################
###############################################################
###                                                         ###
###             NEW ptmvtnorm() FUNCTION                    ###
###                                                         ###
###  This script rebuilds the ptmvtnorm() function from the ###
### tmvtnorm package so that it exposes the algorithm       ###
### argument in the underlying pmvnorm() function.          ###
###                                                         ###
###############################################################
###############################################################

## Define the two tmvtnorm() functions that the function
## requires but for some reason aren't loaded with the
## package. Copied straight from the GitHub page

checkSymmetricPositiveDefinite <- function(x, name="sigma") {
  if (!isSymmetric(x, tol = sqrt(.Machine$double.eps))) {
    stop(sprintf("%s must be a symmetric matrix", name))
  }
  
  if (NROW(x) != NCOL(x)) {
    stop(sprintf("%s must be a square matrix", name))
  }
  
  if (any(diag(x) <= 0)) {
    stop(sprintf("%s all diagonal elements must be positive", name))
  }
  
  if (det(x) <= 0) {
    stop(sprintf("%s must be positive definite", name))
  }
}

# Uses partly checks as in mvtnorm:::checkmvArgs!
checkTmvArgs <- function(mean, sigma, lower, upper)
{
  if (is.null(lower) || any(is.na(lower))) 
    stop(sQuote("lower"), " not specified or contains NA")
  if (is.null(upper) || any(is.na(upper))) 
    stop(sQuote("upper"), " not specified or contains NA")
  if (!is.numeric(mean) || !is.vector(mean)) 
    stop(sQuote("mean"), " is not a numeric vector")
  if (is.null(sigma) || any(is.na(sigma))) 
    stop(sQuote("sigma"), " not specified or contains NA")
  
  if (!is.matrix(sigma)) {
    sigma <- as.matrix(sigma)
  }
  
  if (NCOL(lower) != NCOL(upper)) {
    stop("lower and upper have non-conforming size")
  }
  
  checkSymmetricPositiveDefinite(sigma)
  
  if (length(mean) != NROW(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  
  if (length(lower) != length(mean) || length(upper) != length(mean)) {
    stop("mean, lower and upper must have the same length")
  }
  
  if (any(lower>=upper)) {
    stop("lower must be smaller than or equal to upper (lower<=upper)")
  }
  
  # checked arguments
  cargs <- list(mean=mean, sigma=sigma, lower=lower, upper=upper)
  return(cargs)
}

## Define the ptmvtnorm.new() function

ptmvtnorm.new <- function (lowerx,
                           upperx,
                           mean = rep(0, length(lowerx)),
                           sigma,
                           lower = rep(-Inf, length = length(mean)),
                           upper = rep(Inf, length = length(mean)),
                           maxpts = 25000,
                           abseps = 0.001,
                           releps = 0,
                           algorithm = "GenzBretz"){
  
  cargs <- checkTmvArgs(mean, sigma, lower, upper)
  
  mean <- cargs$mean
  
  sigma <- cargs$sigma
  
  lower <- cargs$lower
  
  upper <- cargs$upper
  
  if(is.null(lowerx) || any(is.na(lowerx))){ 
    stop(sQuote("lowerx"), " not specified or contains NA")
  }
  
  if(is.null(upperx) || any(is.na(upperx))){ 
    stop(sQuote("upperx"), " not specified or contains NA")
  }
  
  if(!is.numeric(lowerx) || !is.vector(lowerx)){ 
    stop(sQuote("lowerx"), " is not a numeric vector")
  }
  
  if(!is.numeric(upperx) || !is.vector(upperx)){ 
    stop(sQuote("upperx"), " is not a numeric vector")
  }
  
  if(length(lowerx) != length(lower) || length(lower) != length(upperx)){ 
    stop("lowerx an upperx must have the same length as lower and upper!")
  }
  
  if(any(lowerx >= upperx)){ 
    stop("lowerx must be smaller than or equal to upperx (lowerx<=upperx)")
  }
  
  f <- pmvnorm(lower = pmax(lowerx, lower),
               upper = pmin(upperx, upper),
               mean = mean,
               sigma = sigma,
               maxpts = maxpts, 
               abseps = abseps,
               releps = releps,
               algorithm = algorithm) / pmvnorm(lower = lower,
                                                upper = upper, 
                                                mean = mean, 
                                                sigma = sigma, 
                                                maxpts = maxpts, 
                                                abseps = abseps, 
                                                releps = releps, 
                                                algorithm = algorithm)
  
  return(f)
  
}


species_relations <- function(correlation_list,species){
  list_species <- lapply(correlation_list,"[",species,)
  names(list_species) <- paste0("temp",seq(1:length(correlation_list)),collapse=NULL)
  species_rels <- as.matrix(data.frame(list_species))
  species_rels <- species_rels[!rownames(species_rels) %in% species, ]
  return(species_rels)
}



window_greta <- function (x, start, finish) {
  mi <- attr(x, "model_info")
  mi <- rlang:::env_clone(mi)
  x <- coda:::window.mcmc.list(x, start, finish)
  mi$raw_draws <- coda:::window.mcmc.list(mi$raw_draws, start, finish)
  attr(x, "model_info") <- mi
  x
}


