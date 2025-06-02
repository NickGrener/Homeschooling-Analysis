setwd("C:/Users/grene/OneDrive/Desktop/6990 Creative Component")

library(units)
library(tidyr)
library(sf)
library(mapview)
library(tmap)
library(leafsync)
library(ggspatial)
library(remotes)
library(readr)
library(car)
library(MASS)
library(Matrix)
library(fdrtool)
library(mvtnorm)
library(dplyr)
library(truncnorm)
library(ggplot2)
library(spdep)
library(here)
library(epinet)

set.seed(300)

OH.acs.df <- readRDS("OH.acs.df")
master.df <- readRDS("master.df")

### Run OLS regression on the 467 rows with complete data and 6 predictors
small.df <- data.frame(scale(master.df[,c("Pupil_Density", "Median_Income", "White", "Experienced_Teachers", "Rep_Vote_2020",
                                          "Grad_Rate_4Yr")]))
cor(small.df)
small.df$Delta_Ratio <- master.df$Delta_Ratio
small.OLS <- lm(Delta_Ratio ~ Pupil_Density + Median_Income + White + Experienced_Teachers + Rep_Vote_2020
                + Grad_Rate_4Yr, data = small.df)
summary(small.OLS)
vif(small.OLS)
#All VIFs are below 2.5; Pupil Density and Experienced Teachers are significant negative predictors;
#White and Tax Effort are significant positive predictors

# Extract residuals and fitted values
residuals <- resid(small.OLS)
fitted_values <- fitted(small.OLS)


# Attach residuals to the original data frame
master.df$OLS.residuals <- NA
master.df$OLS.fitted_values <- NA
master.df$OLS.residuals[!is.na(master.df$Delta_Ratio)] <- residuals
master.df$OLS.fitted_values[!is.na(master.df$Delta_Ratio)] <- fitted_values


### Attach response and predictor data base to map base:
OH.acs.df <- left_join(OH.acs.df, master.df, by = "GEOID")
# Make sure everything checks out in OH.acs.df in terms of district names and county names and then drop redundancies:
OH.acs.df <- OH.acs.df %>% relocate(District_name, County, .after = NAME)
str(OH.acs.df)
# drop empty rows from sf
OH.acs.df <- OH.acs.df[-which(OH.acs.df$invalid), ]

### Plot & save raw response with missing/censored data showing as white
Response.map <- ggplot(OH.acs.df) +
    geom_sf(aes(fill = Delta_Ratio)) +
    scale_fill_gradient2(low = "#1a9850", mid = "#ffffbf", 
                        high = "#d73027", midpoint = 0, na.value = "white") +
    theme_minimal() +
    labs(fill = "Percentage Point Change")
Response.map
ggsave("response_map.png", plot = Response.map, width = 10, height = 8, dpi = 300)

### Plot & save pupil density predictor
Density.map <- ggplot(OH.acs.df) +
    geom_sf(aes(fill = Pupil_Density)) +
    scale_fill_gradient(low = "#feebe2", high = "#7a0177") +
    theme_minimal() +
    labs(fill = "K-12 students / sq. mile")
Density.map
ggsave("density_map.png", plot = Density.map, width = 10, height = 8, dpi = 300)




### Bayesian approach

### First prepare data frame to match main function below:
master.2.df <- master.df %>% 
    dplyr::select(Pupil_Density, Median_Income, White, Experienced_Teachers,
           Rep_Vote_2020, Grad_Rate_4Yr, Delta_Ratio, miss.ind, cen.ind, 
           min.Delta, max.Delta)


### Helper function- MVN prior on beta (log scale)
prior.beta <- function(beta, mu = 0, var = 1000){
    dmvnorm(beta, mean = rep(0, length(beta)), 
            sigma = diag(var, length(beta)), log = TRUE)
}


### Half-normal prior (corresponding to zero-mean normal with var = 100) on sigma^2 (log scale)
prior.s2 <- function(s2, theta = 0.1 * sqrt(pi/2)) {
    dhalfnorm(s2, theta, log = TRUE)
}


### Half-normal prior (corresponding to zero-mean normal with var = 100) on tau^2 (log scale)
prior.tau2 <- function(tau2, theta = 0.1 * sqrt(pi/2)) {
    dhalfnorm(tau2, theta, log = TRUE)
}


### Beta prior for rho (log scale)
prior.rho <- function(rho, shape1 = 18, shape2 = 2) {
    dbeta(rho, shape1, shape2, log = TRUE)
}


### Calculate acceptance rate; restrict samples to a single column (assumed that block updates were used)
acceptance <- function(samps){
    if(is.matrix(samps)){
        s.vec <- samps[,1]
    } else{
        s.vec <- samps
    }
    1 - sum(s.vec[-1] == s.vec[-length(s.vec)]) / (length(s.vec) - 1)
}

### Generate a sample from MVN(mean = Q^{-1} * b, Var = Q^{-1}).
sampleMVN <- function(b, Q, sum.to.zero = FALSE){
    # Input
    # b: transformed mean vector, i.e., mu = Q^{-1} %*% b
    # Q: precision matrix
    # sum.to.zero: Boolean indicating if a sum-to-zero constraint has
    #   been imposed on the random vector.
    # Output
    #   Samples from MVN using Algorithm 2.5 or 2.6 (sum-to-zero) of 
    #   Rue and Held (2005).
    
    if (sum.to.zero){
      # Algorithm 2.6 of Rue and Held
      L.t <- chol(Q); L <- t(L.t); one.vec <- rep(1, nrow(Q))
      w <- solve(L, b)
      mu <- solve(L.t, w)
      z <- rnorm(n = nrow(Q))
      v <- solve(L.t, z)
      x <- mu[,1] + v
      
      V.vec <- solve(L.t, solve(L, one.vec))
      W.val <- c(crossprod(one.vec, V.vec))
      U.vec <- t(V.vec) / W.val 
      
      c.val <- c(crossprod(one.vec, x))
      x.star <- x - (U.vec[1,] *c.val)
      return(x.star)
      
    } else {
      # Algorithm 2.5 of Rue and Held
      L.t <- chol(Q); L <- t(L.t)
      w <- solve(L, b)
      mu <- solve(L.t, w)
      z <- rnorm(n = nrow(Q))
      v <- solve(L.t, z)  
      return(c(mu[,1] + v))
    }
}

### MVN likelihood evaluation. 
dMVN <- function(x, mu, Q, log = TRUE){
    # Input
    #   x: random vector
    #   mu: mean vector
    #   Q: precision matrix
    #   log: Boolean indicating if log-likelihood should be returned.
    
    n <- nrow(Q)
    L.t <- chol(Q)
    x.m.mu <- x - mu
    log.dens <- as.numeric(-n * log(2 * pi) / 2 + 
        sum(log(diag(L.t))) - crossprod(x.m.mu, (Q %*% x.m.mu)) / 2)
    
    if (log){
      return(log.dens)
    } else {
      return(exp(log.dens))
    }
}

### MVN likelihood evaluation to be used when updating tau^2. 
dMVN.tau2 <- function(x, mu, Q, tau2, log = TRUE){

    # Note that this avoids the evaluation of determ(Q), which is useful
    # when using an intrinsic CAR prior for gamma.
    # Input
    #   x: random vector
    #   mu: mean vector
    #   Q: precision matrix
    #   tau2: value of tau^2 (note: Q = Q.tilde / tau^2)
    #   log: Boolean indicating if log-likelihood should be returned.
    
    n <- nrow(Q)
    x.m.mu <- x - mu

    log.dens <- as.numeric(-n * log(2 * pi) - (n * log(tau2)) -
        crossprod(x.m.mu, (Q %*% x.m.mu))) / 2
    
    if (log){
      return(log.dens)    
    } else{
      return(exp(log.dens))
    }
}

#----------------------------------------------------------#
#----------------------------------------------------------#
#----------------------------------------------------------#


### Metropolis-Hastings sampler for Bayesian linear regression- no spatial random effect
# Expectation is that data frame passed to function has as last 5 columns (from the right)
# max for censored, min for censored, censored y indicator, missing y indicator, and y
lmBayes.OH <- function(df, n.post, beta.prop, s2.prop, save.cov = TRUE){
    
    # Split out data frame into relevant parts
    max.cen <- df[,ncol(df)]
    min.cen <- df[,ncol(df) - 1]
    ycen <- df[,ncol(df) - 2]
    ymiss <- df[,ncol(df) - 3]
    y.k <- df[,ncol(df) - 4]
    X <- df[,1:(ncol(df) - 5)]
    
    # Getting error message that the columns are actually lists, so using this fix
    y.k <- as.numeric(y.k[[1]])
    ymiss <- as.logical(ymiss[[1]])
    ycen <- as.logical(ycen[[1]])
    min.cen <- as.numeric(min.cen[[1]])
    max.cen <- as.numeric(max.cen[[1]])
    
    # Dimension of beta vector
    n.p <- ncol(X) + 1
    # Number of observations
    n.y <- nrow(df)
    # Number of censored data
    n.cen <- sum(ycen)
    # Number of missing data
    n.miss <- sum(ymiss)
    
    
    # Standardize columns of predictor matrix and add on column of ones
    X.new <- apply(X, 2, function(column){(column - mean(column)) / sd(column)})
    X.full <- cbind(rep(1, n.y), X.new)
    
    # Generate initial parameter estimates from OLS estimates
    generating.df <- data.frame(cbind(y.k, X.new))
    col_names <- colnames(generating.df)
    formula <- as.formula(paste(col_names[1], paste(col_names[2:(ncol(generating.df))], collapse = " + "), sep = " ~ "))
    lm.res <- lm(formula, data = generating.df)
    beta.k <- unname(coef(lm.res))
    s2.k <- summary(lm.res)$sigma^2
    mu.k <- X.full %*% beta.k
    
    # Create matrix to hold posterior samples and y values
    post.samps <- matrix(NA_real_, nrow = n.post, ncol = n.p + 2)
    y.vectors <- matrix(NA_real_, nrow = n.post / 100, ncol = n.y)   
    
    for (k in 1:n.post){
        
        # 1. Imputation step
        
        y.k[ycen] <- rtruncnorm(n = n.cen, a = min.cen[ycen], b = max.cen[ycen], mean = mu.k[ycen], sd = sqrt(s2.k))
        y.k[ymiss] <- rnorm(n = n.miss, mean = mu.k[ymiss], sd = sqrt(s2.k))
        
        # 2. Posterior sample of beta
        
        # Proposal from MVN
        beta.star <- mvrnorm(n = 1, mu = beta.k, Sigma = beta.prop)
        mu.star <- X.full %*% beta.star
        
        # Acceptance ratio (log scale)
        accept.num <- sum(dnorm(y.k, mean = mu.star, sd = sqrt(s2.k), log = TRUE)) + 
            prior.beta(beta.star)
        
        accept.denom <- sum(dnorm(y.k, mean = mu.k, sd = sqrt(s2.k), log = TRUE)) + 
            prior.beta(beta.k)
        
        accept.beta <- accept.num - accept.denom
        
        
        # Decide whether to accept beta.star or not
        if (log(runif(1)) < accept.beta){
            beta.k <- beta.star
            mu.k <- mu.star
        } 
        
        # 3. Posterior sample of s^2
        
        # Proposal from normal
        s2.star <- rnorm(n = 1, mean = s2.k, sd = s2.prop)
        
        # Acceptance ratio (log scale), but check first if s2.star is valid
        if (s2.star > 0){
            
            accept.num <- sum(dnorm(y.k, mean = mu.k, sd = sqrt(s2.star), log = TRUE)) +
                prior.s2(s2.star)
            
            accept.denom <- sum(dnorm(y.k, mean = mu.k, sd = sqrt(s2.k), log = TRUE)) + 
                prior.s2(s2.k) 
            
            accept.s2 <- accept.num - accept.denom
            
            # Decide whether to accept s2.star or not  
            if (log(runif(1)) < accept.s2){
                s2.k <- s2.star
            } 
        }  
        
        # 4. Likelihood calculation for DIC computation
        D.k = -2 * sum(dnorm(y.k, mean = mu.k, sd = sqrt(s2.k), log = TRUE))
        
        # 5. Save results- all beta and variance values and every 100th y vector
        post.samps[k,] <- c(beta.k, s2.k, D.k)
        if (k %% 100 == 0) {
            y.vectors[k/100,] <- y.k
        }
    }
    
    # Prepare function output
    s2.post <- post.samps[, n.p + 1]
    D.post <- post.samps[, n.p + 2]
    
    # Prepare function output: sample covariance between (transformed) betas
    if (save.cov){
        beta.cov <- cov(post.samps[, c(1:n.p)])
    } else {
        beta.cov <- NA
    }
    
    # Prepare function output: transform beta back to original scale
    beta.post <- t(apply(post.samps[, c(1:n.p)], 1, function(beta, X.m){
        beta.new <- beta[-1] / apply(X, 2, sd)
        beta0.new <- beta[1] - beta.new %*% colMeans(X) 
        c(beta0.new, beta.new)
    }, X.m = X))
    
    
    # Return posterior samples and likelihood evaluation
    list(beta = beta.post, s2 = s2.post, beta.cov = beta.cov, y = y.vectors, D = D.post)
}


### Metropolis-Hastings sampler for Bayesian linear regression with spatial random effect 

lmBayes.Spatial <- function(df, n.post, W.plus, M.plus, var.init = NULL, save.iter = 100,
                       s2.prop, tau2.prop, rho.prop){
    # Metropolis-Hastings sampler for Bayesian linear regression.
    #   The expectation is that data frame passed to function has 
    #   as last 5 columns (from the right):
    #     - max for censored, 
    #     - min for censored, 
    #     - censored y indicator
    #     - missing y indicator, and 
    #     - y (the outcome).
    # Input
    #   df: data frame (see above for expected formatting)
    #   n.post: number of posterior samples
    #   W.plus: sparse row-standardized adjacency matrix 
    #   M.plus: sparse diagonal matrix (M_ii = 1 / W_{i,\cdot})
    #   var.init: list with initial values for (i) beta, (ii) s2, 
    #       (iii) tau2, and (iv) rho.
    #     Notes: 
    #       - var.init = NULL -> initial values are assigned to
    #           pre-specified defaults.       
    #       - If rho is initialized to 1, implements algorithm with ICAR prior.
    #   save.iter: how often to save y, gamma, and resid vectors;
    #       default = every 100 iterations
    #   s2.prop: standard deviation of s2's proposal distribution
    #   tau2.prop: standard deviation of tau2's proposal distribution
    #   rho.prop: standard deviation of rho's proposal distribution
    # Output
    #   Return posterior samples as a list, containing:
    #       (i) beta, (ii) s2, (iii) tau2, (iv) rho, (v) beta.cov,
    #       (vi) y, (vii) gamma, and (viii) resids.
    
    # Split out data frame into relevant parts
    max.cen <- df[,ncol(df)]
    min.cen <- df[,ncol(df) - 1]
    ycen <- df[,ncol(df) - 2]
    ymiss <- df[,ncol(df) - 3]
    y.k <- df[,ncol(df) - 4]
    X <- df[,1:(ncol(df) - 5)]
    
    # Getting error message that the columns are actually lists, so using this fix
    y.k <- as.numeric(y.k[[1]])
    ymiss <- as.logical(ymiss[[1]])
    ycen <- as.logical(ycen[[1]])
    min.cen <- as.numeric(min.cen[[1]])
    max.cen <- as.numeric(max.cen[[1]])
    
    # Dimension of beta vector
    n.p <- ncol(X) + 1
    # Number of observations
    n.y <- nrow(df)
    # Number of censored data
    n.cen <- sum(ycen)
    # Number of missing data
    n.miss <- sum(ymiss)
    
    # Standardize columns of predictor matrix and add on column of ones
    X.new <- apply(X, 2, function(column){(column - mean(column)) / sd(column)})
    X.full <- cbind(rep(1, n.y), X.new)
    
    if (is.null(var.init)){
      # Generate initial parameter estimates (from OLS estimates where appropriate)
      generating.df <- data.frame(cbind(y.k, X.new))
      col_names <- colnames(generating.df)
      formula <- as.formula(paste(col_names[1], paste(col_names[2:(ncol(generating.df))], collapse = " + "), sep = " ~ "))
      lm.res <- lm(formula, data = generating.df)
      beta.k <- unname(coef(lm.res))
      s2.k <- summary(lm.res)$sigma^2
      tau2.k <- .05
      rho.k <- .75
    } else {
      # initial parameter values supplied by var.init
      beta.k <- var.init$beta
      s2.k <- var.init$s2
      tau2.k <- var.init$tau2
      rho.k <- var.init$rho
    }

    
    # Spatial RE matrix created using definition on pg. 42 of Ver Hoef et al. paper
    # Create sparse precision matrix, Q, for CAR prior
    I.mat <- Diagonal(n = nrow(W.plus))
    Q.k <- solve(M.plus, (I.mat - rho.k * W.plus)) / tau2.k
    # Note to self: don't ever solve for the covariance matrix of the prior on gamma- this
    # creates a dense matrix and slows things down by multiple orders of magnitude. 
    
    # initialize gamma to 0 and then initialize mean
    gamma.k <- rep(0.01, n.y)
    mu.k <- as.vector(X.full %*% beta.k + gamma.k)
        
    # Create matrix to hold posterior samples and y values and residuals
    post.samps <- matrix(NA_real_, nrow = n.post, ncol = n.p + 3)
    gamma.vectors <- matrix(NA_real_, nrow = n.post, ncol = n.y) 
    y.vectors <- matrix(NA_real_, nrow = n.post / save.iter, ncol = n.y) 
    gamma.vectors <- matrix(NA_real_, nrow = n.post, ncol = n.y) 
    resid.vectors <- matrix(NA_real_, nrow = n.post / save.iter, ncol = n.y) 
    
    for (k in 1:n.post){
        
        # 1. Imputation step

        y.k[ycen] <- rtruncnorm(n = n.cen, a = min.cen[ycen], b = max.cen[ycen], 
                                mean = mu.k[ycen], sd = sqrt(s2.k))
        y.k[ymiss] <- rnorm(n = n.miss, mean = mu.k[ymiss], sd = sqrt(s2.k))
    
        # 2. (Direct) Posterior sample of beta
        # Based on pg. 5 of Nathan's notes
        
        # sample beta.k
        Q.beta <- (crossprod(X.full) + s2.k / 1000 * diag(ncol(X.full))) / s2.k
        b.beta <- crossprod(X.full, y.k - gamma.k) / s2.k
        beta.k <- sampleMVN(b = b.beta, Q = Q.beta)
        
        # adjust mean for new beta.k
        mu.k <- as.vector(X.full %*% beta.k + gamma.k)

        # 3. Posterior sample of s^2

        # Proposal from normal
        s2.star <- rnorm(n = 1, mean = s2.k, sd = s2.prop)

        # Acceptance ratio (log scale), but check first if s2.star is valid
        if (s2.star > 0){

            accept.num <- sum(dnorm(y.k, mean = mu.k, sd = sqrt(s2.star), log = TRUE)) +
                prior.s2(s2.star)

            accept.denom <- sum(dnorm(y.k, mean = mu.k, sd = sqrt(s2.k), log = TRUE)) +
                prior.s2(s2.k)

            accept.s2 <- accept.num - accept.denom

            # Decide whether to accept s2.star or not
            if (log(runif(1)) < accept.s2){
                s2.k <- s2.star
            }
        }

        
        # 4. Block update of tau2, rho, and gamma
        # Based on pg. 6 of Nathan's notes and algorithm 2.5 on pg. 46 of Rue & Held
        
        # (a) Proposal from normal for tau^2
        tau2.star <- rnorm(n = 1, mean = tau2.k, sd = tau2.prop)
        
        # check that tau2 is in region of appropriate support
        if (tau2.star > 0){
        
          # (b) Propose rho (if rho < 1; otherwise, hold constant)
          if (rho.k < 1){
            # Proposal from normal
            rho.star <- rtruncnorm(n = 1, a = 0, b = 1, mean = rho.k, sd = rho.prop)
            # Update Q
            Q.star <- solve(M.plus, (I.mat - rho.star * W.plus)) / tau2.star
          } else {
            # Update Q
            Q.star <- (Q.k * tau2.k) / tau2.star
          }

          # Used to sample gamma.star and evaluate full conditional density
          Q.gamma.star <- Q.star + (Diagonal(nrow(Q.star)) / s2.k)
          b.gamma.star <- (y.k - (X.full %*% beta.k)) / s2.k
          mu.full.star <- solve(Q.gamma.star, b.gamma.star) # mean of gamma full conditional

          if (rho.k < 1){
              # CAR prior
              gamma.star <- sampleMVN(b = b.gamma.star, Q = Q.gamma.star, sum.to.zero = FALSE)
          } else {
              # ICAR prior with sum-to-zero constraint
              gamma.star <- sampleMVN(b = b.gamma.star, Q = Q.gamma.star, sum.to.zero = TRUE)
          }

          # Update model mean
          mu.star <- as.vector(X.full %*% beta.k + gamma.star)

          # used in evaluation of gamma full conditional
          Q.gamma.k <- Q.k + (Diagonal(nrow(Q.k)) / s2.k)
          b.gamma.k <- (y.k - (X.full %*% beta.k)) / s2.k
          mu.full.k <- solve(Q.gamma.k, b.gamma.k)

          # acceptance rate
          if (rho.k < 1){
            
            # includes prior for rho
            accept.num <- sum(dnorm(y.k, mean = mu.star, sd = sqrt(s2.k), log = TRUE)) +
              dMVN(gamma.star, mu = rep(0, n.y), Q.star) + prior.rho(rho.star) + 
              prior.tau2(tau2.star) + 
              dMVN(gamma.k, mu.full.k, Q.gamma.k, log = TRUE) +
              log(dtruncnorm(rho.k, a = 0, b = 1, mean = rho.star, sd = rho.prop))
            
           accept.denom <- sum(dnorm(y.k, mean = mu.k, sd = sqrt(s2.k), log = TRUE)) +
             dMVN(gamma.k, mu = rep(0, n.y), Q.k) + prior.rho(rho.k) +
             prior.tau2(tau2.k) +
             dMVN(gamma.star, mu.full.star, Q.gamma.star, log = TRUE) + 
             log(dtruncnorm(rho.star, a = 0, b = 1, mean = rho.k, sd = rho.prop))
           
          } else {
            
            # no need to update rho
            accept.num <- sum(dnorm(y.k, mean = mu.star, sd = sqrt(s2.k), log = TRUE)) +
              dMVN.tau2(gamma.star, mu = rep(0, n.y), Q = Q.star, tau2 = tau2.star) +
              prior.tau2(tau2.star) +
              dMVN(gamma.k, mu.full.k, Q.gamma.k, log = TRUE)

            accept.denom <- sum(dnorm(y.k, mean = mu.k, sd = sqrt(s2.k), log = TRUE)) +
              dMVN.tau2(gamma.k, mu = rep(0, n.y), Q = Q.k, tau2 = tau2.k) +
              prior.tau2(tau2.k) +
              dMVN(gamma.star, mu.full.star, Q.gamma.star, log = TRUE)
          }
          
          # log(acceptance ratio)
          accept.tau2 <- accept.num - accept.denom

          # Decide whether to accept proposals or not
          if (log(runif(1)) < accept.tau2){
              
            tau2.k <- tau2.star
            Q.k <- Q.star
            mu.k <- mu.star
            gamma.k <- gamma.star
              
            if (rho.k < 1){
              rho.k <- rho.star
            }
          }
        }
        
        # 5. Save results- all parameter values, likelihood calculation, and every 100th y vector
        post.samps[k,] <- c(beta.k, s2.k, tau2.k, rho.k)
        gamma.vectors[k,] <- gamma.k
        if (k %% save.iter == 0) {
            y.vectors[k/save.iter,] <- y.k
            resid.vectors[k/save.iter,] <- y.k - (X.full %*% beta.k) - gamma.k
        }
    }
    
    # Prepare function output: hyperparameters & sample covariance
    s2.post <- post.samps[, n.p + 1]
    tau2.post <- post.samps[, n.p + 2]
    rho.post <- post.samps[, n.p + 3]
    beta.cov <- cov(post.samps[, c(1:n.p)])
    
    
    # Prepare function output: transform beta back to original scale
    beta.post <- t(apply(post.samps[, c(1:n.p)], 1, function(beta, X.m){
        beta.new <- beta[-1] / apply(X, 2, sd)
        beta0.new <- beta[1] - beta.new %*% colMeans(X) 
        c(beta0.new, beta.new)
    }, X.m = X))
    
    # Return posterior samples
    list(beta = beta.post, s2 = s2.post, tau2 = tau2.post, rho = rho.post, beta.cov = beta.cov, 
         y = y.vectors, gamma = gamma.vectors, resids = resid.vectors)
}

#----------------------------------------------------------#
#--- Fitting the models ------------------------------------#
#----------------------------------------------------------#

### Initial Fit with 10000 samples
post.samps.OH <- lmBayes.OH(master.2.df, n.post = 10000, beta.prop = .0000001 * diag(7), s2.prop = .0001)

acceptance(post.samps.OH$beta)
acceptance(post.samps.OH$s2)

### Secondary fit of 10000 posterior samples with empirical covariance passed in as proposal
post.v2.OH <- lmBayes.OH(master.2.df, n.post = 50000, beta.prop = post.samps.OH$beta.cov, s2.prop = .00005)


# Create sparse W and M matrices (as specified in Ver Hoef et al.)
OH.nb <- poly2nb(OH.acs.df, queen = TRUE)
nblags <- spdep::nblag(neighbours = OH.nb, maxlag = 2)
OH.nb.second <- spdep::nblag_cumul(nblags)
W <- as(nb2mat(OH.nb.second, glist=NULL, style="B", zero.policy= FALSE), "sparseMatrix")
W.plus <- as(nb2mat(OH.nb.second, glist=NULL, style="W", zero.policy= FALSE), "sparseMatrix")
M.plus <- Diagonal(x = 1 / rowSums(W))

#Original working code that just has first-order neighbors
#OH.nb <- poly2nb(OH.acs.df, queen = TRUE)
#W <- as(nb2mat(OH.nb, glist=NULL, style="B", zero.policy= FALSE), "sparseMatrix")
#W.plus <- as(nb2mat(OH.nb, glist=NULL, style="W", zero.policy= FALSE), "sparseMatrix")
#M.plus <- Diagonal(x = 1 / rowSums(W))


### 20000 samples using sparse matrix with second-order neighbors ~ 30 minutes
### 100000 samples using sparse matrix with second-order neighbors ~ 390 minutes
system.time(
CAR.chain <- lmBayes.Spatial(
  master.2.df, n.post = 50000, W.plus = W.plus, M.plus = M.plus,
  var.init = list(
      beta = coefficients(small.OLS),
      s2 = 0.0001,
      tau2 = 0.0001,
      rho = 0.85
  ), save.iter = 25,
  s2.prop = .00004, tau2.prop = .0001, rho.prop = .1
)
)

system.time(
    ICAR.chain <- lmBayes.Spatial(
        master.2.df, n.post = 50000, W.plus = W.plus, M.plus = M.plus,
        var.init = list(
            beta = coefficients(small.OLS),
            s2 = 0.0001,
            tau2 = 0.0001,
            rho = 1
        ), save.iter = 25,
        s2.prop = .00004, tau2.prop = .0001, rho.prop = .1
    )
)


### Acceptance rates
acceptance(post.v2.OH$beta)
acceptance(post.v2.OH$s2)
acceptance(CAR.chain$s2)
acceptance(CAR.chain$tau2)
acceptance(ICAR.chain$s2)
acceptance(ICAR.chain$tau2)

# Posterior Predictive Loss Calculations 
# Based on "Model Choice" article- Gelfand and Ghosh, 1998
# make new vectors/matrices, throwing out burn-in for all models

LIN.betas <- post.v2.OH$beta[10001:50000, ]
LIN.s2 <- post.v2.OH$s2[10001:50000]
CAR.betas <- CAR.chain$beta[10001:50000, ]
CAR.s2 <- CAR.chain$s2[10001:50000]
CAR.tau2 <- CAR.chain$tau2[10001:50000]
CAR.gammas <- CAR.chain$gamma[10001:50000, ]
ICAR.betas <- ICAR.chain$beta[10001:50000, ]
ICAR.s2 <- ICAR.chain$s2[10001:50000]
ICAR.gammas <- ICAR.chain$gamma[10001:50000, ]

# Set design matrix from data frame
X <- master.2.df[,1:(ncol(master.2.df) - 5)]
n.y <- nrow(master.2.df)
X.full <- as.matrix(cbind(rep(1, n.y), X))

# Get predictions for each set of posterior values
M <- 40000
y.new.LIN <- matrix(rep(NA, M * n.y), nrow = M)
y.new.CAR <- matrix(rep(NA, M * n.y), nrow = M)
y.new.ICAR <- matrix(rep(NA, M * n.y), nrow = M)
for (m in 1:M) {
    LIN.mu <- X.full %*% LIN.betas[m, ]
    CAR.mu <- X.full %*% CAR.betas[m, ] + CAR.gammas[m, ]
    ICAR.mu <- X.full %*% ICAR.betas[m, ] + ICAR.gammas[m, ]
    for (j in 1:n.y) {
        y.new.LIN[m, j] <- rnorm(n = 1, mean = LIN.mu[j], sd = LIN.s2[m])
        y.new.CAR[m, j] <- rnorm(n = 1, mean = CAR.mu[j], sd = CAR.s2[m])
        y.new.ICAR[m, j] <- rnorm(n = 1, mean = ICAR.mu[j], sd = ICAR.s2[m])
    }
}

# Calculate penalty term
LIN.penalty <- sum(apply(y.new.LIN, 2, var))
CAR.penalty <- sum(apply(y.new.CAR, 2, var))
ICAR.penalty <- sum(apply(y.new.ICAR, 2, var))

# Calculate goodness-of-fit term: squared error loss for observed y, nothing for missing y, 
# and squared error from nearest endpoint of censored interval for censored data that falls outside
sq.error.LIN <- rep(NA, n.y)
sq.error.CAR <- rep(NA, n.y)
sq.error.ICAR <- rep(NA, n.y)
for (k in 1:n.y) {
    if (master.2.df$miss.ind[k]) {
        sq.error.LIN[k] <- 0 
        sq.error.CAR[k] <- 0 
        sq.error.ICAR[k] <- 0 
    } else if (master.2.df$cen.ind[k]) {
        if (mean(y.new.LIN[,k]) < master.2.df$min.Delta[k]) {
            sq.error.LIN[k] <- (master.2.df$min.Delta[k] - mean(y.new.LIN[,k]))^2 
        } else if (mean(y.new.LIN[,k]) > master.2.df$max.Delta[k]) {
            sq.error.LIN[k] <- (master.2.df$max.Delta[k] - mean(y.new.LIN[,k]))^2 
        } else {
            sq.error.LIN[k] <- 0
        }
        if (mean(y.new.CAR[,k]) < master.2.df$min.Delta[k]) {
            sq.error.CAR[k] <- (master.2.df$min.Delta[k] - mean(y.new.CAR[,k]))^2 
        } else if (mean(y.new.CAR[,k]) > master.2.df$max.Delta[k]) {
            sq.error.CAR[k] <- (master.2.df$max.Delta[k] - mean(y.new.CAR[,k]))^2 
        } else {
            sq.error.CAR[k] <- 0
        }
        if (mean(y.new.ICAR[,k]) < master.2.df$min.Delta[k]) {
            sq.error.ICAR[k] <- (master.2.df$min.Delta[k] - mean(y.new.ICAR[,k]))^2 
        } else if (mean(y.new.ICAR[,k]) > master.2.df$max.Delta[k]) {
            sq.error.ICAR[k] <- (master.2.df$max.Delta[k] - mean(y.new.ICAR[,k]))^2 
        } else {sq.error.ICAR[k] <- 0}
    } else {
        sq.error.LIN[k] <- (master.2.df$Delta_Ratio[k] - mean(y.new.LIN[,k]))^2
        sq.error.CAR[k] <- (master.2.df$Delta_Ratio[k] - mean(y.new.CAR[,k]))^2
        sq.error.ICAR[k] <- (master.2.df$Delta_Ratio[k] - mean(y.new.ICAR[,k]))^2
    }
}

(LIN.PPL <- LIN.penalty + sum(sq.error.LIN))
(CAR.PPL <- CAR.penalty + sum(sq.error.CAR))
(ICAR.PPL <- ICAR.penalty + sum(sq.error.ICAR))


### Trace plots

# beta
plot(CAR.chain$beta[,1], type = "l")
plot(CAR.chain$beta[,2], type = "l")
plot(CAR.chain$beta[,3], type = "l")

white_quants <- quantile(CAR.betas[,4], probs = c(.025, .975))
trace.df <- as.data.frame(CAR.chain$beta)
trace.df$Iteration <- 1:nrow(trace.df)  
white.trace <- ggplot(trace.df) +
    geom_line(aes(x = Iteration, y = White), color = "lightgray") +
    geom_hline(yintercept = white_quants, color = "red", linetype = "dashed") + 
    theme_minimal() +
    labs(y = "White (non-Hispanic) Coefficient", x = "MCMC Iteration")
white.trace
ggsave("white_trace.png", plot = white.trace, width = 10, height = 8, dpi = 300)

plot(CAR.chain$beta[,5], type = "l")
plot(CAR.chain$beta[,6], type = "l")
plot(CAR.chain$beta[,7], type = "l")

# variance components
plot(CAR.chain$s2, type = "l")
plot(CAR.chain$tau2, type = "l")
plot(CAR.chain$rho, type = "l")

# some random effects
plot(CAR.chain$gamma[,1], type = "l")
plot(CAR.chain$gamma[,300], type = "l")
plot(CAR.chain$gamma[,400], type = "l")
plot(CAR.chain$gamma[,500], type = "l")


### Plot residuals from OLS model
ggplot(OH.acs.df) +
    geom_sf(aes(fill = OLS.residuals)) +
    scale_fill_gradient(low = "#ffffcc", high = "#253494", na.value = "white") +
    theme_minimal() +
    labs(title = "Residuals from OLS Model", fill = "Residuals")

### Plot average of residuals (calculated using imputed y values 
### where needed) from selected Bayesian samples
OH.acs.df$MCMC.residuals <- NA
Bayes_residual_means <- apply(CAR.chain$resids, 2, mean)
OH.acs.df$MCMC.residuals[!OH.acs.df$invalid] <- Bayes_residual_means
resid.plot <- ggplot(OH.acs.df) +
    geom_sf(aes(fill = MCMC.residuals)) +
    scale_fill_gradient(low = "#ffffcc", high = "#253494", na.value = "white") +
    theme_minimal() +
    labs(title = "Average of Residuals from Bayesian CAR Model", fill = "Residuals")
resid.plot


# Plot and save mean value map
OH.acs.df$gamma.mean <- colMeans(CAR.gammas)
OH.acs.df$mu <- X.full %*% apply(CAR.betas, 2, mean) + OH.acs.df$gamma.mean
mu.plot <- ggplot(OH.acs.df) +
    geom_sf(aes(fill = mu)) +
    scale_fill_gradient2(low = "#1a9850", mid = "#ffffbf", 
                         high = "#d73027", midpoint = 0) +
    theme_minimal() +
    labs(fill = "Expected Value")
mu.plot
ggsave("mu_map.png", plot = mu.plot, width = 10, height = 8, dpi = 300)


# Plot and save spatial RE map
gamma.plot <- ggplot(OH.acs.df) +
    geom_sf(aes(fill = gamma.mean)) +
    scale_fill_gradient2(low = "#1b7837", mid = "white", high = "#762a83", midpoint = 0) +
    theme_minimal() +
    labs(fill = "Spatial RE Value")
gamma.plot
ggsave("gamma_map.png", plot = gamma.plot, width = 10, height = 8, dpi = 300)

# Look at spatial random effects and try to figure out what is going on with St. Bernard,
# which is an extreme outlier
summary(OH.acs.df$gamma.mean)
gamma.hist <- hist(OH.acs.df$gamma.mean, nclass = 20)
OH.acs.df[which(OH.acs.df$gamma.mean == max(OH.acs.df$gamma.mean)),]

saveRDS(resid.plot, file = "CAR_resid_plot.rds")
saveRDS(gamma.hist, file = "CAR_gamma_hist.rds")
saveRDS(gamma.plot, file = "CAR_gamma_plot.rds")
saveRDS(CAR.chain, file = "CAR_chain.rds")
saveRDS(ICAR.chain, file = "ICAR_chain.rds")
saveRDS(post.v2.OH, file = "LIN_chain.rds")

### Compare residuals from two approaches:
summary(OH.acs.df$OLS.residuals)
summary(OH.acs.df$MCMC.residuals)


### Moran's I analysis
# Assign equal weights to each neighboring district using nb list object created above
OH.lw <- nb2listw(OH.nb, style="W", zero.policy= FALSE)
Morans.I <- moran(OH.acs.df$MCMC.residuals, listw = OH.lw, n = length(OH.nb), Szero(OH.lw))[1]
Morans.I
### Calculate p-values
# residuals (no evidence of remaining spatial autocorrelation)
moran.test(OH.acs.df$MCMC.residuals, OH.lw, alternative="greater")
# spatial random effects (gamma is clearly spatially dependent)
moran.test(OH.acs.df$gamma.mean, OH.lw, alternative="greater")


### Compare summaries of inference on parameters from the three models
summary(small.OLS)
round(t(apply(CAR.betas, 2, quantile, probs = c(.025, .5, .975))), 8)
round(t(apply(CAR.betas, 2, mean)), 8)
quantile(CAR.s2, probs = c(.025, .5, .975))
quantile(CAR.tau2, probs = c(.025, .5, .975))


### Calculate some effective sample sizes
ess(CAR.chain$beta[,1])
ess(CAR.chain$beta[,2])
ess(CAR.chain$beta[,3])
ess(CAR.chain$beta[,4])
ess(CAR.chain$beta[,5])
ess(CAR.chain$beta[,6])
ess(CAR.chain$beta[,7])
ess(CAR.chain$tau2)
ess(CAR.chain$rho)
ess(CAR.chain$s2)

all.sigmas <- matrix(CAR.chain$s2, nrow = 5, byrow = TRUE)
all.beta.ones <- matrix(CAR.chain$beta[,1], nrow = 5, byrow = TRUE)

### Gelman-Rubin diagnostic function
# Reference is Advanced Statistical Computing by Roger Peng, Section 7.4
# as well as pg. 284 of BDA text by Gelmn
# Value of GelmanRubin should be close to one (say 1.1 or 1.2) if stationarity is reached
# Takes matrix of samples for one parameter where each row is a chain 
GelmanRubin <- function(x, burninProportion = 0.1){
    D <- floor(burninProportion * ncol(x))
    x <- x[,(D + 1) : ncol(x)]
    J <- nrow(x)
    L <- ncol(x)
    chain.means <- apply(x, 1, mean)
    chain.vars <- apply(x, 1, var)
    grand.mean <- mean(chain.means)
    within.var <- mean(chain.vars)
    between.var <- sum((chain.means - grand.mean)^2) * L / (J - 1)
    return(((L - 1) * within.var + between.var) / (L * within.var))
}

GelmanRubin(all.sigmas)
GelmanRubin(all.beta.ones)


par(mfrow = c(1, 1))


### Posterior histograms
hist(CAR.chain$beta[,1], nclass = 20, main = "Intercept w/ OLS 95% confidence bounds", xlab = "Intercept")
abline(v = confint(small.OLS)[1,], col="purple")

hist(CAR.chain$beta[,2], nclass = 20, main = "Pupil Density")
abline(v = confint(small.OLS)[2,], col="purple")

hist(CAR.chain$beta[,3], nclass = 20, main = "Median Household Income")
abline(v = confint(small.OLS)[3,], col="purple")

post.df <- as.data.frame(CAR.betas)
white.hist <- ggplot(post.df) +
    geom_histogram(aes(x = White), fill = "skyblue", color = "black") +
    geom_vline(xintercept = white_quants, color = "red", linetype = "dashed", size = 0.5) + 
    theme_minimal() +
    labs(x = "White (non-Hispanic) Coefficient")
white.hist
ggsave("white_hist.png", plot = white.hist, width = 10, height = 8, dpi = 300)

hist(CAR.chain$beta[,4], nclass = 20, main = "% White (non-Hispanic) Students")
abline(v = confint(small.OLS)[4,], col="purple")

hist(CAR.chain$beta[,5], nclass = 20, main = "Teachers with 10+ Years of Experience")
abline(v = confint(small.OLS)[5,], col="purple")

hist(CAR.chain$beta[,6], nclass = 20, main = "Republican Vote Share, 2020 Presidential")
abline(v = confint(small.OLS)[6,], col="purple")

hist(CAR.chain$beta[,7], nclass = 20, main = "4-Year HS Graduation Rate")
abline(v = confint(small.OLS)[7,], col="purple")

hist(CAR.chain$s2, nclass = 20, main = "Error variance w/ OLS estimate")
abline(v = summary(small.OLS)$sigma^2, col="purple")

