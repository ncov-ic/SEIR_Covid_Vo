# Run the SEIR model specified in the paper
wrapper_model <- function (iter_main, data,
                           fixed_parameters, all_looped_parameters,
                           limits, random_walk_rate,
                           mcmc_iterations, sample_spacing, dir_output, id_chain) {

  # extract parameter combination for this iteration of the loop
  looped_parameters <- all_looped_parameters[iter_main, ]

  # timesteps to solve the SEIR model
  times <- seq(from = fixed_parameters$tSeed,
               to = fixed_parameters$time2,
               by = 1L)

  # initialise loop through chains
  loop.acceptance_rate  <- matrix(nrow = 0, ncol = ncol(limits) + 1)
  loop.random_walk_rate <- matrix(nrow = 0, ncol = ncol(limits) + 2)
  loop.posterior        <- matrix(nrow = 0, ncol = ncol(limits) + 3)
  for (idc in id_chain) {

    # set seed
    set.seed(idc)
    print(paste("## Loop number", iter_main, "chain number", idc))

    # initialise fitted parameters within MCMC limits, but make sure there's space for gamma!
    # N.B. inv_nu + inv_delta + inv_gamma = generation_time
    repeat {
      fitted_parameters <- apply(limits, 2, function(lim) runif(1, lim[1], lim[2]))
      if (sum(fitted_parameters[c("inv_nu", "inv_delta")]) < 0.99 * fixed_parameters$generation_time)
        break
    }

    # run mcmc
    results <- run_mcmc(data              = data,
                        times             = times,
                        fixed_parameters  = fixed_parameters,
                        looped_parameters = looped_parameters,
                        fitted_parameters = fitted_parameters,
                        limits            = limits,
                        random_walk_rate  = random_walk_rate,
                        mcmc_iterations   = mcmc_iterations,
                        sample_spacing    = sample_spacing)

    # clean results
    loop.acceptance_rate <- rbind(loop.acceptance_rate,
                                  c(results$acceptance_rate,
                                    chain = idc))
    loop.random_walk_rate <- rbind(loop.random_walk_rate,
                                   cbind(results$random_walk_rate,
                                         chain = idc,
                                         iter_mcmc = seq.int(nrow(results$fitted_parameters))))
    loop.posterior <- rbind(loop.posterior,
                            cbind(ll = results$log_likelihood,
                                  results$fitted_parameters,
                                  chain = idc,
                                  iter_mcmc = seq.int(nrow(results$fitted_parameters))))

  }

  # write to file in R format
  res <- list(acceptance_rate   = loop.acceptance_rate,
              random_walk_rate  = loop.random_walk_rate,
              posterior         = loop.posterior,
              looped_parameters = unlist(looped_parameters),
              fixed_parameters  = unlist(fixed_parameters))
  saveRDS(res, file.path(dir_output, paste0(
    "results_", sprintf("%02d", iter_main), ".rds")))

}

# define mcmc routine
run_mcmc <- function (data, times,
                      fixed_parameters, looped_parameters, fitted_parameters,
                      limits, random_walk_rate,
                      mcmc_iterations, sample_spacing) {

  # initialise model parameter vectors
  fitted_params <- unlist(fitted_parameters)
  looped_params <- unlist(looped_parameters)
  fixed_params  <- unlist(fixed_parameters)

  # compute "old" MCMC values
  params <- c(fitted_params, looped_params, fixed_params)
  sol <- model_gen(user = as.list(params))$run(times)
  old_llike <- compute_loglikelihood(sol, data, params)

  # initialise MCMC parameters
  nr_fitted <- length(fitted_params)
  nr_stored <- mcmc_iterations%/%sample_spacing

  # initialise output variables
  stored_fit <- matrix(0, ncol = nr_fitted, nrow = nr_stored,
                       dimnames = list(NULL, names(fitted_params)))
  stored_rw  <- matrix(0, ncol = nr_fitted, nrow = nr_stored,
                       dimnames = list(NULL, names(fitted_params)))
  stored_llike <- rep(0, length.out = nr_stored)
  nr_accepted  <- rep(0, length.out = nr_fitted) # number of accepted iterations

  # insert into output variables
  stored_fit[1, ] <- fitted_params
  stored_rw[1, ]  <- random_walk_rate
  stored_llike[1] <- old_llike

  # MCMC ----------------------------------------------------------------------
  for (iter in 2:mcmc_iterations) {
    for (parID in 1:nr_fitted) {

      # sample new value
      nameID <- names(fitted_params[parID])
      old_value <- fitted_params[parID]
      if (nameID %in% c("inv_nu", "inv_delta")) {
        nameotherID <- setdiff(c("inv_nu", "inv_delta"), nameID)
        repeat {
          new_value <- old_value * exp(random_walk_rate[parID] * rnorm(1))   # change old value slightly
          if(new_value > limits[1, parID] &
             new_value + fitted_params[nameotherID] < 0.99 * fixed_parameters$generation_time)
            break
        }
      } else {
        repeat {
          new_value <- old_value * exp(random_walk_rate[parID] * rnorm(1))   # change old value slightly
          if(new_value > limits[1, parID] & new_value < limits[2, parID])
            break
        }
      }
      fitted_params[parID] <- new_value

      # compute log-likelihood using new value
      params <- c(fitted_params, looped_params, fixed_params)
      sol <- model_gen(user = as.list(params))$run(times)
      new_llike <- compute_loglikelihood(sol, data, params)

      # decide whether to accept new value or not
      Nu <- new_llike - old_llike
      if (log(runif(1)) < Nu) {
        # accept: keep new values
        old_llike <- new_llike
        nr_accepted[parID] <- nr_accepted[parID] + 1
      } else {
        # reject: revert to old values
        fitted_params[parID] <- old_value
      }

      # adaptive MCMC: tuning sample proposal variance
      Nu0 <- 0.234    #ideal acceptance probability
      if (iter < mcmc_iterations / (sample_spacing * 2)) {
        temp <- random_walk_rate[parID] * exp(0.4 * (exp(min(Nu, 0)) - Nu0) / (35 * iter / mcmc_iterations + 1)) #tuning random_walk_rate
        random_walk_rate[parID] <- ifelse(temp > 1e-10, ifelse(temp < 10, temp, 10), 1e-10) #bounded RW between 0 - 10
      }

      # store subsample of MCMC chain
      if (iter%%sample_spacing == 0 & parID == nr_fitted) {
        print(iter)
        temp <- iter%/%sample_spacing
        stored_fit[temp, ] <- fitted_params
        stored_llike[temp] <- old_llike
        stored_rw[temp, ]  <- random_walk_rate
      }

    }
  }

  # store output in a list
  acceptance_rate <- nr_accepted / mcmc_iterations
  names(acceptance_rate) <- colnames(stored_fit)
  results <- list(acceptance_rate   = acceptance_rate,
                  random_walk_rate  = stored_rw,
                  log_likelihood    = stored_llike,
                  fitted_parameters = stored_fit)

  return(results)

}

model_gen <- odin::odin({

  # define (time-)dependent parameters
  gamma <- 1 / (generation_time - inv_nu - inv_delta)
  nu    <- 1 / inv_nu
  delta <- 1 / inv_delta
  beta  <- if (t < tQ) R0_1 * gamma else w * R0_1 * gamma

  # ODE system
  deriv(S)    <- - beta * (q * TPp + I_A + I_S) * S/N
  deriv(E)    <- beta * (q * TPp + I_A + I_S) * S/N - nu * E
  deriv(TPp)  <- nu * E - delta * TPp
  deriv(I_A)  <- p * delta * TPp - gamma * I_A
  deriv(I_S)  <- (1 - p) * delta * TPp - gamma * I_S
  deriv(TP_A) <- gamma * I_A - sigma * TP_A
  deriv(TP_S) <- gamma * I_S - sigma * TP_S
  deriv(TN)   <- sigma * (TP_A + TP_S)

  # initial conditions
  initial(S)    <- N - seed
  initial(E)    <- 0
  initial(TPp)  <- 0
  initial(I_A)  <- p * seed
  initial(I_S)  <- (1 - p) * seed
  initial(TP_A) <- 0
  initial(TP_S) <- 0
  initial(TN)   <- 0

  # input parameters
  seed      <- user()
  p         <- user()
  w         <- user()
  inv_nu    <- user()
  inv_delta <- user()
  R0_1      <- user()
  sigma     <- user()
  tSeed     <- user()
  time1     <- user()
  time2     <- user()
  tQ        <- user()
  N         <- user()
  generation_time <- user()
  q         <- user()

})

# log-likelihood formula
ll <- function (prop, numerator, denominator) {
  lgamma(denominator + 1) - lgamma(numerator + 1) - lgamma(denominator - numerator + 1)+
    numerator * log(prop) + (denominator - numerator) * log(1-prop)
}

# compute log-likelihood
compute_loglikelihood <- function(sol, data, params) {

  time1 <- params["time1"]
  time2 <- params["time2"]
  p     <- params["p"]
  N     <- params["N"]

  p_presymp1 <- (1 - p) * sol[time1, "TPp"] / N
  p_presymp2 <- (1 - p) * sol[time2, "TPp"] / N
  p_symp1    <- (sol[time1, "I_S"] + sol[time1, "TP_S"]) / N
  p_symp2    <- (sol[time2, "I_S"] + sol[time2, "TP_S"]) / N
  p_asymp1   <- (p * sol[time1, "TPp"] + sol[time1, "I_A"] + sol[time1, "TP_A"]) / N
  p_asymp2   <- (p * sol[time2, "TPp"] + sol[time2, "I_A"] + sol[time2, "TP_A"]) / N

  ll_presymp1 <- ll(p_presymp1, data$Presymptomatic[1], data$Tested[1])
  ll_presymp2 <- ll(p_presymp2, data$Presymptomatic[2], data$Tested[2])
  ll_symp1    <- ll(p_symp1,    data$Symptomatic[1],    data$Tested[1])
  ll_symp2    <- ll(p_symp2,    data$Symptomatic[2],    data$Tested[2])
  ll_asymp1   <- ll(p_asymp1,   data$Asymptomatic[1],   data$Tested[1])
  ll_asymp2   <- ll(p_asymp2,   data$Asymptomatic[2],   data$Tested[2])

  res <- ll_presymp1 + ll_presymp2 + ll_symp1 + ll_symp2 + ll_asymp1 + ll_asymp2

  return(unname(res))
}
