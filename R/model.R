# Run the SEIR model specified in the paper
wrapper_model <- function (iter_main, solve_seir, data,
                           fixed_parameters, all_looped_parameters,
                           limits, random_walk_rate,
                           mcmc_iterations, sample_spacing, dir_output) {

  # extract parameter combination for this iteration of the loop
  looped_parameters <- all_looped_parameters[iter_main, ]

  # initialise fitted parameters: make sure they're within MCMC limits
  fitted_parameters <- sapply(limits, mean)

  # timesteps to solve the SEIR model
  times <- seq(from = fixed_parameters$tSeed, to = fixed_parameters$time2, by = 1)

  # run mcmc
  set.seed(1)
  results <- run_mcmc(solve_seir        = solve_seir,
                      data              = data,
                      times             = times,
                      fixed_parameters  = fixed_parameters,
                      looped_parameters = looped_parameters,
                      fitted_parameters = fitted_parameters,
                      limits            = limits,
                      random_walk_rate  = random_walk_rate,
                      mcmc_iterations   = mcmc_iterations,
                      sample_spacing    = sample_spacing)

  # write to file in R format
  saveRDS(results, file.path(dir_output, paste0("results_", sprintf("%02d", iter_main), ".rds")))

}

# define mcmc routine
run_mcmc <- function (solve_seir, data, times,
                      fixed_parameters, looped_parameters, fitted_parameters,
                      limits, random_walk_rate,
                      mcmc_iterations, sample_spacing) {

  # initialise model parameter vectors
  fitted_params <- unlist(fitted_parameters)
  looped_params <- unlist(looped_parameters)
  fixed_params  <- unlist(fixed_parameters)

  # compute "old" MCMC values
  params <- c(fitted_params, looped_params, fixed_params)
  sol <- solve_seir(params, times)
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
      old_value <- fitted_params[parID]
      repeat {
        new_value <- old_value * exp(random_walk_rate[parID] * rnorm(1))   # change old value slightly
        if(new_value > limits[1, parID] & new_value < limits[2, parID])
          break
      }

      # compute log-likelihood using new value
      fitted_params[parID] <- new_value
      params <- c(fitted_params, looped_params, fixed_params)
      sol <- solve_seir(params, times)
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
                  fitted_parameters = stored_fit,
                  fixed_parameters  = fixed_params,
                  looped_parameters = looped_params)

  return(results)

}

# define model SEIR 1: S, E, I_A, I_S, TP_A, TP_S, TN
seir_model <- function(time, state, parms) {
  with(as.list(c(state, parms)), {

    if (time < tQ)
      w <- 1

    beta <- w * R0_1 * gamma

    # ODE system
    dS    <- - beta * (I_A + I_S) * S/N
    dE    <- beta * (I_A + I_S) * S/N - delta * E
    dI_A  <- p * delta * E - gamma * I_A
    dI_S  <- (1 - p) * delta * E - gamma * I_S
    dTP_A <- gamma * I_A - sigma * TP_A
    dTP_S <- gamma * I_S - sigma * TP_S
    dTN   <- sigma * (TP_A + TP_S)

    return(list(c(dS, dE, dI_A, dI_S, dTP_A, dTP_S, dTN)))

  })
}

# solve model SEIR 1
solve_seir <- function(parms, times){
  with(as.list(c(parms)), {

    # set initial conditions
    initial_values <- c(S    = N - seed,
                        E    = 0,
                        I_A  = p * seed,
                        I_S  = (1 - p) * seed,
                        TP_A = 0,
                        TP_S = 0,
                        TN   = 0)

    out <- as.data.frame(ode(y     = initial_values,
                             times = times,
                             func  = seir_model,
                             parms = parms))
    return(out)

  })
}

# compute log-likelihood
compute_loglikelihood <- function(sol, data, params) {

  time1 <- params["time1"]
  time2 <- params["time2"]
  p     <- params["p"]
  N     <- params["N"]

  p_asymp1 <- (sol$I_A[time1] + sol$TP_A[time1]) / N
  p_symp1  <- (sol$I_S[time1] + sol$TP_S[time1]) / N
  p_asymp2 <- (sol$I_A[time2] + sol$TP_A[time2]) / N
  p_symp2  <- (sol$I_S[time2] + sol$TP_S[time2]) / N

  if (0 < p & p < 1) { # model includes both asymptomatics and symptomatics
    ll_asymp1 <- ll(p_asymp1, data$Asymptomatic[1], data$Tested[1])
    ll_symp1  <- ll(p_symp1,  data$Symptomatic[1],  data$Tested[1])
    ll_asymp2 <- ll(p_asymp2, data$Asymptomatic[2], data$Tested[2])
    ll_symp2  <- ll(p_symp2,  data$Symptomatic[2],  data$Tested[2])
    res <- ll_asymp1 + ll_symp1 + ll_asymp2 + ll_symp2
  } else { # model only includes one out of asymptomatics and symptomatics
    ll_1  <- ll(p_asymp1 + p_symp1, data$Symptomatic[1] + data$Asymptomatic[1],  data$Tested[1])
    ll_2  <- ll(p_asymp2 + p_symp2, data$Symptomatic[2] + data$Asymptomatic[2],  data$Tested[2])
    res <- ll_1 + ll_2
  }

  return(res)
}

# log-likelihood formula for binomial distribution
ll <- function (prop, numerator, denominator) {
  lgamma(denominator+1) - lgamma(numerator+1) - lgamma(denominator-numerator+1)+
    numerator * log(prop) + (denominator - numerator) * log(1-prop)
}
