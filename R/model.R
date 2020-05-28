# Run the SEIR model specified in the paper
wrapper_model <- function (iter_main, data,
                           fixed_parameters, looped_parameters,
                           lims, limits, random_walk_rate,
                           mcmc_iterations, sample_spacing, dir_output, id_chain) {

  # create output folder
  dir.create(dir_output, recursive = TRUE, showWarnings = FALSE)

  # prepare
  nr_iter_main <- nrow(looped_parameters)
  looped_parameters <- looped_parameters[iter_main, ]
  times <- seq(from = fixed_parameters$tSeed, to = fixed_parameters$time2, by = 1L)
  utils::write.table(
    looped_parameters,
    file.path(dir_output, paste0("looped_parameters-iter_main_", iter_main, ".csv")),
    append = FALSE, quote = FALSE,
    col.names = TRUE, row.names = FALSE, sep = ",")
  utils::write.table(
    fixed_parameters,
    file.path(dir_output, paste0("fixed_parameters-iter_main_", iter_main, ".csv")),
    append = FALSE, quote = FALSE,
    col.names = TRUE, row.names = FALSE, sep = ",")

  # loop over chains
  for (idc in id_chain) {

    set.seed(idc)
    print(paste0("## Loop number ", iter_main, "/", nr_iter_main,
                 ", chain number ", idc, "/", length(id_chain)))

    # initialise fitted parameters within MCMC limits,
    # but make sure there's space for the dependent variable 1/gamma
    repeat {
      fitted_parameters <- apply(lims, 2, function(lim) runif(1, lim[1], lim[2]))
      if (sum(fitted_parameters[c("inv_nu", "inv_delta")]) < 0.99 * fixed_parameters$generation_time)
        break
    }

    # run mcmc
    run_mcmc(dir_output        = dir_output,
             iter_main         = iter_main,
             idc               = idc,
             data              = data,
             times             = times,
             fixed_parameters  = fixed_parameters,
             looped_parameters = looped_parameters,
             fitted_parameters = fitted_parameters,
             limits            = limits,
             random_walk_rate  = random_walk_rate,
             mcmc_iterations   = mcmc_iterations,
             sample_spacing    = sample_spacing)

  }

}

# define mcmc routine
run_mcmc <- function (dir_output, iter_main, idc,
                      data, times,
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

  # print to file
  names(random_walk_rate) <- paste0("rwr_", names(random_walk_rate))
  utils::write.table(
    t(c(iter = 0, ll = old_llike, fitted_params)),
    file.path(dir_output, paste0("posterior-iter_main_", iter_main, "-chain_", idc, ".csv")),
    append = FALSE, quote = FALSE,
    col.names = TRUE, row.names = FALSE, sep = ",")
  utils::write.table(
    t(c(iter = 0, random_walk_rate)),
    file.path(dir_output, paste0("random_walk_rate-iter_main_", iter_main, "-chain_", idc, ".csv")),
    append = FALSE, quote = FALSE,
    col.names = TRUE, row.names = FALSE, sep = ",")

  # initialise MCMC parameters
  nr_fitted <- length(fitted_params)
  nr_stored <- mcmc_iterations%/%sample_spacing
  nr_accepted  <- rep(0, length.out = nr_fitted) # number of accepted iterations


  # MCMC ----------------------------------------------------------------------
  for (iter in 2:mcmc_iterations) {
    for (parID in 1:nr_fitted) {

      # sample new value
      nameID <- names(fitted_params[parID])
      old_value <- fitted_params[parID]
      if (nameID %in% c("inv_nu", "inv_delta")) {
        nameotherID <- setdiff(c("inv_nu", "inv_delta"), nameID)
        repeat {
          new_value <- old_value * exp(random_walk_rate[[parID]] * rnorm(1))   # change old value slightly
          if(new_value > limits[1, parID] &
             new_value + fitted_params[nameotherID] < 0.99 * fixed_parameters$generation_time)
            break
        }
      } else {
        repeat {
          new_value <- old_value * exp(random_walk_rate[[parID]] * rnorm(1))   # change old value slightly
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
        temp <- random_walk_rate[[parID]] * exp(0.4 * (exp(min(Nu, 0)) - Nu0) / (35 * iter / mcmc_iterations + 1)) #tuning random_walk_rate
        random_walk_rate[[parID]] <- ifelse(temp > 1e-10, ifelse(temp < 10, temp, 10), 1e-10) #bounded RW between 0 - 10
      }

    }

    # print subsample of MCMC chain
    if (iter%%sample_spacing == 0) {
      print(iter)

      utils::write.table(
        t(c(iter = iter, ll = old_llike, fitted_params)),
        file.path(dir_output, paste0("posterior-iter_main_", iter_main, "-chain_", idc, ".csv")),
        append = TRUE,
        col.names = FALSE, row.names = FALSE, sep = ",")

      utils::write.table(
        t(c(iter = iter, random_walk_rate)),
        file.path(dir_output, paste0("random_walk_rate-iter_main_", iter_main, "-chain_", idc, ".csv")),
        append = TRUE,
        col.names = FALSE, row.names = FALSE, sep = ",")
    }

  }

  # print acceptance rate
  acceptance_rate <- nr_accepted / mcmc_iterations
  names(acceptance_rate) <- names(fitted_params)
  utils::write.table(
    t(acceptance_rate),
    file.path(dir_output, paste0("acceptance_rate-iter_main_", iter_main, "-chain_", idc, ".csv")),
    append = FALSE, quote = FALSE,
    col.names = TRUE, row.names = FALSE, sep = ",")

  return(0)

}

model_gen <- odin::odin({

  # define dependent parameters

  gamma <- 1 / (generation_time - inv_nu - inv_delta)
  nu    <- 1 / inv_nu
  delta <- 1 / inv_delta

  # define time-dependent parameters
  beta  <- if (t < tQ) R0_1 * gamma else w * R0_1 * gamma

  # ODE system
  deriv(S)    <- - beta * (q_TPp * TPp + q_A * I_A + I_S) * S/N
  deriv(E)    <- beta * (q_TPp * TPp + q_A * I_A + I_S) * S/N - nu * E
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
  q_TPp     <- user()
  q_A       <- user()
  generation_time <- user()

})

# log-likelihood formula
ll <- function (prop, numerator, denominator) {
  lgamma(denominator + 1) - lgamma(numerator + 1) - lgamma(denominator - numerator + 1)+
    numerator * log(prop) + (denominator - numerator) * log(1 - prop)
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
