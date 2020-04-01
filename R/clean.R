# Clean the output of the model run
clean_mcmc_results <- function (dir_output, nr_burnin, nr_sample) {

  set.seed(1)

  # initialise lists
  SEIR             <- list()
  acceptance_rate  <- list()
  random_walk_rate <- list()

  # loop over MCMC output files
  files <- list.files(path = dir_output, pattern = "^results", full.names = TRUE)
  for (iter in seq_along(files)) {

    results <- readRDS(files[iter])
    fixed_parameters  <- results$fixed_parameters
    looped_parameters <- results$looped_parameters

    # compute posterior -------------------------------------------------------#
    # remove burnin
    posterior_log <- results$log_likelihood[-seq(1, nr_burnin, 1)]
    posterior_fit <- results$fitted_parameters[-seq(1, nr_burnin, 1), ]

    # clean and write to file
    posterior <- list(posterior = cbind("log-likelihood" = posterior_log,
                                        posterior_fit,
                                        "1-w" = 1 - posterior_fit[, "w"]),
                      fixed_parameters  = fixed_parameters,
                      looped_parameters = looped_parameters)
    saveRDS(posterior, file.path(dir_output, paste0("posterior_", sprintf("%02d", iter), ".rds")))

    # compute SEIR models -----------------------------------------------------#
    # solve nr_sample consecutive SEIR models
    sample_fit <- apply(posterior_fit, 2, sample, nr_sample)
    times <- seq(fixed_parameters["tSeed"], fixed_parameters["time2"], 1)
    SEIR[[iter]] <- lapply(seq_len(nr_sample), function(i) {
      params <- c(sample_fit[i, ], looped_parameters, fixed_parameters)
      sol <- solve_seir(params, times)
      sol <- data.frame(sol, id_sample = i)
    })
    SEIR[[iter]] <- do.call(bind_rows, SEIR[[iter]])
    SEIR[[iter]] <- cbind(SEIR[[iter]], as.list(looped_parameters))

    # save acceptance rate and random walk rate -------------------------------#
    acceptance_rate[[iter]] <- melt(results$acceptance_rate)
    acceptance_rate[[iter]]  <- c(results$acceptance_rate, results$looped_parameters)
    random_walk_rate[[iter]] <- melt(results$random_walk_rate)
    random_walk_rate[[iter]] <- cbind(random_walk_rate[[iter]][, -1],
                                      as.list(results$looped_parameters))

  }

  # clean lists and save to file
  SEIR             <- do.call(bind_rows, SEIR)
  SEIR             <- list(SEIR = SEIR, fixed_parameters = fixed_parameters)
  acceptance_rate  <- do.call(bind_rows, acceptance_rate)
  random_walk_rate <- do.call(bind_rows, random_walk_rate)

  saveRDS(SEIR,             file.path(dir_output, "SEIR.rds"))
  saveRDS(acceptance_rate,  file.path(dir_output, "acceptance_rate.rds"))
  saveRDS(random_walk_rate, file.path(dir_output, "random_walk_rate.rds"))

}
