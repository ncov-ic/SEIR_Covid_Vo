get_SEIR <- function (nr_sample, times, params, looped_parameters, iter) {
  res <- lapply(seq_len(nr_sample), function(i) {
    sol <- suppressWarnings(model_gen(user = as.list(params[i, ]))$run(times))
    sol <- data.frame(sol, p = params$p[i], id_sample = i)
  })
  res <- do.call(bind_rows, res)
  res <- cbind(res, as.list(looped_parameters))
  
  # check disease has died out
  disease_ongoing <- res %>%
    select(S, t, id_sample) %>%
    filter(t %in% c(max(t), max(t) - 1)) %>%
    group_by(id_sample) %>%
    summarise(diff_S = -diff(S)) %>%
    ungroup() %>%
    summarise(upper_quantile = unname(quantile(diff_S, 0.975)) > 0.5) %>%
    pull()
  
  if (disease_ongoing)
    stop("Disease is ongoing", iter)
  
  return(res)
}

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
    
    # read in -----------------------------------------------------------------#
    results <- readRDS(files[iter])
    fixed_parameters  <- results$fixed_parameters
    looped_parameters <- results$looped_parameters
    posterior         <- data.frame(results$posterior)

    # save all posteriors -----------------------------------------------------#
    # remove burnin
    posterior <- posterior[posterior[, "iter_mcmc"] > nr_burnin, ]
    
    posterior_all_chains <- list(posterior         = posterior,
                                 fixed_parameters  = fixed_parameters,
                                 looped_parameters = looped_parameters)
    saveRDS(posterior_all_chains, file.path(dir_output, paste0(
      "posterior_all_chains_", sprintf("%02d", iter), ".rds")))
    
    # save acceptance rate ----------------------------------------------------#
    acceptance_rate[[iter]] <- cbind(results$acceptance_rate,
                                     data.frame(as.list(looped_parameters)))
    
    # save random walk rate ---------------------------------------------------#
    random_walk_rate[[iter]] <- melt(results$random_walk_rate)
    random_walk_rate[[iter]] <- cbind(random_walk_rate[[iter]][, -1],
                                      as.list(looped_parameters))
    
    # select best chain -------------------------------------------------------#
    # best chain
    best_chain <- posterior %>%
      group_by(chain) %>%
      summarise(mean_ll = mean(ll)) %>%
      ungroup() %>%
      filter(mean_ll == max(mean_ll)) %>%
      select(chain) %>%
      pull()
    posterior <- posterior %>% filter(chain == best_chain) %>% select(-chain)
    
    # save to file
    posterior_best_chain <- list(posterior         = posterior,
                                 fixed_parameters  = fixed_parameters,
                                 looped_parameters = looped_parameters)
    saveRDS(posterior_best_chain, file.path(dir_output, paste0(
      "posterior_best_chain_", sprintf("%02d", iter), ".rds")))
    
    # compute SEIR models -----------------------------------------------------#
    # solve nr_sample consecutive SEIR models
    times <- seq(fixed_parameters["tSeed"], 200, 1)
    id_sample <- sample(seq_len(nrow(posterior)), nr_sample)
    
    params <- cbind(posterior[id_sample, ], as.list(looped_parameters), as.list(fixed_parameters))
    SEIR_with_lockdown <- get_SEIR(nr_sample, times, params, looped_parameters, iter)
    SEIR_with_lockdown$Lockdown <- "Yes"
    
    params$w <- 1
    SEIR_without_lockdown <- get_SEIR(nr_sample, times, params, looped_parameters, iter)
    SEIR_without_lockdown$Lockdown <- "No"
    
    SEIR[[iter]] <- rbind(SEIR_with_lockdown, SEIR_without_lockdown)
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
