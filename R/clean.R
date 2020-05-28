get_paths <- function (dir = dir_output, pattern, x) {
  
  list.files(path = dir,
             pattern = paste0(x, "-", pattern),
             full.names = TRUE)
}

get_tables <- function (dir = dir_output, pattern, x) {
  
  files <- get_paths(dir = dir, pattern = pattern, x = x)
  
  df <- lapply(seq_along(files), function(i) {
    read.csv(files[i]) %>%
      mutate(iter  = row_number() - 1,
             chain = i)
  })
  
  return(do.call(bind_rows, df))
  
}

get_SEIR <- function (nr_sample, times,
                      Lockdown, iter,
                      post, looped_parameters, fixed_parameters) {
  
  res <- 
    lapply(seq_len(nr_sample), function(i) {
      params <- c(post[i, ], looped_parameters, fixed_parameters)
      sol <- suppressWarnings(model_gen(user = as.list(params))$run(times))
      sol <- data.frame(sol, p = params$p, id_sample = i)
    }) %>%
    do.call(bind_rows, .) %>%
    tibble::add_column(looped_parameters, fixed_parameters, iter, Lockdown)
  
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
    stop("Disease is ongoing ", iter)
  
  return(res)
}

# Clean the output of the model run
clean_posteriors <- function (dir_output, dir_clean, nr_burnin, nr_sample) {
  
  set.seed(1)
  
  # prepare for loop read
  SEIR             <- list()
  acceptance_rate  <- list()
  random_walk_rate <- list()
  
  iter_string <- 
    list.files(path = dir_output,
               pattern = "^posterior-iter_main",
               full.names = TRUE) %>%
    str_split("-", simplify = TRUE) %>%
    grep("iter_main", ., value = TRUE) %>%
    unique()
  
  # loop over results of separate jobs
  for (iter in seq_along(iter_string)) {
    
    i_string <- iter_string[iter]
    
    # read in -----------------------------------------------------------------#
    fixed_parameters  <- read.csv(get_paths(dir_output, paste0(i_string, ".csv"), "fixed_parameters"))
    looped_parameters <- read.csv(get_paths(dir_output, paste0(i_string, ".csv"), "looped_parameters")) %>%
      mutate(`1/sigma` = 1/sigma)
    posterior        <- get_tables(pattern = paste0(i_string, "-"), x = "posterior")
    acceptance_rate[[iter]] <- get_tables(pattern = paste0(i_string, "-"), x = "acceptance_rate") %>%
      tibble::add_column(looped_parameters)
    random_walk_rate[[iter]] <- get_tables(pattern = paste0(i_string, "-"), x = "random_walk_rate") %>%
      tibble::add_column(looped_parameters)
    
    
    # clean posterior ---------------------------------------------------------#
    
    # remove burnin
    posterior <- posterior %>%
      filter(iter > nr_burnin + 0.5)
    
    # save all chains
    posterior_all_chains <- list(posterior         = posterior,
                                 fixed_parameters  = fixed_parameters,
                                 looped_parameters = looped_parameters)
    saveRDS(posterior_all_chains, file.path(dir_clean, paste0(
      "posterior_all_chains_", iter, ".rds")))
    
    # find best chain and sample
    posterior <- posterior %>%
      group_by(chain) %>%
      mutate(mean_ll = mean(ll)) %>%
      ungroup() %>%
      filter(mean_ll == max(mean_ll)) %>%
      select(-iter, -chain, -mean_ll)
    
    # save best chain
    posterior_best_chain <- list(posterior         = posterior,
                                 fixed_parameters  = fixed_parameters,
                                 looped_parameters = looped_parameters)
    saveRDS(posterior_best_chain, file.path(dir_clean, paste0(
      "posterior_best_chain_", iter, ".rds")))
    
    # sample best posterior and compute SEIR models ---------------------------#
    
    times <- seq(fixed_parameters$tSeed, 240, 1)
    id_sample <- sample(seq_len(nrow(posterior)), nr_sample)
    
    post <- posterior[id_sample, ]
    SEIR_with_lockdown <- get_SEIR(nr_sample, times,
                                   Lockdown = "Yes", iter,
                                   post, looped_parameters, fixed_parameters)
    
    post$w <- 1
    SEIR_without_lockdown <- get_SEIR(nr_sample, times,
                                      Lockdown = "No", iter,
                                      post, looped_parameters, fixed_parameters)
    
    SEIR[[iter]] <- rbind(SEIR_with_lockdown, SEIR_without_lockdown)

  }
  
  # clean lists and save to file
  SEIR             <- do.call(bind_rows, SEIR)
  acceptance_rate  <- do.call(bind_rows, acceptance_rate)
  random_walk_rate <- do.call(bind_rows, random_walk_rate)
  
  saveRDS(SEIR,             file.path(dir_clean, "SEIR.rds"))
  saveRDS(acceptance_rate,  file.path(dir_clean, "acceptance_rate.rds"))
  saveRDS(random_walk_rate, file.path(dir_clean, "random_walk_rate.rds"))
  
}