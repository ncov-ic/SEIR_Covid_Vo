# compute CrI in print format
compute_CrI <- function(X) {
  paste0(round(mean(X), 2), " (",
         round(quantile(X, 0.025), 2), ", ",
         round(quantile(X, 0.975), 2), ")")
}

table_fitted <- function (dir_output, dir_figures) {
  
  # initialise table
  table <- list()
  
  # loop over posterior files
  files <- list.files(path = dir_output, pattern = "^posterior_best_chain_",
                      full.names = TRUE)
  for (iter in seq_along(files)) {
    
    # extract data
    dt <- readRDS(files[iter])
    posterior         <- dt$posterior %>% select(-iter_mcmc)
    looped_parameters <- dt$looped_parameters
    fixed_parameters  <- dt$fixed_parameters
    
    # compute DIC
    times <- seq(fixed_parameters["tSeed"], fixed_parameters["time2"], 1)
    
    params <- c(apply(posterior, 2, mean), looped_parameters, fixed_parameters)
    if (sum(params[c("inv_nu", "inv_delta")]) > 0.99 + params["generation_time"]) {
      DIC.mean <- NA
    } else {
      sol <- suppressWarnings(model_gen(user = as.list(params))$run(times))
      llike_of_statistic <- compute_loglikelihood(sol, data, params)
      DIC.mean <- - 4 * mean(posterior$ll) + 2 * unname(llike_of_statistic)
    }
    
    my.mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    params <- c(apply(posterior, 2, my.mode), looped_parameters, fixed_parameters)
    if (sum(params[c("inv_nu", "inv_delta")]) > 0.99 * params["generation_time"]) {
      DIC.mode <- NA
    } else {
      sol <- suppressWarnings(model_gen(user = as.list(params))$run(times))
      llike_of_statistic <- compute_loglikelihood(sol, data, params)
      DIC.mode <- - 4 * mean(posterior$ll) + 2 * unname(llike_of_statistic)
    }
    
    params <- c(apply(posterior, 2, median), looped_parameters, fixed_parameters)
    if (sum(params[c("inv_nu", "inv_delta")]) > 0.99 + params["generation_time"]) {
      DIC.median <- NA
    } else {
      sol <- suppressWarnings(model_gen(user = as.list(params))$run(times))
      llike_of_statistic <- compute_loglikelihood(sol, data, params)
      DIC.median <- - 4 * mean(posterior$ll) + 2 * unname(llike_of_statistic)
    }
    
    # compute BIC (needs to compute number of fitted parameters)
    BIC <- (ncol(posterior) - 2) * log(nrow(posterior)) - 2 * max(posterior$ll)
    
    # compute dependent parameters
    posterior["inv_gamma"] <- with(as.list(c(fixed_parameters, posterior)),
                                   (generation_time - inv_nu - inv_delta))
    posterior["1-w"] <- 1 - posterior["w"]
    
    # compute 95% credible intervals
    CrI <- apply(posterior, 2, compute_CrI)
    
    # better presentation
    names(posterior) <- sub("^inv_", "1/", names(posterior))
    looped_parameters["1/sigma"] <- 1/looped_parameters["sigma"]
    looped_parameters <- looped_parameters[-which(names(looped_parameters) == "sigma")]
    
    # combine
    table[[iter]] <- c(looped_parameters, CrI,
                       DIC.mean   = round(DIC.mean, 2),
                       DIC.mode   = round(DIC.mode, 2),
                       DIC.median = round(DIC.median, 2),
                       BIC        = round(BIC, 2))
    
  }
  
  # clean lists and save to file
  table <- do.call(bind_rows, table)
  table[[1]] <- as.numeric(table[[1]])
  table[[2]] <- as.numeric(table[[2]])
  table <- arrange(table, R0_1, `1/sigma`)
  saveRDS(table, file.path(dir_output, "table.rds"))
  
}

table_TN <- function (dir_output, dir_figures, nr_sample) {
  
  # get compartment counts
  SEIR <- readRDS(file.path(dir_output, "SEIR.rds"))
  fixed_parameters <- SEIR$fixed_parameters
  dt <- SEIR$SEIR
  
  # compute final size
  dt <- dt %>%
    filter(t == max(t) & Lockdown == "Yes") %>%
    group_by(R0_1, sigma, id_sample) %>%
    summarise(final_size = fixed_parameters["N"] - S) %>%
    ungroup() %>%
    group_by(R0_1, sigma) %>%
    summarise(final_size = compute_CrI(final_size)) %>%
    ungroup() %>%
    mutate(sigma = 1/sigma) %>%
    rename(`1/sigma` = sigma) %>%
    arrange(R0_1, `1/sigma`)
  
  saveRDS(dt, file.path(dir_output, "table_TN.rds"))
  
}
