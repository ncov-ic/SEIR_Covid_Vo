# compute CrI in print format
compute_CrI <- function(X) {
  paste0(round(mean(X), 2), " (",
         round(quantile(X, 0.025), 2), ", ",
         round(quantile(X, 0.975), 2), ")")
}

table_fitted <- function (dir_clean, data) {

  # initialise table
  table <- list()

  # loop over posterior files
  files <- list.files(path = dir_clean, pattern = "^posterior_best_chain_",
                      full.names = TRUE)
  for (iter in seq_along(files)) {

    # extract data
    dt <- readRDS(files[iter])
    posterior         <- dt$posterior
    looped_parameters <- dt$looped_parameters
    fixed_parameters  <- dt$fixed_parameters

    # compute DIC
    times <- seq(fixed_parameters$tSeed, fixed_parameters$time2, 1)
    params <- c(apply(posterior, 2, mean), unlist(looped_parameters), unlist(fixed_parameters))
    sol <- suppressWarnings(model_gen(user = as.list(params))$run(times))
    llike_of_mean <- compute_loglikelihood(sol, data, params)
    DIC <- - 4 * mean(posterior$ll) + 2 * unname(llike_of_mean)

    # compute dependent parameters
    posterior["inv_gamma"] <- with(as.list(c(fixed_parameters, posterior)),
                                   (generation_time - inv_nu - inv_delta))
    posterior["1-w"] <- 1 - posterior["w"]
    posterior["incubation_period"] <- posterior["inv_nu"] + posterior["inv_delta"]

    # compute 95% credible intervals
    CrI <- apply(posterior, 2, compute_CrI)

    # better presentation
    names(posterior) <- sub("^inv_", "1/", names(posterior))

    # combine
    table[[iter]] <- c(looped_parameters, CrI, DIC = round(DIC, 2))

  }

  # clean lists and save to file
  table <- do.call(bind_rows, table) %>%
    arrange(R0_1, `1/sigma`) %>%
    select(-sigma, -ll, -w)

  table

  saveRDS(table, file.path(dir_clean, "table.rds"))

}


table_final_size <- function (dir_clean) {

  # compute final size (%)
  dt <- readRDS(file.path(dir_clean, "SEIR.rds")) %>%
    filter(t == max(t)) %>%
    group_by(R0_1, sigma, Lockdown, id_sample) %>%
    summarise(final_size = 100 * (N - S) / N) %>%
    ungroup() %>%
    group_by(R0_1, sigma, Lockdown) %>%
    summarise(final_size = compute_CrI(final_size)) %>%
    ungroup() %>%
    mutate(Lockdown = paste(Lockdown, "lockdown")) %>%
    tidyr::spread(Lockdown, final_size) %>%
    mutate(sigma = 1/sigma) %>%
    rename(`1/sigma` = sigma) %>%
    arrange(R0_1, `1/sigma`)

  dt

  saveRDS(dt, file.path(dir_clean, "table_final_size.rds"))

}
