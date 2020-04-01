# Create figures and a table of credible intervals for the fitted parameters

fig_chains <- function (dir_output, dir_figures) {

  # loop over posterior files
  files <- list.files(path = dir_output, pattern = "^posterior", full.names = TRUE)
  for (iter in seq_along(files)) {

    # extract data
    dt <- readRDS(files[iter])
    posterior         <- as.data.frame(dt$posterior)
    looped_parameters <- dt$looped_parameters

    # melt data
    posterior$step <- seq_along(posterior[, 1])
    posterior <- melt(posterior, id.vars = "step")

    # title
    title <- paste("R0_1 =",        looped_parameters["R0_1"],
                   ", 1/sigma =", 1/looped_parameters["sigma"])

    p <- ggplot(data = posterior, aes(x = step, y = value)) +
      geom_line() +
      facet_wrap(vars(variable), scales = "free_y", ncol = 2) +
      ggtitle(title)

    p

    ggsave(filename = file.path(dir_figures, paste0("mcmc_chains_", sprintf("%02d", iter), ".png")),
           plot = p, device = "png",
           width = 17, height = 20, units = "cm", dpi = 300, limitsize = TRUE)

  }

}

fig_acceptance_rates <- function (dir_output, dir_figures) {

  # read in
  dt <- readRDS(file.path(dir_output, "acceptance_rate.rds"))

  dt <- melt(dt, id.vars = c("R0_1", "sigma"))

  p <- ggplot(data = dt, aes(x = value)) +
    geom_histogram() +
    facet_wrap(vars(variable))

  p

  ggsave(filename = file.path(dir_figures, "acceptance_rate.png"),
         plot = p, device = "png",
         width = 17, height = 20, units = "cm", dpi = 300, limitsize = TRUE)

}

fig_SEIR <- function (dir_output, dir_figures, data) {

  # SEIR ----------------------------------------------------------------------#
  # get compartment counts
  dt <- readRDS(file.path(dir_output, "SEIR.rds"))
  fixed_parameters <- dt$fixed_parameters
  dt               <- dt$SEIR

  # compute prevalence by asymptomatic and symptomatic
  dt <- dt %>%
    mutate(Asymptomatic = I_A + TP_A,
           Symptomatic  = I_S + TP_S)
  if (sum(dt$Asymptomatic) < 0.01) {  # Asymptomatic compartment has not been implemented
    do.asymp <- FALSE
    do.symp  <- TRUE
    dt$Asymptomatic <- NULL
  } else if (sum(dt$Symptomatic) < 0.01) { # Symptomatic compartment has not been implemented
    do.symp  <- FALSE
    do.asymp <- TRUE
    dt$Symptomatic <- NULL
  } else {
    do.symp  <- TRUE
    do.asymp <- TRUE
  }
  dt <- dt %>%
    melt(id.vars = c("time", "R0_1", "sigma"),
         variable.name = "Compartment") %>%
    group_by(time, R0_1, sigma, Compartment) %>%
    mutate(value = 100 * value / fixed_parameters["N"]) %>%
    summarise(mean = mean(value),
              low  = quantile(value, probs = 0.025),
              high = quantile(value, probs = 0.975)) %>%
    ungroup()

  # parameters I tested different values for: R0_1, sigma
  # adjust columns for plotting
  dt$R0_1  <- paste("R0 =",      dt$R0_1)
  dt$sigma <- paste("1/sigma =", 1/dt$sigma)
  dt$sigma <- factor(dt$sigma, paste("1/sigma =", c(2,4,6,8,10,12)))

  # screening data ------------------------------------------------------------#
  # transform screening data for plot
  data$Time <- fixed_parameters[c("time1", "time2")]
  my.data <- data %>%
    melt(id.vars = c("Time", "Tested"), variable.name = "Compartment")
  if (!do.asymp) { # Asymptomatic compartment has not been implemented
    my.data$Compartment <- "Symptomatic"
  } else if (!do.symp) { # Symptomatic compartment has not been implemented
    my.data$Compartment <- "Asymptomatic"
  }
  my.data <- my.data %>%
    group_by(Time, Tested, Compartment) %>%
    summarise(value = sum(value)) %>%
    ungroup() %>%                           # don't remove this and following line
    group_by(Time, Tested, Compartment) %>% # else result won't be correct!
    mutate(lower = 100 * binom.test(value, Tested)$conf.int[1],
           mean  = 100 * binom.test(value, Tested)$estimate,
           upper = 100 * binom.test(value, Tested)$conf.int[2]) %>%
    ungroup()

  # Dates for x axis ----------------------------------------------------------#
  my.times <- unlist(fixed_parameters[c("tSeed", "tQ", "time2")])
  full.dates <- c(paste0(sprintf("%02d", 4:29), "/02"),
                  paste0(sprintf("%02d", 1:7), "/03"))
  my.dates <- character(length = length(full.dates))
  my.dates[my.times + 1] <- full.dates[my.times + 1]
  my.ticks <- which(my.dates != "")

  # Plot ----------------------------------------------------------------------#
  for (ylab in c("Prevalence", "Compartments")) {

    if (ylab == "Prevalence") {
      selected <- c("Symptomatic", "Asymptomatic")
    } else {
      selected <- c("S", "E", "I_A", "I_S", "TP_A", "TP_S", "TN")
    }

    p <- ggplot(data = filter(dt, Compartment %in% selected),
                aes(x = time, colour = Compartment, fill = Compartment)) +
      geom_line(aes(y = mean)) +
      geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.4, colour = NA) +
      facet_grid(rows = vars(R0_1), cols = vars(sigma)) +
      geom_vline(xintercept = as.numeric(my.times["tQ"]),
                 linetype = "dashed", size = .3) +
      labs(y = paste(ylab, "%")) +
      scale_x_continuous(name = "Date",
                         breaks = unique(dt$time),
                         labels = my.dates) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.ticks.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            legend.position = "bottom",
            legend.margin = margin(0, 0, 0, 0),
            legend.box.margin = margin(-8, 0, 0, 0))

    # colour palette and legend
    if ((!do.asymp | !do.symp) & ylab == "Prevalence") {
      p <- p +
        scale_colour_manual(values = "dodgerblue3") +
        scale_fill_manual(values = "dodgerblue3") +
        theme(legend.position = "none")
    } else {
      p <- p +
        scale_colour_brewer("", palette = "Dark2") +
        scale_fill_brewer("", palette = "Dark2")
    }

    # Screening data
    if (ylab == "Prevalence") {
      p <- p +
        geom_errorbar(data = my.data, aes(x = Time, ymin = lower, ymax = upper),
                      width = 1.2) +
        geom_point(data = my.data, aes(x = Time, y = mean))
    }

    p
    
    ggsave(filename = file.path(dir_figures, paste0(ylab, ".png")),
           plot = p, device = "png",
           width = 17, height = 8, units = "cm", dpi = 300, limitsize = TRUE)

  }

}

table_fitted <- function (dir_output) {

  # initialise table
  table <- list()

  # loop over posterior files
  files <- list.files(path = dir_output, pattern = "^posterior", full.names = TRUE)
  for (iter in seq_along(files)) {

    # extract data
    dt <- readRDS(files[iter])
    posterior         <- dt$posterior
    looped_parameters <- dt$looped_parameters
    fixed_parameters  <- dt$fixed_parameters

    # compute 95% credible intervals
    CrI <- apply(posterior, 2, function(X) {
      paste0(round(mean(X), 2), " (",
             round(quantile(X, 0.025), 2), ", ",
             round(quantile(X, 0.975), 2), ")")
    })

    # compute DIC
    params <- c(apply(posterior, 2, mean), looped_parameters, fixed_parameters)
    times <- seq(fixed_parameters["tSeed"], fixed_parameters["time2"], 1)
    sol <- solve_seir(params, times)
    llike_of_means <- compute_loglikelihood(sol, data, params)
    DIC <- - 4 * mean(posterior[, "log-likelihood"]) + 2 * unname(llike_of_means)

    # looped parameters
    looped_parameters["1/sigma"] <- 1/looped_parameters["sigma"]
    looped_parameters <- looped_parameters[-which(names(looped_parameters) == "sigma")]

    # combine
    table[[iter]] <- c(looped_parameters, CrI, DIC = round(DIC, 2))

  }

  # clean lists and save to file
  table <- do.call(bind_rows, table)
  table[[1]] <- as.numeric(table[[1]])
  table[[2]] <- as.numeric(table[[2]])
  table <- arrange(table, R0_1, `1/sigma`)
  saveRDS(table, file.path(dir_output, "table.rds"))

}

table_TN <- function (dir_output, nr_sample) {
  
  # initialise table
  table <- list()

  # loop over posterior files
  files <- list.files(path = dir_output, pattern = "^posterior", full.names = TRUE)
  for (iter in seq_along(files)) {
    
    # extract data
    dt <- readRDS(files[iter])
    posterior_fit     <- dt$posterior
    looped_parameters <- dt$looped_parameters
    fixed_parameters  <- dt$fixed_parameters
    
    # solve nr_sample consecutive SEIR models: wait for the disease to die out!
    sample_fit <- apply(posterior_fit, 2, sample, nr_sample)
    times <- seq(fixed_parameters["tSeed"], 200, 1)
    SEIR <- lapply(seq_len(nr_sample), function(i) {
      params <- c(sample_fit[i, ], looped_parameters, fixed_parameters)
      sol <- solve_seir(params, times)
      data.frame(sol[nrow(sol), ], id_sample = i) # only save the last row!
    })
    SEIR <- do.call(bind_rows, SEIR)
    SEIR <- SEIR %>% mutate(disease_compartments = E + I_A + I_S + TP_A + TP_S)
    if (quantile(SEIR$disease_compartments, 0.95) > 0.5)
      stop("Disease is ongoing", iter)

    # compute CrI of TN and add to table
    my.fun <- function(X) {
      paste0(round(mean(X), 2), " (",
             round(quantile(X, 0.025), 2), ", ",
             round(quantile(X, 0.975), 2), ")")
    }
    
    # looped parameters
    looped_parameters["1/sigma"] <- 1/looped_parameters["sigma"]
    looped_parameters <- looped_parameters[-which(names(looped_parameters) == "sigma")]
    
    table[[iter]] <- c(looped_parameters, TN = my.fun(SEIR$TN))
    
  }
  
  # clean lists and save to file
  table <- do.call(bind_rows, table)
  table[[1]] <- as.numeric(table[[1]])
  table[[2]] <- as.numeric(table[[2]])
  table <- arrange(table, R0_1, `1/sigma`)
  saveRDS(table, file.path(dir_output, "table_TN.rds"))
  
}
