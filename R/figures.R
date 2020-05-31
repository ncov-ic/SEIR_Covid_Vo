get_dates <- function (my.times) {

  # day 0 is the 4th of Fabruary 2020
  full.dates <- c(paste0(sprintf("%02d", 4:29), "/02"),
                  paste0(sprintf("%02d", 1:31), "/03"),
                  paste0(sprintf("%02d", 1:30), "/04"),
                  paste0(sprintf("%02d", 1:31), "/05"),
                  paste0(sprintf("%02d", 1:30), "/06"),
                  paste0(sprintf("%02d", 1:31), "/07"),
                  paste0(sprintf("%02d", 1:31), "/08"),
                  paste0(sprintf("%02d", 1:31), "/09"),
                  paste0(sprintf("%02d", 1:31), "/10"))
  my.dates <- full.dates[my.times + 1]

  return(my.dates)

}

fig_chains <- function (dir_clean, dir_figures) {

  # loop over posterior files
  files <- list.files(path = dir_clean, pattern = "^posterior_all_chains_",
                      full.names = TRUE)
  for (iter in seq_along(files)) {

    # extract data
    dt <- readRDS(files[iter])
    posterior         <- dt$posterior
    looped_parameters <- dt$looped_parameters

    # melt data
    posterior <- gather(posterior, "variable", "value", -iter, -chain) %>%
      mutate(chain    = as.factor(chain),
             variable = as.factor(variable) %>% relevel("ll"))

    # title
    title <- paste("R0_1 =",        looped_parameters["R0_1"],
                   ", 1/sigma =", 1/looped_parameters["sigma"])

    p <- ggplot(data = posterior, aes(x = iter, y = value,
                                      group = chain, colour = chain)) +
      geom_line() +
      facet_wrap(vars(variable), scales = "free_y", ncol = 2) +
      ggtitle(title)

    p

    ggsave(filename = file.path(dir_figures, paste0("mcmc_chains_", iter, ".png")),
           plot = p, device = "png",
           width = 17, height = 20, units = "cm", dpi = 300, limitsize = TRUE)

  }

}

fig_acceptance_rates <- function (dir_clean, dir_figures) {

  # read in
  dt <- readRDS(file.path(dir_clean, "acceptance_rate.rds"))

  # melt
  dt <- gather(dt, "variable", "value", -R0_1, -sigma, -chain)
  dt$chain <- as.factor(as.integer(dt$chain))

  p <- ggplot(data = dt, aes(x = value, fill = chain)) +
    geom_histogram(alpha = 0.5, position = "identity") +
    facet_wrap(vars(variable))

  p

  ggsave(filename = file.path(dir_figures, "acceptance_rate.png"),
         plot = p, device = "png",
         width = 17, height = 20, units = "cm", dpi = 300, limitsize = TRUE)

}

get_best_fit <- function (dir_clean) {

  # get best in terms of DIC
  readRDS(file.path(dir_clean, "table.rds")) %>%
    filter(R0_1 < 2.8 & `1/sigma` < 13) %>%
    mutate(sigma = 1/as.numeric(`1/sigma`)) %>%
    filter(DIC == min(DIC, na.rm = TRUE)) %>%
    select(R0_1, sigma)

}

fig_prevalence <- function (dir_clean, dir_figures, data, do) {

  # get compartment counts
  dt <- readRDS(file.path(dir_clean, "SEIR.rds"))
  fixed_parameters <- dt %>%
    select(time1, time2, tQ, tSeed, N) %>%
    unique() %>%
    unlist()

  # select best model in terms of log-likelihood
  if (do == "best_DIC") {
    dt <- dt %>% right_join(get_best_fit(dir_clean))
  } else if (do == "paper_estimate") {
    dt <- dt %>% filter(R0_1 == 2.4 & 1/sigma == 4)
  }

  # Compute prevalence by pre-symptomatic, symptomatic and asymptomatic
  dt <- dt %>%
    filter(Lockdown == "Yes" &
             t < fixed_parameters["time2"] + 0.5) %>%
    # summarise compartments of interest
    transmute(t, R0_1, sigma,
              Presymptomatic = (1 - p) * TPp,
              Symptomatic    = I_S + TP_S,
              Asymptomatic   = p * TPp + I_A + TP_A) %>%
    gather("Compartment", "value", -t, -R0_1, -sigma) %>%
    # compute mean and 95% CrI
    group_by(t, R0_1, sigma, Compartment) %>%
    mutate(value = 100 * value / fixed_parameters["N"]) %>%
    summarise(mean = mean(value),
              low  = quantile(value, probs = 0.025),
              high = quantile(value, probs = 0.975)) %>%
    ungroup() %>%
    # column names for plotting
    mutate(R0_1  = paste("R0 =", R0_1),
           sigma = paste("1/sigma =", round(1/sigma)) %>%
             factor(paste("1/sigma =", 1:30)),
           Compartment = sub("Presym", "Pre-sym", Compartment) %>%
             factor(c("Pre-symptomatic", "Symptomatic", "Asymptomatic")))

  # Screening data
  data$Time <- fixed_parameters[c("time1", "time2")]
  my.data <- data %>%
    gather("Compartment", "value", -Time, -Tested) %>%
    group_by(Time, Tested, Compartment) %>%
    summarise(value = sum(value)) %>%
    ungroup() %>%                           # don't remove this and following line
    group_by(Time, Tested, Compartment) %>% # else result won't be correct!
    mutate(mean  = 100 * binom.test(value, Tested)$estimate,
           lower = 100 * binom.test(value, Tested)$conf.int[1],
           upper = 100 * binom.test(value, Tested)$conf.int[2]) %>%
    ungroup() %>%
    # column names for plotting
    mutate(Compartment = sub("Presym", "Pre-sym", Compartment) %>%
             factor(c("Pre-symptomatic", "Symptomatic", "Asymptomatic")))

  # Dates for x axis
  my.times <- fixed_parameters[c("tSeed", "tQ", "time2")]
  my.dates <- get_dates(my.times)

  # Plot
  if (do != "all") {
    dotsize <- 1
    errorsize <- 1.2
    dashsize <- .3
  } else {
    dotsize <- 0.8
    errorsize <- 2
    dashsize <- .2
  }

  p <- ggplot(data = dt,
              aes(x = t, colour = Compartment, fill = Compartment)) +
    # model mean and 95% CrI
    geom_line(aes(y = mean)) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.4, colour = NA) +
    # screening data
    geom_errorbar(data = my.data, aes(x = Time, ymin = lower, ymax = upper),
                  width = errorsize) +
    geom_point(data = my.data, aes(x = Time, y = mean), size = dotsize) +
    # lockdown date
    geom_vline(xintercept = as.numeric(my.times["tQ"]),
               linetype = "dashed", size = dashsize) +
    # axes
    labs(y = "Prevalence (%)") +
    # colour palette and legend
    scale_colour_brewer("", palette = "Dark2") +
    scale_fill_brewer("", palette = "Dark2") +
    scale_x_continuous(name = "Date",
                       breaks = my.times,
                       minor_breaks = seq(min(my.times), max(my.times), 7),
                       labels = my.dates[which(my.dates != "")]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-8, 0, 0, 0),
          # Nature requirement:
          text = element_text(size = 7, family = "sans"))

  # faceting
  if (do != "all") {
    p <- p +
      theme(legend.position = "right",
            plot.tag = element_text(size = 8, face = "bold", family = "sans")) +
      labs(tag = "a")
  } else {
    p <- p +
      facet_grid(rows = vars(R0_1), cols = vars(sigma)) +
      theme(legend.position = "bottom")
  }

  p

  return(p)

}

fig_incidence <- function (dir_clean, do) {

  # get compartment counts
  dt <- readRDS(file.path(dir_clean, "SEIR.rds"))
  fixed_parameters <- dt %>%
    select(time1, time2, tQ, tSeed, N) %>%
    unique() %>%
    unlist()

  # select best model in terms of log-likelihood
  if (do == "best_DIC") {
    dt <- dt %>% right_join(get_best_fit(dir_clean))
  } else if (do == "paper_estimate") {
    dt <- dt %>% filter(R0_1 == 2.4 & 1/sigma == 4)
  }

  # clean data
  dt <- dt %>%
    arrange(t) %>%
    group_by(id_sample, Lockdown) %>%
    # summarise compartments of interest
    transmute(t, incid = -(S - lag(S, default = S[1]))) %>%
    ungroup() %>%
    # compute mean and 95% CrI
    group_by(Lockdown, t) %>%
    mutate(incid = 100 * incid / fixed_parameters["N"]) %>%
    summarise(mean = mean(incid),
              low  = quantile(incid, probs = 0.025),
              high = quantile(incid, probs = 0.975)) %>%
    ungroup()

  # find time of disease termination
  time_final <- dt %>%
    filter(mean > 0.01) %>%
    summarise(max(t)) %>%
    pull()
  dt <- dt %>%
    filter(t < time_final + 0.5)

  # Dates for x axis
  my.times <- c(fixed_parameters[c("tSeed", "tQ", "time2")], time_final)
  ## add a date in the last long interval
  my.times <- sort(c(my.times, as.integer(mean(rev(my.times)[1:2]))))
  my.dates <- get_dates(unname(my.times))

  # Plot
  p <- ggplot(data = dt,
              aes(x = t, colour = Lockdown, fill = Lockdown)) +
    # model mean and 95% CrI
    geom_line(aes(y = mean)) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.4, colour = NA) +
    # lockdown date
    geom_vline(xintercept = as.numeric(my.times["tQ"]),
               linetype = "dashed", size = .3) +
    # axes
    labs(y = "Incidence (%)") +
    # colour palette and legend
    scale_colour_brewer("", palette = "Set1") +
    scale_fill_brewer("", palette = "Set1") +
    scale_y_continuous(breaks = 0:5) +
    scale_x_continuous(name = "Date",
                       breaks = my.times,
                       minor_breaks = seq(min(my.times), max(my.times), 7),
                       labels = my.dates[which(my.dates != "")]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-8, 0, 0, 0),
          # Nature requirement:
          text = element_text(size = 7, family = "sans"),
          plot.tag = element_text(size = 8, face = "bold", family = "sans")) +
    labs(tag = "b")

  p

  return(p)

}

fig_final_size <- function (dir_clean, do) {

  # get compartment counts
  dt <- readRDS(file.path(dir_clean, "SEIR.rds"))
  fixed_parameters <- dt %>%
    select(time1, time2, tQ, tSeed, N) %>%
    unique() %>%
    unlist()
  
  # select best model in terms of log-likelihood
  if (do == "best_DIC") {
    dt <- dt %>% right_join(get_best_fit(dir_clean))
  } else if (do == "paper_estimate") {
    dt <- dt %>% filter(R0_1 == 2.4 & 1/sigma == 4)
  }

  # clean data
  dt <- dt %>%
    # only keep last time step
    filter(t == max(t)) %>%
    # summarise compartments of interest
    mutate(final_size = 100 * (fixed_parameters["N"] - S) / fixed_parameters["N"]) %>%
    # compute mean and 95% CrI
    group_by(Lockdown) %>%
    summarise(mean = mean(final_size),
              low  = quantile(final_size, probs = 0.025),
              high = quantile(final_size, probs = 0.975)) %>%
    ungroup()

  # Plot
  p <- ggplot(data = dt, aes(x = Lockdown, y = mean, fill = Lockdown)) +
    # model mean and 95% CrI
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.2) +
    # axes
    scale_y_continuous(name = "Epidemic final size (%)",
                       limits = c(0, 98),
                       breaks = seq(0L, 100L, by = 25L)) +
    # colour palette and legend
    scale_colour_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          # Nature requirement:
          text = element_text(size = 7, family = "sans"),
          plot.tag = element_text(size = 8, face = "bold", family = "sans")) +
    labs(tag = "c")

  p

  return(p)

}

fig_timeline <- function () {

  dt.dates <- data.table::fread(
    "date,       event
     21/02/2020, Detection of first death and first case
     07/03/2020, Second survey",
    sep = ",") %>%
    mutate(date = as.Date(date, format = "%d/%m"))

  dt.ranges <- data.table::fread(
    "start,      end,        event
     21/02/2020, 29/02/2020, First survey
     24/02/2020, 08/03/2020, Lockdown",
    sep = ",") %>%
    mutate(start = as.Date(start, format = "%d/%m"),
           end   = as.Date(end,   format = "%d/%m")) %>%
    mutate(mid   = start + as.integer(end - start)/2)

  all.dates <- c(dt.dates$date, dt.ranges$start, dt.ranges$end)

  p <- ggplot(data = dt.dates) +
    # labels
    geom_label_repel(data = dt.dates[1, ],
                     aes(x = date, label = event),
                     y = 0, ylim = c(0.35, NA), size = 3) +
    geom_label_repel(data = dt.ranges[1, ],
                     aes(x = mid, label = event, fill = as.factor(start)),
                     y = 0, ylim = c(NA, -0.3), size = 3) +
    geom_label_repel(data = dt.ranges[2, ],
                     aes(x = mid, label = event, fill = as.factor(start)),
                     y = 0, ylim = c(0.15, NA), size = 3) +
    geom_label_repel(data = dt.dates[2, ],
                     aes(x = date, label = event),
                     y = 0, ylim = c(NA, -0.3), size = 3) +
    # coloured segments
    geom_segment(data = dt.ranges[1, ],
                 aes(x = start, xend = end, colour = as.factor(start)),
                 y = -0.06, yend = -0.06, size = 5) +
    geom_segment(data = dt.ranges[2, ],
                 aes(x = start, xend = end, colour = as.factor(start)),
                 y = 0.06, yend = 0.06, size = 5) +
    # time line
    geom_segment(aes(x = min(all.dates)-1, xend = max(all.dates)+1,
                     y = 0, yend = 0)) +
    geom_point(aes(x = date),
               y = 0) +
    # axes
    ylim(-0.45, 0.5) +
    scale_x_date(date_labels = "%d/%m",
                 breaks = all.dates,
                 labels = all.dates) +
    theme(axis.title = element_blank(),
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          text = element_text(size = 10, family = "sans"),
          plot.tag = element_text(size = 10, face = "bold", family = "sans")) +
    labs(tag = "c")

  p

  return(p)

}
