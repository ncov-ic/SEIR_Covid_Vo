central_serial_interval <- function(linelist, 
                                    contacts, 
                                    onset_delay, 
                                    iterations) {
  
  
  message("---------- Calculating serial interval ----------")
  
  
  ## make an empty df to fill with the estimates of parameters
  empty <- rep(0, 4 * iterations)
  empty <- matrix(empty, ncol = 4)
  empty <- as.data.frame(empty)
  
  estimate_df <- empty
  names(estimate_df) <- c("mean", "shape", "scale", "variance")
  
  for(i in 1:iterations) {
    ## ensure that sampling varies each run 
    set.seed(i)
    
    linelist$id %<>% as.character()
    
    ## impute onset dates
    linelist_known <- linelist %>% dplyr::filter(!is.na(onset) == TRUE)
    linelist_unknown <- linelist %>% dplyr::filter(!is.na(onset) == FALSE) 
    
    linelist_unknown %<>% mutate(onset = first_positive - sample(onset_delay,
                                                                 size = n(),
                                                                 replace = TRUE))
    linelist_final <- rbind(linelist_known, linelist_unknown)
    linelist_final$id %<>% as.character()
    
    ## need to match infectors and infectees and add in their onset dats
    df <- contacts
    
    ## suppress messages so that it doesn't crash your R session
    df <- dplyr::left_join(df, linelist_final[,c("id", "onset")], by = c("from" = "id"))
    names(df)[3] <- "onset_from"
    
    df <- dplyr::left_join(df, linelist_final[,c("id", "onset")], by = c("to" = "id"))
    names(df)[4] <- "onset_to"
    
    ## calculate the serial intervals based on imputed onsets -- correct for negatives
    df %<>% 
      dplyr::mutate(serial_interval = abs(as.numeric(difftime(onset_from, onset_to, 
                                                              units = "days"))))
    
    df %<>% dplyr::mutate(correction = ifelse(onset_to > onset_from, "correct", "incorrect"))
    
    df_correct <- df %>% dplyr::filter(correction == "correct")
    
    df_incorrect <- df %>% dplyr::filter(correction == "incorrect") %>%
      dplyr::mutate(correct_to = from, 
                    correct_from = to,
                    correct_onset_to = onset_from,
                    correct_onset_from = onset_to)
    
    ## make sure the orderings match the correct entries
    df_incorrect <- df_incorrect[,c("correct_to", "correct_from", 
                                    "correct_onset_from", "correct_onset_to",
                                    "serial_interval", "correction")]
    
    names(df_incorrect) <- names(df_correct)
    
    df <- rbind(df_correct, df_incorrect)
    
    get_a_row <- function(df){
      if(nrow(df) == 1){
        sample_ind = 1
      } else{
        sample_ind <- sample(1:nrow(df),1)
      }
      df <- df[sample_ind,]
      return(df)
    }
    
    df_out <- NULL
    
    for(t in unique(df$to)) {
      df_tmp <- df %>% filter(to == t)
      df_out %<>% bind_rows(get_a_row(df_tmp))
      df_out
    }
    
    # ## fit the gamma distributions
    parm <- epitrix::fit_disc_gamma(df_out$serial_interval, mu_ini = 7)
    
    ## fill it in the dataframe
    estimate_df[i,1] <- parm$mu
    estimate_df[i,2] <- parm$distribution$parameters$shape
    estimate_df[i,3] <- parm$distribution$parameters$scale
    estimate_df[i,4] <- parm$distribution$parameters$shape*(parm$distribution$parameters$scale)^2
    
  }
  
  estimate_df
  
}

confidence_serial_interval <- function(linelist, 
                                       contacts, 
                                       onset_delay, 
                                       iterations) {
  
  
  message("---------- Calculating CI of serial interval ----------")
  
  
  ## make an empty df to fill with the estimates of parameters
  empty <- rep(0, 4 * iterations)
  empty <- matrix(empty, ncol = 4)
  empty <- as.data.frame(empty)
  
  estimate_df <- empty
  names(estimate_df) <- c("mean", "shape", "scale", "variance")
  
  for(i in 1:iterations) {
    ## ensure that sampling varies each run 
    set.seed(i)
    
    ## prevents the messages from dplyr which crash R session with large iterations
    linelist$id %<>% as.character()
    
    ## randomly sample the infectees
    infectees <- sample(contacts$to, size = 41, replace = TRUE)
    
    ## want 3 data sets - contacts, infectees with their repetition and imputed onsets, 
    ## and the infectors and onsets
    
    infectee_df <- as.data.frame(matrix(rep(0,41), ncol = 1))
    infectee_df[,1] <- infectees
    names(infectee_df) <- "infectees"
    infectee_df <- dplyr::left_join(infectee_df, linelist, by = c("infectees" = "id"))
    
    ## impute onset dates as before
    
    infectee_known <- infectee_df %>%
      dplyr::filter(!is.na(onset) == TRUE)
    infectee_unknown <- infectee_df %>%
      dplyr::filter(!is.na(onset) == FALSE)
    
    infectee_unknown %<>% mutate(onset = first_positive - sample(onset_delay,
                                                                 size = n(),
                                                                 replace = TRUE))
    
    infectee_df <- rbind(infectee_known, infectee_unknown)
    
    ## now do the same for the infectors - want repeated infectors where necessary
    transmission_pairs <- dplyr::left_join(infectee_df, contacts, by = c("infectees" = "to")) 
    transmission_pairs <- transmission_pairs[,c("infectees", "from")]
    
    infector_df <- as.data.frame(matrix(rep(0,41), ncol = 1))
    infector_df[,1] <- transmission_pairs$from
    names(infector_df) <- "infectors"
    
    infector_df <- dplyr::left_join(infector_df, linelist, by = c("infectors" = "id"))
    
    infector_known <- infector_df %>%
      dplyr::filter(!is.na(onset) == TRUE)
    infector_unknown <- infector_df %>%
      dplyr::filter(!is.na(onset) == FALSE)
    
    infector_unknown %<>% mutate(onset = first_positive - sample(onset_delay,
                                                                 size = n(),
                                                                 replace = TRUE))
    
    infector_df <- rbind(infector_known, infector_unknown)
    
    ## join these together 
    
    final <- cbind(transmission_pairs, infectee_df$onset, infector_df$onset)
    names(final)[3:4] <- c("onset_to", "onset_from")
    
    final %<>% dplyr::mutate(serial_interval = abs(as.numeric(difftime(onset_from, onset_to, units = "days"))))
    
    df <- final
    names(df)[1] <- "to"
    
    df %<>% dplyr::mutate(correction = ifelse(onset_to > onset_from, "correct", "incorrect"))
    
    df_correct <- df %>% dplyr::filter(correction == "correct")
    
    df_incorrect <- df %>% dplyr::filter(correction == "incorrect") %>%
      dplyr::mutate(correct_to = from, 
                    correct_from = to,
                    correct_onset_to = onset_from,
                    correct_onset_from = onset_to)
    
    ## make sure the orderings match the correct entries
    df_incorrect <- df_incorrect[,c("correct_to", "correct_from", 
                                    "correct_onset_from", "correct_onset_to",
                                    "serial_interval", "correction")]
    
    names(df_incorrect) <- names(df_correct)
    
    df <- rbind(df_correct, df_incorrect)
    
    ## ensuring that there is only one infector per infectee
    get_a_row <- function(df){
      if(nrow(df) == 1){
        sample_ind = 1
      } else{
        sample_ind <- sample(1:nrow(df),1)
      }
      df <- df[sample_ind,]
      return(df)
    }
    
    df_out <- NULL
    
    for(t in unique(df$to)) {
      df_tmp <- df %>% filter(to == t)
      df_out %<>% bind_rows(get_a_row(df_tmp))
      df_out
    }
    
    # ## fit the gamma distributions
    parm <- epitrix::fit_disc_gamma(df_out$serial_interval, mu_ini = 7)
    
    ## fill it in the dataframe
    estimate_df[i,1] <- parm$mu
    estimate_df[i,2] <- parm$distribution$parameters$shape
    estimate_df[i,3] <- parm$distribution$parameters$scale
    estimate_df[i,4] <- parm$distribution$parameters$shape*(parm$distribution$parameters$scale)^2
    
  }
  
  estimate_df
  
}

assign_lockdown_status <- function(onset_infector, onset_infectee) {
  if(onset_infectee < as.Date("2020-02-24") & onset_infector < as.Date("2020-02-24")) {
    status <- "pre-lockdown"
  } else if(onset_infectee >= as.Date("2020-02-24") & onset_infector >= as.Date("2020-02-24")) {
    status <- "post-lockdown"
  } else {
    ## just in case the direction of transmission is wrong
    first_date <- min(onset_infectee, onset_infector)
    last_date <- max(onset_infectee, onset_infector)
    
    ## if more of the serial interval is pre-lockdown
    if(difftime(as.Date("2020-02-24"), first_date) >= difftime(last_date, as.Date("2020-02-24"))) {
      status <- "pre-lockdown"
    } else if (difftime(as.Date("2020-02-24"), first_date) < difftime(last_date, as.Date("2020-02-24"))) {
      status <- "post-lockdown"
    } else {
      stop("This is an edge-case -- all possible date combinations should be defined")
    }
  }
  status
}

## central estimates pre and post lockdown
central_sensitivity <- function(linelist, 
                                contacts, 
                                onset_delay, 
                                iterations) {
  
  
  message("---------- Calculating central estimates pre and post lockdown ----------")
  
  
  ## make an empty df to fill with the estimates of parameters
  empty <- rep(0, 4 * iterations)
  empty <- matrix(empty, ncol = 4)
  empty <- as.data.frame(empty)
  
  pre_estimate_df <- empty
  names(pre_estimate_df) <- c("mean", "shape", "scale", "variance")
  
  post_estimate_df <- empty
  names(post_estimate_df) <- c("mean", "shape", "scale", "variance")
  
  for(i in 1:iterations) {
    ## ensure that sampling varies each run 
    set.seed(i)
    
    ## impute onset dates
    linelist_known <- linelist %>% dplyr::filter(!is.na(onset) == TRUE)
    linelist_unknown <- linelist %>% dplyr::filter(!is.na(onset) == FALSE) 
    
    linelist_unknown %<>% mutate(onset = first_positive - sample(onset_delay,
                                                                 size = n(),
                                                                 replace = TRUE))
    linelist_final <- rbind(linelist_known, linelist_unknown)
    linelist_final$id %<>% as.character()
    
    ## need to match infectors and infectees and add in their onset dats
    
    df <- contacts
    
    ## suppress messages so that it doesn't crash your R session
    df <- suppressMessages(dplyr::left_join(df, linelist_final[,c("id", "onset")], by = c("from" = "id")))
    names(df)[3] <- "onset_from"
    
    df <- suppressMessages(dplyr::left_join(df, linelist_final[,c("id", "onset")], by = c("to" = "id")))
    names(df)[4] <- "onset_to"
    
    ## calculate the serial intervals based on imputed onsets -- correct for negatives
    df %<>% dplyr::mutate(serial_interval = abs(as.numeric(difftime(onset_from, onset_to, units = "days"))))
    
    df %<>% dplyr::mutate(correction = ifelse(onset_to > onset_from, "correct", "incorrect"))
    
    df_correct <- df %>% dplyr::filter(correction == "correct")
    
    df_incorrect <- df %>% dplyr::filter(correction == "incorrect") %>%
      dplyr::mutate(correct_to = from, 
                    correct_from = to,
                    correct_onset_to = onset_from,
                    correct_onset_from = onset_to)
    
    ## make sure the orderings match the correct entries
    df_incorrect <- df_incorrect[,c("correct_to", "correct_from", 
                                    "correct_onset_from", "correct_onset_to",
                                    "serial_interval", "correction")]
    
    names(df_incorrect) <- names(df_correct)
    
    df <- rbind(df_correct, df_incorrect)
    
    ## ensures there is only one infector per infectee
    get_a_row <- function(df){
      if(nrow(df) == 1){
        sample_ind = 1
      } else{
        sample_ind <- sample(1:nrow(df),1)
      }
      df <- df[sample_ind,]
      return(df)
    }
    
    df_out <- NULL
    
    for(t in unique(df$to)) {
      df_tmp <- df %>% filter(to == t)
      df_out %<>% bind_rows(get_a_row(df_tmp))
      df_out
    }
    
    df_out %<>% 
      rowwise() %>% 
      dplyr::mutate(status = assign_lockdown_status(onset_infector = onset_from,
                                                    onset_infectee = onset_to))
    
    pre_df <- dplyr::filter(df_out, status == "pre-lockdown")
    post_df <- dplyr::filter(df_out, status == "post-lockdown")
    
    # ## fit the gamma distributions
    parm_pre <- epitrix::fit_disc_gamma(pre_df$serial_interval, mu_ini = 7)
    parm_post <- epitrix::fit_disc_gamma(post_df$serial_interval, mu_ini = 7)    
    
    ## fill it in the dataframe
    pre_estimate_df[i,1] <- parm_pre$mu
    pre_estimate_df[i,2] <- parm_pre$distribution$parameters$shape
    pre_estimate_df[i,3] <- parm_pre$distribution$parameters$scale
    pre_estimate_df[i,4] <- parm_pre$distribution$parameters$shape*(parm_pre$distribution$parameters$scale)^2
    
    post_estimate_df[i,1] <- parm_post$mu
    post_estimate_df[i,2] <- parm_post$distribution$parameters$shape
    post_estimate_df[i,3] <- parm_post$distribution$parameters$scale
    post_estimate_df[i,4] <- parm_post$distribution$parameters$shape*(parm_post$distribution$parameters$scale)^2
    
  }
  
  pre_estimate_df$status <- "pre-lockdown"
  post_estimate_df$status <- "post-lockdown"
  
  estimate_df <- rbind(pre_estimate_df, post_estimate_df)
  estimate_df
  
}

## 95% confidence intervals on the central estimates for pre and post lockdown 
confidence_sensitivity <- function(linelist, 
                                   contacts, 
                                   onset_delay, 
                                   iterations) {
  
  
  message("---------- Calculating CI of central estimates pre and post lockdown ----------")
  
  
  ## make an empty df to fill with the estimates of parameters
  empty <- rep(0, 4 * iterations)
  empty <- matrix(empty, ncol = 4)
  empty <- as.data.frame(empty)
  
  pre_estimate_df <- empty
  names(pre_estimate_df) <- c("mean", "shape", "scale", "variance")
  
  post_estimate_df <- empty
  names(post_estimate_df) <- c("mean", "shape", "scale", "variance")
  
  for(i in 1:iterations) {
    ## ensure that sampling varies each run 
    set.seed(i)
    
    ## prevents the messages from dplyr which crash R session with large iterations
    linelist$id %<>% as.character()
    
    ## randomly sample the infectees
    infectees <- sample(contacts$to, size = 41, replace = TRUE)
    
    ## want 3 data sets - contacts, infectees with their repetition and imputed onsets, 
    ## and the infectors and onsets
    
    infectee_df <- as.data.frame(matrix(rep(0,41), ncol = 1))
    infectee_df[,1] <- infectees
    names(infectee_df) <- "infectees"
    infectee_df <- dplyr::left_join(infectee_df, linelist, by = c("infectees" = "id"))
    
    ## impute onset dates as before
    
    infectee_known <- infectee_df %>%
      dplyr::filter(!is.na(onset) == TRUE)
    infectee_unknown <- infectee_df %>%
      dplyr::filter(!is.na(onset) == FALSE)
    
    infectee_unknown %<>% mutate(onset = first_positive - sample(onset_delay,
                                                                 size = n(),
                                                                 replace = TRUE))
    
    infectee_df <- rbind(infectee_known, infectee_unknown)
    
    ## now do the same for the infectors - want repeated infectors where necessary
    transmission_pairs <- dplyr::left_join(infectee_df, contacts, by = c("infectees" = "to")) 
    transmission_pairs <- transmission_pairs[,c("infectees", "from")]
    
    infector_df <- as.data.frame(matrix(rep(0,41), ncol = 1))
    infector_df[,1] <- transmission_pairs$from
    names(infector_df) <- "infectors"
    
    infector_df <- dplyr::left_join(infector_df, linelist, by = c("infectors" = "id"))
    
    infector_known <- infector_df %>%
      dplyr::filter(!is.na(onset) == TRUE)
    infector_unknown <- infector_df %>%
      dplyr::filter(!is.na(onset) == FALSE)
    
    infector_unknown %<>% mutate(onset = first_positive - sample(onset_delay,
                                                                 size = n(),
                                                                 replace = TRUE))
    
    infector_df <- rbind(infector_known, infector_unknown)
    
    ## join these together 
    
    final <- cbind(transmission_pairs, infectee_df$onset, infector_df$onset)
    names(final)[3:4] <- c("onset_to", "onset_from")
    
    final %<>% dplyr::mutate(serial_interval = abs(as.numeric(difftime(onset_from, onset_to, units = "days"))))
    
    df <- final
    names(df)[1] <- "to"
    
    df %<>% dplyr::mutate(correction = ifelse(onset_to > onset_from, "correct", "incorrect"))
    
    df_correct <- df %>% dplyr::filter(correction == "correct")
    
    df_incorrect <- df %>% dplyr::filter(correction == "incorrect") %>%
      dplyr::mutate(correct_to = from, 
                    correct_from = to,
                    correct_onset_to = onset_from,
                    correct_onset_from = onset_to)
    
    ## make sure the orderings match the correct entries
    df_incorrect <- df_incorrect[,c("correct_to", "correct_from", 
                                    "correct_onset_from", "correct_onset_to",
                                    "serial_interval", "correction")]
    
    names(df_incorrect) <- names(df_correct)
    
    df <- rbind(df_correct, df_incorrect)
    
    get_a_row <- function(df){
      if(nrow(df) == 1){
        sample_ind = 1
      } else{
        sample_ind <- sample(1:nrow(df),1)
      }
      df <- df[sample_ind,]
      return(df)
    }
    
    df_out <- NULL
    
    for(t in unique(df$to)) {
      df_tmp <- df %>% filter(to == t)
      df_out %<>% bind_rows(get_a_row(df_tmp))
      df_out
    }
    
    df_out %<>% 
      rowwise() %>% 
      dplyr::mutate(status = assign_lockdown_status(onset_infector = onset_from,
                                                    onset_infectee = onset_to))
    
    pre_df <- dplyr::filter(df_out, status == "pre-lockdown")
    post_df <- dplyr::filter(df_out, status == "post-lockdown")
    # ## fit the gamma distributions
    parm_pre <- epitrix::fit_disc_gamma(pre_df$serial_interval, mu_ini = 7)
    parm_post <- epitrix::fit_disc_gamma(post_df$serial_interval, mu_ini = 7)    
    
    ## fill it in the dataframe
    pre_estimate_df[i,1] <- parm_pre$mu
    pre_estimate_df[i,2] <- parm_pre$distribution$parameters$shape
    pre_estimate_df[i,3] <- parm_pre$distribution$parameters$scale
    pre_estimate_df[i,4] <- parm_pre$distribution$parameters$shape*(parm_pre$distribution$parameters$scale)^2
    
    post_estimate_df[i,1] <- parm_post$mu
    post_estimate_df[i,2] <- parm_post$distribution$parameters$shape
    post_estimate_df[i,3] <- parm_post$distribution$parameters$scale
    post_estimate_df[i,4] <- parm_post$distribution$parameters$shape*(parm_post$distribution$parameters$scale)^2
    
  }
  
  pre_estimate_df$status <- "pre-lockdown"
  post_estimate_df$status <- "post-lockdown"
  
  estimate_df <- rbind(pre_estimate_df, post_estimate_df)
  estimate_df
  
}

