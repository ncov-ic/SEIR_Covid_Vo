central_rep_num <- function(linelist = linelist_villager, contacts = village_contacts, 
                            onset_delay = onset_to_confirmation, 
                            iterations = 1000, 
                            pre_shape = pre_parm$shape, 
                            pre_scale = pre_parm$scale,
                            post_shape = post_parm$shape, 
                            post_scale = post_parm$scale) {
  empty <- rep(0, 1 * iterations)
  empty <- matrix(empty, ncol = 1)
  empty <- as.data.frame(empty)
  
  names(empty) <- c("R")
  
  cohort_1 <- empty
  cohort_2 <- empty
  
  for(i in 1:iterations) {
    
    set.seed(i)
    
    ## impute onset dates
    linelist_known <- linelist %>% dplyr::filter(!is.na(onset) == TRUE)
    linelist_unknown <- linelist %>% dplyr::filter(!is.na(onset) == FALSE) 
    
    linelist_unknown %<>% mutate(onset = first_positive - sample(onset_delay,
                                                                 size = n(),
                                                                 replace = TRUE))
    linelist_final <- rbind(linelist_known, linelist_unknown)
    linelist_final$id %<>% as.character()
    
    first_inf <- min(linelist_final$onset)
    index_case <- linelist$id[linelist_final$onset == first_inf]
    
    ## need to match infectors and infectees and add in their onset dats
    
    df <- contacts
    
    ## suppress messages so that it doesn't crash your R session
    df <- dplyr::left_join(df, linelist_final[,c("id", "onset")], by = c("from" = "id"))
    names(df)[3] <- "onset_from"
    
    df <- dplyr::left_join(df, linelist_final[,c("id", "onset")], by = c("to" = "id"))
    names(df)[4] <- "onset_to"
    
    ## correct the directions of transmission where necessary
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
                                    "correction")]
    
    names(df_incorrect) <- names(df_correct)
    
    df <- rbind(df_correct, df_incorrect)
    
    ## ensure infectors are unique
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
    
    known_infector <- c(unique(df_out$to), index_case)
    
    missing_infector <- linelist_final %>%
      dplyr::filter((id %in% known_infector) == FALSE) %>%
      dplyr::mutate(infector = ".", exposure = as.Date(NA))
    
    ### for each member of the missing_infector, we want to impute their infector 
    ### using the serial interval
    for(j in 1:nrow(missing_infector)) {
      ## sample the date of exposure using the gamma distribution
      
      ## ensure that you use the correct distribution
      if(missing_infector$onset[j] >= as.Date("2020-02-24")) {
        missing_infector$exposure[j] <- as.Date(missing_infector$onset[j]) - 
          as.difftime(rgamma(n = 1, shape = post_shape,
                             scale = post_scale), units = "days")
      } else {
        missing_infector$exposure[j] <- as.Date(missing_infector$onset[j]) - 
          as.difftime(rgamma(n = 1, shape = pre_shape,
                             scale = pre_scale), units = "days")
      }
      
      ## interval of dates where infector must have had symptom onset
      int <- seq.Date(missing_infector$exposure[j] - as.difftime(4, units = "days"), 
                      missing_infector$exposure[j], by = "days")
      
      infectors <- linelist_final %>%
        dplyr::filter(onset %in% int)
      
      
      infectors <- infectors$id
      
      while(length(infectors) == 0) {
        if(missing_infector$onset[j] >= as.Date("2020-02-24")) {
          missing_infector$exposure[j] <- as.Date(missing_infector$onset[j]) - 
            as.difftime(rgamma(n = 1, shape = post_shape,
                               scale = post_scale), units = "days")
        } else {
          missing_infector$exposure[j] <- as.Date(missing_infector$onset[j]) - 
            as.difftime(rgamma(n = 1, shape = pre_shape, 
                               scale = pre_scale), units = "days")
        }
        
        int <- seq.Date(missing_infector$exposure[j] - as.difftime(4, units = "days"), 
                        missing_infector$exposure[j], by = "days")
        
        infectors <- linelist_final %>%
          dplyr::filter(onset %in% int)
        
        infectors <- infectors$id
        
      }
      
      ## impute an infector for all individuals
      missing_infector$infector[j] <- sample(infectors, size = 1)
    }
    
    ## make a contacts df of the imputed infectors and append to the existing contacts df 
    imputed_contacts <- missing_infector[,c("id", "infector")]
    names(imputed_contacts) <- c("to", "from")
    
    df_missing <- dplyr::left_join(imputed_contacts, linelist_final[,c("id", "onset")], by = c("from" = "id"))
    names(df_missing)[3] <- "onset_from"
    
    df_missing <- dplyr::left_join(df_missing, linelist_final[,c("id", "onset")], by = c("to" = "id"))
    names(df_missing)[4] <- "onset_to"
    
    df_out <- df_out[,names(df_missing)]
    
    ## this is the full list of infectors and infectees, in the correct orders
    df <- rbind(df_out, df_missing)
    
    ## people infected by each infector
    df %<>% dplyr::group_by(from) %>%
      dplyr::summarise(secondary_cases = n())
    
    ## split all individuals into 3 cohorts - before 20th, 21st - 28th, 29th+
    onset_to_cohort <- function(date_of_onset) {
      if(date_of_onset <= as.Date("2020-02-20")) {
        cohort = "1"
      } else {
        cohort = "2"
      }
      cohort
    }
    
    linelist_final %<>% 
      dplyr::rowwise() %>%
      dplyr::mutate(cohort = onset_to_cohort(onset))
    
    data <- dplyr::left_join(linelist_final, df, by = c("id" = "from")) %>%
      dplyr::mutate(R = ifelse(is.na(secondary_cases) == TRUE, 0, secondary_cases))
    
    ## data for the 3 cohorts
    data_1 <- data %>% dplyr::filter(cohort == "1") 
    data_2 <- data %>% dplyr::filter(cohort == "2")
    
    ## for each cohort - we want to output the mean of R (R_eff), variance of R and the se
    ### mean R
    cohort_1[i,1] <- mean(data_1$R)
    cohort_2[i,1] <- mean(data_2$R)
    
  }
  
  cohort_1$cohort <- "1"
  cohort_2$cohort <- "2"
  r_estimates <- rbind(cohort_1, cohort_2)
  
}


confidence_rep_num <- function(linelist = linelist_villager, contacts = village_contacts, 
                               onset_delay = onset_to_confirmation, 
                               iterations = 1000, 
                               pre_shape = pre_parm$shape, 
                               pre_scale = pre_parm$scale,
                               post_shape = post_parm$shape, 
                               post_scale = post_parm$scale) {
  empty <- rep(0, 1 * iterations)
  empty <- matrix(empty, ncol = 1)
  empty <- as.data.frame(empty)
  
  names(empty) <- c("R")
  
  cohort_1 <- empty
  cohort_2 <- empty
  
  for(i in 1:iterations) {
    
    set.seed(i)
    
    ## impute onset dates
    linelist_known <- linelist %>% dplyr::filter(!is.na(onset) == TRUE)
    linelist_unknown <- linelist %>% dplyr::filter(!is.na(onset) == FALSE) 
    
    linelist_unknown %<>% mutate(onset = first_positive - sample(onset_delay,
                                                                 size = n(),
                                                                 replace = TRUE))
    linelist_final <- rbind(linelist_known, linelist_unknown)
    linelist_final$id %<>% as.character()
    
    first_inf <- min(linelist_final$onset)
    index_case <- linelist$id[linelist_final$onset == first_inf]
    
    ## need to match infectors and infectees and add in their onset dats
    
    df <- contacts
    
    ## suppress messages so that it doesn't crash your R session
    df <- dplyr::left_join(df, linelist_final[,c("id", "onset")], by = c("from" = "id"))
    names(df)[3] <- "onset_from"
    
    df <- dplyr::left_join(df, linelist_final[,c("id", "onset")], by = c("to" = "id"))
    names(df)[4] <- "onset_to"
    
    ## correct the directions of transmission where necessary
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
                                    "correction")]
    
    names(df_incorrect) <- names(df_correct)
    
    df <- rbind(df_correct, df_incorrect)
    
    ## ensure infectors are unique
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
    
    known_infector <- c(unique(df_out$to), index_case)
    
    missing_infector <- linelist_final %>%
      dplyr::filter((id %in% known_infector) == FALSE) %>%
      dplyr::mutate(infector = ".", exposure = as.Date(NA))
    
    ### for each member of the missing_infector, we want to impute their infector using the serial interval
    for(j in 1:nrow(missing_infector)) {
      ## sample the date of exposure using the gamma distribution
      
      ## ensure that you use the correct distribution
      if(missing_infector$onset[j] >= as.Date("2020-02-24")) {
        missing_infector$exposure[j] <- as.Date(missing_infector$onset[j]) - 
          as.difftime(rgamma(n = 1, shape = post_shape,
                             scale = post_scale), units = "days")
      } else {
        missing_infector$exposure[j] <- as.Date(missing_infector$onset[j]) - 
          as.difftime(rgamma(n = 1, shape = pre_shape,
                             scale = pre_scale), units = "days")
      }
      
      int <- seq.Date(missing_infector$exposure[j] - as.difftime(4, units = "days"), 
                      missing_infector$exposure[j], by = "days")
      
      infectors <- linelist_final %>%
        dplyr::filter(onset %in% int)
      
      
      infectors <- infectors$id
      
      while(length(infectors) == 0) {
        if(missing_infector$onset[j] >= as.Date("2020-02-24")) {
          missing_infector$exposure[j] <- as.Date(missing_infector$onset[j]) - 
            as.difftime(rgamma(n = 1, shape = post_shape,
                               scale = post_scale), units = "days")
        } else {
          missing_infector$exposure[j] <- as.Date(missing_infector$onset[j]) - 
            as.difftime(rgamma(n = 1, shape = pre_shape, 
                               scale = pre_scale), units = "days")
        }
        
        int <- seq.Date(missing_infector$exposure[j] - as.difftime(4, units = "days"), 
                        missing_infector$exposure[j], by = "days")
        
        infectors <- linelist_final %>%
          dplyr::filter(onset %in% int)
        
        infectors <- infectors$id
        
      }
      
      ## impute an infector for all individuals
      missing_infector$infector[j] <- sample(infectors, size = 1)
    }
    
    ## make a contacts df of the imputed infectors and append to the existing contacts df 
    imputed_contacts <- missing_infector[,c("id", "infector")]
    names(imputed_contacts) <- c("to", "from")
    
    df_missing <- dplyr::left_join(imputed_contacts, linelist_final[,c("id", "onset")], by = c("from" = "id"))
    names(df_missing)[3] <- "onset_from"
    
    df_missing <- dplyr::left_join(df_missing, linelist_final[,c("id", "onset")], by = c("to" = "id"))
    names(df_missing)[4] <- "onset_to"
    
    df_out <- df_out[,names(df_missing)]
    
    ## this is the full list of infectors and infectees, in the correct orders
    df <- rbind(df_out, df_missing)
    
    ## people infected by each infector
    df %<>% dplyr::group_by(from) %>%
      dplyr::summarise(secondary_cases = n())
    
    ## split all individuals into 3 cohorts - before 20th, 21st - 28th, 29th+
    onset_to_cohort <- function(date_of_onset) {
      if(date_of_onset <= as.Date("2020-02-20")) {
        cohort = "1"
      } else {
        cohort = "2"
      }
      cohort
    }
    
    linelist_final %<>% 
      dplyr::rowwise() %>%
      dplyr::mutate(cohort = onset_to_cohort(onset))
    
    data <- dplyr::left_join(linelist_final, df, by = c("id" = "from")) %>%
      dplyr::mutate(R = ifelse(is.na(secondary_cases) == TRUE, 0, secondary_cases))
    
    ## data for the 3 cohorts
    data_1 <- data %>% dplyr::filter(cohort == "1") 
    data_2 <- data %>% dplyr::filter(cohort == "2")
    
    ## bootstrapping by sampling infectors with replacement
    data_1 <- data_1[sample(nrow(data_1), nrow(data_1), replace = TRUE), ]
    data_2 <- data_2[sample(nrow(data_2), nrow(data_2), replace = TRUE), ]
    
    ## for each cohort - we want to output the mean of R (R_eff), variance of R and the se
    ### mean R
    cohort_1[i,1] <- mean(data_1$R)
    cohort_2[i,1] <- mean(data_2$R)
  }
  
  cohort_1$cohort <- "1"
  cohort_2$cohort <- "2"
  r_estimates <- rbind(cohort_1, cohort_2)
  
}
