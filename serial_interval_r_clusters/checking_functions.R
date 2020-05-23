# checking functions

#### ----------------------------------------------------------------------------------- ####
### check contact links are unique and offer a warning ###
check_unique_contact_links = function(df){
  #remove duplicates
  df_out = df[!duplicated(t(apply(df[c("from", "to")], 1, sort))),]
  
  if(nrow(df_out)<nrow(df)){
    warning(paste0("There were contact links defined twice."))
  }
}
#### ----------------------------------------------------------------------------------- ####
### check the dates are in the right order ###
check_date_order = function(linelist){
  #check whether onset reported before death
  linelist %<>%
    mutate(dates_in_correct_order = dmy(reported_onset_date) < dmy(death_date) )
  
  linelist %<>% mutate(dates_in_correct_order = ifelse(is.na(dates_in_correct_order), TRUE, dates_in_correct_order))
  
  if(sum(!linelist$dates_in_correct_order, na.rm = T)>0){
    stop(paste0("There is a date out of order, please check ", linelist$id[!linelist$dates_in_correct_order]))
  }

}

#### ----------------------------------------------------------------------------------- ####
### check contact links are in exposure windows ###
check_exposure_timeline = function(linelist, contacts){
  
  #add extra column to contacts to check if the link is feasible wrt exposure windows
  contacts %<>% mutate(INCONSISTENT = FALSE)
  
  #check each contact
  for(i in 1:nrow(contacts)){
    
    #if the linelist id is in the contacts (ie. that the contacts are specified correctly)
    if(sum(linelist$id %in% contacts$from[i])>0
       & sum(linelist$id %in% contacts$to[i])>0){
      
      #get the indices of each in the linelist
      linelist_index_from = which(linelist$id %in% contacts$from[i])
      linelist_index_to = which(linelist$id %in% contacts$to[i])
      
      #check if the exposure occurred before onset
      if(!is.na(linelist$reported_onset_date[linelist_index_from]) & 
         !is.na(linelist$reported_onset_date[linelist_index_to])){
        
        if(ymd(linelist$reported_onset_date[linelist_index_from]) >=
           ymd(linelist$reported_onset_date[linelist_index_to] )){
          
          contacts$INCONSISTENT[i] = TRUE
          
        } 
      }
      
     
    } 
  }
  
  return(contacts)
}
