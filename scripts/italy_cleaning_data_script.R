## Load the magrittr library to define the %>% and %<>% functions
library(magrittr)

set.seed((1))

## fix encoding of characters
full_data <- readxl::read_excel("data/anonymised_data_public_final.xlsx", sheet = 1,
                                guess_max = 2900)

head(full_data)
unique(full_data$"43882")

## find a way of doing this that doesn't require the specification of columns --> patient names
contact_data_long <- tidyr::gather(full_data, contact_id, contact_type, 
                                   "UQfjwKik":"db2cfcea_place", factor_key = TRUE) %>%
  dplyr::filter(!is.na(contact_type) == TRUE) 

contact_data_long <- contact_data_long[,c("id", "household_id", "contact_id", "contact_type")]
contact_ids <- unique(contact_data_long$contact_id)

missed_ids <- setdiff(contact_ids, full_data$id)

### check none of these are in fact individual entries in linelist
if(sum(missed_ids %in% full_data$id) != 0) {
  stop("Some of the additional ids are in the full data already")
}

missing_households <- data.frame(id = missed_ids, 
                                 household_id = c(rep("out_of_vo", length(missed_ids) - 3), 
                                                  rep("location_in_vo", 3)))
household_df <- full_data[,c("id", "household_id")]
household_complete <- rbind(household_df, missing_households)
names(household_complete) <- c("contact_id", "contact_household")

contact_data_long$contact_id <- as.character(contact_data_long$contact_id)

contact_data <- dplyr::full_join(contact_data_long, household_complete, by = "contact_id") %>%
  dplyr::filter(contact_type != "NA")

length(unique(household_df$household_id))

## if two individuals live togeher, add them as contacts
household_contacts <- household_df %>%
  dplyr::group_by(household_id) %>%
  dplyr::mutate(individuals = list(id)) %>%
  dplyr::filter(length(individuals) > 1) 

## testing getting out permutations
combn(c("e5249500", "ca2448ad", "2ebf8b2e", "c5c409c4"), m = 2)
combn(c("e5249500", "ca2448ad"), m = 2)
class(household_contacts$individuals)

permute <- function(x) {
  combn(x, m = 2)
}

test <- c(lapply(household_contacts$individuals, FUN = permute))
head(test)

test <- as.data.frame(test)

household_contacts_long <- t(test)
household_contacts %<>% as.data.frame()
household_contacts_long %<>% as.data.frame()

names(household_contacts_long)[1] <- "id"

household_contacts_long$id <- as.character(household_contacts_long$id)

test <- dplyr::left_join(household_contacts_long, household_contacts, by = "id")
names(test)[c(1:2, 4)] <- c("to", "from", "to_household")

test %<>% dplyr::mutate(from_household = to_household) %>%
  dplyr::mutate(contact_type = 4)

## keep the columns we want
household_data_final <- test[,c("to", "from", "household_id", "contact_type")]

## check that this matches the rest of the contact dataset
head(contact_data)
names(household_data_final) <- c("id", "contact_id", "household_id", "contact_type")
household_data_final$contact_household <- household_data_final$household_id

contact_data <- rbind(contact_data, household_data_final)

# remove duplicates (to --> from = from --> to) -- we are not interested in the type of contact
df_sort <- t(apply(contact_data[,1:2], 1, sort))
contact_final <- contact_data[!duplicated(df_sort),]

## checking that duplicates have been removed
if(nrow(contact_data) == nrow(contact_final)) {
  stop("No duplicates removed")
}

## contact type 3 = indirect contacts
direct_contact <- contact_final %>% dplyr::filter(contact_type != 3)

full_data$id <- as.factor(full_data$id)

## filter for only cases who tested positive
data <- dplyr::filter(full_data, positive == TRUE)

## starts with 4 == the date of testing (this is an artefact of reading in an excel file)
long_data <- tidyr::pivot_longer(data = data, cols = starts_with("4"), names_to = "date_confirmed",
                                 values_to = "test_result") %>%
  ## another artefact of reading in from xlsx -- usually would use !is.na() 
  dplyr::filter(test_result != "NA")

## break clause to check we keep all the entries
if (nrow(data) != length(unique(long_data$id))) {
  stop("Entries have been lost when making the data frame long")
}

length(unique(long_data$id))
nrow(data)
nrow(long_data)
covid_testing_data <- long_data

## fixing the dates of testing -- origin is again an artefact of xlsx
covid_testing_data %<>% dplyr::mutate(date_confirmed = as.Date(as.numeric(date_confirmed), 
                                                               origin = "1899-12-30"))


class(covid_testing_data$date_confirmed)

if(max(covid_testing_data$date_confirmed) != "2020-03-10") {
  stop("Error in formatting date -- max should be 10 March 2020")
}

## tests marked from 3 - 5 March were in fact not undertaken until at least the 6th -- correct to 7th
sort(unique(covid_testing_data$date_confirmed))

## want to only filter through individuals who have at least one positive test
### from looking there are slight issues with household ids -- don't use these for now.
positive_ids <- (full_data[full_data$positive == TRUE,])$id

## covid data has some missing ids off positive people 
positives <- covid_testing_data %>% dplyr::filter(id %in% positive_ids)
positives <- positives[,c("id", "household_id", "date_confirmed", "test_result", "first_symptoms_date")]

## creformatting to get the kind of data frame we want - positive individuals and when they are first and last +ve/-ve 
covid_testing_data_pos <- positives %>% dplyr::group_by(id) %>%
  dplyr::filter(test_result == "Pos") %>%
  dplyr::mutate(first_positive = min(date_confirmed)) %>%
  dplyr::mutate(last_positive = max(date_confirmed))

covid_testing_data_neg <- positives %>% dplyr::group_by(id) %>%
  dplyr::filter(test_result == "Neg") %>%
  dplyr::mutate(first_negative = min(date_confirmed)) %>%
  dplyr::mutate(last_negative = max(date_confirmed))

covid_final <- dplyr::full_join(covid_testing_data_pos, covid_testing_data_neg, 
                                by = c("id","household_id")) %>%
  dplyr::left_join(covid_testing_data, by = c("id", "household_id"))

covid_final <- covid_final[,c("id", "household_id", "first_positive", "last_positive", 
                              "first_negative", "last_negative", "first_symptoms_date")] %>%
  dplyr::distinct() %>%
  dplyr::mutate(onset = first_symptoms_date) %>% 
  dplyr::mutate(onset = as.Date(as.numeric(onset), origin = "1899-12-30"))

## NAs introduced by coercion due to as.Date(as.numeric("NA"), origin = "1899-12-30") -- another excel consequence

linelist <- covid_final
## true means no duplicates
nrow(linelist) == length(unique(linelist$id))

linelist_extras <- data.frame(id = missing_households$id, household_id = missing_households$household_id,
                              first_positive = rep(NA, length(missed_ids)), 
                              last_positive = rep(NA, length(missed_ids)), 
                              first_negative = rep(NA, length(missed_ids)),
                              last_negative = rep(NA, length(missed_ids)), 
                              first_symptoms_date = rep(NA, length(missed_ids)), 
                              onset= rep(NA, length(missed_ids)))

linelist %<>% as.data.frame()
linelist_extras %<>% as.data.frame()

names(linelist) == names(linelist_extras)

linelist <- rbind(linelist, linelist_extras)
linelist <- linelist[,c("id", "household_id", "first_positive", "last_positive", 
                        "first_negative", "last_negative", "onset")]


# save ------------------------------------------------------------------------
  
  
saveRDS(full_data, "full_data.rds")
saveRDS(linelist, "linelist.rds")
saveRDS(contact_final, "all_contact.rds")
saveRDS(direct_contact, "direct_contact.rds")

# write.table(linelist, "linelist.tsv", quote=FALSE, sep='\t', row.names = FALSE)
# write.table(contact_final,"all_contacts.tsv", quote=FALSE, sep='\t', row.names = FALSE)
# write.table(direct_contact,"direct_contacts.tsv", quote=FALSE, sep='\t', row.names = FALSE)
