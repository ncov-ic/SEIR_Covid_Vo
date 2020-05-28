##########################################################################
# Suppression of COVID-19 outbreak in the municipality of Vo', Italy
# Logistic regression on linelist data
# Testing association between positivity and age/gender
##########################################################################

library("readxl")

data <- read_excel("anonymised_data_public_final.xlsx")

####### create variables for logistic regression ####### 

# POS - positive at first or second survey
data$POS <- NA
data$POS[which(data$positive == "FALSE")]  <- 0
data$POS[which(data$positive == "TRUE")]  <- 1
table(data$POS)

# POS1 - positive at first survey
data$POS1 <- NA
data$POS1[which(data$first_sampling == "Positive")]  <- 1
data$POS1[which(data$first_sampling == "Negativized")]  <- 1
data$POS1[which(data$first_sampling == "Negative")]  <- 0
table(data$POS1)

# POS2 - positive at second survey
data$POS2 <- NA
data$POS2[which(data$second_sampling == "Positive")]  <- 1
data$POS2[which(data$second_sampling == "Positivized")]  <- 1
data$POS2[which(data$second_sampling == "Negative")]  <- 0
table(data$POS2)

# aggregate 91+ with 81-90
data$age_group[which(data$age_group == "91+")] <- "81-90"
data <- as.data.frame(unclass(data))
#data$age_group <- as.factor(data$age_group)
data$age_group <- relevel(data$age_group, ref = "21-30")

####### logisitic regression analysis ########

# predict outcome using age and gender  
fit  <- glm(POS ~ age_group + gender, family= binomial(link='logit'), data = data)
summary(fit)

fit  <- glm(POS1 ~ age_group + gender, family= binomial(link='logit'), data = data)
summary(fit)

fit  <- glm(POS2 ~ age_group + gender, family= binomial(link='logit'), data = data)
summary(fit)


# predict outcome using age  
fit  <- glm(POS ~ age_group, family= binomial(link='logit'), data = data)
summary(fit)

fit  <- glm(POS1 ~ age_group, family= binomial(link='logit'), data = data)
summary(fit)

fit  <- glm(POS2 ~ age_group, family= binomial(link='logit'), data = data)
summary(fit)


# predict outcome using gender  
fit  <- glm(POS ~ gender, family= binomial(link='logit'), data = data)
summary(fit)

fit  <- glm(POS1 ~ gender, family= binomial(link='logit'), data = data)
summary(fit)

fit  <- glm(POS2 ~ gender, family= binomial(link='logit'), data = data)
summary(fit)


