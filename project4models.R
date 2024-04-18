#BIOSTAT699 Project 4 
#Model Building

#packages ----
library(survival)
library(dplyr)
library(readxl)
library(lubridate)
library(naniar)
library(UpSetR)

#read data set ----
surv_data <- read_xlsx("final_target_RP.xlsx")
#View(surv_data)


#missing values ----
sum(is.na(surv_data) == T) #35243 NAs
surv_data %>% select(-c(T1,T21, T22, T23, Censoring, asdate, trttype)) %>% gg_miss_upset() #visualize missing data, very small proportion of missing data

#check number of adjuvant events
sum(surv_data$delt == 1) #469 adjuvant events


#Question 1 ----
#Has the use of adjuvant therapy in MUSIC changed as a result of the Randomized Trial Data published in 2020?

#turn timeX into days since time = 0 (RP date)
surv_data <- surv_data %>% mutate(X = difftime(timeX, rrpdate, units = "days"))

#create indicator of RP before or after 2020
surv_data <- surv_data %>% mutate(i_2020 = ifelse(year(rrpdate) < 2020, 0, 1))

mod1 <- coxph(Surv(X, delt) ~ i_2020, data = surv_data)
summary(mod1)
#HR = 0.33, p-value < 0.001, significant decrease in risk of adjuvant therapy after 2020, compared to before 2020


#Question 2 ----
#Do the findings from Objective 1 differ when considering patient and disease characteristics?
#patient characteristics: age, race, family history,comorbidities
#disease characteristics: s_gg (grade group), margin, epe, svi, pre-op PSA
#need to create pre-op PSA and comorbidities sum variables!!

#recoding categorical variables
#fix race categories!!
surv_data <- surv_data %>% mutate(race_bin = case_when(
  race == 4 ~ "white",
  race < 4 | race == 5 ~"non-white",
  race == 6 ~ NA))
surv_data$race_bin <- as.factor(surv_data$race_bin)
surv_data$race_bin <- relevel(surv_data$race_bin, ref = "white") #use white as reference group, since it has the highest sample size

surv_data$s_gg <- as.factor(surv_data$s_gg)
surv_data$s_gg <- relevel(surv_data$s_gg, ref = 1) #use low risk group as reference?

surv_data$margin <- as.factor(surv_data$margin)
surv_data$epe <- as.factor(surv_data$epe)
surv_data$svi <- as.factor(surv_data$svi)

#family history
surv_data <- surv_data %>% mutate(famhx_bin = case_when(
  famhx == 6 ~ "No",
  famhx <= 4 ~ "Yes",
  famhx == 5 ~ "Unknown"
))
table(surv_data$famhx_bin)
surv_data$famhx_bin <- relevel(as.factor(surv_data$famhx_bin), ref = "No")

#count of comorbidities for each patient
comorbid <- surv_data %>% select(patientid, mi, chf, pvd, cvd, dementia, cpd, ctd, ulcer, mildliver, diabnocomp, diaborgdam, hemiplegia, cri, tumor, leukemia, lymphoma, modsevliver, mettumor, aids)
comorbid <- comorbid %>% mutate(com_count = rowSums(comorbid %>% select(-patientid))) %>% select(patientid, com_count)
surv_data <- surv_data %>% merge(comorbid, by = "patientid")

#first version of this model
mod2 <- coxph(Surv(X,delt) ~ i_2020 + age + race_bin + famhx_bin + com_count + s_gg + margin + epe + svi + labvalue, data = surv_data)
summary(mod2)


#Question 3 ----
#Are there practices that may need a quality improvement 
#intervention to address the lack of de-implementation of adjuvant therapy?


