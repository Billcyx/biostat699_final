#BIOSTAT699 Project 4 
#Model Building

#packages ----
library(survival)
library(dplyr)
library(readxl)
library(lubridate)
library(naniar)
library(UpSetR)
library(gtsummary)
library(cmprsk)
library(ggplot2)
library(survminer)
library(xtable)

#read data set ----
surv_data <- read_xlsx("final_target_RP.xlsx")
#View(surv_data)


#missing values ----
sum(is.na(surv_data) == T) #41801 NAs
surv_data %>% select(-c(T1,T21, T22, T23, Censoring, asdate, trttype)) %>% gg_miss_upset() #visualize missing data, very small proportion of missing data

#check number of adjuvant events
sum(surv_data$delt == 1) #474 adjuvant events


#Question 1 ----
#Has the use of adjuvant therapy in MUSIC changed as a result of the Randomized Trial Data published in 2020?

#turn timeX into days since time = 0 (RP date)
surv_data <- surv_data %>% mutate(X = difftime(timeX, rrpdate, units = "days"))

#create indicator of RP before or after 2020
surv_data <- surv_data %>% mutate(i_2020 = ifelse(year(rrpdate) < 2020, 0, 1))

mod1 <- coxph(Surv(X, delt) ~ i_2020, data = surv_data)
summary(mod1)
#HR = 0.3393, p-value < 0.001, significant decrease in hazard of adjuvant therapy after 2020, compared to before 2020

table_mod1 <- tbl_regression(mod1, exponentiate = T, label = list(i_2020 ~ "Post-2020"))

#cumulative incidence function
#define competing risk indicator
surv_data <- surv_data %>% mutate(comprisk = case_when(
  timeX == T1 ~ 1,
  timeX == T2 ~ 2,
  timeX == Censoring ~ 0
))
#run CIF
CIF <- cuminc(surv_data$X, surv_data$comprisk, surv_data$i_2020, cencode = 0)
#plot CIF, adjuvant in red, competing risk in blue
plot(CIF, 
     ylim = c(0,0.15), 
     xlim = c(0,325), 
     lwd = 1.5,
     color = c("red", "red", "blue", "blue"),
     main = "Cumulative Incidence Function for Pre- vs Post-2020",
     xlab = "Days")

#better CIF plot
CIFplot <- ggcompetingrisks(fit = CIF,
                 xlab = "Days",
                 title = "Cumulative Incidence Function for Pre- vs Post-2020",
                 ylim = c(0,0.15),
                 xlim = c(0, 350),
                 lwd = 2,
                 multiple_panels = F) +
  scale_linetype_manual(name="Pre/Post 2020",values=c(1,2),labels=c("Pre-2020","Post-2020")) +
  scale_color_discrete(labels = c(1 ~ "ART", 2 ~ "Competing Event"), type = c("red", "blue"))

CIFplot$mapping <- aes(x = time, y = est, colour = event, linetype = group)

CIFplot + labs(linetype = "Pre/Post 2020", colour = "Event")

#Question 2 ----
#Do the findings from Objective 1 differ when considering patient and disease characteristics?
#patient characteristics: age, race, family history,comorbidities
#disease characteristics: s_gg (grade group), margin, epe, svi, pre-op PSA
#need to create pre-op PSA and comorbidities sum variables!!
surv_data <- surv_data %>% mutate(age_c = age - mean(age, na.rm = T))

#recoding categorical variables
#race
surv_data <- surv_data %>% mutate(race_bin = case_when(
  race == 4 ~ "white",
  race < 4 | race == 5 ~"non-white",
  race == 6 ~ NA))
surv_data$race_bin <- as.factor(surv_data$race_bin)
surv_data$race_bin <- relevel(surv_data$race_bin, ref = "white") #use white as reference group, since it has the highest sample size

#grade group to risk group
surv_data <- surv_data %>% mutate(risk_group = case_when(
  s_gg == 1 ~ "low",
  s_gg == 2 | s_gg == 3 ~ "intermediate",
  s_gg == 4 | s_gg == 5 ~ "high"
))
surv_data$risk_group <- as.factor(surv_data$risk_group)
surv_data$risk_group <- relevel(surv_data$risk_group, ref = "low") #use low risk group as reference?

#margin, epe, svi
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

#pre-op PSA
surv_data <- surv_data %>% rename("preopPSA" = "labvalue")
hist(surv_data$preopPSA) #highly right skewed!
#log transformation of pre-op PSA, add 1 to account for zeros
surv_data <- surv_data %>%  mutate(log_preopPSA = log(preopPSA + 1))
hist(surv_data$log_preopPSA) #looks more normal in distribution :)

#first version of this model
mod2 <- coxph(Surv(X,delt) ~ i_2020 + age_c + race_bin + famhx_bin + com_count + risk_group + margin + epe + svi + log_preopPSA, data = surv_data)
summary(mod2)
#HR for 2020: 0.2913, p-value < 0.001, similar result to unadjusted 
#high risk group has significantly increased hazard of adjuvant therapy, compared to low risk group, makes sense

mod2_v2 <-coxph(Surv(X,delt) ~ i_2020 + age_c + race_bin + famhx_bin + com_count + risk_group + margin + log_preopPSA, data = surv_data)
summary(mod2_v2)


library(finalfit)
explanatory <- c("age", "race_bin", "famhx_bin", "", "heart", "cancer", "pre_psa", "after_2020")
dependent_os  <- "Surv(time, delta1)"
final %>%
  coxphmulti(dependent_os, explanatory) %>%
  cox.zph() %>%
  {zph_result <<- .} %>%
  plot(var=5)
zph_result


table_mod2 <- tbl_regression(mod2, exponentiate = T, label = list(i_2020 ~ "Post-2020", age_c ~ "Age", race_bin ~ "Race", famhx_bin ~ "Family History", com_count ~ "Number of Comorbidities", risk_group ~ "Risk Group", margin ~ "Surgical Margin Status", epe ~ "Extraprostatic extension", svi ~ "Seminal vesicle invasion", log_preopPSA ~ "log(Pre-Operative PSA)"))
table_mod12 <- tbl_merge(list(table_mod1, table_mod2), tab_spanner = c("Model 1", "Model 2")) %>% bold_labels()
table_mod12
#xtable(as.data.frame(table_mod12))

#Question 3 ----
#Are there practices that may need a quality improvement 
#intervention to address the lack of de-implementation of adjuvant therapy?

#identify practices with large sample size
temp <- surv_data %>% filter(delt == 1) 
table(temp$providerid, temp$i_2020) #check how many patients for each provider (pre- vs post-2020)
providers <- surv_data %>% group_by(providerid) %>% summarise(n = n()) %>% filter(n > 200)
surv_data_provider <- surv_data %>% filter(providerid %in% providers$providerid)
surv_data_provider$providerid <- as.factor(surv_data_provider$providerid)

mod3 <- coxph(Surv(X,delt) ~ i_2020 + age_c + race_bin + famhx_bin + com_count + risk_group + margin + epe + svi + log_preopPSA + i_2020:providerid + strata(providerid) , data = surv_data_provider)
summary(mod3)

#another attempt
events <- surv_data %>% filter(delt == 1) %>% group_by(providerid, i_2020) %>% summarise(n = n()) %>% ungroup()
events <- events %>% mutate(event = ifelse( n > 1, 1, 0)) #only want practices with more than one event
events <- events %>% group_by(providerid) %>% summarise(event_sum = sum(event)) %>% mutate(events_prepost = ifelse(event_sum == 2, 1, 0)) #check for events both pre and post 2020
surv_data_3 <- merge(surv_data, events, by = "providerid")
surv_data_3 <- surv_data_3 %>% mutate(provider = case_when(
  (events_prepost == 1) ~ as.character(providerid),
  (events_prepost == 0 & practice_type == 1) ~ "smallAcademic",
  (events_prepost == 0 & practice_type == 2) ~ "smallCommunity",
  (events_prepost == 0 & practice_type == 3) ~ "smallHybrid"
))
surv_data_3 <- surv_data_3 %>% filter(provider != "smallAcademic")
surv_data_3$provider <- relevel(as.factor(surv_data_3$provider), ref = "smallCommunity")
mod3v2 <- coxph(Surv(X,delt) ~ i_2020 + age_c + race_bin + famhx_bin + com_count + risk_group + margin + epe + svi + log_preopPSA + i_2020:strata(provider) , data = surv_data_3)
summary(mod3v2)

table_mod3 <- tbl_regression(mod3v2, exponentiate = T, include = c("i_2020", "i_2020:strata(provider)"), label = list(i_2020 ~ "Post-2020", 'i_2020:strata(provider)' ~ "Post-2020*Provider"))
table_mod3


#table 1 ----
surv_data <- surv_data %>% mutate(i_2020_label = ifelse(i_2020 == 0, "Pre-2020", "Post-2020"))
surv_data$i_2020_label <- as.factor(surv_data$i_2020_label)
surv_data$i_2020_label <- relevel(surv_data$i_2020_label, ref = "Pre-2020")
table_1 <- surv_data %>% select(i_2020_label, age, race_bin, famhx_bin, com_count, risk_group, margin, epe, svi, preopPSA) %>% 
  tbl_summary(by = "i_2020_label", 
              label = list(age ~ "Age", race_bin ~ "Race", famhx_bin ~ "Family History", com_count ~ "Number of Comorbidities", risk_group ~ "Risk Group", margin ~ "Surgical Margin Status", epe ~ "Extraprostatic extension", svi ~ "Seminal vesicle invasion", preopPSA ~ "Pre-Operative PSA"),
              type = list(com_count ~ "continuous"),
              missing_text = "Missing")
table_1
xtable(as.data.frame(table_1))

