#BIOSTAT699 Project 4 
#Model Building

#packages
library(survival)
library(dplyr)
library(readxl)
library(lubridate)

#read data set 
surv_data <- read_xlsx("final_target_RP.xlsx")
#View(surv_data)

#check number of adjuvant events
sum(surv_data$delt == 1) #469 adjuvant events

#turn timeX into days since time = 0 (RP date)
surv_data <- surv_data %>% mutate(X = difftime(timeX, rrpdate, units = "days"))

#Question 1 ----
#Has the use of adjuvant therapy in MUSIC changed as a result of the Randomized Trial Data published in 2020?

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
surv_data$race <- as.factor(surv_data$race)
surv_data$race <- relevel(surv_data$race, ref = 4) #use white as reference group, since it has the highest sample size

#still working on this model!
mod2 <- coxph(Surv(X,delt) ~ i_2020 + age + race + s_gg + margin + epe + svi, data = surv_data)




#Question 3 ----
#Are there practices that may need a quality improvement 
#intervention to address the lack of de-implementation of adjuvant therapy?


