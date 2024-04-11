#install.packages("readxl")
library(readxl)

RP = read_excel("/Users/yuxichen/biostat699/project4/aRT.xlsx", sheet = "RP")
PSA = read_excel("/Users/yuxichen/biostat699/project4/aRT.xlsx", sheet = "PSA")
Treatments = read_excel("/Users/yuxichen/biostat699/project4/aRT.xlsx", sheet = "Other Treatments")

#####first thing, we need to find our target samples

###total sample size, comment:
#### many people recovered so they dont need treatment anymore 
length(unique(RP$patientid)) #9746
length(unique(PSA$patientid)) #9735
length(unique(Treatments$patientid)) #1955

### Include in the analysis only patients who either have their first PSA >28d undetectable (<= 0.1) 
###or have any PSA 1-28d post-prostatectomy undetectable (<= 0.1)

#first step, move surgery data from RP data to PSA data
#comment: 9746 - 9735 = 11 people who dont have PSA measurement available after surgery, they are removed
library(dplyr)
PSA <- PSA %>%
  left_join(RP %>% select(patientid, rrpdate), by = "patientid")

#second step, filter out the target individuals with the criteria
#may need to create a indicator variable to indicate whether it is the first PSA measurement after 28 days or not 

##first step, creating a variable representing the difference between surgery date and PSA measure date
PSA$time_change = as.numeric(PSA$labdate - PSA$rrpdate, units = "days")

##based on the criteria, we need 
##first, filter out individuals who have any PSA 1-28d post-prostatectomy undetectable (<= 0.1)
target =  subset(PSA, time_change >= 0  & time_change <= 28 & labvalue <= 0.1 ) ## only 365 rows of data
length(unique(target$patientid)) # 365 individuals

##second, Include in the analysis only patients who have their first PSA >28d undetectable (<= 0.1) 
##more specifically, we need the first PSA measure to be after 28 days but it needs to be before the secondary treatment 
#so that we know whether it is persistently positive or not. If it is persistently positive, we remove this indiivdual 
#from our analysis.

target2 = subset( PSA, !(time_change >= 0  & time_change <= 28 & labvalue <= 0.1 )) ## 92696 - 365 = 92331 rows of data
target2 = subset(target2, !(patientid %in% unique(target$patientid))) #remove individual in target 1 from target2


##let us remove days of change that is shorter than 28 days.
target2 = subset( target2, time_change > 28)

##now, we need a indicator to indicate whether it is the first PSA measurement after 28 days or not for each individudal
##first PSA measurement -> the time-change is the smallest
target2 <- target2 %>%
  group_by(patientid) %>%
  mutate(first_PSA = as.integer(time_change == min(time_change))) %>%
  ungroup()


##but, we need to first PSA measurement after 28 days, not only undetactable, but also before the secondary treatment 
##to determine whether they are persistently positive or not 
##therefore, we need information from treatment dataset
#now we merge treatment date from Treatments dataset to our target2
#but, for Treatment dataset, one individual may have multiple treatments
#therefore, we need to identify the first treatment individual experience

#we find that an individual can have multilpe treatment at the same time, for example; patientID 102
#but it does not matter because what we care about is the first treatment time. we need the 
#first treatment time to be after the first PSA measurement after 28 days of surgery.
first_treatment <- Treatments %>%
  group_by(patientid) %>%
  filter(asdate == min(asdate))

target2 <- target2 %>%
  left_join(first_treatment %>% select(patientid, asdate), by = "patientid")
####################################wrong
##now, we need to keep rows that 1. psa value is the first psa value after 28 days and less than 0.1
##2. the date for the first psa measurements that satisfy the condition(after 28 days) has to be earlier than the treatment.
##for the second criteria, it removes those individuals who may have first undetectable psa measurements after treatment.
##these individuals will not be eligible for adjuvent RT in the first place.
#target2 = subset(target2, first_PSA == 1 & labvalue <= 0.1 & labdate < asdate)
target2 = subset(target2, first_PSA == 1 & labvalue <= 0.1)

##the individual left in target2 is the patient have their first PSA >28d undetectable (<= 0.1) 
#length(unique(target2$patientid)) #1217 individuals

##next, include individuals who does not receive anything. In other words,
##include patients who do not need treatment after surgery
##this may overlap with target1 that has undetectable PSA within 28 days, (but it may not matter)
#recoverd_IDs =  setdiff(unique(PSA$patientid), Treatments$patientid)
#length(unique(recoverd_IDs)) #7780 individuals just recovered
#target3 = PSA[PSA$patientid %in% recoverd_IDs,]

#########
##now, we know the id of all individuals who satstify the criteria
##merge all interested individuals together
all_target_IDs = c(unique(target$patientid), unique(target2$patientid))
final_target_RP = subset(RP, patientid %in% all_target_IDs) ## final interested populaiton is 9063.

#########because it is survival analysis, now we need to determine follow-up time
####

#i need to remove patient who got surgery on 2020
#install.packages("lubridate")
library(lubridate)
final_target_RP = subset(final_target_RP, year(rrpdate) !=  2020)

##now determine the events, competing events, censoring time
# time0: time of surgery
# time of events: adjuvent radiation therapy, Adt
# time of competing events: salvage, and other treatments.
# 
