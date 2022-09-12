###Project: HARVI Genomic Analysis
###Purpose: Collapse inpatient encounters to single summary row
###Author: Josh Petrie
###Date: 9/12/2022

# =========== Load libraries ==================
library(dplyr)
library(data.table)

# =========== Define Function ==================
#Raw encounter data have many rows for each admission, and 
#not every admission has a summary row.
#This function collapses rows with overlapping start and end times into
#a single summary encounter. Readmission within 24 hours of a discharge
#is considered a continuation of the same admission.

collapse_encounters <- function(df){
  df <- df %>%
    #change dates from character to POSIX
    mutate(AdmitDate = as.POSIXct(AdmitDate,'%m/%d/%Y %H:%M', tz = "utc"),
           DischargeDate = as.POSIXct(DischargeDate,'%m/%d/%Y %H:%M', tz = "utc"),
           InpatientAdmitDate = as.POSIXct(InpatientAdmitDate,'%m/%d/%Y %H:%M', tz = "utc"),
           InpatientDischargeDate = as.POSIXct(InpatientDischargeDate,'%m/%d/%Y %H:%M', tz = "utc"),
           ActivityDate = as.POSIXct(ActivityDate,'%m/%d/%Y %H:%M', tz = "utc"),
           #if discharge date is before admission, recode to be same as admission
           DischargeDate = if_else(DischargeDate<AdmitDate,AdmitDate,DischargeDate)) %>%
    #remove rows for encounters that didn't happen
    filter(is.na(PatientClassCode) | PatientClassCode != "Cancelled",
           is.na(DischargeTypeName) | !grepl("never arrived",DischargeTypeName,ignore.case = TRUE),
           is.na(DischargeTypeName) | !grepl("diverted elsewhere",DischargeTypeName,ignore.case = TRUE)) 
  
  #create age tables - 1 row per individual
  age.yrs <- df %>%
    filter(!is.na(AgeInYears)) %>%
    filter(!duplicated(PatientID)) %>%
    select(PatientID,AgeInYears)
  
  age.mos <- df %>%
    filter(!is.na(AgeInMonths)) %>%
    filter(!duplicated(PatientID)) %>%
    select(PatientID,AgeInMonths)
  
  age.wks <- df %>%
    filter(!is.na(AgeInWeeks)) %>%
    filter(!duplicated(PatientID)) %>%
    select(PatientID,AgeInWeeks)
  
  #this is the main step to collapse encounters falling within a single admission
  df <- df %>%
    #sort by individual patient and admission date
    arrange(PatientID,AdmitDate) %>%
    #for each patient...
    group_by(PatientID) %>%
    #create an index identifying separate admissions
    #if the next (lead) adm date is greater than the max disch date so far (cummax), the index increases
    #a 24 hour allowance (-86400) is given so that disch and re-adm on same day are treated as single enc
    mutate(indx = c(0, cumsum(as.numeric(lead(AdmitDate))-86400 >
                                cummax(as.numeric(DischargeDate)))[-n()])) %>%
    #for each patient and each encounter index...
    group_by(PatientID,indx) %>%
    #the adm date is the min adm date and the disch date is the max disch date w/in the index
    #keep all discharge codes and encounter ids from each row being collapsed
    summarise(AdmitDate = min(AdmitDate), DischargeDate = max(DischargeDate),
              DischargeTypeCode = paste(DischargeTypeCode,collapse = ","),
              DischargeTypeName = paste(DischargeTypeName,collapse = ","),
              EncounterIDlist = paste(EncounterID,collapse = ","))
  
  #remove NAs from discharge codes
  df <- as.data.frame(df) %>%
    mutate(DischargeTypeCode = gsub(",NA","",DischargeTypeCode),
           DischargeTypeName = gsub(",NA","",DischargeTypeName),
           DischargeTypeCode = gsub("NA,","",DischargeTypeCode),
           DischargeTypeName = gsub("NA,","",DischargeTypeName))
  
  #merge in age variables
  df <- left_join(df,age.yrs,by="PatientID")
  df <- left_join(df,age.mos,by="PatientID")
  df <- left_join(df,age.wks,by="PatientID")
  
  return(df)
}

# =========== Set data path and read data ==================
filepath_1718 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 17-18/"
filepath_1920 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 19-20/"
data_1718 <- fread(paste0(filepath_1718, "all_hosp_1718_enc.csv"))
data_1920 <- fread(paste0(filepath_1920, "all_hosp_1920_enc.csv"))

# =========== Apply function ==================
collapsed_1718 <- collapse_encounters(data_1718)
collapsed_1920 <- collapse_encounters(data_1920)

# =========== Write data ==================
fwrite(collapsed_1718,
       paste0(filepath_1718, "derived/collapsed_encounters_1718.csv"), 
       row.names = FALSE)
fwrite(collapsed_1920,
       paste0(filepath_1920, "derived/collapsed_encounters_1920.csv"), 
       row.names = FALSE)