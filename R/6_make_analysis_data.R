###Project: HARVI Genomic Analysis
###Purpose: Make final analysis dataset
###Author: Tiffany Wan, Josh Petrie
###Date: 2022-09-12

# =========== Load Libraries ==================
library(tidyverse)
library(data.table)
library(sqldf)
library(lubridate)

# =========== Define Data Management Functions ==================

merge_data <- function(flu,loc,clusters,demo){
  # only need gender from demographics
  demo <- demo %>%
    select(PatientID, GenderName)
  
  # reformat dates
  flu <- flu %>%
    mutate(
      ORDER_DATE = as.POSIXct(ORDER_DATE,'%Y-%m-%d %H:%M:%S', tz = "utc"),
      AdmitDate = as.POSIXct(AdmitDate,'%Y-%m-%d %H:%M:%S', tz = "utc"),
      DischargeDate = as.POSIXct(DischargeDate,'%Y-%m-%d %H:%M:%S', tz = "utc")
    )
  
  # reformat dates
  loc <- loc %>% 
    mutate(
      StartDate = as.character(StartDate),
      EndDate = as.character(EndDate)
    )
  
  # make transfer data wide
  loc <- loc %>%
    group_by(PatientID,AdmitDate) %>%
    mutate(Transfer = 1:n()) %>%
    gather("StartDate","EndDate","MiChartDeptName","Unit", 
           key = variable, value = value) %>%
    unite(combined, variable, Transfer) %>%
    spread(combined, value) %>%
    mutate_at(
      vars(starts_with("EndDate"),starts_with("StartDate")),
      as.POSIXct,
      format = '%Y-%m-%d %H:%M:%S', 
      tz = "utc"
    ) %>%
    select(
      PatientID, AdmitDate, StartDate_1, EndDate_1, MiChartDeptName_1, Unit_1,
      StartDate_2, EndDate_2, MiChartDeptName_2, Unit_2,
      StartDate_3, EndDate_3, MiChartDeptName_3, Unit_3,
      StartDate_4, EndDate_4, MiChartDeptName_4, Unit_4,
      StartDate_5, EndDate_5, MiChartDeptName_5, Unit_5,
      StartDate_6, EndDate_6, MiChartDeptName_6, Unit_6,
      StartDate_7, EndDate_7, MiChartDeptName_7, Unit_7,
      StartDate_8, EndDate_8, MiChartDeptName_8, Unit_8,
      StartDate_9, EndDate_9, MiChartDeptName_9, Unit_9,
      StartDate_10, EndDate_10, MiChartDeptName_10, Unit_10,
      StartDate_11, EndDate_11, MiChartDeptName_11, Unit_11,
      StartDate_12, EndDate_12, MiChartDeptName_12, Unit_12,
      StartDate_13, EndDate_13, MiChartDeptName_13, Unit_13,
      StartDate_14, EndDate_14, MiChartDeptName_14, Unit_14,
      StartDate_15, EndDate_15, MiChartDeptName_15, Unit_15,
      StartDate_16, EndDate_16, MiChartDeptName_16, Unit_16,
      StartDate_17, EndDate_17, MiChartDeptName_17, Unit_17,
      StartDate_18, EndDate_18, MiChartDeptName_18, Unit_18,
      StartDate_19, EndDate_19, MiChartDeptName_19, Unit_19,
      StartDate_20, EndDate_20, MiChartDeptName_20, Unit_20
    ) %>%
    ungroup()
  
  #Make clusters 1 row per person/admission
  clusters <- clusters %>%
    filter(!duplicated(clusters[,c("PatientID","AdmitDate")])) %>%
    select(PatientID,AdmitDate,cluster)
  
  #merge flu and demo
  merged <- left_join(flu,demo,by="PatientID")
  
  #merge in transfer locations
  merged <- left_join(merged,loc,by=c("PatientID","AdmitDate"))
  
  #merge in cluster ids
  merged <- left_join(merged,clusters,by=c("PatientID","AdmitDate"))
  
  #make short ID
  merged <- merged %>%
    ungroup() %>%
    arrange(desc(analysis_set), ORDER_DATE) %>%
    mutate(ID = ifelse(analysis_set==1,match(PatientID,unique(PatientID)),NA))
  
  return(merged)
}

# =========== Set data path and read data ==================
filepath_1718 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 17-18/"
filepath_1920 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 19-20/"


#patient transfer data
loc_1718 <- fread(
  paste0(filepath_1718, "derived/transfers_long_1718.csv"), header = TRUE
)

loc_1920 <- fread(
  paste0(filepath_1920, "derived/transfers_long_1920.csv"), header = TRUE
)

#inpatient flu cases
flu_1718 <- fread(paste0(filepath_1718,"derived/flu_enc_1718.csv"), 
                  stringsAsFactors = FALSE) 

flu_1920 <- fread(paste0(filepath_1920,"derived/flu_enc_1920.csv"), 
                  stringsAsFactors = FALSE) 

#inpatient flu clusters
clusters_1718 <- fread(paste0(filepath_1718,"derived/clusters_1718.csv"), 
                  stringsAsFactors = FALSE) 

clusters_1920 <- fread(paste0(filepath_1920,"derived/clusters_1920.csv"), 
                  stringsAsFactors = FALSE) 

#demographics
demo_1718 <- fread(paste0(filepath_1718,"all_hosp_1718_demo.csv"), 
                   stringsAsFactors = FALSE,
                   header = TRUE) 

demo_1920 <- fread(paste0(filepath_1920,"all_hosp_1920_demo.csv"), 
                   stringsAsFactors = FALSE,
                   header = TRUE) 

# =========== Apply function ==================
analysis_1718 <- merge_data(flu_1718,loc_1718,clusters_1718,demo_1718)
analysis_1920 <- merge_data(flu_1920,loc_1920,clusters_1920,demo_1920)

# =========== Write data ==================
fwrite(analysis_1718,
       paste0(filepath_1718, "derived/analysis_1718.csv"), 
       row.names = FALSE)
fwrite(analysis_1920,
       paste0(filepath_1920, "derived/analysis_1920.csv"), 
       row.names = FALSE)