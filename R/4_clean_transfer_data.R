###Project: HARVI Genomic Analysis
###Purpose: Clean patient transfer data for those in the analysis set
###Author: Josh Petrie, Tiffany Wan
###Date: 9/12/2022

# =========== Load Packages ==================
library(dplyr)
library(data.table)
library(sqldf)
library(lubridate)

# =========== Define Data Management Function ==================
#The raw transfer data has many duplicate rows with no additional information
#We will also focus on Unit level transfers and eliminate 
#procedure, consult or clinic type transfers.
#The goal of the final data set is to have 1 row per location during each admission.
#Rows with UNKNOWN location will be added as needed so that all time during
#the admission is accounted for i.e. the end date of the current row will be
#equal to the start date of the next row for all time during the admission

collapse_transfers <- function(loc,flu){
  
  # check <- flu[duplicated(flu$PatientID)] %>% 
  #   arrange(PatientID)
  # check2 <- flu[flu$PatientID %in% unique(check$PatientID),] %>% 
  #   arrange(PatientID)
  # 
  # check <- flu[duplicated(flu[,c("PatientID","AdmitDate")])] %>% 
  #   arrange(PatientID)
  # check2 <- flu[flu$PatientID %in% unique(check$PatientID),] %>% 
  #   arrange(PatientID)
  
  # make one row per admission per individual and reformat dates
  flu <- flu[!duplicated(flu[,c("PatientID","AdmitDate")]),] %>%
    select(PatientID, AdmitDate, DischargeDate) %>% 
    mutate(AdmitDate = as.POSIXct(AdmitDate,'%Y-%m-%d %H:%M:%S', tz = "utc"),
           DischargeDate = as.POSIXct(DischargeDate,'%Y-%m-%d %H:%M:%S', tz = "utc"),
           # will add 12 hour buffer to discharge because 
           #some transfers end after discharge
           DischargeDate_plus12 = DischargeDate + hours(12)) 
  
  #reformat dates
  loc <- loc %>% 
    mutate(StartDate = as.POSIXct(StartDate,'%m/%d/%Y %H:%M', tz = "utc"),
           EndDate = as.POSIXct(EndDate,'%m/%d/%Y %H:%M', tz = "utc"))
  
  #merge flu cases and transfer data by Patient ID and relevant dates
  flu_loc <- sqldf("SELECT * FROM flu
                  LEFT JOIN loc
                  ON  flu.PatientID = loc.PatientID
                  AND flu.AdmitDate <= loc.StartDate 
                  AND flu.DischargeDate_plus12 >= loc.EndDate")
  
  # delete duplicated column (PatientID)
  flu_loc <- flu_loc[, !duplicated(colnames(flu_loc))]
  
  # =========== Remove duplicate data rows ==================
  #There are multiple rows, with different facility codes, coding for the same location stay.
  #iterate through the code below to identify which facility codes can be removed 
  
  #List of irrelevant FacilityCode values
  rm_FacilityCode <- c(
    "FADM","UMADM","UBW02","UBW03","UBW09","UBW10","SUHS","UBRL","UECP1","UKER",
    "UMLAB","URBP","UWAA","UEAA","ULHC","UDF","UNHC","PXRA","MPHP","11W","7A1",
    "10E","10W","5B1","5D","6B1","6C","7C","7DN","7E1","8A3","8B2","8DNS","AILI",
    "AMOU","B9BI","BICU","CDR","CEC3","CGCA","CPB3","CVC","CVC4","CVD","DIAL",
    "DVU","ESA","HST","IDC","INU","IRU","M4PS","MBMT","MHIP","MHOR","MHPI","MHRI",
    "MIFB","MMRI","MOTT","MPUA","MRAC","MRI","MUS","NEP","NHGI","NMA","NUC",
    "OADM","OPH","OR","ORC","OSM","PDI","PEL","PFT","PMO","RAAO","RAC","RAO",
    "RAU","RTH","TOP","TSA","TTXP","UADM","UGPL","UH","UHDO","UROP","USSS","XMBD",
    "XREC","XTB3","UECPH","UCHC","UBHC","UCHE","UDHC"
  )
  
  # arrange by admitdate, startdate and enddate, remove irrelevant facility codes, 
  # and subset to minimal columns
  flu_loc <- flu_loc %>%
    arrange(PatientID, AdmitDate, StartDate, EndDate) %>%
    filter(!(FacilityCode %in% rm_FacilityCode)) %>%
    select(PatientID, AdmitDate, DischargeDate, StartDate,EndDate, FacilityCode, 
           FacilityName, MiChartDeptName,Unit, Room, Bed)
  
  # function to create test datasets containing all transfers for individuals 
  # with specified values.
  # create_test_set <- function(df,var,val){
  #   var <- enquo(var)
  #   id_list <- df %>% filter(UQ(var)==val) %>% select(PatientID)
  #   id_list <- id_list <- id_list[['PatientID']]
  #   tmp <- df %>% filter(PatientID %in% id_list)
  #   return(tmp)
  # }
  
  # check FacilityCodes
  
  # table(test$FacilityCode)
  # 
  # test <- arrange(create_test_set(flu_loc,FacilityCode,"NSH"), 
  #                 PatientID,AdmitDate,StartDate,EndDate) #keep
  # test <- arrange(create_test_set(flu_loc,FacilityCode,"SUHS"), 
  #                 PatientID,AdmitDate,StartDate,EndDate)
  # test <- arrange(create_test_set(flu_loc,FacilityCode,"UBRL"),
  #                 PatientID,AdmitDate,StartDate,EndDate)
  # test <- create_test_set(flu_loc,FacilityCode,"UBW10") # outpatient  clinic, remove
  # test <- create_test_set(flu_loc,FacilityCode,"UCC")#keep
  # test <- create_test_set(flu_loc,FacilityCode,"UCVC")#keep
  # test <- create_test_set(flu_loc,FacilityCode,"UECP1")
  # test <- create_test_set(flu_loc,FacilityCode,"UKER")
  # test <- create_test_set(flu_loc,FacilityCode,"UMH")#keep
  # test <- create_test_set(flu_loc,FacilityCode,"UMI")
  # test <- create_test_set(flu_loc,FacilityCode,"UMLAB")
  # test <- create_test_set(flu_loc,FacilityCode,"URBP")
  # test <- create_test_set(flu_loc,FacilityCode,"UTC")#keep
  # test <- create_test_set(flu_loc,FacilityCode,"UUH")#keep
  # test <- create_test_set(flu_loc,FacilityCode,"11W") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"7A1") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"10E") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"10W") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"5B1") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"5D") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"6B1") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"6C") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"7C") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"7DN") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"7E1") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"8A3") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"8B2") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"8DNS") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"AILI") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"AMOU") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"B9BI") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"BICU") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"CDR") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"CEC3") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"CGCA") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"CPB3") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"CVC") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"CVC4") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"CVD") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"DIAL") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"DVU") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"ESA") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"HST") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"IDC") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"INU") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"IRU") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"M4PS") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"MBMT") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"MHIP") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"MHOR") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"MHPI") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"MHRI") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"MIFB") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"MMRI") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"MOTT") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"MPUA") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"MRAC") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"MRI") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"MUS") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"NEP") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"NHGI") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"NMA") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"NSH") 
    #KEEP
  # test <- create_test_set(flu_loc,FacilityCode,"NUC") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"OADM") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"OPH") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"OR") 
    #Captured by other rows - add to rm_FacilityCode and re-run
  # test <- create_test_set(flu_loc,FacilityCode,"ORC") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"OSM") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"PDI") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"PEL") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"PFT") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"PMO") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"RAAO") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"RAC") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"RAO") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"RAU") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"RTH") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"TOP") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"TSA") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"TTXP") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"UADM") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"UCC") #KEEP
  # test <- create_test_set(flu_loc,FacilityCode,"UCVC") #KEEP
  # test <- create_test_set(flu_loc,FacilityCode,"UGPL") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"UH") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"UHDO") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"UMH") #KEEP
  # test <- create_test_set(flu_loc,FacilityCode,"UROP") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"USSS") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"UTC") #KEEP
  # test <- create_test_set(flu_loc,FacilityCode,"UUH") #KEEP
  # test <- create_test_set(flu_loc,FacilityCode,"XMBD") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"XREC") 
    #Captured by other rows - add to rm_FacilityCode
  # test <- create_test_set(flu_loc,FacilityCode,"XTB3") 
    #Captured by other rows - add to rm_FacilityCode
  
  #List of irrelevant MiChartDeptName values
  rm_MiChartDeptName <- c("ADM BED COORD CNTR",
                          "Congenital Heart Center",
                          "Electroencephalogram Lab",
                          "Hepatology Clinic",
                          "Mott 8E Child Psych",
                          "Occupational Therapy",
                          "Pathology Blood Draw",
                          "Ped Physical Medicine Rehab",
                          "Ped Pulmonary Lab",
                          "Ped Transplant Clinic",
                          "Physical Medicine Rehab",
                          "Physical Therapy",
                          "PsychOncology Clinic",
                          "UH 9C ADULT PSYCH",
                          "U of M Adult Dialysis",
                          "U of M Apheresis Procedure Unit",
                          "U of M Breast Imaging",
                          "U of M Cath Lab",
                          "U of M Cardiac Surgery Clinic",
                          "U of M Cardiology Clinic",
                          "U of M Cross Sectional Radiology",
                          "U of M CT Radiology",
                          "U of M Diagnostic Vascular Unit",
                          "U of M Echocardiogram",
                          "U of M Electrophysiology Lab",
                          "U of M Epilepsy Lab",
                          "U of M Family Medicine Clinic",
                          "U of M General Medicine Clinic",
                          "U of M GIGU Radiology",
                          "U of M Hematology Clinic",
                          "U of M Infectious Disease",
                          "U of M Infusion Area",
                          "U of M Interventional Nephrology",
                          "U of M Interventional Radiology",
                          "U of M Medical Procedures Unit",
                          "U of M MRI",
                          "U of M Nephrology Clinic",
                          "U of M Neurology Clinic",
                          "U of M Nuc Med Adult Clinic",
                          "U of M Nuclear Medicine Pet Scan",
                          "U of M Oral Surgery Clinic",
                          "U of M Orthotic Prosthetics Clinic",
                          "U of M Pathology Blood Draw",
                          "U of M Patient Financial Counseling",
                          "U of M Ped Dialysis Unit",
                          "U of M Ped Hematology Oncology",
                          "U of M Ped Pulmonary Lab",
                          "U of M Pediatric Epilepsy Lab",
                          "U of M Pharmacy",
                          "U of M Photophoresis",
                          "UH Psychiatric Emergency Department",
                          "U of M Pulmonary Function Lab",
                          "U of M Pulmonary Clinic",
                          "U of M Radiation Oncology",
                          "U of M Radiology",
                          "U of M Ultrasound",
                          "U of M Urology Clinic",
                          "U of M Diagnostic Vascular Unit",
                          "U of M Bone Marrow Transplant and Leukemia",
                          "Radiation Oncology",
                          "U of M Neurology Oncology Clinic")
  
  #Remove duplicate rows, irrelevant MiChartDeptNames and Units, create test variables
  flu_loc <- flu_loc %>%
    filter(
      !duplicated(
        flu_loc[c("PatientID","StartDate","EndDate",
                  "MiChartDeptName","Unit","Room","Bed")]
      )
    ) %>%
    filter(!(MiChartDeptName %in% rm_MiChartDeptName) & Unit !="UHDO" & Unit != "MPUL") %>%
    group_by(PatientID,AdmitDate) %>%
    mutate(#Make flags for checks
      #Transfer start date is before admission
      before_admit_flag = ifelse(StartDate<AdmitDate,1,0),
      #Transfer end date is after discharge 
      after_disch_flag = ifelse(EndDate>DischargeDate,1,0),
      #Transfer end date is before start date
      end_before_start_flag = ifelse(EndDate<StartDate,1,0),
      #Consecutive transfers have overlapping time 
      overlap_flag = case_when(StartDate < lag(EndDate) ~ 1,
                                    lead(StartDate) < EndDate~ 1,
                                    TRUE ~ 0),
      #There is unaccounted for time in between consecutive transfers
      discontinuos_flag = ifelse(EndDate != lead(StartDate),1,0),
      #The next transfer end date is before the current end date
      enddate_flag = ifelse(lead(EndDate) <= EndDate,1,0),
      #The transfer start and end dates are equal
      no_time_flag = ifelse(StartDate==EndDate,1,0)) %>%
    ungroup()
  
  # check MiChartDeptNames
  
  # table(flu_loc$MiChartDeptName[flu_loc$overlap_flag==1])
  # table(flu_loc$MiChartDeptName)
  # 
  # test <- create_test_set(flu_loc,MiChartDeptName,"ADM BED COORD CNTR") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"Ped Transplant Clinic") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Pulmonary Clinic") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Family Medicine Clinic")   
    #Outpatient- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Radiology") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M General Medicine Clinic") 
    #Outpatient- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Hematology Clinic") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Ped Hematology Oncology") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"Mott 8E Child Psych") 
    #sensitive- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"UH 9C ADULT PSYCH") 
    #sensitive- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"UH Psychiatric Emergency Department") 
    #sensitive- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,
  #                         "U of M Von Voigtlander Women's Hospital") #keep
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Women's Triage Center") #keep 
  # test <- create_test_set(flu_loc,MiChartDeptName,"Congenital Heart Center") 
    #Testing or Procedure- remove 
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Diagnostic Vascular Unit") 
    #Testing or Procedure- remove 
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Interventional Radiology") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Medical Procedures Unit") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Infusion Area") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M MRI") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Cardiac Surgery Clinic") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Adult Dialysis") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,Unit,"UHDO") #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Cardiology Clinic") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Cath Lab") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Infectious Disease") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Interventional Nephrology") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Nephrology Clinic") 
    #Testing or Procedure- remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Orthotic Prosthetics Clinic") 
    #Testing/Procedure-remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"U of M Ped Dialysis Unit") 
    #Testing/Procedure-remove
  # test <- create_test_set(flu_loc,MiChartDeptName,"UH OPERATING ROOMS") #KEEP
  # test <- create_test_set(flu_loc,MiChartDeptName,"UH 6A PHYS MED REHAB") #KEEP
  # test <- create_test_set(flu_loc,MiChartDeptName,
  #                         "U of M Bone Marrow Transplant and Leukemia") #KEEP
  # test <- create_test_set(flu_loc,MiChartDeptName,"Radiation Oncology")
  
  #remove rows that have same start and end date, recreate test variables, 
  #recode is.na(MiChartDeptName)
  flu_loc <- flu_loc %>%
    filter(no_time_flag==0) %>%
    group_by(PatientID,AdmitDate) %>%
    mutate(#Remake flags for checks
      #Transfer start date is before admission
      before_admit_flag = ifelse(StartDate<AdmitDate,1,0),
      #Transfer end date is after discharge 
      after_disch_flag = ifelse(EndDate>DischargeDate,1,0),
      #Transfer end date is before start date
      end_before_start_flag = ifelse(EndDate<StartDate,1,0),
      #Consecutive transfers have overlapping time 
      overlap_flag = case_when(StartDate < lag(EndDate) ~ 1,
                               lead(StartDate) < EndDate~ 1,
                               TRUE ~ 0),
      #There is unaccounted for time in between consecutive transfers
      discontinuos_flag = ifelse(EndDate != lead(StartDate),1,0),
      #The next transfer end date is before the current end date
      enddate_flag = ifelse(lead(EndDate) <= EndDate,1,0),
      #The transfer start and end dates are equal
      no_time_flag = ifelse(StartDate==EndDate,1,0), 
      MiChartDeptName = case_when(
             is.na(MiChartDeptName) & Unit == "MHOR" ~ "Mott Operating Room",
             is.na(MiChartDeptName) & Unit == "ORC" ~ "CVC Operating Room",
             is.na(MiChartDeptName) & Unit == "CVIC5" ~ "CVC-5 CARDIOVASC ICU",
             is.na(MiChartDeptName) & Unit == "7E1" ~ "7 East - Ped BMT",
             is.na(MiChartDeptName) & Unit == "7W1" ~ "7 West - Adult BMT",
             is.na(MiChartDeptName) & Unit == "M4PS" ~ "Mott Operating Room",
             is.na(MiChartDeptName) & Unit == "MSCD" ~ "MH 12E RICU",
             TRUE ~ MiChartDeptName
             )
    ) %>%
    ungroup() %>%
    filter(Unit!="PMSP" & Unit!="SOLI" & Unit!="UHSE" & Unit!="RCBO" & Unit!="UHSF")
  
  # test <- flu_loc[is.na(flu_loc$MiChartDeptName),]
  # test <- flu_loc[flu_loc$Unit=="UHSE",]
  # test <- flu_loc[flu_loc$Unit=="RCBO",]
  # test <- flu_loc[flu_loc$Unit=="UHSF",]
  
  #Collapse rows that have the same location across continuous time
  flu_loc <- flu_loc %>%
    #sort by individual patient and admission date
    arrange(PatientID,AdmitDate,StartDate,EndDate) %>%
    #for each patient...
    group_by(PatientID,AdmitDate) %>%
    #create an index identifying separate transfers
    #if the next (lead) StartDate is not less than or equal the current EndDate OR
    #if the next (lead) MiChartDeptName is not equal to the current MiChartDeptName, 
    #the index increases
    mutate(row_count = max(row_number()),
           indx = c(0, cumsum(as.numeric(lead(StartDate)) > as.numeric(EndDate) |
                                lead(MiChartDeptName) != MiChartDeptName))[-n()]) %>%
    #tidy up the last row for each patient because lead() results in NA on last row
    mutate(
      indx = case_when(
        #the admission only has 1 row, index = 0
        row_count == 1 ~ 0,
        #if it is the last row of the admission, 
        #and the location is the same as the previous and time is continuous
        #index = the previous index
        is.na(indx) & 
          StartDate == lag(EndDate) & 
          MiChartDeptName == lag(MiChartDeptName) ~ lag(indx),
        #if it is the last row of the admission, 
        #and the location is the not same as the previous or time is not continuous
        #index = the previous index + 1
        is.na(indx) ~ lag(indx+1),
        #otherwise index is as previously calculated
        TRUE ~ indx
      )
    ) %>%
    #for each patient and each encounter index...
    group_by(PatientID,AdmitDate,indx) %>%
    #the StartDate is the min StartDate and the EndDate is the max EndDate w/in the index
    #keep all discharge codes and encounter ids from each row being collapsed
    summarise(StartDate = min(StartDate), 
              EndDate = max(EndDate),
              AdmitDate = min(AdmitDate),
              DischargeDate = min(DischargeDate),
              MiChartDeptName = min(MiChartDeptName)) %>%
    ungroup() 
  
  #Add rows with UNKNOWN location for missing time between admission and discharge
  flu_loc <- flu_loc %>%
    #sort by individual patient and admission date
    arrange(PatientID,AdmitDate,StartDate,EndDate) %>%
    #for each patient admission...
    group_by(PatientID,AdmitDate) %>%
    #create an index identifying discontinuous time
    #if the previous (lag) EndDate is not equal the current StartDate, the index increases
    mutate(
      row_count = max(row_number()),
      indx = c(0, cumsum(as.numeric(lead(StartDate)) != as.numeric(EndDate)))[-n()]
      ) %>%
    #tidy up the last row for each patient
    mutate(indx = case_when(row_count == 1 ~ 0,
                            is.na(indx) & StartDate == lag(EndDate) ~ lag(indx),
                            is.na(indx) & StartDate != lag(EndDate) ~ lag(indx)+1,
                            TRUE ~ indx)) %>%
    #for each patient and each encounter index...
    group_by(PatientID,AdmitDate,indx) %>%
    #Add row before 1st PatientID,indx instance - 
    #accounts for admit to 1st location and missing time between transfers
    do(add_row(., .before=0)) %>%
    ungroup() %>%
    #add data to new row
    mutate(PatientID = ifelse(is.na(PatientID),lead(PatientID),PatientID),
           indx = ifelse(is.na(indx),lead(indx),indx),
           AdmitDate = if_else(is.na(AdmitDate),
                               lead(AdmitDate),
                               AdmitDate),
           DischargeDate = if_else(is.na(DischargeDate),
                                   lead(DischargeDate),
                                   DischargeDate),
           MiChartDeptName = ifelse(is.na(MiChartDeptName),"UNKNOWN",MiChartDeptName),
           StartDate = if_else(is.na(StartDate),
                               if_else(indx==0,AdmitDate,lag(EndDate)),
                               StartDate),
           EndDate = if_else(is.na(EndDate),
                             if_else(is.na(lead(StartDate)),DischargeDate,lead(StartDate)),
                             EndDate)) %>%
    #for each patient...
    group_by(PatientID,AdmitDate)  %>%
    #Add row after last instance of PatientID - 
    #accounts for missing time from last transfer to discharge
    do(add_row(.))  %>%
    ungroup() %>%
    #add data to new row
    mutate(PatientID = ifelse(is.na(PatientID),lag(PatientID),PatientID),
           indx = ifelse(is.na(indx),lag(indx),indx),
           AdmitDate = if_else(is.na(AdmitDate),
                               lag(AdmitDate),
                               AdmitDate),
           DischargeDate = if_else(is.na(DischargeDate),
                                   lag(DischargeDate),
                                   DischargeDate),
           MiChartDeptName = ifelse(is.na(MiChartDeptName),"UNKNOWN",MiChartDeptName),
           StartDate = if_else(is.na(StartDate),
                               lag(EndDate),
                               StartDate),
           EndDate = if_else(is.na(EndDate),
                             DischargeDate,
                             EndDate)) %>%
    #recreate test variables
    group_by(PatientID,AdmitDate) %>%
    mutate(#Remake flags for checks
      #Transfer start date is before admission
      before_admit_flag = ifelse(StartDate<AdmitDate,1,0),
      #Transfer end date is after discharge 
      after_disch_flag = ifelse(EndDate>DischargeDate,1,0),
      #Transfer end date is before start date
      end_before_start_flag = ifelse(EndDate<StartDate,1,0),
      #Consecutive transfers have overlapping time 
      overlap_flag = case_when(StartDate < lag(EndDate) ~ 1,
                               lead(StartDate) < EndDate~ 1,
                               TRUE ~ 0),
      #There is unaccounted for time in between consecutive transfers
      discontinuos_flag = ifelse(EndDate != lead(StartDate),1,0),
      #The next transfer end date is before the current end date
      enddate_flag = ifelse(lead(EndDate) <= EndDate,1,0),
      #The transfer start and end dates are equal
      no_time_flag = ifelse(StartDate==EndDate,1,0)
      ) %>%
    ungroup()
  
  #check test variables
  # sum(flu_loc$before_admit_flag, na.rm=TRUE)
  # sum(flu_loc$after_disch_flag, na.rm=TRUE)
  # sum(flu_loc$end_before_start_flag, na.rm=TRUE)
  # sum(flu_loc$overlap_flag, na.rm=TRUE)
  # sum(flu_loc$discontinuos_flag, na.rm=TRUE)
  # sum(flu_loc$enddate_flag, na.rm=TRUE)
  # sum(flu_loc$no_time_flag, na.rm=TRUE)
  # unique(flu_loc$MiChartDeptName)
  # test <- flu_loc[flu_loc$enddate_flag==1,]
  
  #remove rows that 1) start after discharge, have unknown location and 
  #end date is before start date OR 2) the start and end dates or equal
  flu_loc <- flu_loc %>%
    filter(
      !(StartDate>DischargeDate & MiChartDeptName=="UNKNOWN" & end_before_start_flag==1) &
        no_time_flag != 1
    ) %>%
    group_by(PatientID,AdmitDate) %>%
    mutate(before_admit_flag = ifelse(StartDate<AdmitDate,1,0),
           after_disch_flag = ifelse(EndDate>DischargeDate,1,0),
           end_before_start_flag = ifelse(EndDate<StartDate,1,0),
           overlap_flag = case_when(StartDate < lag(EndDate) ~ 1,
                                    lead(StartDate) < EndDate~ 1,
                                    TRUE ~ 0),
           discontinuos_flag = ifelse(EndDate != lead(StartDate),1,0),
           enddate_flag = ifelse(lead(EndDate) <= EndDate,1,0),
           no_time_flag = ifelse(StartDate==EndDate,1,0)) %>%
    ungroup()
  
  #check test variables
  # sum(flu_loc$before_admit_flag, na.rm=TRUE)
  # sum(flu_loc$after_disch_flag, na.rm=TRUE)
  # sum(flu_loc$end_before_start_flag, na.rm=TRUE)
  # sum(flu_loc$overlap_flag, na.rm=TRUE)
  # sum(flu_loc$discontinuos_flag, na.rm=TRUE)
  # sum(flu_loc$enddate_flag, na.rm=TRUE)
  # sum(flu_loc$no_time_flag, na.rm=TRUE)
  # sort(unique(flu_loc$MiChartDeptName))
  # test <- flu_loc[flu_loc$no_time_flag==1,]
  
  #simplify unit names, and select final variables
  flu_loc <- flu_loc %>%
    mutate(Unit = case_when(
      MiChartDeptName == "UNKNOWN" ~ "UNKNOWN",
      
      MiChartDeptName == "Mott Childrens Emergency Department" ~ "Pediatric ED",
      
      MiChartDeptName == "7 East - Ped BMT" ~ "Pediatric Unit 1",
      MiChartDeptName == "MH 10W PCTU" ~ "Pediatric Unit 2",
      MiChartDeptName == "MH 11W GEN CARE" ~ "Pediatric Unit 3",
      MiChartDeptName == "MH 12E GEN CARE" ~ "Pediatric Unit 4",
      MiChartDeptName == "MH 12W GEN CARE" ~ "Pediatric Unit 5",
      
      MiChartDeptName == "MH 10E PED ICU" ~ "Pediatric ICU",
      
      MiChartDeptName == "Mott Operating Room" ~ "Pediatric OR",
      
      MiChartDeptName == "UH Adult Emergency Department" ~ "Adult ED",
      
      MiChartDeptName == "UH MSSU Maize-Taubman" ~ "Adult Short Stay 1",
      MiChartDeptName == "UH MSSU Blue-South" ~ "Adult Short Stay 2",
      
      MiChartDeptName == "UH 4A" ~ "Adult Unit 1",
      MiChartDeptName == "UH 4B SURG ACUTE CARE" ~ "Adult Unit 2",
      MiChartDeptName == "UH 4C SURG ACUTE CARE" ~ "Adult Unit 3",
      MiChartDeptName == "UH 5A SURG ACUTE CARE"~"Adult Unit 4",
      MiChartDeptName == "UH 5B1 Med Acute Care" ~ "Adult Unit 5",
      MiChartDeptName == "UH 5C SURG ACUTE CARE" ~ "Adult Unit 6",
      MiChartDeptName == "UH 6A PHYS MED REHAB" ~ "Adult Unit 7",
      MiChartDeptName == "UH 6B1 Med Acute Care" ~ "Adult Unit 8",
      MiChartDeptName == "UH 6C MED ACUTE CARE" ~ "Adult Unit 9",
      MiChartDeptName == "UH 7A1 GEN MED" ~ "Adult Unit 10",
      MiChartDeptName == "UH 7B MED ACUTE CARE" ~ "Adult Unit 11",
      MiChartDeptName == "UH 7C MED ACUTE CARE" ~ "Adult Unit 12",
      MiChartDeptName == "UH 8A3 Onc Acute Care" ~ "Adult Unit 13",
      MiChartDeptName == "UH 8B2 Med Acute Care" ~ "Adult Unit 14",
      MiChartDeptName == "UH 8C SURG ACUTE CARE" ~ "Adult Unit 15",
      MiChartDeptName == "MH 8E ACUTE CARE" ~ "Adult Unit 16",
      MiChartDeptName == "UH BAC BURN ACUTE" ~ "Adult Unit 17",
      MiChartDeptName == "7 West - Adult BMT" ~ "Adult Unit 18",
      MiChartDeptName == "U of M Von Voigtlander Women's Hospital" ~ "Adult Unit 19",
      MiChartDeptName == "U of M Women's Triage Center" ~ "Adult Unit 20",
      MiChartDeptName == "CVC 5 CARDIO MOD CRE" ~ "Adult Unit 21",
      MiChartDeptName == "MH 12E RICU" ~ "Adult Unit 22",
      
      MiChartDeptName == "UH 4DNI NEUROSURG ICU" ~ "Adult ICU 1",
      MiChartDeptName == "UH 5D SURG ICU" ~ "Adult ICU 2",
      MiChartDeptName == "UH 6D CRIT CARE ICU" ~ "Adult ICU 3",
      MiChartDeptName == "UH 7DN CORONARY CARE" ~ "Adult ICU 4",
      MiChartDeptName == "UH 8DNS" ~ "Adult ICU 5",
      MiChartDeptName == "UH BICU BURN ACUTE ICU" ~ "Adult ICU 6",
      MiChartDeptName == "CVC 4 CARDIOVASC ICU" ~ "Adult ICU 7",
      MiChartDeptName == "CVC-5 CARDIOVASC ICU" ~ "Adult ICU 8",
      
      MiChartDeptName == "UH OPERATING ROOMS" ~ "Adult OR 1",
      MiChartDeptName == "CVC Operating Room" ~ "Adult OR 2"
    )) %>%
    select(PatientID, AdmitDate, DischargeDate, StartDate, EndDate, MiChartDeptName, Unit)
  
  #sum(is.na(flu_loc$Unit))
  return(flu_loc)
}

# =========== Read Data ==================
filepath_1718 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 17-18/"
filepath_1920 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 19-20/"

#patient transfer data
loc_1718 <- fread(
  paste0(filepath_1718, "all_hosp_1718_encLoc.csv"), header = TRUE
)

loc_1920 <- rbind(
  fread(paste0(filepath_1920, "all_hosp_1920_encLoc1.csv"), header = TRUE), 
  fread(paste0(filepath_1920, "all_hosp_1920_encLoc2.csv"), header = TRUE)) %>%
  select(-V1)

#inpatient flu cases in analysis set
flu_1718 <- fread(paste0(filepath_1718,"derived/flu_enc_1718.csv"), 
                  stringsAsFactors = FALSE) %>%
  arrange(COLLECTION_DATE,AdmitDate) %>%
  filter(analysis_set==1)

flu_1920 <- fread(paste0(filepath_1920,"derived/flu_enc_1920.csv"), 
                  stringsAsFactors = FALSE) %>%
  arrange(COLLECTION_DATE,AdmitDate) %>%
  filter(analysis_set==1)

# =========== Apply function ==================
data_1718 <- collapse_transfers(loc_1718,flu_1718)
data_1920 <- collapse_transfers(loc_1920,flu_1920)

# =========== Write data ==================
fwrite(data_1718,
       paste0(filepath_1718, "derived/transfers_long_1718.csv"), 
       row.names = FALSE)
fwrite(data_1920,
       paste0(filepath_1920, "derived/transfers_long_1920.csv"), 
       row.names = FALSE)
