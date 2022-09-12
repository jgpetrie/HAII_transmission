###Project: HARVI Genomic Analysis
###Purpose: 
#Identify contacts of HA flu cases who 1) were in the hospital at the same time,
#2) were "infected" before the HA flu case, and 3) were in a unit at the same
#time as the HA flu case in the 4 days prior to when the HA flu case was infected
###Author: Tiffany Wan, Josh Petrie
###Date: 2022-09-12

# =========== Load Libraries ==================
library(tidyverse)
library(data.table)
library(sqldf)
library(lubridate)

# =========== Define Data Management Functions ==================
make_clusters <- function(flu,loc){
  ## subset of HA-Flu: only hospital acquired fluA48==1
  HA_Flu <-  flu %>%
    arrange(PatientID, ORDER_DATE) %>%
    filter(fluA48==1 & !duplicated(PatientID)) %>%
    select(PatientID, ORDER_DATE, AdmitDate, DischargeDate)
  
  #merge HA-Flu and loc
  HA_Flu_Units <- left_join(HA_Flu,loc,by=c("PatientID","AdmitDate","DischargeDate")) %>%
    filter(
      #keep only locations in the 4 days prior to test order
      (( ORDER_DATE >= StartDate & ORDER_DATE <= EndDate ) |
         ( ORDER_DATE-days(1) >= StartDate & ORDER_DATE-days(1) <= EndDate ) |
         ( ORDER_DATE-days(2) >= StartDate & ORDER_DATE-days(2) <= EndDate ) |
         ( ORDER_DATE-days(3) >= StartDate & ORDER_DATE-days(3) <= EndDate ) |
         ( ORDER_DATE-days(4) >= StartDate & ORDER_DATE-days(4) <= EndDate )) &
        
        #keep only known locations
        !is.na(Unit) &
        Unit != "UNKNOWN"
    )
  
  # make one row per admission per individual
  flu <- flu[!duplicated(flu[,c("PatientID","AdmitDate")]),] %>%
    # filter to analysis set
    filter(analysis_set==1) %>%
    # select needed columns
    select(
      PatientID, AdmitDate, DischargeDate,ORDER_DATE,fluA48
    ) %>%
    arrange(PatientID, ORDER_DATE) %>%
    #create new date variable representing start of influenza infection period:
    #admission for community acquired, order date for hospital acquired
    mutate(new_date = case_when(fluA48==1 ~ ORDER_DATE, 
                                fluA48==0 ~ AdmitDate),
           new_date = as.POSIXct(new_date,'%Y-%m-%d %H:%M:%S', tz = "utc"))
  
  #merge all flu and loc
  flu <- left_join(flu,loc,by=c("PatientID","AdmitDate","DischargeDate")) %>%
    filter(
      #keep only known locations
      !is.na(Unit) &
        Unit != "UNKNOWN"
    )
  
  
  #empty list to store contacts of ha_flu
  contacts <- vector(mode = "list", length = nrow(HA_Flu_Units)) 
  
  for(i in 1:nrow(HA_Flu_Units)){
    
    c_id <- HA_Flu_Units$PatientID[i] #current ha_flu id
    c_admit <- HA_Flu_Units$AdmitDate[i] #current ha_flu admit date
    c_disch <- HA_Flu_Units$DischargeDate[i] #current ha_flu discharge date
    c_order <- HA_Flu_Units$ORDER_DATE[i] #current ha_flu order date
    c_start <- HA_Flu_Units$StartDate[i] #current ha_flu order date
    c_end <- HA_Flu_Units$EndDate[i] #current ha_flu order date
    c_unit <- HA_Flu_Units$Unit[i] #current ha_flu unit
    
    tmp <- flu %>%
      filter( 
        
        #keep only contacts with overlapping hospital admissions
        (( c_admit >= AdmitDate & c_admit <= DischargeDate ) |
           ( c_admit <= AdmitDate & c_disch >= AdmitDate )) &
          
          #keep only contacts with admit/order dates before ha_flu order date
          ( new_date <= c_order ) &
          
          #keep only locations in the 4 days prior to test order
          (( c_order >= StartDate & c_order <= EndDate ) | 
             ( c_order-days(1) >= StartDate & c_order-days(1) <= EndDate ) |
             ( c_order-days(2) >= StartDate & c_order-days(2) <= EndDate ) |
             ( c_order-days(3) >= StartDate & c_order-days(3) <= EndDate ) |
             ( c_order-days(4) >= StartDate & c_order-days(4) <= EndDate )) &
          
          #keep only overlapping units
          Unit == c_unit &
          
          #keep only overlapping time in same unit
          (( c_start >= StartDate & c_start <= EndDate ) |
             ( c_start <= StartDate & c_end >= StartDate ))
      ) %>%
      mutate(ha_flu = ifelse(PatientID==c_id,1,0)) %>% 
      arrange(desc(ha_flu),PatientID,StartDate)
    
    #store contacts dataframe in list 
    contacts[[i]] <- tmp 
    
    #name list element with ha_flu id to help keep track
    names(contacts)[i] <- c_id
  }
  
  # bind rows of the list that have the same name i.e.
  # combine data of HA flu cases who were in more than one location
  contacts <- tapply(contacts, names(contacts), dplyr::bind_rows)
  
  #remove dataframes that contain only the HA flu case i.e. no contacts
  contacts <- keep(contacts, ~length(unique(.x$PatientID)) > 1)
  
  #Find overlapping clusters
  #stack all dataframes in contacts into single data frame, and create cluster 
  #variable that identifies each dataframe of contacts
  clusters <- bind_rows(contacts, .id = "cluster") %>%
    #for each individual...
    group_by(PatientID) %>%
    #...get their minimum cluster id
    mutate(cluster2 = min(as.numeric(cluster))) %>%
    ungroup %>%
    #for each original cluster...
    group_by(cluster) %>%
    #...make a new cluster id that is equal to their minimum cluster id
    mutate(cluster3 = min(cluster2)) %>%
    ungroup %>%
    #do it all again to capture longer linkages
    #e.g. cluster 3 overlaps clusters 1 and 9, but cluster 9 does not overlap 1.
    #cluster 3 would become 1, and cluster 9 would become 3, but all should be 1.
    #There is probably a better way to do this, but 2 iterations captures all 
    #scenarios in the current data.
    group_by(PatientID) %>%
    mutate(cluster4 = min(as.numeric(cluster3))) %>%
    ungroup %>%
    group_by(cluster3) %>%
    mutate(cluster5 = min(cluster4)) %>%
    ungroup 
  
  clusters$cluster5 <- factor(clusters$cluster5,
                             levels = unique(clusters$cluster5[
                               order(clusters$ORDER_DATE)
                                  ]
                               )
  )
  
  #after checking the above data, clean things up
  clusters <- clusters %>%
    group_by(cluster5) %>%
    #make new cluster ids sequential e.g. 1 thru 10 instead of 1, 4, 5, 13, etc.
    mutate(cluster = cur_group_id()) %>% 
    ungroup %>%
    #remove duplicate rows
    filter(!duplicated(clusters[,c("PatientID","Unit")])) %>%
    #sort by cluster and "infection" date
    arrange(cluster,new_date) %>%
    #remove intermediate columns
    select(-cluster2, -cluster3, -cluster4, -cluster5)
  
  return(clusters)
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

# =========== Apply function ==================
clusters_1718 <- make_clusters(flu_1718,loc_1718)
clusters_1920 <- make_clusters(flu_1920,loc_1920)

# =========== Write data ==================
fwrite(clusters_1718,
       paste0(filepath_1718, "derived/clusters_1718.csv"), 
       row.names = FALSE)
fwrite(clusters_1920,
       paste0(filepath_1920, "derived/clusters_1920.csv"), 
       row.names = FALSE)