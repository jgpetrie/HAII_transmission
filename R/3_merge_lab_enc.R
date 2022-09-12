###Project: HARVI Genomic Analysis
###Purpose: 
#############1. merging lab and encounter data (merged) with location data -> epi data
#############2. linking epi data with sequence data and keeping the ones present in tree
###Author: Tiffany Wan, Josh Petrie
###Date: 2022-09-12

# =========== Load libraries ==================
library(data.table)
library(dplyr)
library(sqldf)
library(lubridate)


# =========== Define Data Management Function ==================
manage_data <- function(year,enc,flu,key,tree,resolved){
  # format date variables
  enc <- enc %>% 
    mutate(AdmitDate = as.POSIXct(AdmitDate,'%Y-%m-%d %H:%M', tz = "utc") ,
           DischargeDate = as.POSIXct(DischargeDate,'%Y-%m-%d %H:%M', tz = "utc"),
           AdmitDateminus1 = AdmitDate-86400,
           DischargeDateplus1 = DischargeDate+86400
    )
  
  flu <- flu %>% 
    mutate(COLLECTION_DATE = as.POSIXct(COLLECTION_DATE,'%m/%d/%Y %H:%M', tz = "utc"),
           ORDER_DATE = as.POSIXct(ORDER_DATE,'%m/%d/%Y %H:%M', tz = "utc")
    )
  
  resolved <- resolved %>% 
    mutate(COLLECTION_DATE = as.POSIXct(collectionDateTime,'%m/%d/%Y %H:%M', tz = "utc"),
           COLLECTION_DATE = ceiling_date(COLLECTION_DATE,"minute"),
           #individuals in this table had unsubtyped results resolved by sequencing
           resolved = 1
    ) %>%
    select(PatientID, COLLECTION_DATE,resolved)
  
  
  tree <- left_join(tree, key, by = c("ID" = "Sample.Number")) %>% 
    mutate(COLLECTION_DATE = as.POSIXct(collectionDateTime,'%m/%d/%Y %H:%M', tz = "utc"),
           COLLECTION_DATE = ceiling_date(COLLECTION_DATE,"minute"),
           #individuals in the tree file had specimens successfully sequenced
           sequenced = 1,
           SpecimenID = ID
    ) %>%
    select(PatientID, SpecimenID, COLLECTION_DATE,sequenced)
  
  key <- key %>% 
    mutate(COLLECTION_DATE = as.POSIXct(collectionDateTime,'%m/%d/%Y %H:%M', tz = "utc"),
           COLLECTION_DATE = ceiling_date(COLLECTION_DATE,"minute"),
           #individuals in this table had specimens retrieved
           retrieved = 1
    ) %>%
    select(PatientID, COLLECTION_DATE,retrieved)
  
  #remove duplicates
  key <- key[!duplicated(key[,c("PatientID","COLLECTION_DATE")]),] 
  
  ### merge encounter data to flu positive specimens with sqldf --left join
  data <- sqldf("SELECT * FROM flu
                  LEFT JOIN enc
                  ON  enc.PatientID = flu.PatientID
                  AND flu.COLLECTION_DATE BETWEEN enc.AdmitDateminus1 AND enc.DischargeDateplus1")
  
  # delete duplicated column (PatientID)
  data <- data[, !duplicated(colnames(data))]
  
  ### merge in indicator of whether specimens were sequenced
  data <- left_join(data,tree,by=c("PatientID","COLLECTION_DATE"))
  
  ### merge in indicator of whether specimens were subtyped by sequencing 
  # as subtype we were not interested in (17-18:H1, 19-20:H3)
  data <- left_join(data,resolved,by=c("PatientID","COLLECTION_DATE"))
  # check - should have be same as number of rows in resolved
  # sum(data$resolved,na.rm = T)
  
  ### merge in indicator of whether specimens were retrieved 
  data <- left_join(data,key,by=c("PatientID","COLLECTION_DATE"))
  # check - should have be same as number of rows in key
  # sum(data$retrieved,na.rm = T)
  
  ### create new flu status variables
  data <- data %>%
    mutate(
      #text variables for each possible flu infection type
      text_untyped = case_when(fluA==1 & is.na(fluAH1) & is.na(fluAH3) ~ "A untyped",
                               fluA==1 & fluAH1==0 & fluAH3==0 ~ "A untyped",
                               TRUE ~ NA_character_),
      text_AH3 = ifelse(fluA==1 & fluAH1==1, "AH1", NA_character_),
      text_AH1 = ifelse(fluA==1 & fluAH3==1, "AH3", NA_character_),
      text_B = ifelse(fluB==1, "B", NA_character_),
      
      #combine text variables in single string to identify possible co-infections
      flu_result = paste(text_untyped,text_AH1,text_AH3,text_B,sep = ","),
      #clean up result string to remove extra commas and NAs
      flu_result = gsub(",NA","",flu_result),
      flu_result = gsub("NA,","",flu_result),
      
      #recode results for A unsubtyped specimens that were resolved by sequencing
      flu_result = case_when(year=="17-18" & resolved==1 ~ "AH1",
                             year=="19-20" & resolved==1 ~ "AH3",
                             year=="17-18" & flu_result=="A untyped" & sequenced==1 ~ "AH3",
                             year=="19-20" & flu_result=="A untyped" & sequenced==1 ~ "AH1",
                             TRUE ~ flu_result),
      
      #flag specimens to be included in the analysis 
      #(17-18: A untyped & H3, 19-20: A untyped & H1)
      analysis_set = case_when(
        year=="17-18" & flu_result %in% c("A untyped","AH3") ~ 1,
        year=="19-20" & flu_result %in% c("A untyped","AH1","AH1,B") ~ 1,
        TRUE ~ 0
      )
    ) %>% 
    #clean up intermediate variables
    select(-starts_with("text"), -resolved)
  
  ## flag first flu positives
  data <- data %>%
    #sort by individual, collection date, and descending analysis flag
    arrange(PatientID, COLLECTION_DATE, desc(analysis_set)) %>%
    #we only really care about HA flu in the analysis set, so 
    #group by individual and analysis flag
    group_by(PatientID, analysis_set) %>%
    mutate(
      #number rows within each individual X analysis flag group
      num = row_number(),
      #first flu is the first row for an individual that is flagged 
      #in the analysis set. This is the earliest because we sorted above.
      firstflu = case_when(num==1 & analysis_set==1 ~ 1,
                           TRUE ~ 0)) %>%
    #clean up
    select(-num) %>%
    ungroup()
  
  #make final HA flu flags
  data <- data %>%
    mutate(
      #calculate differnce between admission and test order date
      hours= case_when(!is.na(AdmitDate) & !is.na(ORDER_DATE) ~ 
                         difftime(ORDER_DATE, AdmitDate, units = "hours")
      ),
      hours = as.numeric(hours),
      # flag first flu positives that were ordered gt 48 and 72 hrs after admit
      fluA48 = case_when(hours > 48 & firstflu==1 ~1,
                         is.na(firstflu)~NA_real_,
                         TRUE ~ 0),
      fluA72 = case_when(hours > 72 & firstflu==1 ~1,
                         is.na(firstflu)~NA_real_,
                         TRUE ~ 0),
    ) 
  
  #remove records with out of range dates
  if(year=="17-18"){
    data <- data %>%
      filter(
        AdmitDate >= as.Date("2017-07-01") &
          AdmitDate <= as.Date("2018-06-30") &
          ORDER_DATE >= as.Date("2017-07-01") &
          ORDER_DATE <= as.Date("2018-06-30")
        ) 
    
  } else if(year=="19-20"){
    data <- data %>%
      filter(
        AdmitDate >= as.Date("2019-07-01") &
          AdmitDate <= as.Date("2020-06-30") &
          ORDER_DATE >= as.Date("2019-07-01") &
          ORDER_DATE <= as.Date("2020-06-30")
      ) 
  }
  
  
  return(data)
}

# =========== Set data path and read data ==================
filepath_1718 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 17-18/"
filepath_1920 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 19-20/"

#collapsed encounter data
enc_1718 <- fread(
  paste0(filepath_1718, "derived/collapsed_encounters_1718.csv")
  ) 

#flu positives
flu_1718 <- fread(
  paste0(filepath_1718, "derived/flu_1718.csv")
  )

# link between sequence SpecimenIDs and Patient ID
key_1718 <- fread(
  paste0(filepath_1718, "specimen_PatientID_link.csv")
  ) 

# IDs from phylogenetic tree
tree_1718 <- fread(
  paste0(filepath_1718, "derived/H3N2_tree_id.csv")
)

# Individuals with unsubtyped influenza A resolved by sequencing
resolved_1718 <- fread(
  paste0(filepath_1718, "derived/1718_unsubtyped_H1N1.csv")
)

#collapsed encounter data
enc_1920 <- fread(
  paste0(filepath_1920, "derived/collapsed_encounters_1920.csv")
) 

#flu positives
flu_1920 <- fread(
  paste0(filepath_1920, "derived/flu_1920.csv")
)

# link between sequence SpecimenIDs and Patient ID
key_1920 <- fread(
  paste0(filepath_1920, "specimen_PatientID_link_1920.csv")
) %>%
  rename("Sample.Number" = "Sample.Name")

# IDs from phylogenetic tree
tree_1920 <- fread(
  paste0(filepath_1920, "derived/H1N1_tree_id.csv")
)

# Individuals with unsubtyped influenza A resolved by sequencing
resolved_1920 <- fread(
  paste0(filepath_1920, "derived/1920_unsubtyped_H3N2.csv")
)

# =========== Apply function ==================
data_1718 <- manage_data("17-18",enc_1718,flu_1718,key_1718,tree_1718,resolved_1718)
data_1920 <- manage_data("19-20",enc_1920,flu_1920,key_1920,tree_1920,resolved_1920)

# =========== Write data ==================
fwrite(data_1718,
       paste0(filepath_1718, "derived/flu_enc_1718.csv"), 
       row.names = FALSE)
fwrite(data_1920,
       paste0(filepath_1920, "derived/flu_enc_1920.csv"), 
       row.names = FALSE)