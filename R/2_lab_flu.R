###Project: HARVI Genomic Analysis 
###Purpose: Subset Respiratory Virus Laboratory Results to Get Flu Positives
###Author: Tiffany Wan, Josh Petrie
###Date: 2022-09-12

# =========== Load Libraries ==================
library(data.table)
library(dplyr)

# =========== Set data path and read data and subset ==================
filepath_1718 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 17-18/"
filepath_1920 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 19-20/"

data_1718 <- fread(paste0(filepath_1718, "respiratory_virus_lab_results_long.csv")) %>%
  filter(fluA == 1 | fluB == 1) %>%
  select(PatientID,ORDER_DATE,COLLECTION_DATE,fluA,fluB,fluAH1,fluAH3)

data_1920 <- fread(paste0(filepath_1920, "respiratory_virus_lab_results_long.csv")) %>%
  filter(fluA == 1 | fluB == 1) %>%
  select(PatientID,ORDER_DATE,COLLECTION_DATE,fluA,fluB,fluAH1,fluAH3)

# =========== Write data ==================
fwrite(data_1718,
       paste0(filepath_1718, "derived/flu_1718.csv"), 
       row.names = FALSE)
fwrite(data_1920,
       paste0(filepath_1920, "derived/flu_1920.csv"), 
       row.names = FALSE)