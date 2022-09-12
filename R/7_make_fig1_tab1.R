###Project: HARVI Genomic Analysis
###Purpose: Make figure 1 epicurve and table 1
###Author: Tiffany Wan, Josh Petrie
###Date: 2022-09-12

# =========== Load Libraries ==================
library(tidyverse)
library(data.table)
library(lubridate)
library(gridExtra)
library(gtsummary)

# =========== Define Functions ==================

plot_epicurve <- function(data,year){
  ###plot epicurve
  epicurve <- data %>%
    mutate(week= ceiling_date(COLLECTION_DATE,"week")) %>%
    group_by(week) %>%
    summarise(
      A_untyped = sum(flu_result=="A untyped"),
      AH1 = sum(flu_result %in% c("AH1","AH1,B")),
      AH3 = sum(flu_result=="AH3"),
      B = sum(flu_result=="B")
    ) %>%
    ungroup()
  
  if(year=="17-18"){
    epicurve_start <- as.Date("2017-07-01")
    epicurve_end <- as.Date("2018-06-30")
  } else if(year=="19-20"){
    epicurve_start <- as.Date("2019-07-01")
    epicurve_end <- as.Date("2020-06-30")
  }
  
  full_weeks <- data.frame(
    week = ceiling_date(
      seq.Date(epicurve_start,epicurve_end,by="week"), "week"
    )
  )
  
  epicurve <- full_join(full_weeks,epicurve) %>%
    replace(is.na(.), 0) %>%
    gather("A_untyped","AH1","AH3","B", 
           key = flu_type, value = count)
  
  p <- ggplot(epicurve, aes( x=as.Date(week), y=count, fill=flu_type)) + 
    geom_bar(stat="identity",position="stack", color = "white", size = 0.25) +
    labs(x="Week Ending", y = "Count") +
    ggtitle(ifelse(year=="17-18","A","B")) +
    scale_y_continuous(breaks = seq(0,40,by=10), limits = c(0,40)) +
    scale_x_date(date_breaks = "2 month",date_labels = "%b %Y") + 
    scale_fill_manual(name = "Influenza", 
                      labels = c("A untyped", "A (H1N1)", "A (H3N2)","B"),
                      values = c("#1F77B4FF","#D62728FF","#2CA02CFF","#17BECFFF")) +
    theme_classic()
  
  return(p)
}

make_table1 <- function(data){
  data <- data %>%
    filter(analysis_set==1 & !duplicated(PatientID)) %>%
    mutate(ageGroup = factor(case_when(!is.na(AgeInWeeks) ~ "<18",
                                       !is.na(AgeInMonths) ~ "<18",
                                       AgeInYears<18 ~ "<18",
                                       AgeInYears<65 ~ "18-64",
                                       AgeInYears>=65 ~ ">=65"),
                             levels = c("<18","18-64",">=65")),
           month = factor(format(ORDER_DATE, "%B"),
                          levels = c(
                            "July","August","September","October","November",
                            "December","January","February","March","April","May",
                            "June"
                          ))
    ) %>%
    group_by(PatientID) %>%
    summarise(n_retrieved = sum(retrieved,na.rm = TRUE),
              n_sequenced = sum(sequenced,na.rm = TRUE),
              n_ha = sum(fluA48,na.rm = TRUE),
              month = month[which.min(ORDER_DATE)],
              ageGroup = ageGroup[which.min(ORDER_DATE)],
              GenderName = GenderName[which.min(ORDER_DATE)]
    ) %>%
    mutate(
      source = factor(
        ifelse(n_ha >= 1, "Hospital-acquired", "Community-acquired"),
        levels = c("Community-acquired", "Hospital-acquired")
      ),
      sequenced_label = factor(
        ifelse(n_sequenced >= 1,"Sequenced","Not Sequenced"),
        levels = c("Sequenced","Not Sequenced")
      ),
      retrieved_label = factor(
        ifelse(n_retrieved >= 1,"Retrieved","Not Retrieved"),
        levels = c("Retrieved","Not Retrieved")
      )
    ) 
  
  t1a <- tbl_summary(
    data %>% select(GenderName, ageGroup, source, month),
    digits = list(all_categorical() ~ c(0, 1))
  )
  t1b <- tbl_summary(
    data %>% 
      filter(retrieved_label=="Not Retrieved") %>% 
      select(GenderName, ageGroup, source, month),
    digits = list(all_categorical() ~ c(0, 1))
  ) 
  t1c <- tbl_summary(
    data %>% 
      filter(retrieved_label=="Retrieved") %>% 
      select(GenderName, ageGroup, source, month, sequenced_label),
    digits = list(all_categorical() ~ c(0, 1)),
    by = sequenced_label) %>%
    add_p()
  
  t1 <- tbl_merge(list(t1a,t1b,t1c)) %>%
    modify_spanning_header(
      list(
        stat_0_1 ~ "**Total**",
        stat_0_2 ~ "**Not Retrieved**",
        c(stat_1_3,stat_2_3,p.value_3) ~ "**Retrieved**"
      )
    )
  return(t1)
}

# =========== Set data path and read data ==================
filepath_1718 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 17-18/"
filepath_1920 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 19-20/"
outpath <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/RESULTS/Josh/TiffanyWan_ILE/"

analysis_1718 <- fread(paste0(filepath_1718,"derived/analysis_1718.csv"), 
                       stringsAsFactors = FALSE) 

analysis_1920 <- fread(paste0(filepath_1920,"derived/analysis_1920.csv"), 
                       stringsAsFactors = FALSE) 

# =========== Figure 1: Epi Curve ==================
fig1_1718 <- plot_epicurve(analysis_1718,"17-18")
fig1_1920 <- plot_epicurve(analysis_1920,"19-20")


pdf(paste0(outpath,"figure1_",Sys.Date(),".pdf"), width = 6, height = 8)
grid.arrange(fig1_1718,fig1_1920,ncol=1)
dev.off()


# =========== Table 1 ==================
t1_1718 <- make_table1(analysis_1718)
t1_1920 <- make_table1(analysis_1920)

t1 <- tbl_merge(list(t1_1718,t1_1920))

t1 <- as_tibble(t1)

fwrite(t1,paste0(outpath,"table1_",Sys.Date(),".csv"))


# =========== Results Paragraph 1 ==================
#375 influenza positive specimens 17-18
nrow(analysis_1718) 
#339 individuals 17-18
nrow(analysis_1718[!duplicated(analysis_1718$PatientID)]) 
#number of flu+ specimens by type/subtype
table(analysis_1718$flu_result) 
#number of flu+ individuals by type/subtype
table(analysis_1718$flu_result[!duplicated(analysis_1718$PatientID)])

table(analysis_1718$fluA48)
table(analysis_1718$fluA48[!duplicated(analysis_1718$PatientID)])
table(analysis_1718$fluA48,analysis_1718$flu_result)
table(analysis_1718$fluA72)

table(
  analysis_1718$retrieved[
    analysis_1718$analysis_set==1 & !duplicated(analysis_1718$PatientID)
    ], useNA = "always"
  ) 

table(
  analysis_1718$retrieved[
    analysis_1718$analysis_set==1 & 
      !duplicated(analysis_1718$PatientID) &
      analysis_1718$fluA48==1
  ]
) 

table(
  analysis_1718$sequenced[
    analysis_1718$analysis_set==1 & !duplicated(analysis_1718$PatientID)
  ]
) 

table(
  analysis_1718$sequenced[
    analysis_1718$analysis_set==1 & 
      !duplicated(analysis_1718$PatientID) &
      analysis_1718$fluA48==1
  ]
) 

#283 influenza positive specimens 19-20
nrow(analysis_1920) 
#263 individuals 19-20
nrow(analysis_1920[!duplicated(analysis_1920$PatientID)]) 
#number of flu+ specimens by type/subtype
table(analysis_1920$flu_result) 
#number of flu+ individuals by type/subtype
table(analysis_1920$flu_result[!duplicated(analysis_1920$PatientID)])

table(analysis_1920$fluA48)
table(analysis_1920$fluA48[!duplicated(analysis_1920$PatientID)])
table(analysis_1920$fluA48,analysis_1920$flu_result)
table(analysis_1920$fluA72)

table(
  analysis_1920$retrieved[
    analysis_1920$analysis_set==1 & !duplicated(analysis_1920$PatientID)
  ]
) 

table(
  analysis_1920$retrieved[
    analysis_1920$analysis_set==1 & 
      !duplicated(analysis_1920$PatientID) &
      analysis_1920$fluA48==1
  ]
) 

table(
  analysis_1920$sequenced[
    analysis_1920$analysis_set==1 & !duplicated(analysis_1920$PatientID)
  ]
) 

table(
  analysis_1920$sequenced[
    analysis_1920$analysis_set==1 & 
      !duplicated(analysis_1920$PatientID) &
      analysis_1920$fluA48==1
  ]
) 


