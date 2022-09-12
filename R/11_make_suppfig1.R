###Project: HARVI Genomic Analysis
###Purpose: Make supplemental figure 1 pairwise distance
###Author: Tiffany Wan, Josh Petrie
###Date: 2022-09-12

# =========== Load Libraries ==================
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)

# =========== Define Functions ==================

make_pairwise_plot <- function(meta,dist,year){
  names(dist) <- c("SpecimenID1", "SpecimenID2", "distance")
  
  dist <- dist %>%
    #remove distances to reference strains
    filter(!grepl("A_",SpecimenID1) & !grepl("A_",SpecimenID2)) %>%
    #make SpecimenIDs numeric
    mutate(SpecimenID1 = as.integer(SpecimenID1),
           SpecimenID2 = as.integer(SpecimenID2))
  
  #merge to meta data
  dist <- left_join(
    dist,
    meta %>%
      mutate(SpecimenID1 = SpecimenID,
             ID1 = ID,
             cluster1 = cluster) %>%
      select(SpecimenID1, ID1, cluster1),
    by = "SpecimenID1"
  )
  
  dist <- left_join(
    dist,
    meta %>%
      mutate(SpecimenID2 = SpecimenID,
             ID2 = ID,
             cluster2 = cluster) %>%
      select(SpecimenID2, ID2, cluster2),
    by = "SpecimenID2"
  )
  
  
  
  dist <- dist %>%
    #remove rows that aren't in meta data
    filter(!is.na(ID1) & !is.na(ID2)) %>%
    #remove rows that are equivalent (e.g. 1-2, 2-1)
    filter(
      !duplicated(paste0(pmax(SpecimenID1, SpecimenID2), 
                         pmin(SpecimenID1, SpecimenID2)))
      ) %>%
    #remove distances from self
    filter(SpecimenID1 != SpecimenID2) %>%
    mutate(linked = case_when(is.na(cluster1) | is.na(cluster2) ~ "No",
                              cluster1 == cluster2 ~ "Yes",
                              cluster1 != cluster2 ~ "No"))
  
  quantile(dist$distance[dist$linked=="No"], prob = 0.001)
  
  
  # plot clusters and pairwise distance
  sf1 <- ggplot(dist, aes(x= distance, y =..ncount..,fill=as.factor(linked))) + 
    geom_histogram(binwidth = 2, boundary = 0,position = "dodge", alpha = 0.6) + 
    geom_vline(xintercept=3, linetype="dashed", color = "black", size=0.5) +
    theme_classic() +
    ggtitle(ifelse(year=="17-18","A","B")) +
    xlab("Pairwise Genetic Distance") +
    ylab("Normalized Count") +
    scale_fill_manual(
      name = 'In Time-\nLocation\nGroup', 
      values = c("gray50", "red")
      )
  
  return(sf1)
}

# =========== Set data path and read data ==================
filepath_1718 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 17-18/"
filepath_1920 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 19-20/"
outpath <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/RESULTS/Josh/TiffanyWan_ILE/"

meta_1718 <- fread(paste0(filepath_1718,"derived/analysis_1718.csv"), 
                   stringsAsFactors = FALSE) 

meta_1920 <- fread(paste0(filepath_1920,"derived/analysis_1920.csv"), 
                   stringsAsFactors = FALSE) 

dist_1718 <- fread(
  paste0(filepath_1718,"derived/H3N2_pairwise.csv")
) 

dist_1920 <- fread(
  paste0(filepath_1920,"derived/H1N1_pairwise.csv")
) 

# =========== Supplemental Figure 1: Pairwise Distance ==================
sfig1_1718 <- make_pairwise_plot(meta_1718,dist_1718,"17-18")
sfig1_1920 <- make_pairwise_plot(meta_1920,dist_1920,"19-20")


pdf(paste0(outpath,"supplemental_figure1_",Sys.Date(),".pdf"), width = 8, height = 12)
grid.arrange(sfig1_1718,sfig1_1920,ncol=1)
dev.off()

