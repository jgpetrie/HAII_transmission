###Project: HARVI Genomic Analysis
###Purpose: Patient movement timeline plot
###Author: Josh Petrie, Tiffany Wan
###Date: 9/12/2022

###Example: https://urldefense.com/v3/__https://github.com/wlhamilton/Patient-ward-movement-timelines__;!!J7HzeEKFbK9hUUY!aKuXX_tb1zJbrFs1Twtelo7iQN4RL7-F5CIfWGXDwgVtFEC0rrQke37ObqLGLWa32Do6zNvbzjvK$ 

# =========== Load libraries ==================
library(data.table)
library(tidyverse)
library(vistime)
library(ggpubr)
library(scales)
library(cowplot)

# =========== Read Data ==================
filepath_1718 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 17-18/"
outpath <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/RESULTS/Josh/TiffanyWan_ILE/"


clusters <- fread(paste0(filepath_1718,"derived/analysis_1718.csv"), 
                       stringsAsFactors = FALSE) 

# =========== Prepare Data ==================
dates <- clusters %>%
  filter(cluster %in% c(1,2)) %>%
  mutate(sequenced = ifelse(is.na(sequenced),0,sequenced)) %>%
  arrange(cluster,ID) %>%
  select(ID,cluster,COLLECTION_DATE,sequenced) 

clusters <- clusters %>%
  filter(cluster %in% c(1,2) &
           !duplicated(clusters[,c("ID","StartDate_1")])) %>%
  mutate(status = case_when(fluA72==1 ~ "**",
                        fluA48==1 ~ "*",
                        TRUE ~ ""),
         ID_new = paste0(ID,status)) %>%
  select(ID,ID_new,cluster,starts_with("StartDate"),
         starts_with("EndDate"),starts_with("Unit")) %>%
  arrange(cluster,ID) 

clusters <- reshape(clusters,
                     direction = "long",
                     varying = list(grep("StartDate", colnames(clusters), value=T),
                                    grep("EndDate", colnames(clusters), value=T),
                                    grep("Unit", colnames(clusters), value=T)),
                     timevar = "row",
                     times = seq(1:20),
                     v.names = c("StartDate","EndDate","Unit"),
                     idvar = "ID") %>%
  filter(!is.na(StartDate) & !is.na(EndDate) & !is.na(Unit)) %>%
  arrange(cluster,ID)



units <- unique(clusters$Unit)

colors <- c("#000000",#black
            "#7F7F7FFF",#gray
            "#1F77B4FF",#blue
            "#AEC7E8FF",#light gray
            "#FF7F0EFF",#orange
            "#2CA02CFF",#green
            "#FF9896FF",#light red
            "#D62728FF",#red
            "#FFBB78FF",#light orange
            "#98DF8AFF",#light green
            "#C5B0D5FF",#light purple
            "#17BECFFF",#teal
            "#9467BDFF"#purple
)

# Create mapping of colors to units
color_mapping <- data.frame(
  Unit=units, 
  color=colors)

# merge in the mapping to the df
clusters <- merge(clusters,
                   color_mapping,
                   by="Unit",
                   all.x=T,all.y=T) %>%
  arrange(cluster,ID) %>%
  mutate(new_unit = case_when(Unit=="UNKNOWN" ~ "UNKNOWN",
                              Unit=="Adult ED" ~ "Adult ED",
                              Unit=="Adult Short Stay 2" ~ "Adult Short Stay",
                              Unit=="Adult Unit 3" ~ "Adult Unit 1",
                              Unit=="Adult Unit 9" ~ "Adult Unit 2",
                              Unit=="Adult Unit 10" ~ "Adult Unit 3",
                              Unit=="Adult Unit 11" ~ "Adult Unit 4",
                              Unit=="Adult Unit 12" ~ "Adult Unit 5",
                              Unit=="Adult Unit 14" ~ "Adult Unit 6",
                              Unit=="Adult ICU 2" ~ "Adult ICU 1",
                              Unit=="Adult ICU 3" ~ "Adult ICU 2",
                              Unit=="Adult ICU 5" ~ "Adult ICU 3",
                              Unit=="Adult ICU 7" ~ "Adult ICU 4"),
         new_unit = factor(new_unit, levels = c(
           "UNKNOWN","Adult ED","Adult Short Stay","Adult Unit 1","Adult Unit 2",
           "Adult Unit 3","Adult Unit 4","Adult Unit 5","Adult Unit 6",
           "Adult ICU 1","Adult ICU 2","Adult ICU 3","Adult ICU 4"
         )))


# =========== Plot ==================
# Produce the basic plot 
fig5 <- gg_vistime(data = clusters,
                   col.group = "ID_new", # Each row will be a patient
                   col.event = "new_unit", # Rows will be coloured by the ward
                   col.start = "StartDate",
                   col.end = "EndDate",
                   show_labels = FALSE, # Remove labels indicating the ward
                   linewidth = 10,
                   background_lines = 0)

# Adjust plot formatting
fig5  <- fig5 + theme_classic() +
  ggplot2::theme(
    plot.title = element_text(size=16),
    axis.text.x = element_text(size = 12, color = "black", 
                               angle = 30, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 12, color = "black")) +
  scale_x_datetime(breaks = breaks_width("7 days"), labels = date_format("%b %d, %Y"))

# Adding date of positive swab
sz1 = 2
sz2 = 1
fig5 <- fig5 +
  geom_text(aes(x=min(clusters$StartDate[clusters$cluster==1]) + 2*(86400),
                y=43,
                label="Group 1"), 
            color="black", 
            size=8, fontface="bold") +
  geom_text(aes(x=min(clusters$StartDate[clusters$cluster==2]) + 2*(86400),
                y=2,
                label="Group 2"), 
            color="black", 
            size=8, fontface="bold") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[1]), 
           y = 43, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[2]), 
           y = 41, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[2]), 
           y = 41, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[3]), 
           y = 41, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[4]), 
           y = 39, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[4]), 
           y = 39, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[5]), 
           y = 37, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[5]), 
           y = 37, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[6]), 
           y = 35, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[7]), 
           y = 33, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[8]), 
           y = 33, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[9]), 
           y = 33, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[10]), 
           y = 31, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[11]), 
           y = 29, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[11]), 
           y = 29, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[12]), 
           y = 29, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[13]), 
           y = 29, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[14]), 
           y = 27, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[14]), 
           y = 27, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[15]), 
           y = 25, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[16]), 
           y = 23, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[16]), 
           y = 23, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[17]), 
           y = 21, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[17]), 
           y = 21, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[18]), 
           y = 19, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[18]), 
           y = 19, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[19]), 
           y = 17, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[19]), 
           y = 17, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[20]), 
           y = 15, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[20]), 
           y = 15, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[21]), 
           y = 13, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[21]), 
           y = 13, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[22]), 
           y = 11, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[22]), 
           y = 11, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[23]), 
           y = 9, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[23]), 
           y = 9, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[24]), 
           y = 7, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[24]), 
           y = 7, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[25]), 
           y = 5, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[25]), 
           y = 5, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[27]), 
           y = 5, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[26]), 
           y = 5, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[26]), 
           y = 5, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[28]), 
           y = 3, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[28]), 
           y = 3, size = sz2, colour = "black") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[29]), 
           y = 1, size = sz1, colour = "white") +
  annotate("point", x = as.POSIXct(dates$COLLECTION_DATE[29]), 
           y = 1, size = sz2, colour = "black") +
  geom_rect(xmin = min(clusters$StartDate[clusters$cluster==1]) - 86400,
            xmax = max(clusters$EndDate[clusters$cluster==1]) + 86400,
            ymin = 36, ymax = 44,
            fill = NA, color = "black", size = 2) +
  geom_rect(xmin = min(clusters$StartDate[clusters$cluster==2]) - 86400,
            xmax = max(clusters$EndDate[clusters$cluster==2]) + 86400,
            ymin = 0, ymax = 36,
            fill = NA, color = "black", size = 2)

### Create a legend
#create dummy plot; specify your x and y axes and the group used to create the fill
#scale_fill_manual allows you to customize the legend
dat <- select(clusters, StartDate, ID, new_unit) 
dummy <- ggplot(data=dat, aes(ID,StartDate,fill = new_unit)) + 
  geom_bar(stat="identity") +
  theme(legend.position = "bottom") + 
  scale_fill_manual(
    name="",
    values=c("#000000",#black, UNKNOWN
             "#7F7F7FFF",#gray, Adult ED
             "#AEC7E8FF",#light blue, Adult Short Stay
             "#2CA02CFF",#green, Adult Unit 1
             "#98DF8AFF",#light green, Adult Unit 2
             "#9467BDFF",#purple, Adult Unit 3
             "#C5B0D5FF",#light purple, Adult Unit 4
             "#1F77B4FF",#blue, Adult Unit 5
             "#17BECFFF",#teal, Adult Unit 6
             "#FF9896FF",#light red, Adult ICU 1
             "#D62728FF",#red, Adult ICU 2
             "#FFBB78FF",#light orange, Adult ICU 3
             "#FF7F0EFF" #orange, Adult ICU 4
             )
    )
#extract the legend from the dummy plot
leg <- get_legend(dummy)
#arrange them together
# it defaults to legend on top, so specify legend position if you want
fig5 <- ggarrange(fig5, legend.grob = leg, legend = "bottom" ) 


pdf(paste0(outpath,"figure5_",Sys.Date(),".pdf"), width = 15, height = 10)
fig5
dev.off()







