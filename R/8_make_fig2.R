###Project: HARVI Genomic Analysis
###Purpose: Make figure 2 description of clusters
###Author: Tiffany Wan, Josh Petrie
###Date: 2022-09-12

# =========== Load Libraries ==================
library(tidyverse)
library(data.table)
library(vistime)
library(plotly)
library(ggsci)
library(scales)
library(gridExtra)

# =========== Define Functions ==================
#show_col(pal_lancet("lanonc")(9))
make_figure2 <- function(clusters, year, labelsize, dotsize){
  if(year=="17-18"){
    minDate <- as.POSIXct("2017-12-15")
    maxDate <- as.POSIXct("2018-04-30")
  } else if(year=="19-20"){
    minDate <- as.POSIXct("2019-12-15")
    maxDate <- as.POSIXct("2020-04-30")
  } 
  
  clusters <- clusters %>%
    filter(!is.na(cluster)) %>%
    select(ID,ORDER_DATE,cluster,fluA48,sequenced) %>%
    group_by(ID) %>%
    summarise(group = cluster[which.min(ORDER_DATE)],
              start = as.Date(min(ORDER_DATE)),
              end = start,
              ha_flu = sum(fluA48, na.rm = TRUE),
              sequenced = sum(sequenced, na.rm=TRUE)) %>%
    mutate(color = case_when(ha_flu==1 & sequenced >= 1 ~ "#D62728FF",
                              ha_flu==0 & sequenced >= 1 ~ "black",
                              ha_flu==1 & sequenced == 0 ~ "#FF9896FF",
                              ha_flu==0 & sequenced == 0 ~ "gray",)) %>%
    rename(event = ID) %>%
    select(event,group,start,end,color)
  

  f2 <- gg_vistime(clusters,
                   show_labels = FALSE)
  
  f2 <- f2 + 
    ggrepel::geom_text_repel(
      data=f2$layers[[3]]$data, 
      label=f2$layers[[3]]$data$label, 
      size = labelsize, 
      color="black",
      min.segment.length = Inf,
      max.overlaps = Inf) +
    theme_classic() +
    ggtitle(ifelse(year=="17-18","A","B")) +
    ggplot2::theme(
      axis.text.x = element_text(
        size = 12, color = "black", angle = 30, vjust = 1, hjust = 1
        ),
      axis.text.y = element_text(size = 12, color = "black")
    ) +
    scale_x_datetime(breaks = breaks_width("month"),
                     labels = date_format("%b %Y"),
                     limits = c(minDate,maxDate))
  
  f2$layers[[3]]$aes_params$size <- dotsize
  
  return(f2)
  # f2 <- vistime(clusters,
  #               title = ifelse(year=="17-18","A","B"))

  # transform into a list
  # pf2 <- plotly::plotly_build(f2)
  # 
  # pf2$x$layout$xaxis$tickfont <- list(size = labelsize)
  # pf2$x$layout$yaxis$tickfont <- list(size = labelsize)
  # 
  # pf2$x$data[[3]]$textfont$size <- labelsize
  # 
  # pf2$x$data[[2]]$marker$size <- dotsize
  # 
  # pf2 <- pf2 %>%
  #   layout(xaxis = list(dtick = "M1", 
  #                       tickformat="%b<br>%Y",
  #                       range = c(minDate,maxDate)))
  # 
  # return(pf2)
}

# =========== Set data path and read data ==================
filepath_1718 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 17-18/"
filepath_1920 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 19-20/"
outpath <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/RESULTS/Josh/TiffanyWan_ILE/"

clusters_1718 <- fread(paste0(filepath_1718,"derived/analysis_1718.csv"), 
                       stringsAsFactors = FALSE) 

clusters_1920 <- fread(paste0(filepath_1920,"derived/analysis_1920.csv"), 
                       stringsAsFactors = FALSE)

# =========== Figure 2: Cluster Summary ==================
fig2_1718 <- make_figure2(clusters_1718,"17-18",5,4)
fig2_1920 <- make_figure2(clusters_1920,"19-20",5,4)


pdf(paste0(outpath,"figure2_",Sys.Date(),".pdf"), width = 12, height = 16)
grid.arrange(fig2_1718,fig2_1920,ncol=1)
dev.off()