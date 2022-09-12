### Project: HARVI Genomic Analysis
### Purpose: Plot trees based on consensus sequences and make table 2.
###Author: Tiffany Wan, Josh Petrie
###Date: 2022-09-12

# =========== Load Libraries ==================
library(dplyr)
library(data.table)
library(ggtree)
library(treeio)
library(gtsummary)

# =========== Define Functions ==================
plot_tree <- function(data,tree,year){
  #this saves just the tip labels from your tree into a new dataframe
  #they are in the same order as in the tree file
  tip_labels <- data.frame("id" = tree$tip.label) %>%
    #number the rows so you can get back to this same sort order
    mutate(sort = row_number())
  
  data <- data %>%
    arrange(ID, COLLECTION_DATE)
  
  #find tips that aren't in metadata e.g. they weren't in the inpatient emr data
  #also exclude 133743 - second specimen from patient 136 with read errors
  exclusions <- tip_labels %>%
    filter(
      (!(id %in% unique(data$SpecimenID)) & !grepl("A_|Sing|Hawaii",id)) |
        id == 133743 
      ) 
  
  #drop tip
  tree <- drop.tip(tree, exclusions$id)
  
  #remake tip_labels after exclusions
  tip_labels <- data.frame("id" = tree$tip.label) %>%
    #number the rows so you can get back to this same sort order
    mutate(sort = row_number())
  
  # annotate tree
  data <- data %>%
    group_by(ID) %>%
    mutate(id = as.character(SpecimenID),
           fluA48 = max(fluA48, na.rm = TRUE),
           fluA72 = max(fluA72, na.rm = TRUE),
           status = ifelse(fluA72==1, "**", 
                           ifelse(fluA48==1, "*","")),
           label = paste0(ID," | ",as.Date(COLLECTION_DATE), status )) %>% 
    ungroup() %>%
    select(id, label, cluster)
  
  #merge your new labels into the original labels
  data <- left_join(tip_labels, data, by = "id") %>%
    #copy reference strain labels over to the new label variable
    mutate(label = case_when(
      is.na(label)&id=="A_Texas"~ "A/Texas/50/2012",
      is.na(label)&id=="A_Maryland" ~ "A/Maryland/11/2018",
      is.na(label)&id=="A_Washington" ~ "A/Washington/16/2017",
      is.na(label)&id=="A_Singapore" ~ "A/Singapore/infihm-16-0019/2016",
      is.na(label)&id=="A_California" ~ "A/California/7/2009",
      is.na(label)&id=="A_Wisconsin" ~ "A/Wisconsin/588/2019",
      is.na(label)&id=="A_Michigan" ~ "A/Michigan/45/2015",
      is.na(label)&id=="A_Hawaii" ~ "A/Hawaii/70/2019",
      is.na(cluster) ~ "",
      !is.na(label)~label
    )) %>%
    #arrange to get back to tree order
    arrange(sort)
  
  #save as new dataset for coloring tree plot
  data <- data %>%
    mutate(cluster = factor(
      case_when(
        label == "A/Texas/50/2012" ~ "Reference",
        label == "A/Maryland/11/2018" ~ "Reference",
        label == "A/Washington/16/2017" ~ "Reference",
        label == "A/Singapore/infihm-16-0019/2016" ~ "Reference",
        label == "A/California/7/2009" ~ "Reference",
        label == "A/Wisconsin/588/2019" ~ "Reference",
        label == "A/Michigan/45/2015" ~ "Reference",
        label == "A/Hawaii/70/2019" ~ "Reference",
        is.na(cluster) ~ "None",
        TRUE ~ as.character(cluster)),
      levels = c("None","Reference",sort(unique(data$cluster)))),
      tip.label = label) %>%
    select(tip.label, cluster)
  
  #convert back to vector since this is how it is stored in your tree file list
  new_tip_labels <- c(data$tip.label)
  
  #overwrite the old labels in your tree file
  tree$tip.label <- new_tip_labels
  
  ### plot tree 
  n_groups <- length(levels(data$cluster))
  
  colors <- c("#000000",#black
              "#7F7F7FFF",#gray
              "#1F77B4FF",#blue
              "#FF7F0EFF",#orange
              "#2CA02CFF",#green
              "#D62728FF",#red
              "#9467BDFF",#purple
              "#8C564BFF",#brown
              "#E377C2FF",#pink
              "#BCBD22FF",#yellow green
              "#17BECFFF",#teal
              "#6b1414"#dark red
  )
  if(year=="17-18"){
    x_scale <- 0.012
    y_scale <- 50
  } else if(year=="19-20"){
    x_scale <- 0.0178
    y_scale <- 17
  }
  
  message(levels(data$cluster))
  
  tree.plot <- ggtree(tree) + geom_treescale(x=x_scale, y=y_scale)
  tree.plot <- tree.plot %<+% data + geom_tiplab(aes(color = cluster),size = 5) +
    geom_tippoint(aes(color=cluster), size = 2) +
    scale_color_manual(name = "Group",
                       values = c(colors[1:n_groups]), na.value = "pink") +
    theme(legend.position = c(.92,.4))


  # bs_tibble <- tibble(
  #   node=1:Nnode(tree) + Ntip(tree),
  #   bootstrap = ifelse(tree$node.label>70,tree$node.label,NA))
  # 
  # tree.plot <- tree.plot %<+% bs_tibble +
  #   geom_label(aes(label=bootstrap))
  #   #geom_text(aes(label=bootstrap), hjust=-.25, size = 3)

  print(tree.plot)

}

make_table2 <- function(data,dist){
  dist <- dist %>%
    filter(V3<3) %>%
    filter(V1 != V2)
  
  all <- data %>%
    filter(analysis_set==1 & !is.na(cluster) & !is.na(SpecimenID)) %>%
    mutate(SpecimenID=as.character(SpecimenID)) %>%
    select(SpecimenID,ID,cluster,sequenced) %>%
    rename(ID2 = ID,
           cluster2 = cluster)
  
  tmp <- data %>%
    filter(fluA48==1) %>%
    mutate(SpecimenID=as.character(SpecimenID)) %>%
    select(ID, SpecimenID, fluA72, cluster, sequenced)
  
  tmp <- left_join(tmp, dist, by = c("SpecimenID"="V1")) %>%
    mutate(dist_lte2 = ifelse(is.na(V3),0,1),
           dist_lte2 = ifelse(
             V2 %in% unique(data$SpecimenID[data$analysis_set==1]), dist_lte2, 0
           ),
           in_cluster = ifelse(is.na(cluster),0,1)) %>%
    select(ID, fluA72, cluster, sequenced, in_cluster, dist_lte2, V2, SpecimenID)
  
  # unique(
  #   tmp$V2[!is.na(tmp$V2)]
  #   ) %in%
  #   unique(
  #     data$SpecimenID[!is.na(data$SpecimenID)]
  #     )
  
  tmp <- left_join(tmp,all,by=c("V2"="SpecimenID")) %>%
    mutate(linked = ifelse(ID != ID2 & 
                             !is.na(cluster) & 
                             !is.na(cluster2) & 
                             cluster==cluster2 &
                             dist_lte2==1, 1, 0),
           sequenced = sequenced.x) %>%
    arrange(desc(linked)) %>%
    filter(!duplicated(ID)) %>%
    select(ID,cluster,fluA72,in_cluster,sequenced,dist_lte2,linked, V2, SpecimenID)
  
  cluster_seq_count <- all %>%
    group_by(cluster2) %>%
    summarise(n_seq_in_cluster = sum(sequenced, na.rm = TRUE)) %>%
    ungroup()
  
  tmp <- left_join(tmp,cluster_seq_count,by=c("cluster"="cluster2")) %>%
    mutate(
      gte2_sequenced_in_cluster = ifelse(
        !is.na(sequenced) & 
          sequenced==1 & 
          !is.na(n_seq_in_cluster) & 
          n_seq_in_cluster>1,
        1,
        0
      ),
      sequenced = ifelse(is.na(sequenced),0,sequenced)
      #in_cluster_sequenced =ifelse(!is.na(sequenced) & in_cluster==1 & sequenced==1,1,0)
    ) %>%
    select(ID,fluA72,in_cluster,sequenced,dist_lte2,
           gte2_sequenced_in_cluster,linked)
  
  
  
  t2 <- tbl_summary(
    tmp %>% select(-ID),
    digits = list(all_categorical() ~ c(0, 1)),
    by = fluA72) %>%
    add_overall()
  
  return(t2)
}
# =========== Set data path and read data ==================
filepath_1718 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 17-18/"
filepath_1920 <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/DATA/umich_data/HARVI Data 19-20/"
outpath <- "H:/RPublic/Project Mgtmt/MERC/1947_HARVI/RESULTS/Josh/TiffanyWan_ILE/"

meta_1718 <- fread(paste0(filepath_1718,"derived/analysis_1718.csv"), 
                       stringsAsFactors = FALSE) 

meta_1920 <- fread(paste0(filepath_1920,"derived/analysis_1920.csv"), 
                   stringsAsFactors = FALSE) 

tree_1718 <- read.tree(
  paste0(filepath_1718,"derived/H3N2_tree.fasta.contree")
  )

tree_1920 <- read.tree(
  paste0(filepath_1920,"derived/H1N1_tree.fasta.contree")
)

dist_1718 <- fread(
  paste0(filepath_1718,"derived/H3N2_pairwise.csv")
) 

dist_1920 <- fread(
  paste0(filepath_1920,"derived/H1N1_pairwise.csv")
) 
# ======================== Figures 3 & 4: trees ==========================
pdf(paste0(outpath,"figure3_",Sys.Date(),".pdf"), width = 15, height = 15)
plot_tree(meta_1718,tree_1718,"17-18")
dev.off()

pdf(paste0(outpath,"figure4_",Sys.Date(),".pdf"), width = 15, height = 15)
plot_tree(meta_1920,tree_1920,"19-20")
dev.off()

# =========== Table 2 ==================
t2_1718 <- make_table2(meta_1718,dist_1718)
t2_1920 <- make_table2(meta_1920,dist_1920)

t2 <- tbl_merge(list(t2_1718,t2_1920))

t2 <- as_tibble(t2)

fwrite(t2,paste0(outpath,"table2_",Sys.Date(),".csv"))

# ======================== Text ==========================
dist_1718 <- left_join(
  dist_1718,
  meta_1718 %>%
    mutate(V1 = as.character(SpecimenID),
           ID1 = ID) %>%
    select(V1,ID1),
  by = "V1"
                       )

dist_1718 <- left_join(
  dist_1718,
  meta_1718 %>%
    mutate(V2 = as.character(SpecimenID),
           ID2 = ID) %>%
    select(V2,ID2),
  by = "V2"
)

dist_1718 <- dist_1718 %>%
  filter(!is.na(ID1) & !is.na(ID2)) %>%
  arrange(ID1) %>%
  filter(!duplicated(paste0(pmax(V1, V2), pmin(V1, V2)))) %>%
  filter(V1 != V2)

dist_1718[dist_1718$ID1 %in% c(32,39,49) & dist_1718$ID2 %in% c(32,39,49),]

dist_1718[dist_1718$ID1 %in% c(130, 135, 136, 180) & dist_1718$ID2 %in% c(130, 135, 136, 180),]

dist_1718_lt3 <- dist_1718 %>%
  filter(V3<3)


dist_1920 <- left_join(
  dist_1920,
  meta_1920 %>%
    mutate(V1 = as.character(SpecimenID),
           ID1 = ID) %>%
    select(V1,ID1),
  by = "V1"
)

dist_1920 <- left_join(
  dist_1920,
  meta_1920 %>%
    mutate(V2 = as.character(SpecimenID),
           ID2 = ID) %>%
    select(V2,ID2),
  by = "V2"
)

dist_1920 <- dist_1920 %>%
  filter(!is.na(ID1) & !is.na(ID2)) %>%
  arrange(ID1) %>%
  filter(!duplicated(paste0(pmax(V1, V2), pmin(V1, V2)))) %>%
  filter(V1 != V2)

dist_1920[dist_1920$ID1 %in% c(58,59) & dist_1920$ID2 %in% c(58,59),]

dist_1920_lt3 <- dist_1920 %>%
  filter(V3<3)
