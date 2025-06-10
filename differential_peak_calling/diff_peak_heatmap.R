library(genomation)
library(tidyverse)
library(patchwork)
library(GenomicRanges)
library(DiffBind)
library(viridis)

dir.create("figures")
set.seed(123)

peaks <- readRDS("normCounts_KO_WT_peaks.rds")
filter <- peaks$Fold>=0
counts <- read_tsv("normCounts_KO_WT_counts.tsv",col_names = T) |>
  column_to_rownames("name") |>
  as.matrix()

peaks <- peaks[filter,]
counts <- counts[filter,]

clusters <- readGeneric("../Rif1_analysis/TADs_cluster_nc.bed",chr = 1,start = 2,end = 3,keep.all.metadata = T)
names(mcols(clusters))<-"Cluster"

cnames <- c("LB1_high-RIF1_high","LB1_mid-RIF1_high","LB1_low-RIF1_high",
            "LB1_high-RIF1_low","LB1_mid-RIF1_low","LB1_low-RIF1_low")
clist <- cnames |> map(~subset(clusters,clusters$Cluster==.x))
names(clist) <- cnames

assignCluster <- function(a){
  if(intersect(a,clist$`LB1_high-RIF1_high`) |> length() >0){
    return("LB1_high-RIF1_high")
  }
  if(intersect(a,clist$`LB1_mid-RIF1_high`) |> length() >0){
    return("LB1_mid-RIF1_high")
  }
  if(intersect(a,clist$`LB1_low-RIF1_high`) |> length() >0){
    return("LB1_low-RIF1_high")
  }
  if(intersect(a,clist$`LB1_high-RIF1_low`) |> length() >0){
    return("LB1_high-RIF1_low")
  }
  if(intersect(a,clist$`LB1_mid-RIF1_low`) |> length() >0){
    return("LB1_mid-RIF1_low")
  }
  if(intersect(a,clist$`LB1_low-RIF1_low`) |> length() >0){
    return("LB1_low-RIF1_low")
  }
  else{
    return("NA")
  }
}

assclu <- vector()
for(i in 1:length(peaks)){
  clust <- assignCluster(peaks[i])
  assclu <- append(assclu,clust)
}
pc <- data.frame(name = names(peaks),cluster=factor(assclu,levels=c(cnames,"NA")))       

## HMap

cluster_within_group <- function(mat, group_labels, scaled=T) {
  # Scale rows (z-score): subtract mean, divide by SD across columns
  if(scaled == T){
    mat_scaled <- t(scale(t(mat)))
  }
  else{
    mat_scaled <- mat
  }
  # Split and cluster within group
  groups <- levels(group_labels)
  clustered_rows <- c()
  
  for (grp in groups) {
    rows_in_group <- which(group_labels == grp)
    submat <- mat_scaled[rows_in_group, , drop=FALSE]
    
    # Hierarchical clustering
    hc <- hclust(dist(submat))
    clustered_rows <- c(clustered_rows, rows_in_group[hc$order])
  }
  
  mat_scaled[clustered_rows, , drop=FALSE]
}

annotation_col <- data.frame(
  Condition = colnames(counts) |> str_split_i("_",1)
)
rownames(annotation_col) <- colnames(counts)

matrix <- cluster_within_group(counts,pc$cluster,scaled=F)
rm = pc |> column_to_rownames("name")

annotation_row <- data.frame(
  Cluster = rm[rownames(matrix),"cluster"]
)
rownames(annotation_row) <- rownames(matrix)

cols = c("WT"="#212529","KO"="#6c757d","DPP1"="#AE5B87","HEMI"="#65AFFF")
ann_colours <- list(
  Condition = cols,
  Cluster = c("LB1_high-RIF1_high"="#225981","LB1_mid-RIF1_high"="#3E92CC","LB1_low-RIF1_high"="#8EBFE1",
              "LB1_high-RIF1_low"="#891A36","LB1_mid-RIF1_low"="#D8315B","LB1_low-RIF1_low"="#E88795","NA"="grey20")
)


hm <- pheatmap(matrix,
         cluster_cols = F,
         cluster_rows = F,
         show_rownames = F,
         color = viridis(100),
         scale = "row",
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = ann_colours,border_color = "black",main = "H3K9me3 Differential peaks KO v WT")

pdf("figures/Differential_Peaks_KO_WT_heatmap.down.pdf",height=8,width=8)
hm
dev.off()
