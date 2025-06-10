library(genomation)
library(tidyverse)
library(GenomicRanges)
library(pheatmap)
library(patchwork)

dir.create("figures")

rep_families = c("ERV1","ERVK","ERVL")
ERVs <- readRDS("ERV.rds")

snlevels <- c("WT_B","WT_E","WT_F","KO_KO_2","KO_KO18","KO_KO20","KO_KO24",
             "Hemi_HEMI_1","Hemi_G1","Hemi_H6","DPP1_DPP1_1","DPP1_G11","DPP1_H2")
slevels <- c("B","E","F","KO_2","KO18","KO20","KO24",
              "HEMI_1","G1","H6","DPP1_1","G11","H2")

t <- readRDS("FC_RM_repeat_counts_normalised.rds")
terv <- t |> filter(repFamily %in% rep_families) |>
  mutate(SampleName = factor(SampleName,levels = snlevels),
         Sample = factor(Sample,levels = slevels)) |>
  arrange(SampleName)

filter_t <- function(t,count=10,enrichment=0){
  cf <- t |> filter(Group=="WT") |> group_by(Geneid) |> summarise(K9n = mean(K9n)) |>
    filter(K9n > count) |> pull(Geneid)
  enr <- t |> filter(Group=="WT") |> group_by(Geneid) |> summarise(K9n = mean(K9n),Inputn=mean(Inputn)) |>
    filter(log2(K9n/Inputn) > enrichment) |> pull(Geneid)
  filtered = c(cf[cf %in% enr])
  t |> filter(Geneid %in% filtered)
}

terv |> filter_t(count = 10,enrichment = 0.1) |> write_tsv("Enriched_ERVs.tsv",col_names = T)

cols = c("WT"="#212529","KO"="#6c757d","DPP1"="#AE5B87","Hemi"="#65AFFF")

ann_colours <- list(
  Condition = cols,
  Family = c("ERV1"="#ef476f","ERVK"="#06d6a0","ERVL"="#073b4c")
)

cluster_within_group <- function(mat, group_labels, scaled=T) {
  # Scale rows (z-score): subtract mean, divide by SD across columns
  if(scaled == T){
    mat_scaled <- t(scale(t(mat)))
  }
  else{
    mat_scaled <- mat
  }
  # Split and cluster within group
  groups <- unique(group_labels)
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

hmplot <- function(t,names="SampleName",rows=FALSE,title,scaled=T){
  mdf <- t |>
    pivot_wider(id_cols = c("Geneid","repFamily"),names_from = all_of(names),values_from = "log2n") |>
    arrange(repFamily,Geneid)
  
  matrix <- mdf |> select(-repFamily) |> 
    column_to_rownames("Geneid") |>
    as.matrix()
  
  repMap <- mdf |> select(Geneid,repFamily)
  
  annotation_col <- data.frame(
    Condition = colnames(matrix) |> str_split_i("_",1)
  )
  rownames(annotation_col) <- colnames(matrix)
  
  matrix2 <- cluster_within_group(matrix,repMap$repFamily,scaled=scaled)
  rm = repMap |> column_to_rownames("Geneid")
  
  annotation_row <- data.frame(
    Family = rm[rownames(matrix2),"repFamily"]
  )
  rownames(annotation_row) <- rownames(matrix2)
  
  pheatmap(
    matrix2,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    #scale = "row",         # Normalize across repeat types (optional)
    #color = colorRampPalette(c("white", "lightblue", "navy"))(100),
    main = paste("Repeat Enrichment",title),
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = ann_colours,
    show_rownames = F
  )  
}


pdf("figures/ERV_heatmap.pdf",height=6,width=8)
hmplot(terv |> filter_t(count = 10,enrichment = -100),title="ERV Repeats - Filtered by minimum count",scaled=T)
dev.off()
pdf("figures/ERV_heatmap.enriched.pdf",height=6,width=8)
hmplot(terv |> filter_t(count = 10,enrichment = 0.1),title="ERV Repeats - Filtered by minimum count and WT log2(IP/Input) > 0.1")
dev.off()

## Heatmap grouped

t_g <- terv |> filter_t(count = 10,enrichment = -100) |>
  group_by(Group,Geneid,repFamily) |>
  summarise(K9n = mean(K9n),Inputn=mean(Inputn)) |>
  mutate(log2n = log2(K9n/Inputn))

t_gf <- terv |> filter_t(count = 10,enrichment = 0.1) |>
  group_by(Group,Geneid,repFamily) |>
  summarise(K9n = mean(K9n),Inputn=mean(Inputn)) |>
  mutate(log2n = log2(K9n/Inputn))

pdf("figures/ERV_heatmap.grouped.pdf",height=6,width=8)
hmplot(t_g,names="Group",title="ERV Repeats - Filtered by minimum count")
dev.off()
pdf("figures/ERV_heatmap.enriched.grouped.pdf",height=6,width=8)
hmplot(t_gf,names="Group",title="ERV Repeats - Filtered by minimum count and WT log2(IP/Input) > 0.1")
dev.off()


## Box plots - Filtered by count only

pa <- terv |> filter_t(count = 10,enrichment = -100) |>
  ggplot(aes(Sample,log2n,colour=Group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2,size=0.04,alpha=0.3) +
  facet_wrap(~repFamily,scales = "free") +
  theme_bw() +
  labs(y="log2(IP/Input)",title = "ERV Repeats - Filtered by minimum count") +
  scale_color_manual(values=cols) +
  theme(axis.text.x = element_text(angle=90))

ggsave(pa, filename="figures/ERV_boxplots.individual.pdf",height=5,width=8)

pb <- terv |> filter_t(count = 10,enrichment = 0.1) |>
  ggplot(aes(Sample,log2n,colour=Group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2,size=0.04,alpha=0.3) +
  facet_wrap(~repFamily,scales = "free") +
  theme_bw() +
  labs(y="log2(IP/Input)",title = "ERV Repeats - Filtered by minimum count and WT log2(IP/Input) > 0.1") +
  scale_color_manual(values=cols) +
  theme(axis.text.x = element_text(angle=90))

ggsave(pb, filename="figures/ERV_boxplots.individual.enriched.pdf",height=5,width=8)

t_g <- terv |> filter_t(count = 10,enrichment = -100) |>
    group_by(Group,Geneid,repFamily) |>
    summarise(K9n = mean(K9n),Inputn=mean(Inputn)) |>
    mutate(log2n = log2(K9n/Inputn))

p1 <- ggplot(t_g,aes(Group,log2n,colour=Group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2,size=0.1) +
  facet_wrap(~repFamily,scales = "free") +
  theme_bw() +
  labs(y="log2(IP/Input)",title = "ERV Repeats - Filtered by minimum count") +
  scale_color_manual(values=cols)
  
trel <- t_g |> pivot_wider(id_cols = c(Geneid,repFamily),names_from = "Group",values_from = "log2n")
  
trell <- trel |> mutate(`KO-WT` = KO - WT, `DPP1-Hemi` = DPP1 - Hemi) |>
    select(-c(WT,Hemi,KO,DPP1)) |>
    pivot_longer(-c(Geneid,repFamily),names_to = "Sample",values_to = "Score") |>
    mutate(Sample = factor(Sample,levels = c("KO-WT","DPP1-Hemi")))
  
p2 <- ggplot(trell,aes(Sample,Score,colour=Sample)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2,size=0.1) +
    facet_wrap(~repFamily) +
    theme_bw() +
    labs(y="log2(IP/Input) - log2(IP/Input)",title = "") +
    scale_color_manual(values=c("#6c757d","#AE5B87"))+
    coord_cartesian(ylim = c(-1,1))

ggsave(p1/p2,filename="figures/ERV_boxplots.pdf",height=8,width=8)

## Box plots - Filtered by count and enrichment

t_g <- terv |> filter_t(count = 10,enrichment = 0.1) |>
  group_by(Group,Geneid,repFamily) |>
  summarise(K9n = mean(K9n),Inputn=mean(Inputn)) |>
  mutate(log2n = log2(K9n/Inputn))

p1 <- ggplot(t_g,aes(Group,log2n,colour=Group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2,size=0.1) +
  facet_wrap(~repFamily,scales = "free") +
  theme_bw() +
  labs(y="log2(IP/Input)",title = "ERV Repeats - Filtered by minimum count and WT log2(IP/Input) > 0.1") +
  scale_color_manual(values=cols)

trel <- t_g |> pivot_wider(id_cols = c(Geneid,repFamily),names_from = "Group",values_from = "log2n")

trell <- trel |> mutate(`KO-WT` = KO - WT, `DPP1-Hemi` = DPP1 - Hemi) |>
  select(-c(WT,Hemi,KO,DPP1)) |>
  pivot_longer(-c(Geneid,repFamily),names_to = "Sample",values_to = "Score") |>
  mutate(Sample = factor(Sample,levels = c("KO-WT","DPP1-Hemi")))

p2 <- ggplot(trell,aes(Sample,Score,colour=Sample)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2,size=0.1) +
  facet_wrap(~repFamily) +
  theme_bw() +
  labs(y="log2(IP/Input) - log2(IP/Input)",title = "") +
  scale_color_manual(values=c("#6c757d","#AE5B87"))+
  coord_cartesian(ylim = c(-1,1))

ggsave(p1/p2,filename="figures/ERV_boxplots.enriched.pdf",height=8,width=8)

