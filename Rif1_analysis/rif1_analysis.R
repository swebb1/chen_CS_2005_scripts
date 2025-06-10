library(genomation)
library(tidyverse)
library(regioneR)
library(AnnotationHub)
library(patchwork)
library(GenomicRanges)


set.seed(123)
dir.create("figures")

p11 <- readBroadPeak("11_ESC_peaks.broadPeak")
p16 <- readBroadPeak("16_ESC_peaks.broadPeak")
rif1_peak <- c(p11,p16) |> GenomicRanges::reduce()

rif1_domain <- readBed("domain_intersect_homer_format.bed")

clusters <- readGeneric("TADs_cluster_nc.bed",chr = 1,start = 2,end = 3,keep.all.metadata = T)
names(mcols(clusters))<-"Cluster"

ERVs <- readRDS("../Repeat_analysis/ERV.rds")
ERVK <- subset(ERVs,ERVs$repFamily =="ERVK")
ERV1 <- subset(ERVs,ERVs$repFamily =="ERV1")
ERVL <- subset(ERVs,ERVs$repFamily =="ERVL")

enr <- read_tsv("../Repeat_analysis/Enriched_ERVs.tsv",col_names = T)

ERV_enr <- subset(ERVs,ERVs$repName %in% (enr$Geneid |> unique()))
ERVK_enr <- subset(ERVK,ERVK$repName %in% (enr$Geneid |> unique()))
ERV1_enr <- subset(ERV1,ERV1$repName %in% (enr$Geneid |> unique()))
ERVL_enr <- subset(ERVL,ERVL$repName %in% (enr$Geneid |> unique()))


## Testing overlaps


cumulativeCoverage <- function(A, B, ...) {
  
  # Find overlaps
  hits <- findOverlaps(A, B)
  
  # Get overlapping regions
  overlaps <- pintersect(A[queryHits(hits)], B[subjectHits(hits)])
  
  # Reduce to merge overlapping ranges and avoid double-counting
  reduced_overlaps <- reduce(overlaps)
  
  # Return total base pairs covered
  sum(width(reduced_overlaps))
}

## Get masked regions
ah <- AnnotationHub()
mask <- ah[["AH107322"]]

covPerm <- function(A,B,n=100){
  permTestRes <- permTest(A = A,
                          B = B,
                          ntimes = n,
                          genome = "mm10",
                          randomize.function=randomizeRegions,
                          randomize.params = list(mask = mask),
                          evaluate.function=cumulativeCoverage,
                          alternative = "greater")
}

## Rif1 Domains

erv_perm <- list(ERVK_enr,ERVL_enr,ERV1_enr) |> map(~covPerm(.x,rif1_domain,n = 100))
erv_names <- c("ERVK","ERVL","ERV1")
names(erv_perm) <- erv_names

df_erv_perm <- erv_perm |> imap(function(x,y){
    res <- x
    df <- data.frame(Coverage = c(res$cumulativeCoverage$observed,res$cumulativeCoverage$permuted),
           Score = c("Observed",rep("Permuted",length(res$cumulativeCoverage$permuted))),
           Feature = y,
           Pval = case_when(res$cumulativeCoverage$pval<=0.05~"<0.05",.default = "ns")) 
}) |>
  bind_rows()

pp <- df_erv_perm |> ggplot(mapping=aes(Feature,log(Coverage))) +
  geom_boxplot(data=df_erv_perm |> filter(Score=="Permuted"),colour="#828A95") +
  geom_point(data=df_erv_perm |> filter(Score=="Observed"),colour="#AE5B87",size=3) +
  geom_text(data=df_erv_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 0.5) +
  theme_bw() +
  labs(y="Coverage overlap (log)",title = "Rif1 Domain overlap - Permutation test")

pp2 <- df_erv_perm |>
  group_by(Feature,Score) |>
  summarise(Coverage = mean(Coverage),Pval=Pval[1]) |>
  ggplot(mapping=aes(Feature,Coverage,fill=Score)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values=c("Observed"="#AE5B87","Permuted"="#828A95")) +
  geom_text(data=df_erv_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 1) +
  theme_bw() +
  labs(y="Coverage overlap",title = "Rif1 Domain overlap - Permutation test")

dfo <- list("ERV1"=ERV1,"ERVK"=ERVK,"ERVL"=ERVL) |> imap(function(x,y){
  data.frame(Feature = y,
             Coverage = sum(width(x)),
             Overlap = cumulativeCoverage(x,rif1_domain)) |>
             mutate(No_Overlap = Coverage - Overlap)
}) |>
  bind_rows() |>
  select(-Coverage) |>
  pivot_longer(-Feature,names_to="Overlap",values_to = "Coverage") |>
  mutate(Overlap = factor(Overlap,levels=c("Overlap","No_Overlap")))
po <- dfo |>  ggplot(aes(Feature,Coverage,fill=Overlap)) +
  geom_col(position = "stack") +
  theme_bw() +
  scale_fill_manual(values=c("Overlap"="#4D6D93","No_Overlap"="#828A95")) +
  labs(title="Rif1 Domains")
  
ggsave(po/pp/pp2,filename="figures/Rif1_domain_ERV_overlaps.pdf",height=12,width=6)
write_tsv(dfo,"figures/Rif1_domain_ERV_overlaps.tsv",col_names = T)
write_tsv(df_erv_perm,"figures/Rif1_domain_ERV_permutation.tsv",col_names = T)

## Rif1 Peaks

erv_perm <- list(ERVK_enr,ERVL_enr,ERV1_enr) |> map(~covPerm(.x,rif1_peak),n=100)
erv_names <- c("ERVK","ERVL","ERV1")
names(erv_perm) <- erv_names

df_erv_perm <- erv_perm |> imap(function(x,y){
  res <- x
  df <- data.frame(Coverage = c(res$cumulativeCoverage$observed,res$cumulativeCoverage$permuted),
                   Score = c("Observed",rep("Permuted",length(res$cumulativeCoverage$permuted))),
                   Feature = y,
                   Pval = case_when(res$cumulativeCoverage$pval<=0.05~"<0.05",.default = "ns"))
}) |>
  bind_rows()


pp <- df_erv_perm |> ggplot(mapping=aes(Feature,log(Coverage))) +
  geom_boxplot(data=df_erv_perm |> filter(Score=="Permuted"),colour="#828A95") +
  geom_point(data=df_erv_perm |> filter(Score=="Observed"),colour="#AE5B87",size=3) +
  geom_text(data=df_erv_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 0.5) +
  theme_bw() +
  labs(y="Coverage overlap (log)",title = "Rif1 Peak overlap - Permutation test")

pp2 <- df_erv_perm |>
  group_by(Feature,Score) |>
  summarise(Coverage = mean(Coverage),Pval=Pval[1]) |>
  ggplot(mapping=aes(Feature,Coverage,fill=Score)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values=c("Observed"="#AE5B87","Permuted"="#828A95")) +
  geom_text(data=df_erv_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 1) +
  theme_bw() +
  labs(y="Coverage overlap",title = "Rif1 Peak overlap - Permutation test")

dfo <- list("ERV1"=ERV1,"ERVK"=ERVK,"ERVL"=ERVL) |> imap(function(x,y){
  data.frame(Feature = y,
             Coverage = sum(width(rif1_peak)),
             Overlap = cumulativeCoverage(rif1_peak,x)) |>
    mutate(No_Overlap = Coverage - Overlap)
}) |>
  bind_rows() |>
  select(-Coverage) |>
  pivot_longer(-Feature,names_to="Overlap",values_to = "Coverage") |>
  mutate(Overlap = factor(Overlap,levels=c("Overlap","No_Overlap")))
po <- dfo |>  ggplot(aes(Feature,Coverage,fill=Overlap)) +
  geom_col(position = "stack") +
  theme_bw() +
  scale_fill_manual(values=c("Overlap"="#4D6D93","No_Overlap"="#828A95")) +
  labs(title="Rif1 Peaks")

ggsave(po/pp/pp2,filename="figures/Rif1_peak_ERV_overlaps.pdf",height=8,width=6)
write_tsv(dfo,"figures/Rif1_peak_ERV_overlaps.tsv",col_names = T)
write_tsv(df_erv_perm,"figures/Rif1_peak_ERV_permutation.tsv",col_names = T)


### ERVK at TAD clusters

cnames <- c("LB1_high-RIF1_high","LB1_mid-RIF1_high","LB1_low-RIF1_high",
            "LB1_high-RIF1_low","LB1_mid-RIF1_low","LB1_low-RIF1_low")
clist <- cnames |> map(~subset(clusters,clusters$Cluster==.x))
names(clist) <- cnames

clus_perm <- clist |> map(~covPerm(ERVK_enr,.x,n=100))
names(clus_perm) <- cnames

df_clus_perm <- clus_perm |> imap(function(x,y){
  res <- x
  df <- data.frame(Coverage = c(res$cumulativeCoverage$observed,res$cumulativeCoverage$permuted),
                   Score = c("Observed",rep("Permuted",length(res$cumulativeCoverage$permuted))),
                   Feature = y,
                   Pval = case_when(res$cumulativeCoverage$pval<=0.05~"<0.05",.default = "ns")) 
}) |>
  bind_rows() |>
  mutate(Feature = factor(Feature,levels=cnames))
  

pp <- df_clus_perm |> ggplot(mapping=aes(Feature,log(Coverage))) +
  geom_boxplot(data=df_clus_perm |> filter(Score=="Permuted"),colour="#828A95") +
  geom_point(data=df_clus_perm |> filter(Score=="Observed"),colour="#AE5B87",size=3) +
  geom_text(data=df_clus_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 0.5) +
  theme_bw() +
  labs(y="Coverage overlap (log)",title = "ERVK overlap - Permutation test") +
  theme(axis.text.x = element_text(angle=90))

pp2 <- df_clus_perm |>
  group_by(Feature,Score) |>
  summarise(Coverage = mean(Coverage),Pval=Pval[1]) |>
  ggplot(mapping=aes(Feature,Coverage,fill=Score)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values=c("Observed"="#AE5B87","Permuted"="#828A95")) +
  geom_text(data=df_clus_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 1) +
  theme_bw() +
  labs(y="Coverage overlap",title = "ERVK overlap - Permutation test") +
  theme(axis.text.x = element_text(angle=90))

dfo <- clist |> imap(function(x,y){
  data.frame(Feature = y,
             Coverage = sum(width(x)),
             Overlap = cumulativeCoverage(ERVK_enr,x)) |>
    mutate(No_Overlap = Coverage - Overlap)
}) |>
  bind_rows() |>
  select(-Coverage) |>
  mutate(Feature = factor(Feature,levels=cnames)) |>
  pivot_longer(-Feature,names_to="Overlap",values_to = "Coverage") |>
  mutate(Overlap = factor(Overlap,levels=c("Overlap","No_Overlap")))
po <- dfo |> ggplot(aes(Feature,Coverage,fill=Overlap)) +
  geom_col(position = "stack") +
  theme_bw() +
  scale_fill_manual(values=c("Overlap"="#4D6D93","No_Overlap"="#828A95")) +
  labs(title="ERVK overlap") +
  theme(axis.text.x = element_text(angle=90))

ggsave(po/pp/pp2,filename="figures/ERVK_TAD_cluster_overlaps.pdf",height=12,width=6)
write_tsv(dfo,"figures/ERV_TAD_cluster_overlaps.tsv",col_names = T)
write_tsv(df_clus_perm,"figures/ERV_TAD_cluster_permutation.tsv",col_names = T)

