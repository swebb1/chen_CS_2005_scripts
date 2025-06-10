library(genomation)
library(tidyverse)
library(regioneR)
library(AnnotationHub)
library(patchwork)
library(GenomicRanges)

dir.create("figures")
set.seed(123)

## Differential peaks KO v WT
peaks <- readRDS("normCounts_KO_WT_peaks.rds")
peaksd <- subset(peaks, peaks$Fold > 0) |> makeGRangesFromDataFrame()
peaksdf <- subset(peaks,peaks$Fold >= 0.75) |> makeGRangesFromDataFrame()

## Rif1
p11 <- readBroadPeak("../Rif1_analysis/11_ESC_peaks.broadPeak")
p16 <- readBroadPeak("../Rif1_analysis/16_ESC_peaks.broadPeak")
rif1_peak <- c(p11,p16) |> GenomicRanges::reduce()

rif1_domain <- readBed("../Rif1_analysis/domain_intersect_homer_format.bed")

clusters <- readGeneric("../Rif1_analysis/TADs_cluster_nc.bed",chr = 1,start = 2,end = 3,keep.all.metadata = T)
names(mcols(clusters))<-"Cluster"

cnames <- c("LB1_high-RIF1_high","LB1_mid-RIF1_high","LB1_low-RIF1_high",
            "LB1_high-RIF1_low","LB1_mid-RIF1_low","LB1_low-RIF1_low")
clist <- cnames |> map(~subset(clusters,clusters$Cluster==.x))
names(clist) <- cnames

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

## All Peaks down

rlist <- list(rif1_domain,rif1_peak)
names(rlist) <- c("Rif1_Domains","Rif1_Peaks")

p_perm <- rlist |> map(~covPerm(peaksd,.x))
names(p_perm) <- c("Rif1_Domains","Rif1_Peaks")

df_p_perm <- p_perm |> imap(function(x,y){
  res <- x
  df <- data.frame(Coverage = c(res$cumulativeCoverage$observed,res$cumulativeCoverage$permuted),
                   Score = c("Observed",rep("Permuted",length(res$cumulativeCoverage$permuted))),
                   Feature = y,
                   Pval = case_when(res$cumulativeCoverage$pval<=0.05~"<0.05",.default = "ns")) 
}) |>
  bind_rows() |>
  mutate(Feature = factor(Feature,levels=names(p_perm)))

pp <- df_p_perm |> ggplot(mapping=aes(Feature,log(Coverage))) +
  geom_boxplot(data=df_p_perm |> filter(Score=="Permuted"),colour="#828A95") +
  geom_point(data=df_p_perm |> filter(Score=="Observed"),colour="#AE5B87",size=3) +
  geom_text(data=df_p_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 0.5) +
  theme_bw() +
  labs(y="Coverage overlap (log)",title = "H3K9me3 peaks down in KO v WT overlap with Rif1 - Permutation test") +
  theme(axis.text.x = element_text(angle=90))

pp2 <- df_p_perm |>
  group_by(Feature,Score) |>
  summarise(Coverage = mean(Coverage),Pval=Pval[1]) |>
  ggplot(mapping=aes(Feature,Coverage,fill=Score)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values=c("Observed"="#AE5B87","Permuted"="#828A95")) +
  geom_text(data=df_p_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 1) +
  theme_bw() +
  labs(y="Coverage overlap",title = "H3K9me3 peaks down in KO v WT overlap with Rif1 - Permutation test") +
  theme(axis.text.x = element_text(angle=90))

dfo <- rlist |> imap(function(x,y){
  data.frame(Feature = y,
             Coverage = sum(width(peaksd)),
             Overlap = cumulativeCoverage(peaksd,x)) |>
    mutate(No_Overlap = Coverage - Overlap)
}) |>
  bind_rows() |>
  select(-Coverage) |>
  mutate(Feature = factor(Feature,levels=names(p_perm))) |>
  pivot_longer(-Feature,names_to="Overlap",values_to = "Coverage") |>
  mutate(Overlap = factor(Overlap,levels=c("Overlap","No_Overlap")))

po <- dfo |> ggplot(aes(Feature,Coverage,fill=Overlap)) +
  geom_col(position = "stack") +
  theme_bw() +
  scale_fill_manual(values=c("Overlap"="#4D6D93","No_Overlap"="#828A95")) +
  labs(title="H3K9me3 peaks down in KO v WT overlap with Rif1") +
  theme(axis.text.x = element_text(angle=90))

ggsave(po/pp/pp2,filename="figures/Peaks_down_Rif1_overlaps.pdf",height=12,width=6)
write_tsv(dfo,"figures/Peaks_down_Rif1_overlaps.tsv",col_names = T)
write_tsv(df_p_perm,"figures/Peaks_down_Rif1_permutation.tsv",col_names = T)

## Filtered Peaks down 0.75

p_perm <- rlist |> map(~covPerm(peaksdf,.x))
names(p_perm) <- c("Rif1_Domains","Rif1_Peaks")

df_p_perm <- p_perm |> imap(function(x,y){
  res <- x
  df <- data.frame(Coverage = c(res$cumulativeCoverage$observed,res$cumulativeCoverage$permuted),
                   Score = c("Observed",rep("Permuted",length(res$cumulativeCoverage$permuted))),
                   Feature = y,
                   Pval = case_when(res$cumulativeCoverage$pval<=0.05~"<0.05",.default = "ns")) 
}) |>
  bind_rows() |>
  mutate(Feature = factor(Feature,levels=names(p_perm)))

pp <- df_p_perm |> ggplot(mapping=aes(Feature,log(Coverage))) +
  geom_boxplot(data=df_p_perm |> filter(Score=="Permuted"),colour="#828A95") +
  geom_point(data=df_p_perm |> filter(Score=="Observed"),colour="#AE5B87",size=3) +
  geom_text(data=df_p_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 0.5) +
  theme_bw() +
  labs(y="Coverage overlap (log)",title = "H3K9me3 peaks down in KO v WT overlap with Rif1 - Permutation test") +
  theme(axis.text.x = element_text(angle=90))

pp2 <- df_p_perm |>
  group_by(Feature,Score) |>
  summarise(Coverage = mean(Coverage),Pval=Pval[1]) |>
  ggplot(mapping=aes(Feature,Coverage,fill=Score)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values=c("Observed"="#AE5B87","Permuted"="#828A95")) +
  geom_text(data=df_p_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 1) +
  theme_bw() +
  labs(y="Coverage overlap",title = "H3K9me3 peaks down in KO v WT overlap with Rif1 - Permutation test") +
  theme(axis.text.x = element_text(angle=90))

dfo <- rlist |> imap(function(x,y){
  data.frame(Feature = y,
             Coverage = sum(width(peaksd)),
             Overlap = cumulativeCoverage(peaksd,x)) |>
    mutate(No_Overlap = Coverage - Overlap)
}) |>
  bind_rows() |>
  select(-Coverage) |>
  mutate(Feature = factor(Feature,levels=names(p_perm))) |>
  pivot_longer(-Feature,names_to="Overlap",values_to = "Coverage") |>
  mutate(Overlap = factor(Overlap,levels=c("Overlap","No_Overlap")))

po <- dfo |> ggplot(aes(Feature,Coverage,fill=Overlap)) +
  geom_col(position = "stack") +
  theme_bw() +
  scale_fill_manual(values=c("Overlap"="#4D6D93","No_Overlap"="#828A95")) +
  labs(title="H3K9me3 peaks down in KO v WT overlap with Rif1") +
  theme(axis.text.x = element_text(angle=90))

ggsave(po/pp/pp2,filename="figures/Peaks_down_Rif1_overlaps.filtered.pdf",height=12,width=6)
write_tsv(dfo,"figures/Peaks_down_Rif1_overlaps.filtered.tsv",col_names = T)
write_tsv(df_p_perm,"figures/Peaks_down_Rif1_permutation.filtered.tsv",col_names = T)


## By TAD Cluster - domains

cplist <- clist |> map(~peaksd[findOverlaps(peaksd,.x) |> queryHits() |> unique(),])

cp_perm <- cplist |> map(~covPerm(.x,rif1_domain))
names(cp_perm) <- cnames

df_cp_perm <- cp_perm |> imap(function(x,y){
  res <- x
  df <- data.frame(Coverage = c(res$cumulativeCoverage$observed,res$cumulativeCoverage$permuted),
                   Score = c("Observed",rep("Permuted",length(res$cumulativeCoverage$permuted))),
                   Feature = y,
                   Pval = case_when(res$cumulativeCoverage$pval<=0.05~"<0.05",.default = "ns")) 
}) |>
  bind_rows() |>
  mutate(Feature = factor(Feature,levels=cnames))

pp <- df_cp_perm |> ggplot(mapping=aes(Feature,log(Coverage))) +
  geom_boxplot(data=df_cp_perm |> filter(Score=="Permuted"),colour="#828A95") +
  geom_point(data=df_cp_perm |> filter(Score=="Observed"),colour="#AE5B87",size=3) +
  geom_text(data=df_cp_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 0.5) +
  theme_bw() +
  labs(y="Coverage overlap (log)",title = "H3K9me3 peaks down in KO v WT overlap with Rif1 domains - Permutation test") +
  theme(axis.text.x = element_text(angle=90))

pp2 <- df_cp_perm |>
  group_by(Feature,Score) |>
  summarise(Coverage = mean(Coverage),Pval=Pval[1]) |>
  ggplot(mapping=aes(Feature,Coverage,fill=Score)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values=c("Observed"="#AE5B87","Permuted"="#828A95")) +
  geom_text(data=df_cp_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 1) +
  theme_bw() +
  labs(y="Coverage overlap",title = "H3K9me3 peaks down in KO v WT overlap with Rif1 domains - Permutation test") +
  theme(axis.text.x = element_text(angle=90))

dfo <- cplist |> imap(function(x,y){
  data.frame(Feature = y,
             Coverage = sum(width(x)),
             Overlap = cumulativeCoverage(x,rif1_domain)) |>
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
  labs(title="H3K9me3 peaks down in KO v WT overlap with Rif1 domains") +
  theme(axis.text.x = element_text(angle=90))

ggsave(po/pp/pp2,filename="figures/Peaks_down_Rif1_domain_overlaps.pdf",height=12,width=8)
write_tsv(dfo,"figures/Peaks_down_Rif1_domain_overlaps.tsv",col_names = T)
write_tsv(df_cp_perm,"figures/Peaks_down_Rif1_domain_permutation.tsv",col_names = T)


cplist <- clist |> map(~peaksdf[findOverlaps(peaksdf,.x) |> queryHits() |> unique(),])

cp_perm <- cplist |> map(~covPerm(.x,rif1_domain))
names(cp_perm) <- cnames

df_cp_perm <- cp_perm |> imap(function(x,y){
  res <- x
  df <- data.frame(Coverage = c(res$cumulativeCoverage$observed,res$cumulativeCoverage$permuted),
                   Score = c("Observed",rep("Permuted",length(res$cumulativeCoverage$permuted))),
                   Feature = y,
                   Pval = case_when(res$cumulativeCoverage$pval<=0.05~"<0.05",.default = "ns")) 
}) |>
  bind_rows() |>
  mutate(Feature = factor(Feature,levels=cnames))

pp <- df_cp_perm |> ggplot(mapping=aes(Feature,log(Coverage))) +
  geom_boxplot(data=df_cp_perm |> filter(Score=="Permuted"),colour="#828A95") +
  geom_point(data=df_cp_perm |> filter(Score=="Observed"),colour="#AE5B87",size=3) +
  geom_text(data=df_cp_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 0.5) +
  theme_bw() +
  labs(y="Coverage overlap (log)",title = "H3K9me3 peaks down in KO v WT (<= 0.75) overlap with Rif1 domains - Permutation test") +
  theme(axis.text.x = element_text(angle=90))

pp2 <- df_cp_perm |>
  group_by(Feature,Score) |>
  summarise(Coverage = mean(Coverage),Pval=Pval[1]) |>
  ggplot(mapping=aes(Feature,Coverage,fill=Score)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values=c("Observed"="#AE5B87","Permuted"="#828A95")) +
  geom_text(data=df_cp_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 1) +
  theme_bw() +
  labs(y="Coverage overlap",title = "H3K9me3 peaks down in KO v WT (<= 0.75) overlap with Rif1 domains - Permutation test") +
  theme(axis.text.x = element_text(angle=90))

dfo <- cplist |> imap(function(x,y){
  data.frame(Feature = y,
             Coverage = sum(width(x)),
             Overlap = cumulativeCoverage(x,rif1_domain)) |>
    mutate(No_Overlap = Coverage - Overlap)
}) |>
  bind_rows() |>
  select(-Coverage) |>
  mutate(Feature = factor(Feature,levels=cnames)) |>
  pivot_longer(-Feature,names_to="Overlap",values_to = "Coverage") |>
  mutate(Overlap = factor(Overlap,levels=c("Overlap","No_Overlap")))
po <- dfo |>  ggplot(aes(Feature,Coverage,fill=Overlap)) +
  geom_col(position = "stack") +
  theme_bw() +
  scale_fill_manual(values=c("Overlap"="#4D6D93","No_Overlap"="#828A95")) +
  labs(title="H3K9me3 peaks down in KO v WT (<= 0.75) overlap with Rif1 domains") +
  theme(axis.text.x = element_text(angle=90))

ggsave(po/pp/pp2,filename="figures/Peaks_down_Rif1_domain_overlaps.filtered.pdf",height=12,width=8)
write_tsv(dfo,"figures/Peaks_down_Rif1__domain_overlaps.filtered.tsv",col_names = T)
write_tsv(df_cp_perm,"figures/Peaks_down_Rif1_domain_permutation.filtered.tsv",col_names = T)


## By TAD Cluster - peaks

cplist <- clist |> map(~peaksd[findOverlaps(peaksd,.x) |> queryHits() |> unique(),])

cp_perm <- cplist |> map(~covPerm(.x,rif1_peak))
names(cp_perm) <- cnames

df_cp_perm <- cp_perm |> imap(function(x,y){
  res <- x
  df <- data.frame(Coverage = c(res$cumulativeCoverage$observed,res$cumulativeCoverage$permuted),
                   Score = c("Observed",rep("Permuted",length(res$cumulativeCoverage$permuted))),
                   Feature = y,
                   Pval = case_when(res$cumulativeCoverage$pval<=0.05~"<0.05",.default = "ns")) 
}) |>
  bind_rows() |>
  mutate(Feature = factor(Feature,levels=cnames))

pp <- df_cp_perm |> ggplot(mapping=aes(Feature,log(Coverage))) +
  geom_boxplot(data=df_cp_perm |> filter(Score=="Permuted"),colour="#828A95") +
  geom_point(data=df_cp_perm |> filter(Score=="Observed"),colour="#AE5B87",size=3) +
  geom_text(data=df_cp_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 0.5) +
  theme_bw() +
  labs(y="Coverage overlap (log)",title = "H3K9me3 peaks down in KO v WT overlap with Rif1 peaks - Permutation test") +
  theme(axis.text.x = element_text(angle=90))

pp2 <- df_cp_perm |>
  group_by(Feature,Score) |>
  summarise(Coverage = mean(Coverage),Pval=Pval[1]) |>
  ggplot(mapping=aes(Feature,Coverage,fill=Score)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values=c("Observed"="#AE5B87","Permuted"="#828A95")) +
  geom_text(data=df_cp_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 1) +
  theme_bw() +
  labs(y="Coverage overlap",title = "H3K9me3 peaks down in KO v WT overlap with Rif1 peaks - Permutation test") +
  theme(axis.text.x = element_text(angle=90))

dfo <- cplist |> imap(function(x,y){
  data.frame(Feature = y,
             Coverage = sum(width(x)),
             Overlap = cumulativeCoverage(x,rif1_peak)) |>
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
  labs(title="H3K9me3 peaks down in KO v WT overlap with Rif1 peaks") +
  theme(axis.text.x = element_text(angle=90))

ggsave(po/pp/pp2,filename="figures/Peaks_down_Rif1_peaks_overlaps.pdf",height=12,width=8)
write_tsv(dfo,"figures/Peaks_down_Rif1_peaks_overlaps.tsv",col_names = T)
write_tsv(df_cp_perm,"figures/Peaks_down_Rif1_peaks_permutation.tsv",col_names = T)


cplist <- clist |> map(~peaksdf[findOverlaps(peaksdf,.x) |> queryHits() |> unique(),])

cp_perm <- cplist |> map(~covPerm(.x,rif1_peak))
names(cp_perm) <- cnames

df_cp_perm <- cp_perm |> imap(function(x,y){
  res <- x
  df <- data.frame(Coverage = c(res$cumulativeCoverage$observed,res$cumulativeCoverage$permuted),
                   Score = c("Observed",rep("Permuted",length(res$cumulativeCoverage$permuted))),
                   Feature = y,
                   Pval = case_when(res$cumulativeCoverage$pval<=0.05~"<0.05",.default = "ns")) 
}) |>
  bind_rows() |>
  mutate(Feature = factor(Feature,levels=cnames))

pp <- df_cp_perm |> ggplot(mapping=aes(Feature,log(Coverage))) +
  geom_boxplot(data=df_cp_perm |> filter(Score=="Permuted"),colour="#828A95") +
  geom_point(data=df_cp_perm |> filter(Score=="Observed"),colour="#AE5B87",size=3) +
  geom_text(data=df_cp_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 0.5) +
  theme_bw() +
  labs(y="Coverage overlap (log)",title = "H3K9me3 peaks down in KO v WT (<= 0.75) overlap with Rif1 peaks - Permutation test") +
  theme(axis.text.x = element_text(angle=90))

pp2 <- df_cp_perm |>
  group_by(Feature,Score) |>
  summarise(Coverage = mean(Coverage),Pval=Pval[1]) |>
  ggplot(mapping=aes(Feature,Coverage,fill=Score)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values=c("Observed"="#AE5B87","Permuted"="#828A95")) +
  geom_text(data=df_cp_perm |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 1) +
  theme_bw() +
  labs(y="Coverage overlap",title = "H3K9me3 peaks down in KO v WT (<= 0.75) overlap with Rif1 peaks - Permutation test") +
  theme(axis.text.x = element_text(angle=90))

dfo <- cplist |> imap(function(x,y){
  data.frame(Feature = y,
             Coverage = sum(width(x)),
             Overlap = cumulativeCoverage(x,rif1_peak)) |>
    mutate(No_Overlap = Coverage - Overlap)
}) |>
  bind_rows() |>
  select(-Coverage) |>
  mutate(Feature = factor(Feature,levels=cnames)) |>
  pivot_longer(-Feature,names_to="Overlap",values_to = "Coverage") |>
  mutate(Overlap = factor(Overlap,levels=c("Overlap","No_Overlap")))
po <- dfo |>  ggplot(aes(Feature,Coverage,fill=Overlap)) +
  geom_col(position = "stack") +
  theme_bw() +
  scale_fill_manual(values=c("Overlap"="#4D6D93","No_Overlap"="#828A95")) +
  labs(title="H3K9me3 peaks down in KO v WT (<= 0.75) overlap with Rif1 peaks") +
  theme(axis.text.x = element_text(angle=90))

ggsave(po/pp/pp2,filename="figures/Peaks_down_Rif1_peak_overlaps.filtered.pdf",height=12,width=8)
write_tsv(dfo,"figures/Peaks_down_Rif1_peak_overlaps.filtered.tsv",col_names = T)
write_tsv(df_p_perm,"figures/Peaks_down_Rif1_peak_permutation.filtered.tsv",col_names = T)
