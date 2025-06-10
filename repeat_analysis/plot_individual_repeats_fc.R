library(genomation)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(ggrepel)
library(ggpubr)  

## Load repeat annotations
reps <- readGeneric("/home/genomes/mouse/mm10/annotation/UCSC_mm10_repeat_masker.tsv",chr = 1,start = 2,end = 3,strand=4,keep.all.metadata = T,header = T)
names(mcols(reps))[1] <- "name"

## Selected repeat families
rep_families = c("ERVK", "ERV1", "ERVL","L1", "Alu","B2","B4","Satellite")

repList <- rep_families |> map(function(n){reps[grep(n,reps$repFamily),]}) |> set_names(rep_families)

##Split by list
repf <- subset(reps,repFamily %in% rep_families) |>
  mcols() |> as_tibble() |> dplyr::select(Repeat=name,Family=repFamily) |> unique()

repfCount = repf |> group_by(Family) |> summarise(count=n()) |> pull(count)
names(repfCount) = repf |> group_by(Family) |> summarise(count=n()) |> pull(Family)

t <- readRDS("FC_RM_repeat_counts_normalised.rds")

comparisons <- list(c("KO", "WT"), c("Hemi", "DPP1"))
dir.create("individual_repeats")

cols = c("WT"="#212529","KO"="#6c757d","DPP1"="#274060","Hemi"="#65AFFF","DMSO"="#e36414","IAA"="#fb8b24")

t$Geneid |> unique() |> map(function(x){
  frep <- t |> filter(Geneid == x)  
  if(nrow(frep)>0 & sum(frep$K9n)>0){
    p <- frep |> ggplot(aes(Group,log2n,colour=Group,label=Sample)) + 
      geom_boxplot(outlier.shape = NA) +
      geom_point(alpha=0.8) +
      geom_label_repel() +
      theme_bw() +
      scale_colour_manual(values = cols) +
      theme(axis.text.x = element_text(angle=90,hjust=1)) +
      labs(x="Sample",y="H3K9 log2(IP/Input)") +
      guides(colour="none") +
      ggtitle(paste("Repeat Masker -",x,repf |> filter(Repeat==x) |> pull(Family))) +
      stat_compare_means(comparisons = comparisons, method = "t.test",vjust = 2) #,label = "p.signif")
    
    ggsave(plot = p,filename = paste0("individual_repeats/",x,".log2.pdf"),height = 6,width=6)
  }
})

t$Geneid |> unique() |> map(function(x){
  frep <- t |> filter(Geneid == x) |>
    select(Geneid,K9n,Inputn,Sample,Group) |>
    pivot_longer(-c(Geneid,Sample,Group), names_to = "Type", values_to = "Count_RPM")
  if(nrow(frep)>0 & sum(frep$Count_RPM)>0){
    p <- frep |> ggplot(aes(Group,Count_RPM,colour=Group,label=Sample)) + 
      geom_boxplot(outlier.shape = NA) +
      geom_point(alpha=0.8) +
      geom_label_repel() +
      theme_bw() +
      scale_colour_manual(values = cols) +
      theme(axis.text.x = element_text(angle=90,hjust=1)) +
      labs(x="Sample",y="Normalised Read Counts (RPM)") +
      ggtitle(paste("Repeat Masker -",x,repf |> filter(Repeat==x) |> pull(Family))) +
      guides(colour="none") +
      facet_wrap(~Type,ncol = 1,scales = "free") +
      stat_compare_means(comparisons = comparisons, method = "t.test",vjust = 2) #,label = "p.signif")
    
    ggsave(plot = p,filename = paste0("individual_repeats/",x,".pdf"),height = 8,width=6 )
  }
})
