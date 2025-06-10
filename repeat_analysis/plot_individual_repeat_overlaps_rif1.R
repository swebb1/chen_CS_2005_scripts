 library(genomation)
library(tidyverse)
library(regioneR)
library(AnnotationHub)
library(patchwork)
library(GenomicRanges)

set.seed(123)
dir.create("individual_overlaps")

p11 <- readBroadPeak("../Rif1_analysis/11_ESC_peaks.broadPeak")
p16 <- readBroadPeak("../Rif1_analysis/16_ESC_peaks.broadPeak")
rif1_peak <- c(p11,p16) |> GenomicRanges::reduce()

ERVs <- readRDS("../Repeat_analysis/ERV.rds")

enr <- read_tsv("../Repeat_analysis/Enriched_ERVs.tsv",col_names = T)

ERV_enr <- subset(ERVs,ERVs$repName %in% (enr$Geneid |> unique()))

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


ERV_enr$repName |> unique() |> map(function(x){
  
  frep <- subset(ERV_enr,ERV_enr$repName==x)  
    
  res <- covPerm(frep,rif1_peak,n=100)
    
  df <- data.frame(Coverage = c(res$cumulativeCoverage$observed,res$cumulativeCoverage$permuted),
                   Score = c("Observed",rep("Permuted",length(res$cumulativeCoverage$permuted))),
                   Feature = x,
                   Pval = case_when(res$cumulativeCoverage$pval<=0.05~"<0.05",.default = "ns"))
    
  pp <- df |>
      ggplot(mapping=aes(Feature,Coverage,fill=Score)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values=c("Observed"="#AE5B87","Permuted"="#828A95")) +
      geom_text(data=df |> filter(Score=="Observed"),mapping=aes(label=Pval),nudge_y = 1) +
      theme_classic() +
      labs(y="Coverage overlap",title = "Rif1 Peak overlap - Permutation test")
    
  
  dfo <-data.frame(Feature = x,
                 Coverage = sum(width(frep)),
                 Overlap = cumulativeCoverage(frep,rif1_peak)) |>
        mutate(No_Overlap = Coverage - Overlap) |>
        select(-Coverage) |>
      pivot_longer(-Feature,names_to="Overlap",values_to = "Coverage") |>
      mutate(Overlap = factor(Overlap,levels=c("Overlap","No_Overlap")))
  
  po <- dfo |>  ggplot(aes(Feature,Coverage,fill=Overlap)) +
      geom_col(position = "stack") +
      theme_bw() +
      scale_fill_manual(values=c("Overlap"="#4D6D93","No_Overlap"="#828A95")) +
      labs(title="Rif1 peaks")
    
    ggsave(po/pp,filename=paste0("individual_overlaps/",x,"_Rif1_peak_overlap.pdf"),height=6,width=6)
})

