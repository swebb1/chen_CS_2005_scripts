library(readr)
library(DiffBind)

samples <- read_tsv("samples.tsv",col_names = T)

dba <- dba(sampleSheet=samples)

plot(dba)

dba <- dba.count(dba,minOverlap = 2,summits = TRUE)
saveRDS(dba,"dba.count.rds")

plot(dba)

dba <- dba.normalize(dba)
saveRDS(dba,"dba.normalize.rds")

png("Normalised_peak_counts_heatmap.png")
plot(dba)
dev.off()

dba <- dba.contrast(dba,minMembers = 2,reorderMeta=list(Condition="WT"))
saveRDS(dba,"dba.contrast.rds")

dba <- dba.analyze(dba)
saveRDS(dba,"dba.analyze.rds")

dba.show(dba, bContrasts=TRUE)
plot(dba, contrast=1)

labs<-c("WT_B","WT_E","WT_F","KO_18","KO_20","KO_24","KO_2","HEMI_G1","HEMI_H6","HEMI_1","DPP1_G11","DPP1_1","DPP1_H2")
contrasts = dba.show(dba, bContrasts = T) |> unite(Contrast,Group,Group2,sep=".") |> pull(Contrast)
db_peaks = dba.show(dba, bContrasts = T) |> pull(DB.DESeq2)

rtracklayer::export.bed(dba$peaks |> as.data.frame() |> GenomicRanges::makeGRangesFromDataFrame(),"peaks.consensus.bed")

## Export normalised counts for WT v KO
diffPeaks <- dba.report(dba, th = 0.05,contrast = 1)
indices <- names(diffPeaks)

library(DESeq2)
# Extract DESeq2 object from the first contrast (for all samples)
deseq <- dba$DESeq2$DEdata
# Get normalized counts (DESeq2-normalized, not log-transformed)
norm_counts <- counts(deseq, normalized=TRUE)
norm_counts_KO_WT <- norm_counts[indices,]
write_tsv(norm_counts_KO_WT |> as.data.frame() |> rownames_to_column("name"),"normCounts_KO_WT_counts.tsv",col_names = T)

saveRDS(diffPeaks,"normCounts_KO_WT_peaks.rds")
rtracklayer::export.bed(diffPeaks,"normCounts_KO_WT_peaks.bed")


## Heatmaps
dir.create("diff_hm")

1:5 |> map(function(x){
  
  file = paste0("diff_hm/DB_peaks.",contrasts[x],".",db_peaks[x])
  peaks.DB <- dba.report(dba,contrast=x)
  write_csv(peaks.DB |> as.data.frame(),paste0(file,".tsv"),col_names = T)
  
  rtracklayer::export.bed(peaks.DB,paste0(file,".bed"))
  
})

########

plot(dba, contrast=1)
dba.plotHeatmap(dba)

hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

## Read cores
dba.plotHeatmap(dba, contrast=1, correlations=FALSE,scale="row", colScheme = hmap)

1:5 |> map(function(x){
  profiles <- dba.plotProfile(dba,samples = 1:13,sites = x,merge = F,labels=labs)
  
  pdf(paste0("diff_hm/DB_heatmap.",contrasts[x],".",db_peaks[x],".pdf"))
  dba.plotProfile(profiles)
  dev.off()
  
})

