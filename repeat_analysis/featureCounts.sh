

featureCounts -p -M --fraction -T 10 -a rm_mm10.gtf -o rm.fc.tsv ../../ChIP-seq_dec24/bwa/*/*.primary_NR_bl.bam
featureCounts -p -M --fraction -T 10 -a rm_mm10.gtf -o rm.fc.mar.tsv ../../ChIP-seq_mar25/bwa/*/*.primary_NR_bl.bam


