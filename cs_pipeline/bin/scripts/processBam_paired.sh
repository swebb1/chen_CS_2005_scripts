
##Convert to bam

samtools sort -@ 5 bwa/$1/$1.sam -o bwa/$1/$1.bam -T bwa/$1/$1 -O BAM
samtools view -F 2304 -bh bwa/$1/$1.bam -o bwa/$1/$1.total.bam
samtools index bwa/$1/$1.total.bam
samtools flagstat bwa/$1/$1.total.bam > bwa/$1/$1.total.flagstat
fastqc bwa/$1/$1.total.bam

##Get primary alignments - properly paired
samtools view -bh -F 260 -f 3 bwa/$1/$1.total.bam -o bwa/$1/$1.primary.bam
samtools index bwa/$1/$1.primary.bam
samtools flagstat bwa/$1/$1.primary.bam > bwa/$1/$1.primary.flagstat
fastqc bwa/$1/$1.primary.bam

##filter for uniquely mapped reads (mapping quality > 20)
samtools view -bh -q 20 bwa/$1/$1.primary.bam -o bwa/$1/$1.unique.bam
samtools index bwa/$1/$1.unique.bam
samtools flagstat bwa/$1/$1.unique.bam > bwa/$1/$1.unique.flagstat
fastqc bwa/$1/$1.unique.bam
