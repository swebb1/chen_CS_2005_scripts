
## Convert bigWig to ucsc compatible bigWig by adding chr before each contig
#bash bw_ucsc.sh <bw_file> <chr_map> <genome len file ucsc>

bigWigToWig $1.bw $1.wig
perl bin/scripts/chrMap.pl --in $1.wig --chr $2 --remove > $1.ucsc.wig
wigToBigWig $1.ucsc.wig $3 $1.ucsc.bw
rm $1.wig $1.ucsc.wig
