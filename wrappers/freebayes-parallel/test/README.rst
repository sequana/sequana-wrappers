
Using some data from the Hm2 data set, we ran the variant calling pipeline and
fro the sambamba filtered data file, extracted the first 10000 reads

bamtools random -in data.filter.sorted.bam  -n 1000 > test.bam
bamtools sort -in test.bam -out test.sorted.bam
 

