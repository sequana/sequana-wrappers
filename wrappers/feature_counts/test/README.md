Here, we use a 50-reads paired data mapped onto a measles.fa genome.
We built a dummy GFF file to search for reads in the region 7000-7400 expecting
3 reads

    bwa index measles.fa
    bwa mem measles.fa a_R1_.fastq  a_R2_.fastq > test.bam
    samtools view -hb test.sam  > test.bam
    bamtools sort -in test.bam -out test.sorted.bam
    bamtools index -in test.sorted.bam
    featureCounts -T 1 -a test.gff  -o test.csv test.sorted.bam -s 0 -t gene -g ID

