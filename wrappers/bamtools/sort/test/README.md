We took 50 reads and map them onto a reference (measles) to generate a small bam
file:

    bwa index genome.fasta
    bwa mem genome.fasta a_R1_.fastq  a_R2_.fastq   > test.sam
    samtools view -Sbh test.sam -o test.bam

