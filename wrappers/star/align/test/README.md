Genome is meales K01711. Fasta and GFF3 were downloaded from NCBI.
GTF was built from gffread


Here is an example with legacy code (2 passes and genome generate):


STAR --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles genome.fasta  --genomeSAindexNbases 5

STAR --genomeDir star_index/ --readFilesIn test.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix star_mapping/_init_

STAR --runMode genomeGenerate --genomeDir star_mapping --genomeFastaFiles genome.fasta  --genomeSAindexNbases 5 --sjdbFileChrStartEnd star_mapping/_init_SJ.out.tab 

STAR --genomeDir star_mapping/ --readFilesIn test.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix star_mapping/_ --sjdbFileChrStartEnd star_mapping/_init_SJ.out.tab


