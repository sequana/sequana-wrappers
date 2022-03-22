Data obtained from biomics.pasteur.fr project B4525.
We used sample 1 (IP) and 6 (control) mapped on Vibrio Cholera.
Two peaks were present at postion ~ 544700 and 548288. We selected the
mapped reads in the region 500000-600000 and subsampled selecting 
20% of the reads. Applying macs3 as follows should identify those peaks.::

    macs3 callpeak -B --SPMR -t 1.sorted.bam  -c 6.sorted.bam -g 4000000 -n 1_vs_6 --bw 300 --keep-dup all -q 0.05 -f BAMPE  --outdir macs3/narrow 

Selection of the reads was done as follows::

    samtools view input.sorted.bam "NC_002505.1:535000-553000" -s 0.2 -b  > output.sorted.bam 

