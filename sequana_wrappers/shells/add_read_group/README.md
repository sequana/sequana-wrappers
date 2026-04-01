The old conda wrapper uses whatever samtools was available when the environment was first created (likely an older version). Newer samtools versions are stricter about coordinate sorting requirements for samtools index — they now enforce that:
 1. The BAM must be coordinate-sorted 
 2. Unmapped reads (NO_COOR) must all be at the end in a single block
 
 The input test.bam has no SO:coordinate tag in its @HD header (neither BAM does — both headers are identical and lack an SO field), so picard preserves the original order by default. Older samtools would index this anyway; newer samtools reject it with the NO_COOR reads not in a single block error 

 The container (e.g. sequana_tools_26.1.14.img) ships a newer samtools that enforces this, while the conda environment resolved an older 
 samtools version that was permissive. Adding -SO coordinate to picard makes the output properly sorted before indexing, which is the correct fix regardless of samtools version. 
