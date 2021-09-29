# Documentation

Counts reads mapped to genomic regions using featureCounts (subread)

Reference: http://bioinf.wehi.edu.au/featureCounts/

Required input:

- the input BAM file
- the GFF annotation file

Required output:

- the counts file (tabulated file)

Required Parameters:

- feature: a valid feature to be found in the GFF file
- attribute: a valid attribute to be found in the GFF file
- strandness: a value equal to 0, 1, or 2 
- options: any other parameters accepted by featureCounts

Log:

- a log file

# Configuration

    feature_counts:
        options:
        feature: gene    # could be exon, mRNA, etc
        attribute: ID    # could be ID, gene_id, etc
        strandness: 0    # set to 0,1 or 2


# Example

    rule feature_counts:
        input:
            bam="{sample}/bamfile/{sample}.sorted.bam,
            gff="genome.gff"
        output:
            counts="{sample}/feature_counts/{sample}_feature.out",
            summary="{sample}/feature_counts/{sample}_feature.out.summary"
        params:
            options=config["feature_counts"]["options"],
            feature=config["feature_counts"]["feature"],
            attribute=config["feature_counts"]["attribute"],
            strandness=config["feature_counts"]["strandness"]
        threads: 
            config["feature_counts"]['threads']
        log:
            "{sample}/feature_counts/feature_counts.log"
        wrapper:
            "main/wrappers/feature_counts"

