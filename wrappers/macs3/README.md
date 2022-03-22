# Documentation

This wrapper calls peaks in ChIP-seq and ATAQ-seq like sequencing runs.

**Required input:**

- 1 or several sorted BAM files (IPs)
- 1 or several sorted BAM files (controls)

You must set the names IP and control in the input list (see example below).

**Required output:**

- A set of files, created by macs3

**The output files are processed to extract the prefix to be used by macs3**. 
Their prefixes must be the same. You may provide only one output directory or list
the different output files created by macs3.

**Log:**

- a log file with stdout/stderr

**Optional parameters:**

* qvalue is optional and defaults to 0.05
* options is optional and defaults to '--keep-dup all'

**Required parameters:**

* bandwidth: 300   # default bandwidth option --bw from macs3
* broad_cutoff: 0.05 is used and required only if the mode is set to 'broad'
* genome_size: this is the mappable size of your genome and is a very important parameter 
* mode: 'broad'    # broad or narrow. compulsary parameter
* paired: True  is the data paired or not. Using sequana manager, you can simply set it to manager.paired
* prefix: tag name. All files will be saved in an output directory according to your output files. However, you alos
  need to provide a prefix to all the files that will be generated. Set to e.g. 'macs3' or the name of a comparison.


# Configuration


    ############################################################################################
    # macs3 peak caller
    #
    # bandwidth: 300           # default bandwidth option --bw from macs3
    # broad_cutoff: 0.05
    # genome_size: 4000000  # the mappable size of your genome. 
    # mode: 'broad'    # broad or narrow. compulsary 
    # options: --keep-dup all      ## we keep all duplicates and let macs3 do the job
    # paired: True 
    # prefix: tag name 
    # qvalue: 0.05



# Example

For a set of inputs, use a list if required. The outputs must contain the world
"narrow" or "broad". Based on this word, macs will be called in narrow and/or broad
mode.


    rule macs3:
        input:
            IP = ["1.bam", "2.bam"],
            control = "3.bam",
        output:
            "macs3/narrow/tag_peaks.xls",
        params:
            bandwidth       = 300,
            broad_cutoff    = 0.05,
            genome_size     = 4000000,
            mode            = "narrow",
            options         = " --keep-dup all ",
            paired          = True,
            prefix          = "tag",
            qvalue          = 0.05
        log:
            out.log
        wrapper:
            "main/wrappers/macs3"
