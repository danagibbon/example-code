# example-code

## Bash pipeline:

This set of scripts are  part of the pre-processing steps before running variant calling in GATK

* file 1:  [pre-processing-pipeline.sh](example-code/pre-processing-pipeline.sh)
    + align reads + add read groups
    + mark duplicates
* file 2:  [prepp-realign-indels.sh](example-code/prepp-realign-indels.sh)
    + Realign around INDELs
* file 3:  [HC_BQSR_01.sh](example-code/HC_BQSR_01.sh)
    + BQSR round 1
    + variant calling
    + joint genotyping
    + filter SNPs and Indels
    + BaseRecalibrator 1 & 2
    + AnalyzeCovariates to get plots


## RMarkdown report example:

