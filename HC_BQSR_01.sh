#!/bin/bash
# Date:
## 10/20/17
# Description
## Use surface fish genome
## BQSR: round 1
### variant calling: g.vcf, joint genotyping, separate SNPs and Indels- filter
### Run BQSR round 1 and 2, plot comparison
### *loops were run separtely

# Tools
GATK="/n/local/stage/gatk/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar"
## Variables
# reference genome
genome_dir="/n/core/Genomics/Analysis/Rohner/nicolas_rohner/nro3/genome"
ref_genome="$genome_dir/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.fna"

# MOLNG orders
molng_orders="MOLNG-1690 MOLNG-2003"
fc_MOLNG_1690="HV2W7BCXX HV32LBCXX"
fc_MOLNG_2003="HKN2LBCXY HKTGWBCXY HKTCHBCXY"
bwa_dir="/n/core/Genomics/Analysis/Rohner/nicolas_rohner/nro3/data/pre-processing/BWA"
bqsr_dir="/n/core/Genomics/Analysis/Rohner/nicolas_rohner/nro3/data/pre-processing/BQSR"
bqsr_round="BQSR_01"
## Prep Sample report for Read Groups
# get location of where the script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# GATK: Run HaplotypeCaller to make g.vcf files
for order in  ${molng_orders[*]};
do
    #set sample report to get barcode
    echo $order
    if [ $order == "MOLNG-1690" ]; then
	flow_cells=$fc_MOLNG_1690
    elif [ $order == "MOLNG-2003" ]; then
	flow_cells=$fc_MOLNG_2003
    fi
    # Get barcodes
    fc_arr=($flow_cells) # make array
    sr_01=$DIR/Sample_Report_$order-${fc_arr[0]}.csv
    # get all unique barcodes
    barcodes="$(tail -n +2 $sr_01 | cut -d "," -f7 | sort -u)"
    #base directory
    base_dir=$bqsr_dir
    ## Cluster information
    cpu="8"
    qsubBase="qsub -o $base_dir/$bqsr_round/g.vcf_files/$order -j y -terse -b y -wd $base_dir/$bqsr_round/g.vcf_files/$order -M dag@stowers.org -m ae"
    qsubJob="$qsubBase -l mem_free=10.6G,h_vmem=10.6G"
    qsubPE="$qsubJob -pe by_node $cpu"
    for barcode in  ${barcodes[*]};
    do
	#input files
	inf_01="$bwa_dir/$order/merged/${barcode}_realigned_reads.bam"
	#output file
	outf="$base_dir/$bqsr_round/g.vcf_files/$order/${barcode}.g.vcf"
	# GATK: Run HaplotypeCaller to make g.vcf files
	if [ -f $outf ]; then
	    echo "File '$outf' exists"
	else
	echo "java -Xmx4g -jar $GATK -T HaplotypeCaller -R $ref_genome --emitRefConfidence GVCF -I $inf_01 --genotyping_mode DISCOVERY --variant_index_type LINEAR --variant_index_parameter 128000 -stand_call_conf 30 -o $outf -nct 8"
        callVar="java -Xmx4g -jar $GATK -T HaplotypeCaller -R $ref_genome --emitRefConfidence GVCF -I $inf_01 --genotyping_mode DISCOVERY --variant_index_type LINEAR --variant_index_parameter 128000 -stand_call_conf 30 -o $outf -nct 8"
	varCallId=$(ssh aspen "$qsubPE -N gvcfBQSR01_$barcode '$callVar'")

	fi
    done
done


# joint call
wd="/n/core/Genomics/Analysis/Rohner/nicolas_rohner/nro3/data/pre-processing/BQSR/BQSR_01/g.vcf_files"
cpu="8"
qsubBase="qsub -o $wd -j y -terse -b y -wd $wd -M dag@stowers.org -m ae"
qsubJob="$qsubBase -l mem_free=20.6G,h_vmem=20.6G"
qsubPE="$qsubJob -pe by_node $cpu"

raw_vcf="joint_called_variants_raw.vcf"
if [ -f $wd/$raw_vcf ]; then
    echo "File '$raw_vcf' exists"
else
    joint_call="java -Xmx16g -jar $GATK -T GenotypeGVCFs -R $ref_genome -V gvcfs.list -o $raw_vcf -nt $cpu"
    echo "$joint_call"
    joint_callId=$(ssh aspen "$qsubPE -N jointcall_vars '$joint_call'")

fi

# NOW filter RAW vcf file
# Cluster info
#  Select Indels
#input files
inf_01="$wd/joint_called_variants_raw.vcf"
#output file
outf="$wd/joint_called_INDELs_raw.vcf"
if [ -f $outf ]; then
    echo "File '$outf' exists"
else
    echo "java -Xmx16g -Djava.io.tempdir=temp -jar $GATK -T SelectVariants -R $ref_genome -V $inf_01 -selectType INDEL -o $outf -nt 8"
    selectIn="java -Xmx16g -Djava.io.tempdir=temp  -jar $GATK -T SelectVariants -R $ref_genome -V $inf_01 -selectType INDEL -o $outf -nt 8"
    selectInId=$(ssh aspen "$qsubPE -N rawIn_BQSR01  -hold_jid $joint_callId '$selectIn'")
fi
# output of filtering
outf_2="$wd/joint_called_INDELs_filt.vcf"
if [ -f $outf_2 ]; then
    echo "File '$outf_2' exists"
else
    echo "java -Xmx16g -jar $GATK -T VariantFiltration -R $ref_genome -V $outf --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0 \" --filterName \"low_qual_indels\" -o $outf_2"
    filtIn="java -Xmx16g -jar $GATK -T VariantFiltration -R $ref_genome -V $outf --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0 \" --filterName \"low_qual_indels\" -o $outf_2"
    filtInId=$(ssh aspen "$qsubPE -N filtIn_BQSR01 -hold_jid $joint_callId, $selectInId '$filtIn'")
fi
# Select SNPs
#output file
outf="$wd/joint_called_SNPs_raw.vcf"
# GATK: select SNPs
if [ -f $outf ]; then
    echo "File '$outf' exists"
else
    echo "java -Xmx16g -jar $GATK -T SelectVariants -R $ref_genome -V $inf_01 -selectType SNP -o $outf -nt 8"
    selectSNPs="java -Xmx16g -jar $GATK -T SelectVariants -R $ref_genome -V $inf_01 -selectType SNP -o $outf -nt 8"
    selectSNPsId=$(ssh aspen "$qsubPE -N rawSNP_BQSR01 -hold_jid $joint_callId '$selectSNPs'")
fi
# output of filtering
outf_2="$wd/joint_called_SNPs_filt.vcf"
if [ -f $outf_2 ]; then
    echo "File '$outf_2' exists"
else
    echo "java -Xmx16g -jar $GATK -T VariantFiltration -R $ref_genome -V $outf --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0\" --filterName \"low_qual_snps\" -o $outf_2"
    filtSNPs="java -Xmx16g -jar $GATK -T VariantFiltration -R $ref_genome -V $outf --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0 \" --filterName \"low_qual_snps\" -o $outf_2"
    filtSNPsId=$(ssh aspen "$qsubPE -N filtSNPs_BQSR01 -hold_jid $joint_callId, $selectSNPsId '$filtSNPs'")
fi


# BQSR round 1 (FINALLY)
wd="/n/core/Genomics/Analysis/Rohner/nicolas_rohner/nro3/data/pre-processing/BQSR/BQSR_01"
cpu="8"
qsubBase="qsub -o $wd -j y -terse -b y -wd $wd -M dag@stowers.org -m ae"
qsubJob="$qsubBase -l mem_free=10.6G,h_vmem=10.6G"
qsubPE="$qsubJob -pe by_node $cpu"

for order in  ${molng_orders[*]};
do
    #set sample report to get barcode
    echo $order
    if [ $order == "MOLNG-1690" ]; then
	flow_cells=$fc_MOLNG_1690
    elif [ $order == "MOLNG-2003" ]; then
	flow_cells=$fc_MOLNG_2003
    fi
    # Get barcodes
    fc_arr=($flow_cells) # make array
    sr_01=$DIR/Sample_Report_$order-${fc_arr[0]}.csv
    # get all unique barcodes
    barcodes="$(tail -n +2 $sr_01 | cut -d "," -f7 | sort -u)"
    for barcode in  ${barcodes[*]};
    do
	#input files
	inf_01="$bwa_dir/$order/merged/${barcode}_realigned_reads.bam"
	#output file
	outf="$wd/$order/${barcode}_BQSR_01.table"
	# GATK: Run BaseRecalibrator round 1
	if [ -f $outf ]; then
	    echo "File '$outf' exists"
	else
	    BQSR_01="java -Xmx8g -jar $GATK -T BaseRecalibrator  -R $ref_genome -I $inf_01 -knownSites $wd/g.vcf_files/joint_called_SNPs_filt.vcf  -knownSites $wd/g.vcf_files/joint_called_SNPs_filt.vcf -o $outf -nct 8"
	    echo $BQSR_01
	    BQSR_01Id=$(ssh aspen "$qsubPE -N ${barcode}_BQSR01 '$BQSR_01'") 
	fi
	#output file
	outf_2="$wd/$order/${barcode}_BQSR_02.table"
	# GATK: Run BaseRecalibrator round 2
	if [ -f $outf_2 ]; then
	    echo "File '$outf_2' exists"
	else
	    BQSR_02="java -Xmx8g -jar $GATK -T BaseRecalibrator  -R $ref_genome -I $inf_01 -knownSites $wd/g.vcf_files/joint_called_SNPs_filt.vcf  -knownSites $wd/g.vcf_files/joint_called_SNPs_filt.vcf -BQSR $outf -o $outf_2 -nct 8"
	    echo $BQSR_02
	    BQSR_02Id=$(ssh aspen "$qsubPE -N ${barcode}_BQSR02 -hold_jid $BQSR_01Id '$BQSR_02'")
	fi
	#output file
	outf_3="$wd/$order/${barcode}_BQSR01_plot.pdf"
	# GATK: Run: Analyze Covariates- plot
	if [ -f $outf_3 ]; then
	    echo "File '$outf_3' exists"
	else
	    plot="java -Xmx8g -jar $GATK -T AnalyzeCovariates -R $ref_genome -before $outf -after $outf_2 -plots $outf_3"
	    echo $plot
	    #plotId=$(ssh aspen "$qsubPE -N ${barcode}_plot -hold_jid $BQSR_01Id, $BQSR_02Id '$plot'")
	    plotId=$(ssh aspen "$qsubPE -N ${barcode}_plot '$plot'")
	fi
	
    done
done


