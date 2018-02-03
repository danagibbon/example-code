#!/bin/bash
# Date:
## 10/19/17
# Description
## Use surface fish genome
## Realign around INDELs

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
## Prep Sample report for Read Groups
# get location of where the script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# GATK: RealignerTargetCreator
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
    base_dir=$bwa_dir
    ## Cluster information
    cpu="4"
    qsubBase="qsub -o $base_dir/$order/merged -j y -terse -b y -wd $base_dir/$order -M dag@stowers.org -m ae"
    qsubJob="$qsubBase -l mem_free=10.6G,h_vmem=10.6G"
    qsubPE="$qsubJob -pe by_node $cpu"
    if [ $order == "MOLNG-1690" ]; then
	for barcode in  ${barcodes[*]};
	do
	    #input files
	    inf_01="./${fc_arr[0]}/$barcode-1.picardIntermediate.sorted.bam"
	    inf_02="./${fc_arr[0]}/$barcode-2.picardIntermediate.sorted.bam"
	    inf_03="./${fc_arr[1]}/$barcode-1.picardIntermediate.sorted.bam"
	    inf_04="./${fc_arr[1]}/$barcode-2.picardIntermediate.sorted.bam"
	    #output file
	    outf="$base_dir/$order/merged/${barcode}_target_intervals.list"
	    # GATK: RealignerTargetCreator
	    if [ -f $outf ]; then
		echo "File '$outf' exists"
	    else
		echo "java -jar $GATK -T RealignerTargetCreator -R $ref_genome -I $inf_01 -I $inf_02 -I $inf_03 -I $inf_04 -o $outf -nt 4"
		RealignTC="java -jar $GATK -T RealignerTargetCreator -R $ref_genome -I $inf_01 -I $inf_02 -I $inf_03 -I $inf_04 -o $outf -nt 4"
		RealignTCId=$(ssh aspen "$qsubPE -N RealignTC_$barcode '$RealignTC'")
	    fi
	    # GATK: IndelRealigner
	    outfin="$base_dir/$order/merged/${barcode}_realigned_reads.bam"
	    if [ -f $outfin ]; then
		echo "File '$outfin' exists"
	    else
		echo "java -Xmx4g -jar $GATK -T IndelRealigner -R $ref_genome -I $inf_01 -I $inf_02 -I $inf_03 -I $inf_04 -targetIntervals $outf -o $outfin"
		IndelReal="java -Xmx4g -jar $GATK -T IndelRealigner -R $ref_genome -I $inf_01 -I $inf_02 -I $inf_03 -I $inf_04 -targetIntervals $outf -o $outfin"
		IndelRealId=$(ssh aspen "$qsubPE -N RealignTC_$barcode -hold_jid $RealignTCId  '$IndelReal'")
	    fi
	done
    elif [ $order == "MOLNG-2003" ]; then
	for barcode in  ${barcodes[*]};
	do
	    #input files
	    inf_01="./${fc_arr[0]}/$barcode-1.picardIntermediate.sorted.bam"
	    inf_02="./${fc_arr[0]}/$barcode-2.picardIntermediate.sorted.bam"
	    inf_03="./${fc_arr[1]}/$barcode-1.picardIntermediate.sorted.bam"
	    inf_04="./${fc_arr[1]}/$barcode-2.picardIntermediate.sorted.bam"
	    inf_05="./${fc_arr[2]}/$barcode-1.picardIntermediate.sorted.bam"
	    inf_06="./${fc_arr[2]}/$barcode-2.picardIntermediate.sorted.bam"
	    #output file
	    outf="$base_dir/$order/merged/${barcode}_target_intervals.list"
	    # GATK: RealignerTargetCreator
	    if [ -f $outf ]; then
		echo "File '$outf' exists"
	    else
		echo "java -jar $GATK -T RealignerTargetCreator -R $ref_genome -I $inf_01 -I $inf_02 -I $inf_03 -I $inf_04 -I $inf_05 -I $inf_06 -o $outf -nt 4"
		RealignTC="java -jar $GATK -T RealignerTargetCreator -R $ref_genome -I $inf_01 -I $inf_02 -I $inf_03 -I $inf_04 -I $inf_05 -I $inf_06 -o $outf -nt 4"
		RealignTCId=$(ssh aspen "$qsubPE -N RealignTC_$barcode '$RealignTC'")
	    fi
	    # GATK: IndelRealigner
	    outfin="$base_dir/$order/merged/${barcode}_realigned_reads.bam"
	    if [ -f $outfin ]; then
		echo "File '$outfin' exists"
	    else
		echo "java -Xmx4g -jar $GATK -T IndelRealigner -R $ref_genome -I $inf_01 -I $inf_02 -I $inf_03 -I $inf_04 -I $inf_05 -I $inf_06 -targetIntervals $outf -o $outfin"
		IndelReal="java -Xmx4g -jar $GATK -T IndelRealigner -R $ref_genome -I $inf_01 -I $inf_02 -I $inf_03 -I $inf_04 -I $inf_05 -I $inf_06 -targetIntervals $outf -o $outfin"
		IndelRealId=$(ssh aspen "$qsubPE -N RealignTC_$barcode -hold_jid $RealignTCId  '$IndelReal'")
	    fi
	done
    fi
done

