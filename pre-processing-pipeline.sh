#!/bin/bash
# Date:
## 10/17/17
# Description
# Use surface fish genome
# align reads to genome with BWA, mark duplicates

# Tools
GATK="/n/apps/CentOS7/bin/GenomeAnalysisTK.jar"
picard="/n/apps/CentOS7/bin/picard.jar"

## Variables
# reference genome
genome_dir="/n/core/Genomics/Analysis/Rohner/nicolas_rohner/nro3/genome"
ref_genome="$genome_dir/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.fna"

# Prep for BWA
# make BWA index if it does not exist already
BWA_index="$genome_dir/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.fna.pac"
if [ -f $BWA_index ]; then
    echo "File '$BWA_index'  exists"
else
    bwa index $ref_genome
fi

# make fasta file index if it does not exist
fa_index="$genome_dir/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.fna.fai"
if [ -f $fa_index ]; then
    echo "File '$fa_index'  exists"
else
    samtools faidx $ref_genome
fi

# make sequence dictionary index if it does not exist
ref_dict="$genome_dir/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.dict"
if [ -f $ref_dict ]; then
    echo "File '$ref_dict'  exists"
else
    java -jar $picard CreateSequenceDictionary  REFERENCE=$ref_genome OUTPUT=$ref_dict
fi

# MOLNG orders
molng_orders="MOLNG-1690 MOLNG-2003"
fc_MOLNG_1690="HV2W7BCXX HV32LBCXX"
fc_MOLNG_2003="HKN2LBCXY HKTGWBCXY HKTCHBCXY"    
bwa_dir="/n/core/Genomics/Analysis/Rohner/nicolas_rohner/nro3/data/pre-processing/BWA"
## Prep Sample report for Read Groups
# get location of where the script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# move Sample_Report.csv to DIR
for order in  ${molng_orders[*]};
do
    #set LIMS order and get associated flow cells 
    limsDir="/n/analysis/Rohner/kzg/$order"
    echo $order
    if [ $order == "MOLNG-1690" ]; then
	fc_s=$fc_MOLNG_1690
    elif [ $order == "MOLNG-2003" ]; then
	fc_s=$fc_MOLNG_2003
    fi
    # Loop through flow cells associated with MOLNG order 
    for flow_c in ${fc_s[*]};
    do
	sr_1="$limsDir/$flow_c/Sample_Report.csv"
	sr_01=$DIR/Sample_Report_$order-$flow_c.csv
	if [ -f $sr_01 ]; then
	    echo "File '$sr_01'  exists"
	else
	    cp $sr_1 $DIR
	    mv $DIR/Sample_Report.csv $sr_01
	fi

	# make BWA directories
	BWA_dir="$bwa_dir/$order/$flow_c"
	if [ -d $BWA_dir ]; then
	    echo "Directory '$BWA_dir'  exists"
	else
	    mkdir $BWA_dir
	fi
	# Now run while loop on sample reports 
	while read line; do
	    ## split the lines by comma, store in a array
	    oIFS=$IFS # input file separator
	    IFS=,
	    line=( $line )
	    IFS=$oIFS

	    read_pair=${line[8]}
	    ## Grab the first item, should be a fastq file
	    file=${line[0]}
	    if [ $read_pair == 1 ]; then
		# Grab the first item, should be a fastq file
		file=${line[0]}
		## Create the fastq file path
		fq=$limsDir/$flow_c/$file
		echo $fq

		## Continue only if the fastq exists
		[ ! -e $fq ] && continue

		## Extract some info from the report.csv
		lane=${line[3]}
		BC=${line[6]}
		lib=${line[5]}
		sample=${line[4]}
		id=$flow_c.lane-$lane.$BC
		rg="@RG\tID:$id\tSM:$sample\tPL:illumina\tLB:$lib\tPU:$flow_c-$BC.$lane"

		if [ $file == s_1_1_$BC.fastq.gz ]; then
		    read_2="s_1_2_$BC.fastq.gz"
		else
		    read_2="s_2_2_$BC.fastq.gz"
		fi

		# run BWA
		output="$BWA_dir/$BC-$lane"
		if [ -f $output.sam ]; then
		    echo "File '$output.sam'  exists"
		else
		    bwacmd="bwa mem -t 8 -M -R $rg $ref_genome $limsDir/$flow_c/$file $limsDir/$flow_c/$read_2 > $output.sam"
		    echo $bwacmd
		    bwa mem -t 8 -M -R "$rg" $ref_genome $limsDir/$flow_c/$file $limsDir/$flow_c/$read_2 > $output.sam
		    # samtools
		    samcmd="samtools view -bS $output.sam | samtools sort - -o $output.sorted.bam"
		    echo $samcmd
		    samtools view -bS $output.sam | samtools sort - -o $output.sorted.bam
		fi

		# Picard: MarkDuplicates
		if [ -f $output.metrics.dat ]; then
		    echo "File '$output.metrics.dat' exists"
		else
		    dedupcmd="java -Xmx4g -jar $picard MarkDuplicates I=$output.sorted.bam O=$output.picardIntermediate.sorted.bam METRICS_FILE=$output.metrics.dat MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 CREATE_INDEX=TRUE"
		    echo $dedupcmd
		    $dedupcmd
		fi

	    fi	    

	done < $sr_01	

    done
    
done
    
