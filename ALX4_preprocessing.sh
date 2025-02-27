#!bin/bash

# (( TRACE )) && set -x
# set -eo pipefail
# trap 'Something went wrong. Script error on line #$LINENO."' ERR
# Processes files from GEO to mapping step


## Grab important file paths
BASEDIR="$(pwd)"
CODEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
GENOMES="/data/campbell-gebelein-lab/Genomes/"


## Load modules
module load fastqc/0.11.2
module load sratoolkit/3.0.0
module load trimgalore/0.6.6
module load bowtie2/2.3.4.1
module load samtools/1.9.0
module load homer/4.9
module load picard/2.18.22
module load bedtools/2.27.0

## Grab necessary arguments
pflag=''
while getopts ":a:b:c:d:e:g:" opt; do
	case $opt in
	a)
		REP1_SRR=${OPTARG}
		echo "$REP1_SRR"
		[[ "$REP1_SRR" == SRR* ]] || 
		echo "	[-a] Enter SRR accession for replicate 1" 
		;;
	b)
		REP2_SRR=${OPTARG}
		[[ "$REP2_SRR" == SRR* ]] || 
		echo "	[-b] Enter SRR accession for replicate 2 if applicable"
		;;
	c)
		INPUT1_SRR=${OPTARG}
		[[ "$INPUT1_SRR" == SRR* ]] || 
		echo "	[-c] Enter SRR accession for input replicate 1 if applicable"
		;;
	d)
		INPUT2_SRR=${OPTARG}
		[[ "$INPUT2_SRR" == SRR* ]] || 
		echo "	[-d] Enter SRR accession for input replicate 2 if applicable"
		;;
	e)
		EXP_name=${OPTARG}
		;;
	g)
		genome=${OPTARG}
		[[ "$genome" == hg* || "$genome" == mm* || "$genome" == dm* ]] ||
		echo "	[-g] Enter a valid genome to align to such as hg19, mm10...
		Default is mm10"
		;;
	?)
		echo 	
		"Script usage:
		
		Summary:
		Processes ChIP-seq or CUT&RUN data for dimer motifs analysis. 
		Accepts up to 2 IP and input SRR ID's and processes files up 
		to de novo motif analysis and COSMO analysis.
		
		Pipeline:
		1. Downloads files with fastq-dump
		2. Trims with trim-galore (powered by cut-adapt)
		3. FASTQC check
		4. Reads aligned with Bowtie2
		5. Duplicates removed with Picard
		6. Calls peaks with Homer
		7. Filters out peaks that fall into SRRODE blacklist
		8. Performs de novo motif analysis with Homer
		
		Program file inputs: 
		[-a] Enter SRR accession for replicate 1
		[-b] Enter SRR accession for replicate 2 if applicable
		[-c] Enter SRR accession for input replicate 1 if applicable
		[-d] Enter SRR accession for input replicate 2 if applicable
		[-e] Enter Experiment name. Default is ANALYSIS_FILE
		[-g] Enter a valid genome to align to such as hg19, mm10..."
		exit 1
		;;
	esac
done

shift $((OPTIND-1))
	
## Set up working directory
if [ ! -z ${EXP_name+x} ]
then
	mkdir -p "$EXP_name"; 
else
	mkdir 'ANALYSIS_FILE'; EXP_name='ANALYSIS_FILE'
fi

## Check for index files for genome; prep if required
if [ -z ${genome+x} ]
then
	echo "No genome selected. Will align to mm10. If this is not correct,
	CTRL+C now!"
	genome=mm10
fi

if [ -d "$GENOMES"/"$genome"/BowtieIndex ]
then
	echo "Bowtie index files exist!"
else
	echo "Generating Bowtie index files..."
	cd "$GENOMES"/"$genome" || echo 'Is the genome downloaded?'
	if [ -e "$genome".fa.gz ]
	then
		gzip -d "$genome".fa.gz
	fi
	bowtie2-build -f "$genome".fa "$genome"
	mkdir BowtieIndex
	mv *.bt2 BowtieIndex
fi
	
# Process individual fastq's
for File in "$REP1_SRR" "$REP2_SRR" "$INPUT1_SRR" "$INPUT2_SRR" 
do

	if [[ "$File" == SRR* ]]
	then
		cd "$BASEDIR"/"$EXP_name"
		echo "
	
		Starting to process $File...
		"
		# Download files
		mkdir "$File"; 
		cd "$File"
		fastq-dump --split-files --gzip "$File"
		pe=$(find ./*_2.fastq.gz 2> /dev/null | wc -l) || pe=0
		
		if [ "$pe" -gt 0 ]
		then
			### PAIRED END SEQUSRRING ###
			
			# ## Trim sequSRRes and QC	
			echo "Trimming $REP1_SRR..."
			trim_galore --paired --cores 4 \
				"$File"_1.fastq.gz "$File"_2.fastq.gz > \
				"$File"_trim_galore_run.log
			
			mkdir Fastqc; mv -f *fastqc* Fastqc
			
			# ## Align to genome
			echo "Aligning $REP1_SRR to genome..."
			bowtie2 -x "$GENOMES"/"$genome"/BowtieIndex/"$genome" \
				--local --very-sensitive-local --no-unal --no-mixed \
				--no-discordant --threads 4 \
				-1 "$File"_1_val_1.fq.gz -2 "$File"_2_val_2.fq.gz \
				-S "$File".sam 

		else
			### SINGLE END SEQUSRRING ###
			echo "Trimming $REP1_SRR..."
			trim_galore --cores 4 \
				"$File"_1.fastq.gz > \
				"$File"_1_trim_galore_run.log
			
			mkdir Fastqc; mv -f *fastqc* Fastqc
			
			echo "Aligning $REP1_SRR to genome..."
			bowtie2 -x "$GENOMES"/"$genome"/BowtieIndex/"$genome" \
				--local --very-sensitive-local --no-unal --no-mixed \
				--no-discordant --threads 4 \
				-U "$File"_1_trimmed.fq.gz \
				-S "$File".sam 

		fi
				
		samtools view -b -o "$File".bam "$File".sam && rm "$File".sam
		samtools sort -o "$File".sort.bam "$File".bam && rm "$File".bam
		picard MarkDuplicates I="$File".sort.bam O="$File".marked.bam REMOVE_DUPLICATES=false M="$File"_rmdup.log
		samtools view -h -b -q 30 -F 1024 -o "$File".rmdup.bam "$File".marked.bam
		samtools index "$File".rmdup.bam 
		
	fi
done
 
