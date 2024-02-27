
module load samtools/1.18.0
module load bedtools

SRR="$1"

SRRID="${SRR%.*}"

# ## Run read statistics on file
samtools stats "$SRR" > "$SRR"_alignment_stats.txt

# ## Filter out reads greater than 120 bp (CUT&RUN only)
samtools view -e 'length(seq)<120'  -O BAM -o "$SRRID"_lt120.bam "$SRRID".bam 
  
samtools index "$SRRID"_lt120.bam
 
## Generate 5' coverage files
bedtools genomecov -bga -5 -strand - -ibam "$SRRID"_lt120.bam > "$SRRID"_minus5.cov
bedtools genomecov -bga -5 -strand + -ibam "$SRRID"_lt120.bam > "$SRRID"_plus5.cov

## Grab 5' only coverage around slopped ALX4 CUT&RUN peaks
bedtools intersect -wa -a "$SRRID"_minus5.cov \
-b /data/campbell-gebelein-lab/zz_BC/ALX4_TWIST/TWIST1FV/All_ALX4_CR_peaks_slopped150bp.bed >\
"$SRRID"_minus5_inPeaks.cov; rm -f "$SRRID"_minus5.cov

bedtools intersect -wa -a "$SRRID"_plus5.cov \
-b /data/campbell-gebelein-lab/zz_BC/ALX4_TWIST/TWIST1FV/All_ALX4_CR_peaks_slopped150bp.bed >\
"$SRRID"_plus5_inPeaks.cov; rm -f "$SRRID"_plus5.cov

## Generate coverage files
bedtools genomecov -bga -strand - -ibam "$SRRID".bam > "$SRRID"_minus.cov
bedtools genomecov -bga -strand + -ibam "$SRRID".bam > "$SRRID"_plus.cov

## Grab coverage around slopped ALX4 CUT&RUN peaks
bedtools intersect -wa -a "$SRRID"_minus.cov \
-b /data/campbell-gebelein-lab/zz_BC/ALX4_TWIST/TWIST1FV/All_ALX4_CR_peaks_slopped150bp.bed >\
"$SRRID"_minus_inPeaks.cov; rm -f "$SRRID"_minus.cov

bedtools intersect -wa -a "$SRRID"_plus.cov \
-b /data/campbell-gebelein-lab/zz_BC/ALX4_TWIST/TWIST1FV/All_ALX4_CR_peaks_slopped150bp.bed >\
"$SRRID"_plus_inPeaks.cov; rm -f "$SRRID"_plus.cov