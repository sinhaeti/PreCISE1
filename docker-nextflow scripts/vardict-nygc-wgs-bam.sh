#!/bin/bash

#SBATCH --job-name=PreCISE1_Script
#SBATCH --partition=panda_physbio
#SBATCH --begin=now
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=28G
#SBATCH --time=3-00:00:00
#SBATCH --exclude=node002.panda.pbtech,node007.panda.pbtech,node010.panda.pbtech
source ~/.bashrc

CONFIG_FILE=$1
INPUT_DIR=$2
OUTPUT_DIR=$3
AF_MIN=$4
SAMPLE_PREFIX=$5
FADE_BOOLEAN=$6

######### TRACK JOB #########
#Module Loading and Sourcing
DATE=$(date '+%d/%m/%Y %H:%M:%S');
echo "$DATE"

# CONFIG FILE
source ${CONFIG_FILE}

## CHANGE TO DIFFERENT FOLDER (TMP)
cd $TMPDIR
TEMP=$(basename ${OUTPUT_DIR})
mkdir ${TEMP}; cd ${TEMP}

####### FILE NAME PREFIX INFO  #########
echo "Sample name is ${SAMPLE_PREFIX}"
echo  "=============================="

# Read 1-3
READ1=${INPUT_DIR}/${SAMPLE_PREFIX}_R1_001.fastq.gz
UMI_READ=${INPUT_DIR}/${SAMPLE_PREFIX}_R2_001.fastq.gz
READ2=${INPUT_DIR}/${SAMPLE_PREFIX}_R3_001.fastq.gz

# New files with UMI attached
R1_UMI=${SAMPLE_PREFIX}_R1_UMI.fastq.gz
UMI_OUT=${SAMPLE_PREFIX}_UMI.fastq.gz
R2_UMI=${SAMPLE_PREFIX}_R2_UMI.fastq.gz

# New files with adapter trimming and error correction
R1_TRIM=${OUTPUT_DIR}/${SAMPLE_PREFIX}_R1_trim.fastq.gz
R2_TRIM=${OUTPUT_DIR}/${SAMPLE_PREFIX}_R2_trim.fastq.gz
HTML=${SAMPLE_PREFIX}_fastp.html
JSON=${OUTPUT_DIR}/${SAMPLE_PREFIX}_fastp.json

# Aligned reads (bwa mem and gencore output)
READGROUP="@RG\tID:${SAMPLE_PREFIX}\tLB:${SAMPLE_PREFIX}\tPL:illumina\tSM:${SAMPLE_PREFIX}"
BWA_SAM=${SAMPLE_PREFIX}_precorrected.sam
BWA_BAM=${SAMPLE_PREFIX}_precorrected.bam
GENCORE_BAM=${OUTPUT_DIR}/${SAMPLE_PREFIX}_corrected.bam

# Picard HSMetrics Output File
HSMETRICS=${OUTPUT_DIR}/${SAMPLE_PREFIX}_hs_metrics.txt
INSERTSIZEMETRICS=${OUTPUT_DIR}/${SAMPLE_PREFIX}_insert_size_metrics.txt 

# Fade output file
FADE_BAM=${SAMPLE_PREFIX}_corrected_fadetmp.bam
FADE_BAM_SORT=${SAMPLE_PREFIX}_corrected_fadetmp_sort.bam
FADE_BAM_OUT=${SAMPLE_PREFIX}_corrected_fadetmp2.bam
FADE_BAM_FINAL=${SAMPLE_PREFIX}_corrected_fade.bam	

# Vardict output files
VARDICT_AUX=${SAMPLE_PREFIX}_aux.txt
VARDICT_TMP1_VCF=${SAMPLE_PREFIX}_varcall.tmp1.vcf.gz
VARDICT_TMP2_VCF=${SAMPLE_PREFIX}_varcall.tmp2.vcf.gz
VARDICT_TMP3_VCF=${SAMPLE_PREFIX}_varcall.tmp3.vcf.gz
VARDICT_FINAL_VCF=${SAMPLE_PREFIX}_varcall.marked.sort.vcf.gz

# VEP output file
VEP_OUT=${SAMPLE_PREFIX}_varcall.marked.sort.vep.annot.vcf
VEP_OUT_GZ=${SAMPLE_PREFIX}_varcall.marked.sort.vep.annot.vcf.gz
VEP_COSMIC=${SAMPLE_PREFIX}_varcall.marked.sort.vep.annot.cosmic.vcf
VEP_TCGA=${OUTPUT_DIR}/${SAMPLE_PREFIX}_varcall.marked.sort.vep.annot.cosmic.tcga.vcf.gz

# samtools depth
SAMTOOLS_DEPTH=${OUTPUT_DIR}/${SAMPLE_PREFIX}.regions.bed
SAMTOOLS_DEPTH_IG=${OUTPUT_DIR}/${SAMPLE_PREFIX}.ig.regions.bed
## Check if file has alreadt been processed
#if [ -f ${OUTPUT_DIR}/${FINAL_OUTPUT_VCF} ]
#then
    #echo "Sample ${SAMPLE_PREFIX} has already been processed. The script will move onto the next file."
    ##echo "final file will be overwritten."
    ##rm -rf ${OUTPUT_DIR}/${FINAL_OUTPUT_VCF}
    #continue
#fi

####### ALIGNMENT AND VARIANT CALLING  #########
#STEP 1: variant calling with VarDict
    
    
    echo  "=============================="
    echo "(STEP 1) Call Variants with VarDict java for sample ${SAMPLE_PREFIX}"
    conda activate samtools
    samtools view -b -L ${TARGETS} ${bam_from_nygc}> ${GENCORE_BAM}
    samtools index ${GENCORE_BAM}

    VarDict \
                -th 10 \
                -G $REF_FASTA \
                -f $AF_MIN \
                -N ${SAMPLE_PREFIX} \
                -b ${GENCORE_BAM} \
                -z 1 \
                -k 1 \
                -r 2 \
                -x $EXTEND_BED \
                -I $MAX_INSERT \
                -Q $MIN_MAPPING_QUALITY \
                -c 1 \
                -S 2 \
                -E 3 \
                -g 4 \
                ${TARGETS} > ${VARDICT_AUX}

    echo "(STEP 1B): Filter strand bias using VarDict supplementary tools"
    spack load /5vatd3o
    cat ${VARDICT_AUX} |
        ${vardict_script_dir}/teststrandbias.R | \
        ${vardict_script_dir}/var2vcf_valid.pl \
        -N ${SAMPLE_PREFIX} \
        -A \
        -f $AF_MIN > output_var2vcf_valid_{SAMPLE_PREFIX}.vcf

        conda activate tabix 
        bgzip -c output_var2vcf_valid_{SAMPLE_PREFIX}.vcf > ${VARDICT_TMP1_VCF}
    tabix -f ${VARDICT_TMP1_VCF}



    echo  "=============================="
    echo "(STEP 2): Annotate variants with bcftools"
    conda activate bcftools
    bcftools annotate \
        -o ${VARDICT_TMP2_VCF} \
        -O z ${VARDICT_TMP1_VCF}
    conda deactivate
    tabix -f ${VARDICT_TMP2_VCF}

    echo  "=============================="
    echo "(STEP 3): Filter Variants with vcflib and bcftools"
    conda activate bcftools
    zcat ${VARDICT_TMP2_VCF} | \
        ${vcflib_dir}vcfglxgt | \
        ${vcflib_dir}vcfbreakmulti | \
        ${vcflib_dir}vcffixup - | \
        ${vcflib_dir}vcfentropy -w 10 -f $REF_FASTA | \
        sed 's/=-0/=0/g' | \
        bcftools filter -e "SBF[0] < ${SB_PVAL} | ODDRATIO[0] > ${SB_ODDRATIO} | SBF[1] < ${SB_PVAL} | ODDRATIO[1] > ${SB_ODDRATIO}" -s "strictSB" -m + | \
        bcftools filter -e "INFO/DP < ${MIN_DEPTH}" -s "d${MIN_DEPTH}" -m + | \
        bcftools filter -e "INFO/VD < ${MIN_VAR_READS}" -s "hi.v.${MIN_VAR_READS}" -m + | \
        bcftools filter -e '((EntropyLeft < 1 | EntropyRight < 1 | EntropyCenter < 1)) & TYPE!="SNP" & HIAF < 0.05' -s "IndelEntropy" -m + > ${VARDICT_TMP3_VCF}

        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP\t%AF\t%AD]\n' ${VARDICT_TMP3_VCF} > ${SAMPLE_PREFIX}bcftools_output.txt
    
    echo  "=============================="
    echo "(STEP 4): Annotate Variants with ANNOVAR"
        perl /home/zil4004/annovar/convert2annovar.pl -format vcf4 ${VARDICT_TMP3_VCF} > {SAMPLE_PREFIX}vardict_TMP3.avinput
        perl /home/zil4004/annovar/table_annovar.pl {SAMPLE_PREFIX}vardict_TMP3.avinput /home/zil4004/annovar/humandb/ -buildver hg38 -out {SAMPLE_PREFIX} -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a,gnomad211_genome -operation g,r,f,f,f,f -nastring . -csvout

        echo  "=============================="
        echo "Variant calling already completed"
    fi


    echo " Script $0 is complete for sample ${SAMPLE_PREFIX}"
    DATE=$(date '+%d/%m/%Y %H:%M:%S');
    echo "$DATE"
