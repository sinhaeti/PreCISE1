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
UMI_BOOLEAN=$7

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
GERMLINE_VEP=${OUTPUT_DIR}/${SAMPLE_PREFIX}_annotated_germline_variants.vcf 

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

if [ ! -f "${R1_TRIM}" ]; then
    # STEP 1: add UMI to R1 (reads1) ... use 1 thread only to ensure fastq records stay synchronous. w=1
    echo "(STEP 1) ADD UMI to R1 and R2 sequence identifiers for sample ${SAMPLE_PREFIX}"
    conda activate fastp
        
    if [[ ${UMI_BOOLEAN} == "no_umi" ]]; then
        R1_NO_UMI=${INPUT_DIR}/${SAMPLE_PREFIX}_R1_001.fastq.gz
        R2_NO_UMI=${INPUT_DIR}/${SAMPLE_PREFIX}_R2_001.fastq.gz

        fastp \
        --correction \
        -w ${THREAD_NUM} \
        --in1 ${R1_NO_UMI} \
        --in2 ${R2_NO_UMI} \
        --out1 ${R1_TRIM} \
        --out2 ${R2_TRIM} \
        --html ${HTML} \
        --json ${JSON} \
        --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    else
        #${fastp} \
        fastp \
            --in1 ${READ1} \
            --in2 ${UMI_READ} \
            --out1 ${R1_UMI} \
            --out2 ${UMI_OUT} \
            --umi \
            --umi_loc=read2 \
            --umi_len=${INDEX_LENGTH} \
            -Q \
            -A \
            -L \
            -w 1 \
            -u 100 \
            -n ${INDEX_LENGTH} \
            -Y 100 \
            -G

        #${fastp} \
        fastp \
            --in1 ${READ2} \
            --in2 ${UMI_READ} \
            --out1 ${R2_UMI} \
            --out2 ${UMI_OUT} \
            --umi \
            --umi_loc=read2 \
            --umi_len=${INDEX_LENGTH} \
            -Q \
            -A \
            -L \
            -w 1 \
            -u 100 \
            -n ${INDEX_LENGTH} \
            -Y 100 \
            -G

        # STEP 2: now do adapter trimming and overlap correction (optional) ... use as many threads as you can/want
        echo  "=============================="
        echo "(STEP 2) Trim adapters and do overlap correction with fastp for ${SAMPLE_PREFIX}"
        #${fastp} \
        fastp \
            --correction \
            -w ${THREAD_NUM} \
            --in1 ${R1_UMI} \
            --in2 ${R2_UMI} \
            --out1 ${R1_TRIM} \
            --out2 ${R2_TRIM} \
            --html ${HTML} \
            --json ${JSON} \
            --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
        #transfer R1 and R2 files
        #rsync -avP ${R1_TRIM} ${OUTPUT_DIR}   
    fi
    #rsync -avP ${R2_TRIM} ${OUTPUT_DIR}
    conda deactivate 
else
    echo  "=============================="
    echo "UMI placement and adapter trimming already completed"
fi

if [ ! -f "${GENCORE_BAM}" ]; then
    # STEP 3: align the output from STEP 3 to reference genome using BWA MEM producing a "precorrected.bam" .... duplicate marking is not needed since we'll use gencore
    echo  "=============================="
    echo "(STEP 3) Align reads with bwa mem for sample ${SAMPLE_PREFIX}"
    bwa mem \
        -t 10 \
        -v 2 \
        -R "${READGROUP}" \
        ${REF_FASTA} \
        ${R1_TRIM} \
        ${R2_TRIM} | samtools sort -@ 2 -o ${BWA_BAM} - 
    samtools index -b ${BWA_BAM} 
    #exit 2
    # STEP 4: use gencore to error correct ... duplicates are "collapsed" into consensus reads using start positions of fwd/rev reads and UMI (if UMI is present).  Gencore will also work without UMI (making it a generalizable step)
    echo  "=============================="
    echo "(STEP 4) Error correct and collaspe duplicates with gencore for sample ${SAMPLE_PREFIX}"
    conda activate gencore

    #${gencore} \
    gencore \
        --ref=${REF_FASTA} \
        --in=${BWA_BAM} \
        --out=${GENCORE_BAM} \
        --umi_prefix="" \
        -s $MIN_READS_PER_MOLECULE
    conda deactivate
    samtools index -b ${GENCORE_BAM}

    # transfer gencore outputs	
    #rsync -avP ${HTML} ${OUTPUT_DIR}
    #rsync -avP ${GENCORE_BAM}* ${OUTPUT_DIR}
else
    echo  "=============================="
    echo "Alignment already completed"
fi

# FASTQC and Picard HSMetrics on the BAM file
if [ ! -f "${HSMETRICS}" ]; then
    
    #FASTQC
    fastqc -o ${OUTPUT_DIR} ${READ1}  
    fastqc -o ${OUTPUT_DIR} ${READ2}  

    # Picard Files for GrCH37
    BED_FILE=${TARGETS}
    REF_DICT="/home/es984/oelab_es984/bed_files/37_genome_files/human_g1k_b37.fasta.dict"

    ## Picard Files for Grch38
    #BED_FILE="/home/es984/oelab_es984/PreCISE/CoveredRegions_Cornell-ClonalHematop_NGSTECustom_0001537d_hg38.bed"
    #INT_LIST="/home/es984/oelab_es984/PreCISE/CoveredRegions_Cornell-ClonalHematop_NGSTECustom_0001537d_hg38.interval_list"
    #REF_DICT="/home/es984/oelab_es984/38_genome_files/genome.dict"

    #picard BedToIntervalList \
        #-I ${BED_FILE} \
        #-O ${INT_LIST} \
        #-SD ${REF_DICT}

    picard CollectHsMetrics \
        -I ${GENCORE_BAM} \
        -O ${HSMETRICS} \
        -R ${REF_FASTA} \
        --BAIT_INTERVALS ${INT_LIST} \
        --TARGET_INTERVALS ${INT_LIST} \
        --COVERAGE_CAP 2000

else
    echo "==========================="
    echo "FASTQC and Picard HSMetrics already run for sample ${SAMPLE_PREFIX}"
fi

if [ ! -f "${INSERTSIZEMETRICS}" ]; then
    picard CollectInsertSizeMetrics \
        --INPUT ${GENCORE_BAM} \
        --OUTPUT ${INSERTSIZEMETRICS} \
        --Histogram_FILE ${SAMPLE_PREFIX}_insert_size_histogram.pdf

else
    echo "==========================="
    echo "Picard insertsizemetrics already run for sample ${SAMPLE_PREFIX}"
fi

#continue

# Get basepair depth information for CH Genes
if [ ! -f "${SAMTOOLS_DEPTH}" ]; then
    echo "Computing base pair depth for CH Genes..."
    # samtools depth \
    #     -b "${TARGETS}" \
    #     "${GENCORE_BAM}" | \
    #     awk -v x="${SAMPLE_PREFIX}" '{print $0, x}'  > "${SAMTOOLS_DEPTH}"
    # conda activate bedtools
    # bedtools coverage \
    #     -a "/home/es984/oelab_es984/bed_files/igv-regions-panheme-all-gene-coords-annot-hg19.bed" \
    #     -b ${GENCORE_BAM} | \
    #     awk -v x="${SAMPLE_PREFIX}" '{print $0, x}' > ${SAMTOOLS_DEPTH}
    # conda deactivate
    conda activate mosdepth
    mosdepth \
        -b "/home/es984/oelab_es984/bed_files/igv-regions-panheme-all-gene-coords-annot-grch37.bed" \
        --no-per-base \
        ${SAMPLE_PREFIX} \
        ${GENCORE_BAM}
    gunzip \
        -c ${SAMPLE_PREFIX}.regions.bed.gz | \
        awk -v x="\t ${SAMPLE_PREFIX}" '{print $0, x}' > ${SAMTOOLS_DEPTH}
    conda deactivate
    echo "...base pair depth generated for CH Genes"
else
    echo  "=============================="
    echo "Base pair depth already generated for CH Genes"
fi

# Get basepair depth information for Ig Genes
if [ ! -f "${SAMTOOLS_DEPTH_IG}" ]; then
    echo "Computing base pair depth for Ig Genes..."
    conda activate mosdepth
    mosdepth \
        -b "/home/es984/oelab_es984/bed_files/merged_probe_file_Cornell_Nuria_IG_Reduced_TE-93864814_grch37_210504170209_num.bed" \
        --no-per-base \
        ${SAMPLE_PREFIX} \
        ${GENCORE_BAM}
    gunzip \
        -c ${SAMPLE_PREFIX}.regions.bed.gz | \
        awk -v x="\t ${SAMPLE_PREFIX}" '{print $0, x}' > ${SAMTOOLS_DEPTH_IG}
    conda deactivate
    echo "...base pair depth generated for Ig Genes"
else
    echo  "=============================="
    echo "Base pair depth already generated for Ig Genes"
fi

# rename VCF output file
if [[ ${FADE_BOOLEAN} == "enzymatic" ]]; then
    VEP_TCGA=${OUTPUT_DIR}/${SAMPLE_PREFIX}_fade.varcall.marked.sort.vep.annot.cosmic.tcga.vcf.gz
fi

if [ ! -f "${VEP_TCGA}" ]; then
    
    # TODO: add FADE to filtet out artifacts from enzymatic calls
    # if statement for FADE input ($6 == "enzymatic")
    # first run FADE and generate new BAM
    # change output file name to include FADE
    if [[ ${FADE_BOOLEAN} == "enzymatic" ]]; then
        echo "(OPTIONAL STEP) Remove enzymatic artifacts using FADE "
        spack load htslib@1.9%gcc@6.3.0
        export LD_LIBRARY_PATH=/home/es984/oelab_es984/tools/parasail-2.4.2-manylinux1_x86_64/lib:$LD_LIBRARY_PATH
        FADE=/home/es984/oelab_es984/tools/fade
        THREAD=10

        # Run FADE
        $FADE annotate -t $THREAD $GENCORE_BAM $REF_FASTA > $FADE_BAM 2> /dev/null
        samtools sort -@ 10 -n $FADE_BAM > $FADE_BAM_SORT 
        $FADE out -t $THREAD $FADE_BAM_SORT > $FADE_BAM_OUT 
        samtools sort -@ 10 $FADE_BAM_OUT > $FADE_BAM_FINAL 
        samtools index -b ${FADE_BAM_FINAL}
        spack unload htslib@1.9%gcc@6.3.0	

        # change variable pointer for bam for variant calling
        GENCORE_BAM=${FADE_BAM_FINAL}
        VEP_TCGA=${OUTPUT_DIR}/${SAMPLE_PREFIX}_fade.varcall.marked.sort.vep.annot.cosmic.tcga.vcf.gz
    fi
    
    
    #STEP 5: variant calling with VarDict
    echo  "=============================="
    echo "(STEP 5) Call Variants with VarDict java for sample ${SAMPLE_PREFIX}"
    spack load -r jdk@8u172-b11%gcc@6.3.0
    #spack load -r vardictjava@1.5.1^gcc@6.3.0
    spack load vardictjava
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

    echo "(STEP 5B): Filter strand bias using VarDict supplementary tools"
    cat ${VARDICT_AUX} |
        ${vardict_script_dir}/teststrandbias.R | \
        ${vardict_script_dir}/var2vcf_valid.pl \
        -N ${SAMPLE_PREFIX} \
        -A \
        -f $AF_MIN | \
        bgzip -c > ${VARDICT_TMP1_VCF}
    tabix -f ${VARDICT_TMP1_VCF}

    spack unload jdk@8u172-b11%gcc@6.3.0
    spack unload vardictjava

    # NOTE: removed VCFHEADER_CONTIG variable
    #include vcfheader contigs in the configuration file
    echo  "=============================="
    echo "(STEP 6): Annotate variants with bcftools"
    bcftools annotate \
        -o ${VARDICT_TMP2_VCF} \
        -O z ${VARDICT_TMP1_VCF}
    tabix -f ${VARDICT_TMP2_VCF}

    echo  "=============================="
    echo "(STEP 7): Filter Variants with vcflib and bcftools"

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

    conda activate vawk  
    vawk \
        --header '{split(S$*$DP,DP,"\t");split(S$*$AF,AF,"\t");split(S$*$MQ,MQ,"\t");split(S$*$NM,NM,"\t");QUAL=$6;((AF[1] * DP[1] < 6) && ((MQ[1] < 55.0 && NM[1] > 1.0) || (MQ[1] < 60.0 && NM[1] > 2.0) || (DP[1] < 10) || (QUAL < 45))) ? $7="bcbio" : $7 = $7; print $0}' ${VARDICT_TMP3_VCF} | \
        sed 's/=-0/=0/g' | \
        vcf-sort 2>/dev/null | \
        bgzip -c > ${VARDICT_FINAL_VCF}
    conda deactivate
    tabix -f ${VARDICT_FINAL_VCF}

    #transfer filtered vardict vcf
    #rsync -avP ${VARDICT_FINAL_VCF} ${OUTPUT_DIR}


    echo  "=============================="
    echo "(STEP 8): Annotate Variants with VEP"

    # remove files so we can re-write over files
    rm -rf ${VEP_OUT} ${VEP_OUT_GZ} ${VEP_COSMIC} ${VEP_TCGA} ${FINAL_OUTPUT_VCF} 
    rm -rf *mosdepth*
    rm -rf *UMI*

    #${vep_dir}/vep \
    conda activate vep
    vep \
        -i ${VARDICT_FINAL_VCF} \
        -o ${VEP_OUT} \
        --species homo_sapiens \
        --everything \
        --symbol \
        --per_gene \
        --cache \
        --sift b \
        --dir "/home/es984/oelab_es984/37_genome_files" \
        --vcf \
        --offline \
        --vcf_info_field ANN \
        --force_overwrite \
        --fasta ${REF_FASTA} \
        --hgvs \
        --pick_allele \
        --filter "AF > $AF_MIN" 
    #--filter "SYMBOL in ACTA2, ACTC1, AD000671.6, ANGPTL4, ANKRD26, APC, APOA5, APOB, APOC3, ARID1A, ASXL1, ATM, BAP1, BARD1, BCL2, BCOR, BCORL1, BMPR1A, BRAF, BRCA1, BRCA2, BRIP1, CARD11, CBL, CBX3P3, CDH1, CDK4, CDKN2A, CEBPA, CHD2, CHEK2, COL3A1, DDX41, DNMT3A, DSC2, DSG2, DSP, EPCAM, ETV6, FBN1, FLT3, GATA1, GATA2, GLA, GNAS, GNB1, GREM1, HIST1H1E, IDH1, IDH2, JAK2, KCNH2, KCNQ1, KDM1A, KMT2D, KRAS, LDLR, LMNA, LPA, LPL, MBD4, MITF, MLH1, MPL, MSH2, MSH6, MUTYH, MYBPC3, MYD88, MYH11, MYH7, MYL2, MYL3, NBN, NOTCH1, NPC1L1, NPM1, NRAS, PALB2, PAX5, PCSK9, PIGA, PIK3CA, PIM1, PKP2, PMS2, POLD1, POLE, POT1, POU2AF1, PPM1D, PRKAG2, PTEN, RAD51C, RAD51D, RP11-244F12.3, RP11-34P13.7, RP11-758N13.1, RPL21P4, RPL36A-HNRNPH2, RUNX1, RYR2, SCN5A, SF3B1, SMAD3, SMAD4, SPEN, SRCAP, SRP72, SRSF2, STK11, TET2, TETT, TGFBR1, TGFBR2, TMEM43, TNFRSF14, TNNI3, TNNT2, TOE1, TP53, TPM1, U2AF1, U2AF1L4, XPC, XPO1, ZRSR2" 
    conda deactivate
    bgzip -f ${VEP_OUT}
    tabix -f ${VEP_OUT_GZ}

    # Annotate with SnpSift for Cosmic and tcga
    #java -Xmx4g -jar ${snpeff_dir}SnpSift.jar annotate \
    conda activate snpeff
    SnpSift annotate \
        -v \
        -tabix \
        $cosmicfile \
        ${VEP_OUT_GZ} > ${VEP_COSMIC} 2> /dev/null

    # TODO
    #${snpeff_dir}/scripts/vcfEffOnePerLine.pl | \
    #SnpEff vcfEffOnePerLine | \
    #java -Xmx4g -jar ${snpeff_dir}SnpSift.jar annotate \
        #bcftools view - -T ^${lcrbed} -O v | \
    SnpSift annotate \
        -info TCGA_SAMPLE \
        -noId $tcgafile \
        ${VEP_COSMIC} | \
        ${snpeff_dir}/scripts/vcfEffOnePerLine.pl | \
        bgzip -c > ${VEP_TCGA}  2> /dev/null
    conda deactivate
    tabix -f ${VEP_TCGA}

    #echo  "=============================="
    #echo "Transfer final file from TMPDIR to ${OUTPUT_DIR}"
    #rsync -avP ${VEP_TCGA} ${OUTPUT_DIR}

    echo "Final file can be found at ${VEP_TCGA}"
else
    echo  "=============================="
    echo "Variant calling already completed"
fi


if [ ! -f "${GERMLINE_VEP}" ]; then
    echo "(STEP 8): Germline Variant Calling"

    KNOWN_SITES_VCF="/home/es984/oelab_es984/PreCISE/All_20180423.nochr.vcf.gz"
    conda activate gatk

    # Mark dupplicates - 4 mins
    gatk MarkDuplicates \
        -I ${GENCORE_BAM} \
        -O ${SAMPLE_PREFIX}_marked_duplicates.bam \
        -M ${SAMPLE_PREFIX}_marked_dup_metrics.txt
    
    # Base recalibration - 1.25 hours
    gatk BaseRecalibrator \
        --reference ${REF_FASTA} \
        --input ${SAMPLE_PREFIX}_marked_duplicates.bam \
        --output ${SAMPLE_PREFIX}_recalibration_report.table \
        --known-sites ${KNOWN_SITES_VCF} 

    # 5 mins
    gatk ApplyBQSR \
        -R ${REF_FASTA} \
        -I ${SAMPLE_PREFIX}_marked_duplicates.bam \
        --bqsr-recal-file ${SAMPLE_PREFIX}_recalibration_report.table \
        -O ${SAMPLE_PREFIX}_recalibrated.bam

    # 5 mins
    gatk HaplotypeCaller -R ${REF_FASTA} \
        -I ${SAMPLE_PREFIX}_recalibrated.bam \
        -O ${SAMPLE_PREFIX}_raw_variants.vcf \
        -L ${INT_LIST} 

    # <1 min
    gatk VariantFiltration -R ${REF_FASTA} \
        -V ${SAMPLE_PREFIX}_raw_variants.vcf \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O ${SAMPLE_PREFIX}_filtered_variants.vcf
    conda deactivate

    # 5 mins
    conda activate vep
    vep \
        --fasta ${REF_FASTA} \
        --input_file ${SAMPLE_PREFIX}_filtered_variants.vcf \
        --output_file ${GERMLINE_VEP} \
        --species homo_sapiens \
        --everything \
        --symbol \
        --per_gene \
        --cache \
        --sift b \
        --dir "/home/es984/oelab_es984/37_genome_files" \
        --vcf \
        --offline \
        --vcf_info_field ANN \
        --force_overwrite \
        --fasta ${REF_FASTA} \
        --hgvs \
        --pick_allele 
    conda deactivate
else   
    echo  "=============================="
    echo "Germline variant calling already completed"
fi

#move files from one output folder to a folder for each sample
#cd ${OUTPUT_DIR}
#mkdir ${SAMPLE_PREFIX}
#mv ${SAMPLE_PREFIX}* ${OUTPUT_DIR}/${SAMPLE_PREFIX}

echo " Script $0 is complete for sample ${SAMPLE_PREFIX}"
DATE=$(date '+%d/%m/%Y %H:%M:%S');
echo "$DATE"
