#!/usr/bin/env nextflow
/*
========================================================================================
                         PreCISE
========================================================================================
 eipm/precise
 https://github.com/eipm/precise/
----------------------------------------------------------------------------------------
*/
//define the input samples, will modify later based on the location and the file name
Channel.fromPath("/home/zil4004/*.fastq")
       .map { file -> tuple(file.baseName.replace(".fastq", ""), file) }
       .set { fastq_files }

// Fastp to process fastq input files
process FASTQPREPROCESS {

    tag "$SAMPLE_PREFIX" 
    container 'dockerfile-panheme'

    input:
    tuple val(SAMPLE_PREFIX), path(READ1), path(UMI_READ), path(INDEX_LENGTH), path(THREAD_NUM)

    output:
    tuple val(SAMPLE_PREFIX), path(R1_TRIM), path(R2_TRIM), emit: fastq

    script:
            """
            echo "(STEP 1) ADD UMI to R1 and R2 sequence identifiers for sample ${SAMPLE_PREFIX}"
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

            echo  "=============================="
            echo "(STEP 2) Trim adapters and do overlap correction with fastp for ${SAMPLE_PREFIX}"
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
            else
                echo  "=============================="
                echo "UMI placement and adapter trimming already completed"
            fi

            """
}

//ALIGHMENT
process BWAALIGNMENT {

    tag "$SAMPLE_PREFIX" 
    //replace with varidct container
    container 'dockerfile-panheme'

    input:
    tuple val(SAMPLE_PREFIX), path(R1_TRIM), path(R2_TRIM), path(REF_FASTA)
    output:
    tuple val(SAMPLE_PREFIX), path(BWA_BAM), emit: bwa

    script:
            """
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
            """


    }


//BCFtools annotation
process GENECOREERRORCORRECTION {

    tag "$SampleID" 
    container 'dockerfile-panheme'

    input:
    tuple val(SAMPLE_PREFIX), path(BWA_BAM), path(REF_FASTA), path(TARGETS), path(INT_LIST)

    output:
    tuple val(SAMPLE_PREFIX), path(GENCORE_BAM), path(INSERTSIZEMETRICS) emit: gencore

    script:
            """
            echo  "=============================="
            echo "(STEP 4) Error correct and collaspe duplicates with gencore for sample ${SAMPLE_PREFIX}"
            gencore \
                --ref=${REF_FASTA} \
                --in=${BWA_BAM} \
                --out=${GENCORE_BAM} \
                --umi_prefix="" \
                -s $MIN_READS_PER_MOLECULE
            samtools index -b ${GENCORE_BAM}


            
            // FASTQC
            fastqc -o ${OUTPUT_DIR} ${READ1}  
            fastqc -o ${OUTPUT_DIR} ${READ2}  

            // Picard Files for GrCH37
            BED_FILE=${TARGETS}
            REF_DICT="/home/es984/oelab_es984/bed_files/37_genome_files/human_g1k_b37.fasta.dict"

            picard CollectHsMetrics \
                    -I ${GENCORE_BAM} \
                    -O ${HSMETRICS} \
                    -R ${REF_FASTA} \
                    --BAIT_INTERVALS ${INT_LIST} \
                    --TARGET_INTERVALS ${INT_LIST} \
                    --COVERAGE_CAP 2000



            if [ ! -f "${INSERTSIZEMETRICS}" ]; then
                picard CollectInsertSizeMetrics \
                    --INPUT ${GENCORE_BAM} \
                    --OUTPUT ${INSERTSIZEMETRICS} \
                    --Histogram_FILE ${SAMPLE_PREFIX}_insert_size_histogram.pdf
    
            else
                echo "==========================="
                echo "Picard insertsizemetrics already run for sample ${SAMPLE_PREFIX}"
            fi
            """

}

//ALIGHMENT
process MOSEDPETHINFO {

    tag "$SAMPLE_PREFIX" 
    //replace with varidct container
    container 'dockerfile-panheme'

    input:
    tuple val(SAMPLE_PREFIX), path(GENCORE_BAM), path(FADE_BOOLEAN)
    output:
    tuple val(SAMPLE_PREFIX), path(FADE_BAM_FINAL), path(VEP_TCGA) emit: mosedepth

    script:
        """
        // Get basepair depth information for CH Genes
        if [ ! -f "${SAMTOOLS_DEPTH}" ]; then
            echo "Computing base pair depth for CH Genes..."
           

            mosdepth \
                -b "/home/es984/oelab_es984/bed_files/igv-regions-panheme-all-gene-coords-annot-grch37.bed" \
                --no-per-base \
                ${SAMPLE_PREFIX} \
                ${GENCORE_BAM}
            gunzip \
                -c ${SAMPLE_PREFIX}.regions.bed.gz | \
                awk -v x="\t ${SAMPLE_PREFIX}" '{print $0, x}' > ${SAMTOOLS_DEPTH}
            
            echo "...base pair depth generated for CH Genes"
        else
            echo  "=============================="
            echo "Base pair depth already generated for CH Genes"
        fi

        // Get basepair depth information for Ig Genes
        if [ ! -f "${SAMTOOLS_DEPTH_IG}" ]; then
            echo "Computing base pair depth for Ig Genes..."
            mosdepth \
                -b "/home/es984/oelab_es984/bed_files/merged_probe_file_Cornell_Nuria_IG_Reduced_TE-93864814_grch37_210504170209_num.bed" \
                --no-per-base \
                ${SAMPLE_PREFIX} \
                ${GENCORE_BAM}
            gunzip \
                -c ${SAMPLE_PREFIX}.regions.bed.gz | \
                awk -v x="\t ${SAMPLE_PREFIX}" '{print $0, x}' > ${SAMTOOLS_DEPTH_IG}
            echo "...base pair depth generated for Ig Genes"
        else
            echo  "=============================="
            echo "Base pair depth already generated for Ig Genes"
        fi

        // rename VCF output file
        if [[ ${FADE_BOOLEAN} == "enzymatic" ]]; then
            VEP_TCGA=${OUTPUT_DIR}/${SAMPLE_PREFIX}_fade.varcall.marked.sort.vep.annot.cosmic.tcga.vcf.gz
        fi

        if [ ! -f "${VEP_TCGA}" ]; then

            if [[ ${FADE_BOOLEAN} == "enzymatic" ]]; then

                # Run FADE
                $FADE annotate -t $THREAD $GENCORE_BAM $REF_FASTA > $FADE_BAM 2> /dev/null
                samtools sort -@ 10 -n $FADE_BAM > $FADE_BAM_SORT 
                $FADE out -t $THREAD $FADE_BAM_SORT > $FADE_BAM_OUT 
                samtools sort -@ 10 $FADE_BAM_OUT > $FADE_BAM_FINAL 
                samtools index -b ${FADE_BAM_FINAL}

                // change variable pointer for bam for variant calling
                GENCORE_BAM=${FADE_BAM_FINAL}
                VEP_TCGA=${OUTPUT_DIR}/${SAMPLE_PREFIX}_fade.varcall.marked.sort.vep.annot.cosmic.tcga.vcf.gz
        fi
    
        """


    }

//
process VARDICTVARIANTSCALLING {

    tag "$SAMPLE_PREFIX" 
    //replace with varidct container
    container 'dockerfile-panheme'

    input:
    tuple val(SAMPLE_PREFIX), path(REF_FASTA), path(AF_MIN), path(GENCORE_BAM), path(EXTEND_BED),path(MAX_INSERT), path(MIN_MAPPING_QUALITY), path(TARGETS)
    output:
    tuple val(SAMPLE_PREFIX), path(VARDICT_TMP1_VCF), emit: vardict

    script:
            """
            echo  "=============================="
            echo "(STEP 5) Call Variants with VarDict java for sample ${SAMPLE_PREFIX}"

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


            """



    }
process BCFTOOLSANNOTATION {

    tag "$SAMPLE_PREFIX" 
    //replace with varidct container
    container 'dockerfile-panheme'

    input:
    tuple val(SAMPLE_PREFIX), path(VARDICT_TMP1_VCF), path(vcflib_dir), path(SB_PVAL), path(SB_ODDRATIO), path(MIN_DEPTH),path(MIN_VAR_READS)
    output:
    tuple val(SAMPLE_PREFIX), path(VARDICT_FINAL_VCF), emit: bcftools

    script:
           """
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

            vawk \
                --header '{split(S$*$DP,DP,"\t");split(S$*$AF,AF,"\t");split(S$*$MQ,MQ,"\t");split(S$*$NM,NM,"\t");QUAL=$6;((AF[1] * DP[1] < 6) && ((MQ[1] < 55.0 && NM[1] > 1.0) || (MQ[1] < 60.0 && NM[1] > 2.0) || (DP[1] < 10) || (QUAL < 45))) ? $7="bcbio" : $7 = $7; print $0}' ${VARDICT_TMP3_VCF} | \
                sed 's/=-0/=0/g' | \
                vcf-sort 2>/dev/null | \
                bgzip -c > ${VARDICT_FINAL_VCF}
            tabix -f ${VARDICT_FINAL_VCF}

            """


    }
process VEPANNOTATION {

    tag "$SAMPLE_PREFIX" 
    //replace with varidct container
    container 'dockerfile-panheme'

    input:
    tuple val(SAMPLE_PREFIX), path(VARDICT_FINAL_VCF), path(REF_FASTA), path(snpeff_dir)
    output:
    tuple val(SAMPLE_PREFIX), path(VEP_TCGA), emit: vardict

    script:
            """
            echo  "=============================="
            echo "(STEP 8): Annotate Variants with VEP"

            // remove files so we can re-write over files
            rm -rf ${VEP_OUT} ${VEP_OUT_GZ} ${VEP_COSMIC} ${VEP_TCGA} ${FINAL_OUTPUT_VCF} 
            rm -rf *mosdepth*
            rm -rf *UMI*

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
            bgzip -f ${VEP_OUT}
            tabix -f ${VEP_OUT_GZ}

            // Annotate with SnpSift for Cosmic and tcga
            SnpSift annotate \
                -v \
                -tabix \
                $cosmicfile \
                ${VEP_OUT_GZ} > ${VEP_COSMIC} 2> /dev/null

            SnpSift annotate \
                -info TCGA_SAMPLE \
                -noId $tcgafile \
                ${VEP_COSMIC} | \
                ${snpeff_dir}/scripts/vcfEffOnePerLine.pl | \
                bgzip -c > ${VEP_TCGA}  2> /dev/null

            tabix -f ${VEP_TCGA}

            echo "Final file can be found at ${VEP_TCGA}"
            """



    }

workflow {
    timestamp = "${new Date().format('yyyyMMddHHmmss')}"
    data = channel.fromPath('/some/path/*.txt')

    // Run the process to call variants and annotate the results using annovar
    FASTQPREPROCESS(SAMPLE_PREFIX,READ1,UMI_READ,INDEX_LENGTH,THREAD_NUM)
    BWAALIGNMENT(SAMPLE_PREFIX,R1_TRIM,R2_TRIM,REF_FASTA)
    GENECOREERRORCORRECTION(SAMPLE_PREFIX,BWA_BAM,REF_FASTA,TARGETS,INT_LIST)
    MOSEDPETHINFO(SAMPLE_PREFIX,GENCORE_BAM,FADE_BOOLEAN)
    VARDICTVARIANTSCALLING(SAMPLE_PREFIX,REF_FASTA,AF_MIN,GENCORE_BAM,EXTEND_BED,MAX_INSERT,MIN_MAPPING_QUALITY,TARGETS)
    BCFTOOLSANNOTATION(SAMPLE_PREFIX,VARDICT_TMP1_VCF,vcflib_dir,SB_PVAL,SB_ODDRATIO,MIN_DEPTH,MIN_VAR_READS)
    VEPANNOTATION(SAMPLE_PREFIX,VARDICT_FINAL_VCF,REF_FASTA,snpeff_dir)
}
