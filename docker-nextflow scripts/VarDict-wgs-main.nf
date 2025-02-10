#!/usr/bin/env nextflow
/*
========================================================================================
                         PreCISE
========================================================================================
 eipm/precise
 https://github.com/eipm/precise/
----------------------------------------------------------------------------------------
*/
// WORK IN PROGRESS

Channel.fromPath(params.sampleInfoFile)
        .map{ row -> tuple(row.SampleID, file(row.BAM) )}
       .set { bam_files }
       bam_file_ch = bam_files.groupTuple(by: 0)

if ( !params.sampleInfoFile && !new File(params.sampleInfoFile).exists() ) {
    exit 1, "sampleInfoFile not defined or an invalid path was provided."
}
if ( !params.resultsDir && !new File(params.resultsDir).exists() ) {
    exit 1, "resultsDir not defined or an invalid path was provided."
}

//ALIGHMENT
process SAMTOOLSALIGNMENT {

    tag "$SAMPLE_PREFIX" 
    //replace with varidct container
    container 'dockerfile'

    input:
    tuple val(SAMPLE_PREFIX), path(TARGETS), path(BAM)
    output:
    tuple val(SAMPLE_PREFIX), path(PANHEME_BAM), emit: bwa

    script:
            """
            echo  "=============================="
            echo "(STEP 1) Align reads with samtools for sample ${SAMPLE_PREFIX}"
	    samtools view -b -L $TARGETS ${BAM} > ${PANHEME_BAM}
	    samtools index ${PANHEME_BAM}

            """

    }


//
process VARDICTVARIANTSCALLING {

    tag "$SAMPLE_PREFIX" 
    //replace with varidct container
    container 'dockerfile'

    input:
    tuple val(SAMPLE_PREFIX), path(REF_FASTA), path(AF_MIN), path(PANHEME_BAM), path(EXTEND_BED),path(MAX_INSERT), path(MIN_MAPPING_QUALITY), path(TARGETS)
    output:
    tuple val(SAMPLE_PREFIX), path(VARDICT_TMP1_VCF), emit: vardict

    script:
            """
            echo  "=============================="
            echo "(STEP 2) Call Variants with VarDict java for sample ${SAMPLE_PREFIX}"

            VarDict \
                        -th 10 \
                        -G $REF_FASTA \
                        -f $AF_MIN \
                        -N ${SAMPLE_PREFIX} \
                        -b ${PANHEME_BAM} \
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
    container 'dockerfile'

    input:
    tuple val(SAMPLE_PREFIX), path(VARDICT_TMP1_VCF), path(vcflib_dir), path(SB_PVAL), path(SB_ODDRATIO), path(MIN_DEPTH),path(MIN_VAR_READS)
    output:
    tuple val(SAMPLE_PREFIX), path(VARDICT_TMP3_VCF), emit: bcftools

    script:
           """
            echo  "=============================="
            echo "(STEP 3): Annotate variants with bcftools"
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
	        """

process VAWKFORMATTING{

    tag "$SAMPLE_PREFIX" 
    //replace with varidct container
    container 'dockerfile'

    input:
    tuple val(SAMPLE_PREFIX), path(VARDICT_TMP3_VCF),path(DP),path(AF)
    output:
    tuple val(SAMPLE_PREFIX), path(VARDICT_FINAL_VCF), emit: vawk

    script:
           """
            echo  "=============================="
            echo "(STEP 3): Annotate variants with bcftools"
            vcftools --vcf ${VARDICT_TMP3_VCF} --recode --stdout | \
	    awk 'BEGIN {FS=OFS="\t"} {
   	       if ($0 ~ /^#/) { print; next }
   	       split($8, info, ";")
  	       for (i in info) {
   	          split(info[i], kv, "=")
   	          if (kv[1] == "DP") dp = kv[2]
   	          if (kv[1] == "AF") af = kv[2]
  	          if (kv[1] == "MQ") mq = kv[2]
  	          if (kv[1] == "NM") nm = kv[2]
		}
	        qual = $6
	        if ((af * dp < 6) && ((mq < 55.0 && nm > 1.0) || (mq < 60.0 && nm > 2.0) || (dp < 10) || (qual < 45))) {
	              $7 = "bcbio"
	        }
	        print
	    }' | \
	    vcf-sort | \
	    bgzip -c > ${VARDICT_FINAL_VCF}

            tabix -f ${VARDICT_FINAL_VCF}

            """


    }
process VEPANNOTATION {

    tag "$SAMPLE_PREFIX" 
    //replace with varidct container
    container 'ensembl-vep'

    input:
    tuple val(SAMPLE_PREFIX), path(VARDICT_FINAL_VCF), path(REF_FASTA), path(snpeff_dir)
    output:
    tuple val(SAMPLE_PREFIX), path(VEP_OUT_GZ), emit: vardict

    script:
            """
            echo  "=============================="
            echo "(STEP 8): Annotate Variants with VEP"

            // remove files so we can re-write over files
            //rm -rf ${VEP_OUT} ${VEP_OUT_GZ} ${VEP_COSMIC} ${VEP_TCGA} ${FINAL_OUTPUT_VCF} 
            //rm -rf *mosdepth*
            //rm -rf *UMI*

            vep \
                -i ${VARDICT_FINAL_VCF} \
                -o ${VEP_OUT} \
                --species homo_sapiens \
                --everything \
                --symbol \
                --per_gene \
                --cache \
                --sift b \
                --dir "/eipm-research/lims/home/zil4004/37_genome_files" \
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

process SNPSFTANNOTATION {

    tag "$SAMPLE_PREFIX" 
    container 'dceoy/snpeff'

    input:
    tuple val(SAMPLE_PREFIX), path(VEP_OUT_GZ), path(REF_FASTA), path(cosmicfile), path(tcgafile), path(snpeff_dir)
    output:
    tuple val(SAMPLE_PREFIX), path(VEP_TCGA), emit: vardict

    script:
            """
            echo  "=============================="
            echo "(STEP 8): Annotate Variants with SnpSft”

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
    SAMTOOLSALIGNMENT(SAMPLE_PREFIX,TARGETS,BAM)
    VARDICTVARIANTSCALLING(SAMPLE_PREFIX,REF_FASTA,AF_MIN,GENCORE_BAM,EXTEND_BED,MAX_INSERT,MIN_MAPPING_QUALITY,TARGETS)
    BCFTOOLSANNOTATION(SAMPLE_PREFIX,VARDICT_TMP1_VCF,vcflib_dir,SB_PVAL,SB_ODDRATIO,MIN_DEPTH,MIN_VAR_READS)
    VEPANNOTATION(SAMPLE_PREFIX,VARDICT_FINAL_VCF,REF_FASTA,snpeff_dir)
    SNPSFTANNOTATION(SAMPLE_PREFIX,VEP_OUT_GZ,REF_FASTA,cosmicfile,tcgafile,snpeff_dir)
}
