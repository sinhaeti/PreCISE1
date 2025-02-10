#!/bin/bash

# check if all inputs are specified
if [ -z $1 ] || [ -z $2 ] || [ -z $3 ] || [ -z $4 ] || [ -z $5 ]; then
    echo "run_BWA_mutect [config_file] [fastq_directory] [output_directory] [min_allele_fraction] [enzymatic_or_sonication]"
    echo ""
    echo "config_file: path to configuration file with parameters"
    echo "fastq_directory: path to raw .fastq or .fastq.gz files"
    echo "output_directory: path for BWA and vardict output"
    echo "min_allele_fraction: minimum VAF, (regular=0.005, lowVAF=0.001, noVAF=0.0001)"
    echo "FADE_BOOLEAN: enter 'enzymatic' if sample used enzymatic fragmentation, 'sonication' if sample used sonication"
    echo "UMI_BOOLEAN: enter 'no_umi' if there is no UMI barcodes in the fastq files"
else
    # Grab 4 inputs
    CONFIG_FILE=$1
    INPUT_DIR=$2
    OUTPUT_DIR=$3
    AF_MIN=$4
    FADE_BOOLEAN=$5
    UMI_BOOLEAN=$6
    
    PARENT_DIRECTORY=$(dirname ${INPUT_DIR}) #get parent directory of $fastq_directory
    CODE_DIR="/athena/elementolab/scratch/es984/PreCISE/" #specify location of precise_variant_calling.sh #TODO: update
    RUN_NAME="$(basename ${INPUT_DIR})"
    FASTQ_PATHS="${PARENT_DIRECTORY}/${RUN_NAME}_samples.txt" #give a path to a file to store the paths to the fastq files in $fastq_directory

    ls "${INPUT_DIR}/" | grep "fastq" | sed -e 's/_R1.*$//g' | sed -e 's/_R2.*$//g' | sed -e 's/_R3.*$//g' | sort -u > ${FASTQ_PATHS} #generate list of full paths to fastq files and save to the file in $fastq_list

    ####### ALIGNMENT AND VARIANT CALLING  #########
    FASTQ_PATHS_LIST=$(cat $FASTQ_PATHS)
    for SAMPLE_PREFIX in ${FASTQ_PATHS_LIST}
    do
        echo "Sending job for ${SAMPLE_PREFIX} in sequencing run ${RUN_NAME}"
        sbatch \
        -o "${OUTPUT_DIR}/${SAMPLE_PREFIX}_${FADE_BOOLEAN}_variant_calling_log" \
        -e "${OUTPUT_DIR}/${SAMPLE_PREFIX}_${FADE_BOOLEAN}_variant_calling_log" \
        "${CODE_DIR}/precise_variant_calling.sh" \
        ${CONFIG_FILE} \
        ${INPUT_DIR} \
        ${OUTPUT_DIR} \
        ${AF_MIN} \
        ${SAMPLE_PREFIX} \
        ${FADE_BOOLEAN} \
        ${UMI_BOOLEAN}
    done

    ####### FILTER VARIANTS  #########
    # when last job is complete

    # check for number of output file in output folder (${OUTPUT_DIR}/${SAMPLE_PREFIX}.vep...) and if it matches the # of lines in FASTQ_PATHS
        # if number doesn't match, don't do next step
        # also the string to match depends on enzymatic or sonication ($FADE_BOOLEAN)

    # generate complete report with multiqc in the output folder

    # add Rscript for filtering variants from this sample
        # script combines variants to one file
        # filters out artifacts -- blacklist, 30% common, etc
    # Rscript aggregate_variants_mutect.R ${mutect_directory} ${output_file_name} ${mutect_pattern}
    # rsync to local mount of DropBox when aggregation is complete
fi


echo " Script $0 is complete for sequencing run ${RUN_NAME}"
	DATE=$(date '+%d/%m/%Y %H:%M:%S');
	echo "$DATE"