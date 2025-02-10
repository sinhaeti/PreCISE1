#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process VCF {
    // Define input and output directories
    publishDir "/Users/zichan/Documents/CH/wgs-nygc/output/${sampleName}", mode: 'copy', overwrite: true
    input:
    path vcf_file
    path reference_csv

    output:
    path "${sampleName}_selected_rows.txt"

    script:
    """
    library(vcfR)

    # Load the reference data for CH variants
    CH_variants <- read.csv("$reference_csv")

    # Extract sample name from the input VCF file
    sampleName <- fileBaseName(vcf_file)

    # Read VCF data
    vcf_data <- read.vcfR("$vcf_file")
    df <- vcfR2tidy(vcf_data)
    df_fix <- df\$fix
    df_CSQ <- data.frame(do.call('rbind', strsplit(as.character(df_fix\$CSQ), '|', fixed = TRUE)))
    df_CSQ\$positions <- sub(".*:", "", df_CSQ\$X13)
    df_CSQ\$variants <- paste(df_CSQ\$X15, df_CSQ\$positions, sep = " ")

    # Compare with the reference data
    reference_variants <- CH_variants\$variant
    selected_rows <- df_CSQ[df_CSQ\$variants %in% reference_variants, ]

    # Write selected rows to output file
    write.table(selected_rows, file="${sampleName}_selected_rows.txt", sep="\t", quote=FALSE, row.names=FALSE)
    """
}

workflow {
    // Specify input parameters
    // params.vcfDirectory = "/Users/zichan/Documents/CH/wgs-nygc/test"
    params.referenceCSV = "/Users/zichan/Documents/CH/CH variants.csv"


    vcf_ch = channel.fromPath("/Users/zichan/Documents/CH/wgs-nygc/test/*.vcf.gz") 
    VCF (
        vcf_ch,params.referenceCSV
    )



}
