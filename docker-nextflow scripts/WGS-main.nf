#!/usr/bin/env nextflow

/*
   DSL2 script to use bcftools to read the VCF file and VEP to annotate the results
   
*/

// Define input parameters
params {
    vcf_file = "/home/ipm-research/homes/zichan/pipelines/PreCISE-NYGCWGS/${SampleID_Tumor}.snv.indel.final.v6.annotated.vcf" // Input VCF file
    vep_cache = "/home/ipm-research/homes/zichan/pipelines/PreCISE-NYGCWGS/vep_cache" // Path to VEP cache directory
    output_dir = "/home/ipm-research/homes/zichan/pipelines/PreCISE-NYGCWGS/output_annotated" // Output directory for annotated VCF
}

// Define process to read VCF using bcftools
process READVCF {
    output:
    file "${output_dir}/${vcf_file.baseName}.bcftools.vcf", emit: bcftoolsProcessed


    script:
    """
    bcftools view ${vcf_file} -Ov -o ${output_dir}/${vcf_file.baseName}.bcftools.vcf
    """
}

// Define process to annotate the results using VEP
process ANNOTATEWITHVEP {
    input:
    file vcfFile 

    output:
    file "${output_dir}/${vcfFile.baseName}.vep_annotated.vcf", emit: annotatedVCF

    script:
    """
    vep --cache --dir ${params.vep_cache} --input_file ${vcfFile} --output_file ${output_dir}/${vcfFile.baseName}.vep_annotated.vcf --vcf --force_overwrite
    """
}

// Define the workflow
workflow {

    // Run the annotateWithVEP process to annotate the results using VEP
    READVCF | ANNOTATEWITHVEP
}
