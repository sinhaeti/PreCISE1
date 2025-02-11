# PreCISE

PreCISE1, or Predictive Clonal hematopoeisis, cardiovascular, and cancer Integrated Sequencing Evaluation, is an end-to-end assay and pipeline that takes in UMI-barcoded fastq reads and calls variants for any targeted sequencing panel. PreCISE1 can detect somatic variants at low VAF (<1%) with high sensitivity and positive predictive value (PPV) of 90% sensitivity and 100% PPV, while retaining 96% sensitivity and 100% positive predictive value for variants at higher VAF (>10%). Furthermore, the test also detects germline variants related to many associated diseases to clonal hematopoiesis, such as cardiovascular disease and solid tumor cancers. PreCISE1 has already been used to monitor patients to study basic biology, such as initiation of CH in the bone marrow, and answer translational questions, such as monitoring clonal hematopoiesis clones in cancer cohorts with treatments.


## Pipeline Information
### How to Run the Pipeline
`./run_precise_jobarray.sh [config_file] [fastq_directory] [output_directory] [min_allele_fraction] [enzymatic_or_sonication]`

config_file: path to configuration file with parameters \
fastq_directory: path to raw .fastq or .fastq.gz files \
output_directory: path for BWA and vardict output \
min_allele_fraction: minimum VAF, (regular=0.005, lowVAF=0.001, noVAF=0.0001) \
enzymatic_or_sonication: enter 'enzymatic' if sample used enzymatic fragmentation, 'sonication' if sample used sonication 

### Summary
1. Add UMIs to read 1 (fasptp)
2. Trim adapters and do overlap correction (fastp)
3. Align reads (bwa)
4. Generate WC reports (fastqc, picard hsmetrics, mosdepth)
5. Error correct, collaspe duplicate reads using UMI and start positions of fwd/rev strands (gencore)
5b. (optional) for enzymatic fragmentation, remove articaftual reads in bam file (fade)
6. Call variants (vardict), Filter strand bias (vardict supplementary tools)
7. Annotate (bcftools)
8. Filter variants (vcflib, bcftools)
9. Annotate (VEP), with cosmic and tcga (snpeff)
10. Save as a tsv file (awk, vawk)

Pipeline is run in parallel using SLURM on a TMPDIR. 

### Pipeline Output
1. fastq files with UMI attched, and adapter trimmed: ``${SAMPLE_PREFIX}_R1_trim.fastq.gz` `${SAMPLE_PREFIX}_R2_trim.fastq.gz``
2. HTML report of fastp output: `${SAMPLE_PREFIX}_fastp.json`
3. HTML report of hsmetrics output: `${SAMPLE_PREFIX}_hs_metrics.txt`
4. Error corrected and collasped bam: `${SAMPLE_PREFIX}_corrected.bam`
5. Mosdepth depth calculations: `${SAMPLE_PREFIX}.regions.bed` `${SAMPLE_PREFIX}.ig.regions.bed`
5. Variant Calling File: `${SAMPLE_PREFIX}_varcall.marked.sort.vep.annot.cosmic.tcga.vcf.gz`

### Next Steps
#### Aggregate QC data using multiqc
  `cd [output_directory]` \
  `conda activate multiqc` \
  `multiqc *` 

#### Aggregate variant calls using RScript
`cd ${OUTPUT_DIR}` \
`Rscript ./aggregate_variants_mutect.R [output_directory] [rawcalls_outputfile] [filteredcalls_output]" *varcall.marked.sort.vep.annot.cosmic.tcga.vcf.gz$"`

### Tools Used
__From athena server:__ \
`spack load -r bwa@0.7.15%gcc@6.3.0` \
`spack load -r samtools@1.8^gcc@6.3.0` \
`spack load -r jdk@8u172-b11^gcc@6.3.0` \
`spack load -r vardictjava@1.5.1^gcc@6.3.0` \
`spack load -r vcftools@0.1.14^gcc@6.3.0` \
`spack load -r python@3.6.0%gcc@6.3.0` \
`spack load -r bcftools@1.6` \
`spack load -r r@3.5.0` \
`spack load -r /j3gbsxv #glib@2.53.1` \
`spack load perl-try-tiny` \
`spack load -r perl-dbd-mysql` \
`spack load -r /g7nwnaf #gcc` \
`spack load -r /orql4pe #perl` 

__From conda:__ \
fade: /home/es984/miniconda3/envs/fade \
fastp: /home/es984/miniconda3/envs/fastp \
gatk: /home/es984/miniconda3/envs/gatk \
gencore: /home/es984/miniconda3/envs/gencore \
mosdepth: /home/es984/miniconda3/envs/mosdepth \
multiqc: /home/es984/miniconda3/envs/multiqc \
vawk: /home/es984/miniconda3/envs/vawk 

__More installed packages:__ \
vcflib_dir \
snpEff \
ensembl-vep 
