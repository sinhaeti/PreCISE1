library(readxl)
library(tidyverse)
library(magrittr)

SetNames <- function(col_names, col_values) {
    names(col_values) <- col_names
    col_values
}

command_args <- commandArgs(trailingOnly = TRUE)
mutect_directory <- command_args[1]
output_file_name <- command_args[2]
filter_output_file_name <- command_args[3]
mutect_pattern <- command_args[4]

# Example Inputs
# mutect_directory <- "/Users/etisinhawork/Dropbox (Guzman Lab)/Guzman Lab/Eti Sinha/Elemento-NMT-11919_comb_2022_03_08_output/VCFs"
# output_file_name <- "Elemento-NMT-11919_comb_2022_03_08.vardict.summary.parse.txt"
# filter_output_file_name <- "Elemento-NMT-11919_comb_2022_03_08.vardict.summary.parse.filtered.tx"
# mutect_pattern <- "*_varcall.marked.sort.vep.annot.cosmic.tcga.vcf.gz$"

# mutect_directory <- "/Users/etisinhawork/Dropbox (Guzman Lab)/Guzman Lab/Eti Sinha/Elemento-NMT-11000_2021_08_11/Vardict+Fade/"
# output_file_name <- "Elemento-NMT-11000_2021_08_11_fade.vardict.summary.parse.txt"
# filter_output_file_name <- "Elemento-NMT-11000_2021_08_11_fade.vardict.summary.parse.filtered.txt"
# mutect_pattern <- "*_fade.varcall.marked.sort.vep.annot.cosmic.tcga.vcf.gz$"

#panel_coordinates <- "/home/es984/oelab_es984/precise1.b37.sort.anno.bed"
#mutect_directory <- "/home/es984/oelab_es984/PreCISE/Elemento-NMT-9706_2020_11_20_re_output_lowVAF/precise_pipeline_output/"
#output_file_name <- "Elemento-NMT-9706_2020_11_20_re_output_lowVAF_mutect.vep.summary.txt"
#mutect_pattern <- "*_mutectcall.filter.vep.out.pass.vcf.gz$"
#output_file_name <- "Elemento-NMT-9707_2020_11_20_re_output_lowVAF_mutect.fade.vep.summary.txt"
#mutect_pattern <- "*_mutectcall.fade.filter.vep.out.pass.vcf.gz$"
#mutect_vcf <- "/home/es984/oelab_es984/PreCISE/Elemento-NMT-9706_2020_11_20_re_output_lowVAF/precise_pipeline_output/Sample_00d0002-Control-1-10-125ngEF-A-PreCISE1v3_S32_mutectcall.fade.filter.vep.out.vcf.gz"
# Sample_00d0001-NA12878-125ngEF-M-PreCISE1v3_S31_mutectcall.filter.vep.out.vcf.gz

print(c("Aggregating Variants for ", mutect_directory, "with samples ", mutect_pattern, " and the output is ", output_file_name))

print("(STEP 1) Grab Column Names from Header")
sample_names <- list.files(mutect_directory, pattern = mutect_pattern) %>% str_remove_all("_[^_]+$")

mutect_vcf_files <- list.files(mutect_directory, pattern = mutect_pattern, full.names = TRUE)
mutect_vcf_header <- read_lines(mutect_vcf_files[1])

vep_columns <- str_subset(mutect_vcf_header, "ID=ANN") %>% str_remove_all("^.*: ") %>% str_remove_all('\\".*$') %>% str_split("\\|") %>% extract2(1)

vcf_colnames <- str_subset(mutect_vcf_header, "^#") %>% str_subset("^##", negate = T) %>% str_split("\\t") %>% extract2(1)
vcf_colnames[length(vcf_colnames)] <- "DATA"
mutect_vcf_list <- map(mutect_vcf_files, read_lines) %>% map(str_subset, "^#", negate = TRUE) %>% map(str_split_fixed, "\\t", length(vcf_colnames)) %>% map(set_colnames, vcf_colnames) %>% map(as_tibble) %>% set_names(sample_names)
#mutect_vcf_all <- bind_rows(mutect_vcf_list, .id = "Sample")
#mutect_vcf_list <- map(mutect_vcf_files, read_lines)  %>% map(set_colnames, vcf_colnames) %>% map(as_tibble) %>% set_names(sample_names)
mutect_vcf_all <- bind_rows(mutect_vcf_list, .id = "Sample")

#mutect_vcf_all <- mutect_vcf_header %>% str_subset("^#", negate = TRUE) %>% str_split_fixed( "\\t", length(vcf_colnames)) %>% set_colnames(vcf_colnames) %>% as_tibble()

mutect_info_names <- str_split(mutect_vcf_all$INFO, ";") %>% map(str_remove_all, "=.*$") 
mutect_info <- str_split(mutect_vcf_all$INFO, ";") %>% map(str_remove_all, "^.*\\=") 
mutect_info_df <- map2(mutect_info_names, mutect_info, SetNames) %>% bind_rows

vep_info <- str_remove_all(mutect_info_df$ANN, "^\\[") %>% str_remove_all("\\]$") %>% str_split_fixed("\\|", length(vep_columns))  %>% as_tibble %>% set_colnames(vep_columns)
mutect_info_df_final <- select(mutect_info_df, -ANN) %>%  select(-DP)

print("(STEP 2) Parse VCF Based on Column Names")
mutect_data_names <- str_split(mutect_vcf_all$FORMAT, ":")  
mutect_data <- str_split(mutect_vcf_all$DATA, ":") 
mutect_data_df <- map2(mutect_data_names, mutect_data, SetNames) %>% bind_rows

print("(STEP 3) Format Column Names")
mutect_vcf_bind <- select(mutect_vcf_all, Sample:FILTER) %>% bind_cols(mutect_data_df) %>% bind_cols(vep_info) %>% bind_cols(mutect_info_df_final) 
colnames(mutect_vcf_bind) %<>% str_replace_all("AF...13", "tumor_f") %>% str_replace_all("#CHROM", "Chr")

mutect_vcf_ncol <- str_split(mutect_vcf_bind$AD, ",") %>% map_int(length) %>% max
mutect_vcf_ad <- str_split_fixed(mutect_vcf_bind$AD, ",", mutect_vcf_ncol) %>% as_tibble %>% mutate(across(everything(), as.numeric))
colnames(mutect_vcf_ad) <- c("t_ref_count", str_c("t_alt_count_", seq(1:(mutect_vcf_ncol - 1))))

mutect_vcf_af <- str_split_fixed(mutect_vcf_bind$tumor_f, ",", (mutect_vcf_ncol - 1)) %>% as_tibble %>% mutate(across(everything(), as.numeric))
colnames(mutect_vcf_af) <- str_c("tumor_f", seq(1:(mutect_vcf_ncol - 1)))

mutect_vcf_filter3 <- mutect_vcf_bind %>% select(-AD) %>% bind_cols(mutect_vcf_ad)

#### change column names to match previous pipeline ####
print("(STEP 6)  Filtered VCF File")
mutect_vcf_filter3_ald <- str_split_fixed(mutect_vcf_filter3$ALD, ",", 2) %>% as_tibble %>% mutate(across(everything(), as.numeric))
colnames(mutect_vcf_filter3_ald) <- c("AltFwd", str_c("AltRev", seq(1:(2 - 1))))
mutect_vcf_filter3_rd <- str_split_fixed(mutect_vcf_filter3$RD, ",", 2)  %>% as_tibble %>% mutate(across(everything(), as.numeric))
colnames(mutect_vcf_filter3_rd) <- c("RefFwd", str_c("RefRev", seq(1:(2 - 1))))

mutect_vcf_filter4 <- mutect_vcf_filter3 %>% 
  select(-ALD, -RD, -Gene) %>% 
  bind_cols(mutect_vcf_filter3_ald, mutect_vcf_filter3_rd) %>% 
  subset(!(SYMBOL == ""))

# change column names
colnames(mutect_vcf_filter4) %<>% 
  str_replace_all("^ID$", "CosmicID") %>%
  str_replace_all("Sample$", "SampleID") %>% 
  str_replace_all("#CHROM", "Chr") %>% 
  str_replace_all("POS", "Pos") %>% 
  str_replace_all("REF", "Ref") %>%
  str_replace_all("ALT", "Alt") %>%
  str_replace_all("EXON", "Exon") %>%
  str_replace_all("SYMBOL", "Gene") %>%
  str_replace_all("tumor_f", "VAF") %>%
  str_replace_all("Existing_variation", "ExistingVariation") %>%
  str_replace_all("PolyPhen", "Polyphen") %>%
  str_replace_all("DOMAINS", "Domains") %>%
  str_replace_all("CLIN_SIG", "ClinSig") %>%
  str_replace_all("SOMATIC", "Somatic") %>%
  str_replace_all("TCGA_SAMPLE", "TCGA") %>%
  str_replace_all("QUAL...7", "QUAL") %>%
  str_replace_all("ODDRATIO", "ODDRATIO") %>%
  str_replace_all("HICNT", "AltReads") %>% 
  str_replace_all("HICOV", "Depth") %>% 
  str_replace_all("HIAF", "HIAF") %>% 
  str_replace_all("AltRev1", "AltRev") %>%
  str_replace_all("RefRev1", "RefRev") %>%  
  str_replace_all("QUAL...104", "QMEAN")

# create a RefReads column, based on Depth - AltReads
# remove t_ref_count, t_alt_count_

print("(STEP 5) Save Aggregated VCF File")
mutect_basic_file <- str_c(mutect_directory, "/", output_file_name)
write_tsv(mutect_vcf_filter4, mutect_basic_file)


#### Filtering ####

# Un-comment for troubleshooting only
#mutect_basic_file <- str_c(mutect_directory, output_file_name)
# mutect_vcf_filter4 <- read_tsv(file = mutect_basic_file, col_types = cols(.default = col_character(),
#                                                                    VAF = col_numeric(), 
#                                                                    AltFwd = col_numeric(),
#                                                                    AltRev = col_numeric(),
#                                                                    Depth = col_numeric(),
#                                                                    EntropyLeft = col_numeric(),
#                                                                    EntropyCenter = col_numeric(),
#                                                                    EntropyRight = col_numeric(),
#                                                                    QUAL = col_numeric(), 
#                                                                    PSTD = col_numeric(),
#                                                                    PMEAN = col_numeric(), 
#                                                                    ODDRATIO = col_numeric(),
#                                                                    SBF = col_numeric(),
#                                                                    MSI = col_numeric(),
#                                                                    MSILEN = col_numeric()))

mutect_vcf_filter5 <- mutect_vcf_filter4 %>% filter(FILTER == "PASS")
mutect_vcf_filter5 <- mutect_vcf_filter5 %>%
  mutate_at('VAF', as.numeric) %>%
  mutate_at('AltFwd', as.numeric) %>%
  mutate_at('AltRev', as.numeric) %>%
  mutate_at('Depth', as.numeric) %>%
  mutate_at('EntropyLeft', as.numeric) %>%
  mutate_at('EntropyCenter', as.numeric) %>%
  mutate_at('EntropyRight', as.numeric) %>%
  mutate_at('QUAL', as.numeric) %>%
  mutate_at('PSTD', as.numeric) %>%
  mutate_at('PMEAN', as.numeric) %>%
  mutate_at('ODDRATIO', as.numeric) %>%
  mutate_at('SBF', as.numeric) %>%
  mutate_at('MSI', as.numeric) %>%
  mutate_at('MSILEN', as.numeric)

mutect_vcf_filter5$lesionID <- str_c(mutect_vcf_filter5$Chr,mutect_vcf_filter5$Pos,mutect_vcf_filter5$Ref,mutect_vcf_filter5$Alt,sep="|")
mutect_vcf_filter5$positionID <- str_c(mutect_vcf_filter5$Chr,mutect_vcf_filter5$Pos,mutect_vcf_filter5$Ref,sep="|")
mutect_vcf_filter5$lesionID2 <- str_c(mutect_vcf_filter5$Gene,mutect_vcf_filter5$HGVSc,mutect_vcf_filter5$HGVSp,sep="|")
# 
# filter_variants <- function(rawCalls){
#   rawCalls <- rawCalls[!duplicated(rawCalls),]
#   
#   lesionID <- as.character(paste(rawCalls$Chr,rawCalls$Pos,rawCalls$Ref,rawCalls$Alt,sep="|"))
#   positionID <- as.character(paste(rawCalls$Chr,rawCalls$Pos,rawCalls$Ref,sep="|"))
#   
#   rawCallsL <- data.frame(rawCalls,lesionID)
#   rawCallsL$lesionID2<-with(rawCallsL, paste(Gene,HGVSc,HGVSp, sep="|"))
#   
#   filteredCalls <- subset(rawCallsL,!(lesionID %in% blacklistLS) & !(positionID %in% blacklistPOS))
#   
#   #return filtered calls
#   filteredCalls
# }
# 
# 
# mutect_vcf_filter6 <- filter_variants(mutect_vcf_filter5)

# mutations in over 30% of samples
mutect_vcf_filter5$SpecimenID <- mutect_vcf_filter5$SampleID %>% 
  gsub("^Sample_", "", .) %>% 
  gsub("Cornell-", "",. ) %>% 
  gsub("_S[0-9][0-9]$", "",.) %>% 
  gsub("-06d[0-9][0-9][0-9][0-9]-PreCISE1", "", .) %>% 
  gsub("-PreCISE1v3_S[0-9][0-9]_L00[1-2]", "", .) %>% 
  gsub("-PreCISE1v3_S[0-9]_L00[1-2]", "", .) %>% 
  gsub("CHIP-", "", .) %>% 
  gsub("-06d[0-9][0-9][0-9][0-9]", "", .) %>% 
  gsub("-PB*", "", .) %>% 
  gsub("-BM*", "", .) %>% 
  gsub("_S[0-9]$", "", .) %>% 
  gsub("-[0-9]$", "",.) %>% 
  gsub("-Chicago-TWpilot$", "", .) %>% 
  gsub("Chicago-C", "", .)
mutect_vcf_filter5$SpecimenID <- unlist(mutect_vcf_filter5$SpecimenID)

# get counts for each variant
lesionidtb <- mutect_vcf_filter5 %>% 
  select(SpecimenID, lesionID2) %>% 
  unique() %>% 
  count(lesionID2) 
# calculate num of samples that would lead to a cut off at 35%
freq_cutoff <- mutect_vcf_filter5$SpecimenID %>% 
  unique() %>% 
  length() %>% 
  multiply_by(0.4)
# use cut off to filter variants with large counts
recurrent <- lesionidtb %>% 
  filter(n > freq_cutoff) %>% 
  select(lesionID2)

# Make sure COSMIC and TCGA variants are not filtered out
mutect_vcf_filter5$CosmicCount <- mutect_vcf_filter5$CNT %>% str_split(",") %>% map(1) %>% map(1) %>% as.numeric()
known.somatic <- mutect_vcf_filter5 %>% subset(CosmicCount >= 10 | grepl("TCGA", mutect_vcf_filter5$TCGA))
likely.artifact <- recurrent %>% subset(!(lesionID2 %in% known.somatic$lesionID2))

# remove variants in blacklist
args <- c("/home/es984/oelab_es984/PreCISE/blacklist.November2019.txt", "/home/es984/oelab_es984/PreCISE/blacklist.position.November2019.txt")

#args <- c("/Users/etisinhawork/Dropbox (Guzman Lab)/Guzman Lab/Eti Sinha/blacklist.January2022.txt",
#          "/Users/etisinhawork/Dropbox (Guzman Lab)/Guzman Lab/Eti Sinha/blacklist.position.January2022.txt")

# args <- c("/Users/etisinhawork/Dropbox (Guzman Lab)/Guzman Lab/Nuria Mencia/Precise1_samples/frequently_used_objects/blacklist.November2019.txt",
#           "/Users/etisinhawork/Dropbox (Guzman Lab)/Guzman Lab/Nuria Mencia/Precise1_samples/frequently_used_objects/blacklist.position.November2019.txt")

blacklistLS <- read.delim(args[1],as.is=T,header=F)[,1]
blacklistPOS <- read.delim(args[2],as.is=T,header=F)[,1]

# Filter reoccurring variants and blacklist variants
mutect_vcf_filter6 <- mutect_vcf_filter5 %>% 
  filter(!(lesionID %in% blacklistLS)) %>% 
  filter(!(positionID %in% blacklistPOS)) %>% 
  filter(!(lesionID2 %in% likely.artifact$lesionID2))

# filter on other columns
mutect_vcf_filter7 <- subset(mutect_vcf_filter6,
                               VAF > 0.01 &
                                 AltFwd >= 3 &
                                 AltRev >=3 &
                                 VAF*Depth >= 8 &
                                 Depth > 50 &
                                 VAF < 0.3 &
                                 (EntropyLeft > 1 & EntropyCenter > 1 & EntropyRight > 1) &
                                 !(Gene=="GNAS" & Exon == '1/13') &
                                 grepl('missense|splice_[ad]|stop_|frameshift|inframe_[id]',Consequence) &
                                 QUAL > 45 &
                                 PSTD != 0 &
                                 PMEAN > 25 &
                                 PMEAN < 75 &
                                 !((ODDRATIO > 2.5 & SBF < 0.1) | ODDRATIO == 0) &
                                 QMEAN > 30 &
                                 (MSI < 10 | (MSI*MSILEN > 5 & VAF > 0.1)))
mutect_vcf_filter7 <- mutect_vcf_filter7 %>% unique()

# Remove unnecessary columns
mutect_vcf_filter8 <- mutect_vcf_filter7 %>% select(SampleID:QUAL, GT, DP, VAF:VARIANT_CLASS, ENSP, SIFT:Domains, gnomAD_AF, ClinSig:VAR_SYNONYMS, EntropyAlt:SBF, CosmicCount:SpecimenID)

print("(STEP 7) Save Filtered VCF File")
mutect_basic_file <- str_c(mutect_directory, "/", filter_output_file_name)
write_tsv(mutect_vcf_filter8, mutect_basic_file)

