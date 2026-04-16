library(dplyr)
library(tidyr)
library(readr)
library(stringr)
setwd(".../") #set working directory


# DATASET GENERATION ----

## BA1 Lookup File ----
## Creates a lookup file using gnomAD v4.1 for variants which meet stand-alone BA1 criteria for benignity
## read in gnomAD v4.1 exome+genome files which contain variant data from all five genes (BRCA1, BRCA2, PALB2, ATM, CHEK2)
## vcf files can be downloaded from the gnomAD website (https://gnomad.broadinstitute.org/data#v4-core-dataset)
gnomAD_all <- read_tsv(".../")

# Quality filter and data cleaning
## AC0: Genotype quality too low, allele count is 0 after removing low quality genotypes
## AS VSQR: Failed VSQR filter
## RF: Flag for failing the random forest model quality check
## Exomes filtered = the exome dataset failed one of the above tests (usually AC0). Same for Genomes filtered.
## If both filtered, remove variant, as failed both datasets. 
gnomAD_all %>%
  filter(FILTER!="BOTH_FILTERED") %>%
  mutate(gene = case_when(grepl("17",`#CHROM`) ~ "BRCA1",
                          grepl("13",`#CHROM`) ~ "BRCA2",
                          grepl("16",`#CHROM`) ~ "PALB2",
                          grepl("11",`#CHROM`) ~ "ATM",
                          grepl("22",`#CHROM`) ~ "CHEK2"),
         `#CHROM` = gsub("chr","", `#CHROM`),
         ID = paste("chr",`#CHROM`, POS, REF, ALT, sep="_"),
         `#CHROM` = as.double(`#CHROM`)) -> gnomAD_all_genes

gnomAD_all_genes %>%
  mutate(
    # extract Grpmax FAF (which per gnomAD website does not include bottlenecked populations)
    grpmax_faf_v4 = case_when(FILTER=="PASS" ~ str_extract(INFO, "(?<=faf95_max_joint=).+?(?=;)"),
                              FILTER=="EXOMES_FILTERED" ~ str_extract(INFO, "(?<=faf95_max_genomes=).+?(?=;)"),
                              FILTER=="GENOMES_FILTERED" ~ str_extract(INFO, "(?<=faf95_max_exomes=).+?(?=;)")),
    across(ends_with("v4"), ~ as.numeric(.x))) -> gnomAD_all_genes_FAF

gnomAD_all_genes_FAF %>%
  # Identify variants which are more common than the MTAF for that gene
  mutate(BA1_filter_gnomad = case_when(gene=="BRCA1" & grpmax_faf_v4>BRCA1_MTAF ~ "BA1", 
                                       gene=="BRCA2" & grpmax_faf_v4>BRCA2_MTAF ~ "BA1", 
                                       gene=="PALB2" & grpmax_faf_v4>PALB2_MTAF ~ "BA1", 
                                       gene=="ATM" & grpmax_faf_v4>ATM_MTAF ~ "BA1", 
                                       gene=="CHEK2" & grpmax_faf_v4>CHEK2_MTAF ~ "BA1", 
                                       .default="Rare")) %>%
  filter(BA1_filter_gnomad == "BA1") -> gnomAD_BA1

write_csv(gnomAD_BA1, "gnomAD_BA1_filter_prelift.csv")

# liftover to build 37 using https://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core (Ensembl release v113)
# most of the datasets for this analysis use GrCh37, this stage allows merging with downstream datasets
read_tsv("gnomAD_BA1_filter_postlift.vcf", comment="##") -> gnomAD_BA1

gnomAD_BA1 %>% 
  mutate(ID=paste(`#CHROM`,POS,REF,ALT,sep="-")) %>% 
  select(ID) %>% 
  write_csv("gnomAD_vv_prep.csv")

# run liftover variants in the vv_prep csv file through variant validator https://variantvalidator.org/service/validate/batch/ 
read_tsv("gnomAD_vv_output.txt",comment="#") -> hgvs_annotations #variant validator output file
hgvs_annotations %>%
  # use MANE select transcript and produce column using HGVS nomenclature
  filter(`Select transcript`=="MANE") %>%
  separate_wider_delim("HGVS_transcript",delim=":",names_sep="part") %>%
  select("Input","gene"="Gene_Symbol","hgvs_cdna"="HGVS_transcriptpart2") -> hgvs_annotations

gnomAD_BA1 %>%
  # annotate BA1 lookup file with HGVS nomenclature and chromosome position information in format to match other datasets
  mutate(Input=paste(`#CHROM`,POS,REF,ALT,sep="-")) %>% 
  left_join(hgvs_annotations, by="Input") %>%
  mutate(ID=paste0("chr",`#CHROM`,"_",POS,"_",REF,"_",ALT)) %>% 
  mutate(rare_flag="BA1_eligible") %>%
  select("ID","gene","hgvs_cdna","rare_flag") -> gnomAD_BA1_list

# output final BA1 lookup file
write_csv(gnomAD_BA1_list, "gnomAD_BA1_filter_final.csv")


## UK Biobank ----
## White female breast cancer cases
## UK Biobank data is available on successful application: https://www.ukbiobank.ac.uk/use-our-data/apply-for-access/ 
## read in vcf file of cases here. vcf file should contain a per-variant list of het, hom, and non-carriers, broken down by self-reported patient ancestry.
## This vcf provided count breakdowns for white female individuals with and without the variant, and with and without breast cancer (see Sup Methods)
UKB_case <- read_tsv(".../")
UKB_case %>% 
  # convert allele counts to counts of individuals
  mutate(White_BC_female_carriers = white_female_BC_AC - white_female_BC_hom,
         White_BC_female_total = white_female_BC_AN/2) %>%
  select(-starts_with("white_", ignore.case=FALSE)) -> UKB_case

# Total non-breast cancer controls (male and female)
## UK Biobank data is available on successful application: https://www.ukbiobank.ac.uk/use-our-data/apply-for-access/ 
## read in vcf file. vcf file should contain a per-variant list of het, hom, and non-carriers, broken down by self-reported patient ancestry.
## This vcf provided count breakdown for all (male and female) individuals with BC (breast cancer) vs nonBC (no breast cancer) (see Supplementary Methods)
UKB_control <- read_tsv(".../") 
UKB_control %>%
  # create column of total tested (all non-carriers, all heterozygous carriers, and all homozygous carriers (het_hom = het carriers + hom carriers))
  mutate(White_nonBC_total = White_nonBC_het_hom+White_nonBC_non_carriers) %>%
  select(1:5, White_nonBC_total_carriers=White_nonBC_het_hom, White_nonBC_total) -> UKB_control
UKB_case %>% full_join(UKB_control, by=c("#CHROM","POS","ID","REF","ALT")) -> UKB

## variantvalidator check
# note: chromosomal coordiantes used in case HGVSc is different when using different transcripts
# VariantValidator run in batch with the 5 chosen transcripts (see Sup Methods) as options
UKB %>% mutate(ID = gsub("_",":",ID)) -> UKB_annotated

#this file submitted for batch validation in VariantValidator
write_csv(UKB_annotated, "UKB_prevv.csv") 

#VariantValidator output
read_tsv("UKB_postvv.txt",skip=2) -> ukb_vv 
left_join(UKB_annotated, ukb_vv, by=c("ID"="Input")) -> UKB

# filter out variants which fail variantvalidator check
UKB %>% filter(!is.na(HGVS_transcript)) -> UKB_vv

# excluding PTVs which are excluded from the BRIDGES dataset (in preparation for later combination of datasets)
UKB_vv %>%
  separate_wider_delim(HGVS_transcript, names=c("gene","hgvs"),delim=":") %>%
  filter(White_BC_female_carriers > 0 | White_nonBC_total_carriers > 0) %>%
  ## amending variant validator HGVS for this variant due to alternate allele presence, which should be T>C as reported in UKB.
  mutate(hgvs = str_replace(hgvs,"c.7397=","c.7397T>C")) %>%
  ## exclude PTVs in the last exon
  mutate(final_exon = case_when(SYMBOL=="BRCA1" & EXON=="23/23" ~ 1,
                                SYMBOL=="BRCA2" & EXON=="27/27" ~ 1,
                                SYMBOL=="ATM" & EXON=="63/63" ~ 1,
                                SYMBOL=="CHEK2" & EXON=="15/15" ~ 1,
                                SYMBOL=="PALB2" & EXON=="13/13" ~ 1,
                                .default=0)) %>%
  filter((final_exon==0 & grepl(paste(consequence_filter,collapse="|"),Consequence)) | grepl("missense",Consequence)) %>%
  ## exclude splice variants affecting the penultimate exon for genes where exon skipping not confirmed pathogenic at those sites (per Dorling et al., 2021)
  filter((SYMBOL=="BRCA2" & !grepl(paste("c.9502-1G","c.9502-2A","c.9648+1G","c.9648+2T",sep="|"),HGVSc)) | SYMBOL %in% c("BRCA1","PALB2","ATM","CHEK2")) %>%
  filter((SYMBOL=="CHEK2" & !grepl(paste("c.1462-1G","c.1462-2A","c.1542+1G","c.1542+2T",sep="|"),HGVSc)) | SYMBOL %in% c("BRCA1","PALB2","ATM","BRCA2")) %>%
  ## exclude 7 canonical BRCA1 VUS variants per ENIGMA: c.594-2A>C, c.4096+1G>A, c.4096+2T>C, c.4096+1G>A, c.4186-2A>G, c.4358-1G>C, c.4358-2del
  filter((SYMBOL=="BRCA1" & !grepl(paste("c.594-2A>C","c.4096+1G>A","c.4096+2T>C","c.4096+1G>A","c.4186-2A>G","c.4358-1G>C","c.4358-2del",sep="|"),HGVSc)) | SYMBOL %in% c("CHEK2","PALB2","ATM","BRCA2")) %>%
  select(SYMBOL, hgvs, ID, Consequence, starts_with("White")) -> UKB_filtered

# Split UKB controls across NDRS and UKB breast cancer cases
## separation of UKB controls into proportional segments (segment 1: for use with NDRS cases, segment 2: for use with UKB cases)
## aiming for best use of control data, division of controls should get the case/control ratio as close to 1 as possible for both case sets to reduce impact of imbalance (see Rowlands et al., 2024)
## Sup Figure 4: A case:control ratio of 6.62 is the lowest ratio obtainable across NDRS cases and UKB cases.
## Therefore, 122,232 controls are retained for UKB cases, and 297,141 controls are used for NDRS cases, and carriers are sampled where 29.1% are retained for UKB and 70.8% are 
UKB_filtered %>% mutate(White_nonBC_noncarrier = White_nonBC_total-White_nonBC_total_carriers) -> UKB_filtered

UKB_filtered %>%
  rowwise() %>%
  mutate(White_nonBC_total_noncarrier_UKBpart = round(White_nonBC_total*0.2914637),
         White_nonBC_total_noncarrier_NDRSpart = White_nonBC_total-White_nonBC_total_noncarrier_UKBpart,
         White_nonBC_total_carriers_UKBpart = sample_controls(White_nonBC_total_carriers),
         White_nonBC_total_carriers_NDRSpart = White_nonBC_total_carriers - White_nonBC_total_carriers_UKBpart) -> UKB_filtered_split

write_csv(UKB_filtered_split, "UKB_BRCA_Filtered_split_ptvmiss.csv")


## CARRIERS ----
## Data access for the CARRIERS dataset was attained on request and is subject to application and project approval from the respective team leads
carriers <- read_csv(".../")

# reformat nomenclature columns to align with other datasets
carriers %>%
  mutate(unique_id = paste0(CAVA_GENE, ":", c., ":", p.), .before = CAVA_CSN) %>%
  mutate(ID = paste0("chr",`#CHROM`,"_",GRCh37Location,"_",REF,"_",ALT), .before=unique_id) %>%
  select(-Sample_ID) %>%
  group_by(unique_id, ID, CAVA_SO, `case control status`, `race/ethnicity`) %>%
  summarise(count = n()) %>%
  pivot_wider(id_cols=c(unique_id,ID,CAVA_SO,`race/ethnicity`), names_from=`case control status`, values_from=`count`) %>%
  separate_wider_delim(unique_id, delim = ":", names=c("gene", "c.", "p.")) %>%
  mutate(c. = if_else((grepl("p.", p.) & !grepl("NA", p.)), c., paste(c., p., sep="_")),
         p. = if_else((grepl("p.", p.) & !grepl("NA", p.)), p., NA),
         c. = gsub("_NA","",c.),
         p. = gsub(" \\(p.Phe3065_Gln3066delinsLeuTer\\)","",p.),
         Case = replace_na(Case, 0),
         Control = replace_na(Control, 0)) %>%
  # this is needed because of the BRCA2 c.9196 variant with strange p. nomenclature, which we amended above
  group_by (gene, c., p., ID, CAVA_SO, `race/ethnicity`) %>%
  summarise (Case = sum(Case),
             Control = sum(Control)) -> carriers_per_var

# no VEP transcript annotations needed as carriers dataset already using build 37 versions of MANE select transcripts, and consequence already annotated

## variantvalidator check
carriers_per_var %>%
  mutate(vep_hgvs = case_when(gene == "BRCA1" ~ "NM_007294.3",
                              gene == "BRCA2" ~ "NM_000059.3",
                              gene == "PALB2" ~ "NM_024675.3",
                              gene == "ATM" ~ "NM_000051.3",
                              gene == "CHEK2" ~ "NM_007194.3"),
         vep_hgvs = paste0(vep_hgvs,":",c.)) -> carriers_per_var

#this file submitted for batch validation in VariantValidator
write_csv(carriers_per_var, "carriers_prevv.csv")

#variant validator output file
read_tsv("carriers_postvv.txt",skip=2) %>% distinct() -> carriers_vv 

left_join(carriers_per_var, carriers_vv, by=c("vep_hgvs"="Input")) -> carriers_per_var

# filter out variants which fail variantvalidator check
carriers_per_var %>% filter(!is.na(HGVS_transcript)) -> carriers_per_var_vv

# excluding PTVs which are excluded from the BRIDGES dataset (in preparation for later combination of datasets)
carriers_per_var_vv %>%
  separate_wider_delim(HGVS_transcript, names=c("hgvs_transcript","hgvs"),delim=":") %>%
  mutate (Case_denom = 32247,
          Control_denom = 32544, .after=Control) %>%
  filter(`race/ethnicity` == "Non_Hispanic_White") %>%
  filter (gene %in% c("BRCA1","BRCA2","PALB2","ATM","CHEK2")) %>%
  ## exclude PTVs in the last exon
  mutate(final_exon = case_when(gene=="BRCA1" & GRCh37_POS <= 41197819 ~ 1,
                                gene=="BRCA2" & GRCh37_POS >= 32972299 ~ 1,
                                gene=="ATM" & GRCh37_POS >= 108236052 ~ 1,
                                gene=="CHEK2" & GRCh37_POS <= 29083974 ~ 1,
                                gene=="PALB2" & GRCh37_POS <= 23614990 ~ 1,
                                .default=0)) %>%
  filter((final_exon==0 & grepl(paste(consequence_filter,collapse="|"),CAVA_SO)) | grepl("missense",CAVA_SO)) %>%
  ## exclude splice variants affecting the penultimate exon for genes where exon skipping not confirmed pathogenic at those sites
  ### BRCA2: 9502-1,-2, 9648+1,+2. CHEK2: c.1462-1,-2, 1542+1,+2
  filter((gene=="BRCA2" & !grepl(paste("c.9502-1G","c.9502-2A","c.9648+1G","c.9648+2T",sep="|"),c.)) | gene %in% c("BRCA1","PALB2","ATM","CHEK2")) %>%
  filter((gene=="CHEK2" & !grepl(paste("c.1462-1G","c.1462-2A","c.1542+1G","c.1542+2T",sep="|"),c.)) | gene %in% c("BRCA1","PALB2","ATM","BRCA2")) %>%
  ## exclude 7 canonical BRCA1 VUS variants per ENIGMA: c.594-2A>C, c.4096+1G>A, c.4096+2T>C, c.4096+1G>A, c.4186-2A>G, c.4358-1G>C, c.4358-2del
  filter((gene=="BRCA1" & !grepl(paste("c.594-2A>C","c.4096+1G>A","c.4096+2T>C","c.4096+1G>A","c.4186-2A>G","c.4358-1G>C","c.4358-2del",sep="|"),c.)) | gene %in% c("CHEK2","PALB2","ATM","BRCA2")
  ) -> carriers_per_var_clean

write_csv (carriers_per_var_clean, "carriers_per_var_clean_ptvmiss.csv")


## BRIDGES ----
## Summary data used from the BRIDGES consortium is available at the following link: https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/bridges-summary-results
## the below files contain variant counts for the five genes under examination from the missense and PTV files available at the above link
read_csv(".../") -> BRIDGES_missense
read_csv(".../") -> BRIDGES_ptv

# VEP transcript annotations not needed as bridges dataset already using build 37 versions of MANE select transcripts, and consequence already annotated

# generate list of IDs for variant validator 
## (need hgvs nomenclature for merging with other datasets, genomic coords not suitable because the coords in UKB aren't always accurate/appropriate, like chr13_32911461_TGAACAAATGGGCA_TGGACAAATGGGCA for a missense variant)
# note: chromosomal coordinates used in case HGVSc is different when using different transcripts
# VariantValidator run in batch with the 5 chosen transcripts as options
BRIDGES_missense %>%
  filter(Cases>0 | Controls >0) %>%
  mutate(VV_ID = ID,
         VV_ID = str_replace_all(VV_ID,"_","-"),
         VV_ID = str_replace_all(VV_ID,"chr","")) -> BRIDGES_missense
BRIDGES_missense %>% select(VV_ID) -> VV_IDs
write_csv(VV_IDs, "bridges_missense_vv_ids.csv")

BRIDGES_ptv %>%
  filter(Cases>0 | Controls >0) %>%
  mutate(VV_ID = ID,
         VV_ID = str_replace_all(VV_ID,"_","-"),
         VV_ID = str_replace_all(VV_ID,"chr","")) -> BRIDGES_ptv
BRIDGES_ptv %>% select(VV_ID) -> VV_IDs
write_csv(VV_IDs, "bridges_truncating_vv_ids.csv")

# read in variant validator outputs and clean up BRIDGES dataset for merging
read_tsv("missense_vv_output.txt",comment="#") -> BRIDGES_vv_missense
read_tsv("truncating_vv_output.txt",comment="#") -> BRIDGES_vv_ptv

BRIDGES_missense %>% 
  left_join(BRIDGES_vv_missense, by=c("VV_ID"="Input")) %>%
  # remove variants which fail variant validator check
  filter(!is.na(HGVS_transcript)) %>%
  separate_wider_delim(HGVS_transcript, names=c("gene","hgvs"),delim=":") %>%
  mutate(Case_denom = 42062, Control_denom = 44035, .after=Controls) %>%
  mutate(Gene_Symbol = case_when(grepl("chr11", ID) ~ "ATM",
                                 grepl("chr17", ID) ~ "BRCA1",
                                 grepl("chr13", ID) ~ "BRCA2",
                                 grepl("chr16", ID) ~ "PALB2",
                                 grepl("chr22", ID) ~ "CHEK2"), .before=ID) -> BRIDGES_processed_missense

BRIDGES_ptv %>% 
  left_join(BRIDGES_vv_ptv, by=c("VV_ID"="Input")) %>%
  # remove variants which fail variant validator check
  filter(!is.na(HGVS_transcript)) %>%
  separate_wider_delim(HGVS_transcript, names=c("gene","hgvs"),delim=":") %>%
  mutate(Case_denom = 42062, Control_denom = 44035, .after=Controls) %>%
  mutate(Gene_Symbol = case_when(grepl("chr11", ID) ~ "ATM",
                                 grepl("chr17", ID) ~ "BRCA1",
                                 grepl("chr13", ID) ~ "BRCA2",
                                 grepl("chr16", ID) ~ "PALB2",
                                 grepl("chr22", ID) ~ "CHEK2"), .before=ID) -> BRIDGES_processed_trunc

BRIDGES <- rbind(BRIDGES_processed_trunc,BRIDGES_processed_missense) 
write_csv(BRIDGES,"bridges_ptvmiss.csv")


## NDRS-UK ----
## Data access for NDRS data was attained on request and is subject to application and project approval from the respective team leads
## the below two files are .csv files of BRCA1 and BRCA2 variant counts from white, female individuals with breast cancer. 
read_csv(".../") -> ndrs_indels #annotated, but not pre-filtered for incorrect nomenclature
read_csv(".../") -> ndrs_singles #pre-filtered for incorrect nomenclature

# nomenclature formatting updates to align with other datasets
ndrs_indels %>%
  mutate(hgvs_cdna = gsub("del]", "del", hgvs_cdna),
         hgvs_cdna = gsub("del_","del",hgvs_cdna),
         hgvs_cdna = gsub("bp","",hgvs_cdna),
         hgvs_cdna = gsub("c.3847_c.3848heT_dell","c.3847_3848del",hgvs_cdna),
         hgvs_cdna = gsub("c.n", "c.", hgvs_cdna),
         hgvs_cdna = gsub("heT_", "", hgvs_cdna),
         hgvs_cdna = gsub("heTdel", "del", hgvs_cdna),
         hgvs_cdna = gsub("del-p", "del", hgvs_cdna),
         hgvs_cdna = gsub("dup-p", "del", hgvs_cdna),
         hgvs_cdna = gsub("c./", "c.", hgvs_cdna),
         hgvs_cdna = gsub("1)2", "1_2", hgvs_cdna),
         hgvs_cdna = gsub("7-1", "7_1", hgvs_cdna),
         hgvs_cdna = gsub("\\(\\)", "", hgvs_cdna),
         hgvs_cdna = gsub("del9ins", "delins", hgvs_cdna),
         hgvs_cdna = gsub("c.3680_3681del\\(leu1227Glnfs*5\\)", "c.3680_3681del", hgvs_cdna),
         hgvs_cdna = gsub("c.10095delinsGAATTATATCT\\(p.ser3366fs\\)", "c.10095delinsGAATTATATCT", hgvs_cdna),
         hgvs_cdna = gsub("(r645fsx15)", "", hgvs_cdna),
         hgvs_cdna = gsub("c.2800_01del", "c.2800_2801del", hgvs_cdna),
         hgvs_cdna = gsub("c.3847_48del", "c.3847_3848del", hgvs_cdna),
         hgvs_cdna = gsub("c.1942-1945del","c.1942_1945del",hgvs_cdna),
         hgvs_cdna = gsub("c.4065_68del", "c.4065_4068del", hgvs_cdna),
         hgvs_cdna = gsub("c.4065_4068del(p.Asn1355fs)","c.4065_4068del",hgvs_cdna),
         hgvs_cdna = gsub("c.4065-4068del","c.4065_4068del",hgvs_cdna),
         hgvs_cdna = gsub("c.6280_81del", "c.6280_6281del", hgvs_cdna),
         hgvs_cdna = gsub("c.6763_64insA", "c.6763_6764insA", hgvs_cdna),
         hgvs_cdna = gsub("c.983_86del", "c.983_986del", hgvs_cdna),
         hgvs_cdna = gsub("c.5266dup-p", "c.5266dupC", hgvs_cdna),
         hgvs_cdna = gsub("c.5266dupC,c.288C>T", "c.5266dupC", hgvs_cdna),
         hgvs_cdna = gsub("c.2157del\\(\\)", "c.2157del", hgvs_cdna),
         hgvs_cdna = gsub("c.3847_3847del", "c.3847_3848del", hgvs_cdna),
         hgvs_cdna = gsub("c.6275_6275del", "c.6275_6276del", hgvs_cdna),
         hgvs_cdna = gsub("c.3481_3491del11\\*\\*","c.3481_3491del11",hgvs_cdna),
         hgvs_cdna = gsub("c.891_899delinsGATACTTAGc.891_899delinsGATACTTCAG", "c.891_899delinsGATACTTCAG", hgvs_cdna),
         hgvs_cdna = gsub("c.3036-3039del","c.3036_3039del",hgvs_cdna)
  ) -> ndrs_indels

# Initial formatting of NDRS data to group based on cDNA HGVS
ndrs_singles %>%
  select(gene, hgvs_cdna, WhiteCount) %>%
  mutate(WhiteCount = replace_na(WhiteCount, 0)) %>%
  group_by (gene, hgvs_cdna) %>%
  summarise (across (everything (), ~ sum(.x))) %>%
  mutate(vartype="missense") %>%
  filter(WhiteCount != 0) -> ndrs_singles_grouped

ndrs_indels %>%
  select(gene, hgvs_cdna, WhiteCount=`Full screen count`) %>%
  mutate(WhiteCount = replace_na(WhiteCount, 0)) %>%
  group_by (gene, hgvs_cdna) %>%
  summarise (across (everything (), ~ sum(.x))) %>%
  mutate(vartype="PTV") %>%
  filter(WhiteCount != 0) -> ndrs_indels_grouped

rbind(ndrs_singles_grouped, ndrs_indels_grouped) -> ndrs

# VEP annotation - converts HGVS to VEP format using same transcripts as other datasets (see Supplementary Methods)
ndrs %>%
  mutate(vep_hgvs = if_else(gene == "BRCA1",
                            # BRCA1 transcript
                            paste0("ENST00000357654.3", ":", hgvs_cdna),
                            # BRCA2 transcript
                            paste0("ENST00000544455.1", ":", hgvs_cdna))) -> ndrs

# the below file is run through VEP
write_csv(ndrs, "ndrs_vep_file.csv")
# VEP call results
read_tsv("ndrs_vep_annotation_output_preVEPcleaning.txt") -> vep_annotations

vep_annotations %>%
  filter(grepl("ENST00000357654.3|ENST00000544455.1", Feature)) %>%
  select("#Uploaded_variation", "Location", "Consequence", "SYMBOL", "EXON", "INTRON") %>%
  filter(grepl("BRCA", SYMBOL)) -> vep_annotations

left_join(ndrs, vep_annotations, join_by("vep_hgvs" == "#Uploaded_variation")) -> ndrs_annotated

# variants which were identified by VEP:
filter(ndrs_annotated, !is.na(Consequence)) -> ndrs_annotated

# variant validator file (standardise the hgvs nomenclature to match the other datasets)
ndrs_annotated %>%
  mutate(vep_hgvs = case_when(gene == "BRCA1" ~ "NM_007294.3",
                              gene == "BRCA2" ~ "NM_000059.3",
                              gene == "PALB2" ~ "NM_024675.3",
                              gene == "ATM" ~ "NM_000051.3",
                              gene == "CHEK2" ~ "NM_007194.3"),
         vep_hgvs = paste0(vep_hgvs,":",hgvs_cdna)) -> ndrs_annotated

# below file is submitted for batch validation in VariantValidator
write_csv(ndrs_annotated, "ndrs_prevv.csv")

# read in results from VariantValidator
read_tsv("ndrs_postvv.txt",skip=2) %>% distinct() -> ndrs_vv #variant validator output
left_join(ndrs_annotated, ndrs_vv, by=c("vep_hgvs"="Input")) -> ndrs_annotated

# filter out variants which fail VariantValidator check
ndrs_annotated %>% filter(!is.na(HGVS_transcript)) -> ndrs_annotated_vv

# excluding PTVs which are excluded from the BRIDGES dataset (in preparation for later combination of datasets)
ndrs_annotated_vv %>%
  separate_wider_delim(HGVS_transcript, names=c("hgvs_transcript","hgvs"), delim=":") %>%
  # sum counts across variants which are now identical ie have same nomenclature + consequence
  group_by (gene, hgvs, Consequence, EXON, INTRON) %>%
  summarise (WhiteCount = sum(WhiteCount)) %>%
  mutate(denominator = 44917) %>%
  ## exclude PTVs in the last exon and restrict to just PTVs and missense variants
  mutate(final_exon = case_when(gene=="BRCA1" & EXON=="23/23" ~ 1,
                                gene=="BRCA2" & EXON=="27/27" ~ 1,
                                .default=0)) %>%
  filter((final_exon==0 & grepl(paste(consequence_filter,collapse="|"),Consequence)) | grepl("missense",Consequence)) %>%
  ## exclude splice variants affecting the penultimate exon for genes where exon skipping not confirmed pathogenic at those sites
  ### BRCA2: 9502-1,-2, 9648+1,+2. CHEK2: c.1462-1,-2, 1542+1,+2
  filter((gene=="BRCA2" & !grepl(paste("c.9502-1G","c.9502-2A","c.9648+1G","c.9648+2T",sep="|"),hgvs)) | gene %in% c("BRCA1","PALB2","ATM","CHEK2")) %>%
  ## exclude 7 canonical BRCA1 VUS variants per ENIGMA: c.594-2A>C, c.4096+1G>A, c.4096+2T>C, c.4096+1G>A, c.4186-2A>G, c.4358-1G>C, c.4358-2del
  filter((gene=="BRCA1" & !grepl(paste("c.594-2A>C","c.4096+1G>A","c.4096+2T>C","c.4096+1G>A","c.4186-2A>G","c.4358-1G>C","c.4358-2del",sep="|"),hgvs)) | gene %in% c("CHEK2","PALB2","ATM","BRCA2")
  ) -> ndrs_annotated_cleaned

# merge with partition of UKB controls
read_csv("UKB_BRCA_Filtered_split_ptvmiss.csv") -> UKB_filtered_split

ndrs_annotated_cleaned %>%
  ungroup() %>%
  full_join(UKB_filtered_split, by=c("gene"="SYMBOL","hgvs")) %>%
  filter (gene %in% c("BRCA1", "BRCA2")) %>%
  select(gene,hgvs,WhiteCount,denominator,White_nonBC_total_carriers_NDRSpart,White_nonBC_total_noncarrier_NDRSpart) -> ndrs_ukb_combined

write_csv(ndrs_ukb_combined, "ndrs_ptvmiss_whiteonly.csv")


## AMBRY + gnomAD ----
### Ambry ----
### Data access for Ambry data was attained on request and is subject to application and project approval from the respective team leads.
read_csv(".../") -> Ambry

### filter for the five genes of interest, align nomenclature format, and summarise phenotype information
Ambry %>%
  distinct() %>%
  filter(symbol %in% c("ATM","CHEK2","PALB2","BRCA1","BRCA2")) %>%
  mutate(nucleotide_id = case_when(symbol=="ATM" ~ "NM_000051.4",
                                   symbol=="CHEK2" ~ "NM_007194.4",
                                   symbol=="PALB2" ~ "NM_024675.4",
                                   symbol=="BRCA1" ~ "NM_007294.4",
                                   symbol=="BRCA2" ~ "NM_000059.4")) %>%
  mutate(ID=paste0(nucleotide_id,":",c_variant),.before=gender) %>%
  mutate(gene_ID=paste0(symbol,":",c_variant),.before=gender) %>%
  mutate(ID=gsub("DEL","del",ID)) %>%
  mutate(ID=gsub("DUP","dup",ID)) %>%
  mutate(ID=gsub("INS","ins",ID)) %>%
  mutate(ID=gsub("Ex","EX",ID)) %>%
  mutate(gene_ID=gsub("DEL","del",gene_ID)) %>%
  mutate(gene_ID=gsub("DUP","dup",gene_ID)) %>%
  mutate(gene_ID=gsub("INS","ins",gene_ID)) %>%
  mutate(gene_ID=gsub("Ex","EX",gene_ID)) %>%
  mutate(genotype=gsub("N/A",NA,genotype)) %>%
  mutate(phenotype=case_when(gender=="female" & !is.na(`PHx Breast`) ~ "female_breast",
                             gender=="female" & !is.na(`PHx OvCa`) ~ "female_ovarian",
                             gender=="male" & !is.na(`PHx Breast`) ~ "male_breast",
                             gender=="male" & !is.na(`PHx OvCa`) ~ "male_ovarian",
                             (gender=="unknown" | is.na(gender)) & !is.na(`PHx Breast`) ~ "other_breast",
                             (gender=="unknown" | is.na(gender)) & !is.na(`PHx OvCa`) ~ "other_ovarian")) %>%
  mutate(fhx = case_when(!is.na(`FHx Breast`) &  !is.na(`FHX OvCa`) ~ "fhx_breastovarian",
                         !is.na(`FHx Breast`) ~ "fhx_breast",
                         !is.na(`FHX OvCa`) ~ "fhx_ovarian")) -> ambry_interim

# column of total number of variant observations
ambry_interim %>%
  select(ID,gene_ID,genotype,phenotype,fhx) %>%
  group_by(ID,gene_ID) %>%
  summarise(n=n()) %>%
  ungroup() -> ambry_totalno

# number of each gentoype observed
ambry_interim %>%
  select(ID,gene_ID,genotype,phenotype,fhx) %>%
  group_by(ID,gene_ID,genotype) %>%
  summarise(n=n()) %>%
  mutate(genotype=replace_na(genotype, "PHASE_UNK")) %>%
  ungroup() %>%
  pivot_wider(id_cols=c(ID,gene_ID),names_from=genotype,values_from=n) -> ambry_genotypes

# number of each cancer phenotype observed
ambry_interim %>%
  select(ID,gene_ID,genotype,phenotype,fhx) %>%
  group_by(ID,gene_ID,phenotype) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  pivot_wider(id_cols=c(ID,gene_ID),names_from=phenotype,values_from=n) -> ambry_phenotypes

# combine into one summary table
ambry_totalno %>%
  left_join(ambry_genotypes,by=c("ID","gene_ID")) %>%
  left_join(ambry_phenotypes,by=c("ID","gene_ID")) %>%
  mutate(breast_denom=187642,.after=female_breast) -> ambry_summary

# submit the below file for VEP
write_csv(ambry_summary,"ambry_per_variant.csv")

# read in and add VEP annotations. Use gene_id column to submit data through VEP.
read_tsv("ambry_VEP_build38.txt") -> build38_VEP #VEP output file
build38_VEP %>% 
  filter(Feature %in% c("ENST00000675843.1","ENST00000357654.9","ENST00000404276.6","ENST00000261584.9","ENST00000380152.8")) %>%
  select(`#Uploaded_variation`, Consequence, EXON) -> build38_VEP

ambry_summary %>% 
  full_join(build38_VEP,by=c("gene_ID"="#Uploaded_variation")) %>%
  # filter for variants seen at least once in female breast cancer patients
  filter(!is.na(female_breast)) -> ambry_summary_VEP

# run through variant validator (confirm nomenclature is standard/accurate before merging)
## NOTE: VV run using Build 38 transcripts here using the 'ID' column!!
write_csv(ambry_summary_VEP, "ambry_pre_vv.csv")

# variant validator output file
read_tsv("ambry_postvv.txt",skip=2) %>% distinct() -> ambry_vv #
left_join(ambry_summary_VEP, ambry_vv, by=c("ID"="Input")) -> ambry_summary_VEP

# filter out variants which fail variantvalidator check
ambry_summary_VEP %>% filter(!is.na(HGVS_transcript)) -> ambry_summary_VEP_VV

# note: liftover done as part of variantvalidator check in build 38 (provides coordinates for both build 37 and 38 when providing gene and hgvs)
# re-load build 37 coordinates with correct transcripts to re-generate hgvs in case any difference between builds
ambry_summary_VEP_VV %>%
  mutate(ID_build37 = paste0("chr",GRCh37_CHR,":",GRCh37_POS,":",GRCh37_REF,":",GRCh37_ALT),.before=ID) -> ambry_summary_VEP_VV
write_csv(ambry_summary_VEP_VV, "ambry_per_variant_VEP_VV_prelift.csv")

read_tsv("ambry_per_variant_VEP_VV_postlift.txt", skip=2) -> ambry_vep_build37 #variant validator output file using build 37 genomic coordinates

ambry_summary_VEP_VV %>%
  left_join(ambry_vep_build37, by=c("ID_build37"="Input")) %>%
  # some variants appear duplicated because they were listed in the Ambry file under both old and new HGVS nomenclature (so left join appears many-to-many)
  # nomenclature has been standardized using variantvalidator, therefore can be combined. 
  distinct() %>%
  select(-ID,-gene_ID) %>%
  separate_wider_delim(HGVS_transcript.y, names=c("gene","hgvs"), delim=":") %>%
  group_by (ID_build37, gene, hgvs) %>%
  summarise (female_breast = sum(female_breast),
             breast_denom=breast_denom) %>%
  distinct() -> ambry_counts_fixed

write_csv(ambry_counts_fixed, "ambry_counts_fixed.csv")


### gnomAD v4.1 ----
### Load in exome data files. NOTE: gnomAD v4.1 is in build 38
## Exome files available from the gnomAD website https://gnomad.broadinstitute.org/data#v4-core-dataset 
gnomAD_all <- read_tsv(".../")

## AC0: Genotype quality too low, allele count is 0 after removing low quality genotypes
## AS VSQR: Failed VSQR filter
## Filter to retain only variants which pass
gnomAD_all %>%
  filter(FILTER=="PASS") %>%
  mutate(gene = case_when(grepl("17",`#CHROM`) ~ "BRCA1",
                          grepl("13",`#CHROM`) ~ "BRCA2",
                          grepl("16",`#CHROM`) ~ "PALB2",
                          grepl("11",`#CHROM`) ~ "ATM",
                          grepl("22",`#CHROM`) ~ "CHEK2"),
         `#CHROM` = gsub("chr","", `#CHROM`),
         ID = paste(`#CHROM`, POS, REF, ALT, sep="-"),
         `#CHROM` = as.double(`#CHROM`)) -> gnomAD_exome_all_genes

gnomAD_exome_all_genes %>%
  mutate(
    # extract AC and AN for female NFE non-UKB partition - for Ambry data prep only
    AC_NFE_v4 = str_extract(INFO, "(?<=AC_non_ukb_nfe=)[:digit:]+"),
    AN_NFE_v4 = str_extract(INFO, "(?<=AN_non_ukb_nfe=)[:digit:]+"),
    HOM_NFE_v4 = str_extract(INFO, "(?<=nhomalt_non_ukb_nfe=)[:digit:]+"),
    grpmax_faf_NFE_v4 = str_extract(INFO, "(?<=faf95_non_ukb_nfe=).+?(?=;)"),
    across(ends_with("v4"), ~ as.numeric(.x))) %>%
  
  # remove variants if not seen (AC_NFE_v4 = 0)
  filter(AC_NFE_v4 != 0) %>%
  
  # convert counts to individuals
  mutate(AC_NFE_v4 = AC_NFE_v4-HOM_NFE_v4,
         AN_NFE_v4 = AN_NFE_v4/2) -> gnomAD_prevep

write_csv(gnomAD_prevep, "gnomAD_exomes_prevep.csv")

# annotate gnomAD data with VEP annotations including HGVS nomenclature
read_tsv("gnomAD_vep_annotation_exomes.txt") -> VEP #VEP output file

VEP %>%
  filter(Feature %in% c("ENST00000675843.1","ENST00000357654.9","ENST00000404276.6","ENST00000261584.9","ENST00000380152.8")) %>%
  separate_wider_delim(HGVSc, delim=":", names=c("Transcript","hgvs")) %>%
  select(`#Uploaded_variation`, hgvs, Consequence, EXON) -> VEP

gnomAD_prevep %>%
  left_join(VEP,by=c("ID"="#Uploaded_variation")) -> gnomAD_postvep

write_csv(gnomAD_postvep, "gnomAD_exomes_postvep.csv")

# liftover gnomAD to build 37
gnomAD_postvep %>% select(1:7) %>% mutate(INFO = ".") -> gnomAD_preliftover
write_csv(gnomAD_preliftover, "gnomAD_exomes_preliftover.csv")

# liftover output file
read_tsv("gnomAD_exomes_preliftover.vcf", col_names=FALSE) -> gnomAD_liftover 

gnomAD_postvep %>%
  left_join(gnomAD_liftover, by=c("ID"="X3")) %>%
  rename("CHROM_BUILD37"="X1", "POS_BUILD37"="X2","REF_BUILD37"="X4","ALT_BUILD37"="X5","CHROM_BUILD38"="#CHROM","POS_BUILD38"="POS","ID_BUILD38"="ID","REF_BUILD38"="REF","ALT_BUILD38"="ALT") %>%
  select(-X6,-X7,-X8) %>%
  mutate(gene_ID = paste0(gene,":",hgvs)) -> gnomAD_annotated

write_csv(gnomAD_annotated, "gnomAD_exomes_annotated.csv")

# final variant validator check using build 37 coordinates (in case hgvs format from build 37 VEP is different; ensures all datasets use the same origin for hgvs nomenclature)
gnomAD_annotated %>% 
  mutate(ID_build37 = paste0("chr",CHROM_BUILD37,":",POS_BUILD37,":",REF_BUILD37,":",ALT_BUILD37)) -> gnomAD_annotated_prevv

write_csv(gnomAD_annotated_prevv, "gnomAD_annotated_prevv.csv")

#variant validator output file
read_tsv("gnomAD_annotated_postvv.txt",skip=2) -> gnomAD_annotated_postvv 
left_join(gnomAD_annotated_prevv, gnomAD_annotated_postvv, by=c("ID_build37"="Input")) -> gnomAD_annotated_postvv

# filter out variants which fail variantvalidator check
gnomAD_annotated_postvv %>% filter(!is.na(HGVS_transcript))-> gnomAD_annotated_postvv
write_csv(gnomAD_annotated_postvv,"gnomAD_clean.csv")


### Combine ----
read_csv("ambry_per_variant_VEP_VV_prelift.csv") -> ambry_summary_VEP_VV
read_tsv("ambry_per_variant_VEP_VV_postlift.txt", skip=2) -> ambry_vep_build37
read_csv("ambry_counts_fixed.csv") -> ambry_counts_fixed
read_csv("gnomAD_clean.csv") -> gnomAD_annotated_postvv

# excluding PTVs from Ambry which are excluded from the BRIDGES dataset (in preparation for later combination of datasets)
ambry_summary_VEP_VV %>%
  ## exclude PTVs in the last exon
  mutate(final_exon = case_when(Gene_Symbol=="BRCA1" & EXON=="23/23" ~ 1,
                                Gene_Symbol=="BRCA2" & EXON=="27/27" ~ 1,
                                Gene_Symbol=="ATM" & EXON=="63/63" ~ 1,
                                Gene_Symbol=="CHEK2" & EXON=="15/15" ~ 1,
                                Gene_Symbol=="PALB2" & EXON=="13/13" ~ 1,
                                .default=0)) %>%
  filter((final_exon==0 & grepl(paste(consequence_filter,collapse="|"),Consequence)) | grepl("missense",Consequence)) %>%
  
  # join with counts and VV annotations
  left_join(ambry_vep_build37, by=c("ID_build37"="Input")) %>%
  select(-ID,-gene_ID,-n:-Warnings.x) %>%
  distinct() %>%
  left_join(ambry_counts_fixed,by="ID_build37") %>%
  filter(!is.na(HGVS_transcript.y)) %>% 
  
  ## exclude splice variants affecting the penultimate exon for genes where exon skipping not confirmed pathogenic at those sites
  ### BRCA2: 9502-1,-2, 9648+1,+2. CHEK2: c.1462-1,-2, 1542+1,+2
  filter((gene=="NM_000059.3" & !grepl(paste("c.9502-1G","c.9502-2A","c.9648+1G","c.9648+2T",sep="|"),hgvs)) | gene %in% c("NM_007294.3","NM_024675.3","NM_000051.3","NM_007194.3")) %>%
  filter((gene=="NM_007194.3" & !grepl(paste("c.1462-1G","c.1462-2A","c.1542+1G","c.1542+2T",sep="|"),hgvs)) | gene %in% c("NM_007294.3","NM_024675.3","NM_000051.3","NM_000059.3")) %>%
  ## exclude 7 canonical BRCA1 VUS variants per ENIGMA: c.594-2A>C, c.4096+1G>A, c.4096+2T>C, c.4096+1G>A, c.4186-2A>G, c.4358-1G>C, c.4358-2del
  filter((gene=="NM_007294.3" & !grepl(paste("c.594-2A>C","c.4096+1G>A","c.4096+2T>C","c.4096+1G>A","c.4186-2A>G","c.4358-1G>C","c.4358-2del",sep="|"),hgvs)) | gene %in% c("NM_007194.3","NM_024675.3","NM_000051.3","NM_000059.3")
  ) %>%
  mutate(gene = case_when(gene=="NM_007294.3" ~ "BRCA1",
                          gene=="NM_000059.3" ~ "BRCA2",
                          gene=="NM_024675.3" ~ "PALB2",
                          gene=="NM_000051.3" ~ "ATM",
                          gene=="NM_007194.3" ~ "CHEK2")) -> ambry_summary_VEP

# excluding PTVs from gnomAD which are excluded from the BRIDGES dataset (in preparation for later combination of datasets)
gnomAD_annotated_postvv %>%
  select(-gene,-hgvs) %>%
  separate_wider_delim(HGVS_transcript, names=c("gene","hgvs"), delim=":") %>%
  mutate(gene = case_when(gene=="NM_007294.3" ~ "BRCA1",
                          gene=="NM_000059.3" ~ "BRCA2",
                          gene=="NM_024675.3" ~ "PALB2",
                          gene=="NM_000051.3" ~ "ATM",
                          gene=="NM_007194.3" ~ "CHEK2")) %>%
  ## exclude PTVs in the last exon
  mutate(final_exon = case_when(gene=="BRCA1" & EXON=="23/23" ~ 1,
                                gene=="BRCA2" & EXON=="27/27" ~ 1,
                                gene=="ATM" & EXON=="63/63" ~ 1,
                                gene=="CHEK2" & EXON=="15/15" ~ 1,
                                gene=="PALB2" & EXON=="13/13" ~ 1,
                                .default=0)) %>%
  filter((final_exon==0 & grepl(paste(consequence_filter,collapse="|"),Consequence)) | grepl("missense",Consequence)) %>%
  ## exclude splice variants affecting the penultimate exon for genes where exon skipping not confirmed pathogenic at those sites
  ### BRCA2: 9502-1,-2, 9648+1,+2. CHEK2: c.1462-1,-2, 1542+1,+2
  filter((gene=="BRCA2" & !grepl(paste("c.9502-1G","c.9502-2A","c.9648+1G","c.9648+2T",sep="|"),hgvs)) | gene %in% c("BRCA1","PALB2","ATM","CHEK2")) %>%
  filter((gene=="CHEK2" & !grepl(paste("c.1462-1G","c.1462-2A","c.1542+1G","c.1542+2T",sep="|"),hgvs)) | gene %in% c("BRCA1","PALB2","ATM","BRCA2")) %>%
  ## exclude 7 canonical BRCA1 VUS variants per ENIGMA: c.594-2A>C, c.4096+1G>A, c.4096+2T>C, c.4096+1G>A, c.4186-2A>G, c.4358-1G>C, c.4358-2del
  filter((gene=="BRCA1" & !grepl(paste("c.594-2A>C","c.4096+1G>A","c.4096+2T>C","c.4096+1G>A","c.4186-2A>G","c.4358-1G>C","c.4358-2del",sep="|"),hgvs)) | gene %in% c("CHEK2","PALB2","ATM","BRCA2")
  ) %>%
  
  # merge gnomAD counts with Ambry counts
  full_join(ambry_summary_VEP, by=c("gene","hgvs")) %>%
  # c,7397 variant is the same variant as in the UKB dataset which appears incorrectly in the variant validator nomenclature due to alternative allele presence. Amended here to match with UKB dataset.
  mutate(hgvs = str_replace_all(hgvs,"c.7397=","c.7397T>C")) -> ambry_gnomAD_merged

write_csv(ambry_gnomAD_merged, "ambry_ptvmiss_merged.csv")

## Collate variant data ----
## (the six output files below produced from the 'Integrating_Enriched_Data_Cleaning.R' script) 
UKB <- read_csv("UKB_BRCA_Filtered_split_ptvmiss.csv")
CARRIERS <- read_csv("carriers_per_var_clean_ptvmiss.csv")
BRIDGES <- read_csv("bridges_ptvmiss.csv")
NDRS <- read_csv("ndrs_ptvmiss_whiteonly.csv")
AMBRY <- read_csv("ambry_ptvmiss_merged.csv")
gnomAD_BA1_list <- read_csv("gnomAD_BA1_filter_final.csv")

NDRS %>%
  full_join(AMBRY, join_by("hgvs", "gene")) %>%
  full_join(UKB, join_by("hgvs", "gene"=="SYMBOL")) %>%
  full_join(BRIDGES, join_by("hgvs", "gene"=="Gene_Symbol")) %>%
  full_join(CARRIERS, join_by("hgvs", "gene")) -> joined

# remove very common variants
joined %>%
  left_join(gnomAD_BA1_list, by=c("gene","hgvs"="hgvs_cdna")) %>%
  filter(is.na(rare_flag)) -> joined_rare 

joined_rare %>%
  # set empty values, add or estimate denominators for all variants
  select(gene,hgvs,
         NDRS_case_carrier=WhiteCount,NDRS_case_total=denominator,NDRS_control_carrier=White_nonBC_total_carriers_NDRSpart.x, NDRS_control_total=White_nonBC_total_noncarrier_NDRSpart.x,
         AMBRY_case_carrier=female_breast,AMBRY_case_total=breast_denom,AMBRY_control_carrier=AC_NFE_v4,AMBRY_control_total=AN_NFE_v4,
         UKB_case_carrier=White_BC_female_carriers,UKB_case_total=White_BC_female_total,UKB_control_carrier=White_nonBC_total_carriers_UKBpart,UKB_control_total=White_nonBC_total_noncarrier_UKBpart, # for partitioned controls - UKB split
         BRIDGES_case_carrier=Cases,BRIDGES_case_total=Case_denom.x,BRIDGES_control_carrier=Controls,BRIDGES_control_total=Control_denom.x,
         CARRIERS_case_carrier=Case,CARRIERS_case_total=Case_denom.y,CARRIERS_control_carrier=Control,CARRIERS_control_total=Control_denom.y) %>%
  mutate(across(3:22, ~ replace_na(.x, 0))) %>%
  # Replace 0 values in the denominators with average allele denominator for that dataset
  mutate(BRIDGES_control_total=44035,
         BRIDGES_case_total=42062,
         UKB_control_total = if_else(UKB_control_carrier==0 & UKB_control_total==0,122232,UKB_control_total), #for partitioned controls (UKB partition)
         UKB_case_total = if_else(UKB_case_carrier==0 & UKB_case_total==0,18477,UKB_case_total),
         NDRS_control_total = if_else(NDRS_control_carrier==0 & NDRS_control_total==0,297141,NDRS_control_total), #for partitioned controls (NDRS partition)
         NDRS_case_total = 44917,
         AMBRY_control_total = if_else(AMBRY_control_carrier==0 & AMBRY_control_total==0,175054,AMBRY_control_total),
         AMBRY_case_total = 187642,
         CARRIERS_case_total=32247,
         CARRIERS_control_total=32544) -> joined_names

write_csv(joined_names, "combined_ptvmissvar_rare_cleaned.csv")
