## Load library
library(tidyverse)

## Define data paths
in_path_varscan = 'data/maf/6c93f518-1956-4435-9806-37185266d248/TCGA.BRCA.varscan.6c93f518-1956-4435-9806-37185266d248.DR-10.0.somatic.maf.gz'
in_path_muse = 'data/maf/b8ca5856-9819-459c-87c5-94e91aca4032/TCGA.BRCA.muse.b8ca5856-9819-459c-87c5-94e91aca4032.DR-10.0.somatic.maf.gz'
in_path_ss = 'data/maf/7dd592e3-5950-4438-96d5-3c718aca3f13/TCGA.BRCA.somaticsniper.7dd592e3-5950-4438-96d5-3c718aca3f13.DR-10.0.somatic.maf.gz'
in_path_mutect = 'data/maf/995c0111-d90b-4140-bee7-3845436c3b42/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.gz'


## Define column types
col_types = cols_only(
  Chromosome = 'c',
  Start_Position = 'i',
  End_Position = 'i',
  SYMBOL = 'c',
  Reference_Allele = 'c',
  Allele = 'c', 
  Variant_Classification = 'c', 
  IMPACT = 'c',
  Variant_Type = 'c',
  Tumor_Sample_Barcode = 'c'
)

## Concatenate vertically
in_paths = c(
  in_path_varscan,
  in_path_muse,
  in_path_ss,
  in_path_mutect
)
df = bind_rows(lapply(
  in_paths, 
  read_tsv, 
  col_types = col_types,
  comment = '#'
))

## Rename
col_name_map = list(
  CHR = 'Chromosome',
  START = 'Start_Position',
  END = 'End_Position',
  GENE = 'SYMBOL',
  REF = 'Reference_Allele',
  ALT = 'Allele',
  CLASS = 'Variant_Classification',
  IMPACT = 'IMPACT',
  TYPE = 'Variant_Type',
  BARCODE = 'Tumor_Sample_Barcode'
)
df = rename(df, !!! col_name_map)

## Make new colums
df$MUTATION = paste(
  df$CHR, 
  df$GENE,
  df$START,
  df$END,
  df$REF,
  df$ALT,
  sep = ':'
)
df$SAMPLE = str_sub(df$BARCODE, 1, 12)

## Remove duplicates
df = unique(df)

## Summary
df_summary = df %>%
  group_by(GENE) %>%
  summarise(
    N_MUTATIONS = n(),
    N_UNIQUE_MUTATIONS = length(unique(MUTATION)),
    SAMPLES = length(unique(BARCODE))
)
