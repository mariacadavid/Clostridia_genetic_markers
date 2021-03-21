#Retrieve accession numbers from GTDB for reference genome
#(nov-2020) based on genome_retrieval_tutorial.Rmb from Jacobo delaCuesta

library(tidyverse)

#Load tables
gtdb_bac  = read_tsv("https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_metadata_r95.tar.gz", na = "none")
gtdb_arc =  read_tsv("https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/ar122_metadata_r95.tar.gz", na = "none")

gtdb_bac2 = gtdb_bac %>% 
  rename("GB_Accession" = "bac120_metadata_r95.tsv") %>% 
  mutate(across(matches("ncbi_"), as.character))
gtdb_arc2 = gtdb_arc %>% 
  rename("GB_Accession" = "ar122_metadata_r95.tsv") %>% 
  mutate(across(matches("ncbi_"), as.character))
gtdb_df = bind_rows(gtdb_bac2, gtdb_arc2) %>% 
  filter(GB_Accession != "")

head(gtdb_df)

#Select genomes- filter with desired conditions:
#Belong to the class Ruminococcaceae
#Have high completeness (>= 95% completeness)
#Have low contamination (< 1%)
#Have fewer than 10 ssu copies (often a sign of misassembly)
selected_genomes = gtdb_df %>% 
  filter(str_detect(gtdb_taxonomy, "c__Clostridia"), 
         checkm_completeness >= 95,
         checkm_contamination < 1, 
         ssu_count < 10)

write_tsv(selected_genomes, 
          "/Users/mariacadavid/Desktop/Tesis_Clostridiales/genome_retrieval_GTDB/selected_clostridia_genomes.tsv", 
          col_names = TRUE)





