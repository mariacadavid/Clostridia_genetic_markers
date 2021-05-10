#Compare gen distances (dnak1 vs. gyrb vs. 16s)
#Intention: produce Pairwise distances plots (similar to Gou 2019 preprint graphs)
#Alignment needed 
#By: Maria Cadavid Feb2021


#libraries
library(seqinr) #read.alignment
library(ape) #dist.dna
library (spaa) #dist2list
library(tidyr)#drop_na
library(cowplot) #plot_grid
#plot density
library(MASS)
library(ggplot2)
library(viridis)
theme_set(theme_bw(base_size=12))
#Add r^2 to graph
library(ggpubr)
library(ggpmisc)



#DNAK1
gene_aln  <- read.alignment(file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Alineamientos/Alineamientos_genes_GTDB/genes_dnak1_GTDB_tblastn_aligned_MAFFT_rename_duplicates.fasta", format = "fasta") #Input alignment
#dist.dna() function requires the input alignment to be in a special format known as “DNAbin” format, so we must use the as.DNAbin() function
aln_bin <- as.DNAbin(gene_aln) # Convert the alignment to "DNAbin" format
dist_object <- dist.dna(aln_bin, model= "K80", as.matrix = FALSE, pairwise.deletion = TRUE) #Calculate the genetic distance in dist object
dist_list_dnak1 <- dist2list(dist_object)
colnames(dist_list_dnak1)[which(names(dist_list_dnak1) == "value")] <- "dist_dnak1"
#Create a new dataframe to merge with other genes
dist_list_dnak1$comparison_dnak1 <- paste(dist_list_dnak1$col, dist_list_dnak1$row,  sep="_compared_to_")
compare_dataframe_dnak1 <- dist_list_dnak1
compare_dataframe_dnak1$col <- NULL 
compare_dataframe_dnak1$row <- NULL 
write.table(compare_dataframe_dnak1, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/compare_dataframe_dnak1.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE)

#16S
gene_aln  <- read.alignment(file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Alineamientos/Alineamientos_genes_GTDB/genes_16s_GTDB_blastn_aligned_MAFFT_sin_comas.fasta", format = "fasta") #Input alignment
aln_bin <- as.DNAbin(gene_aln) # Convert the alignment to "DNAbin" format
dist_object <- dist.dna(aln_bin, model= "K80", as.matrix = FALSE, pairwise.deletion = TRUE) # Calculate the genetic distance in dist object
dist_list_16s <- dist2list(dist_object)
colnames(dist_list_16s)[which(names(dist_list_16s) == "value")] <- "dist_16s"
#Create a new dataframe to merge with other genes
dist_list_16s$comparison_16s <- paste(dist_list_16s$col, dist_list_16s$row,  sep="_compared_to_")
compare_dataframe_16s <- dist_list_16s
compare_dataframe_16s$col <- NULL 
compare_dataframe_16s$row <- NULL 
write.table(compare_dataframe_16s, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/compare_dataframe_16s_dos.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE)


#GYRB
gene_aln  <- read.alignment(file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Alineamientos/Alineamientos_genes_GTDB/genes_gyrb_GTDB_tblastn_aligned_MAFFT_cut1281_3955.fasta", format = "fasta") #Input alignment
aln_bin <- as.DNAbin(gene_aln) # Convert the alignment to "DNAbin" format
dist_object <- dist.dna(aln_bin, model= "K80", as.matrix = FALSE, pairwise.deletion = TRUE) #Calculate the genetic distance in dist object
dist_list_gyrb <- dist2list(dist_object)
colnames(dist_list_gyrb)[which(names(dist_list_gyrb) == "value")] <- "dist_gyrb"
#Create a new dataframe to merge with other genes
dist_list_gyrb$comparison_gyrb <- paste(dist_list_gyrb$col, dist_list_gyrb$row,  sep="_compared_to_")
compare_dataframe_gyrb <- dist_list_gyrb
compare_dataframe_gyrb$col <- NULL 
compare_dataframe_gyrb$row <- NULL 
write.table(compare_dataframe_gyrb, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/compare_dataframe_gyrb.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE)

#Important: Here we merge datasets first try (didnt load in R, so I did it in python with script "merge_df_distance.py")


#Read tables merged by python and calculate identity

#3 genes complete
distances <- read.csv(file='/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/distances_complete.csv', sep = "\t", header = TRUE)
common_pairwise <- drop_na(distances)

#16s vs dnak1
distances_16s_dnka1 <- read.csv(file='/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/distances_16s_dnak1.csv', sep = "\t", header = TRUE)
common_pairwise_16s_dnak1 <- drop_na(distances_16s_dnka1)
common_pairwise_16s_dnak1$simil_16s <- ((1 - common_pairwise_16s_dnak1$dist_16s)*100 )
common_pairwise_16s_dnak1$simil_dnak1 <- ((1 - common_pairwise_16s_dnak1$dist_dnak1)*100)
write.table(common_pairwise_16s_dnak1, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_dnak1_common.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE)

#16s vs gyrb
distances_16s_gyrb <- read.csv(file='/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/distances_16s_gyrb.csv', sep = "\t", header = TRUE)
common_pairwise_16s_gyrb <- drop_na(distances_16s_gyrb)
common_pairwise_16s_gyrb$simil_16s <- ((1 - common_pairwise_16s_gyrb$dist_16s)*100 )
common_pairwise_16s_gyrb$simil_gyrb <- ((1 - common_pairwise_16s_gyrb$dist_gyrb)*100)
write.table(common_pairwise_16s_gyrb, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_gyrb_common.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE)

#dnak1 vs gyrb
distances_dnak1_gyrb <- read.csv(file='/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/distances_dnak1_gyrb.csv', sep = "\t", header = TRUE)
common_pairwise_dnak1_gyrb <- drop_na(distances_dnak1_gyrb)
common_pairwise_dnak1_gyrb$simil_dnak1 <- ((1 - common_pairwise_dnak1_gyrb$dist_dnak1)*100 )
common_pairwise_dnak1_gyrb$simil_gyrb <- ((1 - common_pairwise_dnak1_gyrb$dist_gyrb)*100)
write.table(common_pairwise_dnak1_gyrb, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_dnak1_gyrb_common.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE)





####PLOTS 16s vs. DNAK1#####################################################
#ORGANIZE DATA FOR PLOTING
##Open tables of common pairwise comparisons with taxonomy analysis 
tax_analysis_16s_dnak1 <- as.data.frame(read.csv(file="/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_dnak1_common_taxonomy_analysis.tsv", sep = "\t", header = TRUE))

#filtrar "same genomes" para no incluirlas en la grafica
tax_analysis_16s_dnak1_true <-tax_analysis_16s_dnak1[tax_analysis_16s_dnak1$same_specie == 'true',]
tax_analysis_16s_dnak1_false <-tax_analysis_16s_dnak1[tax_analysis_16s_dnak1$same_specie == 'false',]
tax_analysis_16s_dnak1_no_diag <- rbind(tax_analysis_16s_dnak1_true, tax_analysis_16s_dnak1_false)


#Eliminate duplicate records 
#First separate contrast
tax_analysis_16s_dnak1_no_diag$contrast_char <- as.character(tax_analysis_16s_dnak1_no_diag$contrast) 
tax_analysis_16s_dnak1_no_diag$assembly1 <- sapply(strsplit(tax_analysis_16s_dnak1_no_diag$contrast_char, "_compared_to_"),`[`, 1) 
tax_analysis_16s_dnak1_no_diag$assembly2 <- sapply(strsplit(tax_analysis_16s_dnak1_no_diag$contrast_char, "_compared_to_"),`[`, 2)
#Export file with separated cols for assemblies 
write.table(tax_analysis_16s_dnak1_no_diag, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_dnak1_common_taxonomy_analysis_2.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)


#Bring back file with separated cols for assemblies to eliminate duplicates
conn_input <- file("/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_dnak1_common_taxonomy_analysis_2.tsv");
conn_output <- file("/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_dnak1_common_taxonomy_analysis_3.tsv", open="wt");
open(conn_input);
open(conn_output);
## escribir el encabezado al archivo de salida
line <- readLines(conn_input, n=1);
writeLines(line,conn_output, sep='\n')  ;
data <- list()
line <- readLines(conn_input, n=1);
while(length(line) > 0) {
  fields = strsplit(line,"\t");
  key = paste(fields[[1]][17],'-',fields[[1]][18], sep='');
  invkey = paste(fields[[1]][18],'-',fields[[1]][17], sep='');
  if (invkey %in% names(data)) {
    line <- readLines(conn_input, n=1);
    next;
  }
  data[[key]] <- line;
  writeLines(line,conn_output, sep='\n');
  line <- readLines(conn_input, n=1);
}
close(conn_input);
close(conn_output)



#PLOTS 1 (16s vs. dnak1):plot color by density: 
#Alternative option check here: https://slowkow.com/notes/ggplot2-color-by-density/

#Import the table with eliminated duplicates
tax_analysis_16s_dnak1 <- read.csv(file="/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_dnak1_common_taxonomy_analysis_3.tsv", sep = "\t", header = TRUE)

#Calculate density 
x <- densCols(tax_analysis_16s_dnak1$identity_16s, tax_analysis_16s_dnak1$identity_dnak1, colramp=colorRampPalette(c("black", "white")))
tax_analysis_16s_dnak1$dens <- col2rgb(x)[1,] + 1L

#similarity plot 16s vs dnak1
plot_16s_dnak1_similiarity <- ggplot(tax_analysis_16s_dnak1, aes(x=identity_16s, y=identity_dnak1, colour=dens)) +
  geom_point(aes(identity_16s, identity_dnak1), size= 0.5) +
  scale_color_viridis(option = "D") +
  labs(colour = "Number of pairwise \ncomparisons") +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) +
  xlab ("16S rRNA identity (%)") + ylab ("dnaK1 identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_dnak1_similiarity 
#Export as TIFF (1000 x 500) 


#PLOTS 2 (16s vs. dnak1): next we are going to plot coloring by the TAXONOMY ANALYSIS
##Organize tables to add col ("interest_fam") of f__Ruminococcaceae and f__Lachnospiraceae (to add shapes in graphs)
tax_analysis_16s_dnak1$interest_fam_1 <- as.factor(ifelse(tax_analysis_16s_dnak1$family_genome_1 =='f__Ruminococcaceae' ,"Ruminococcaceae", 
                                                          ifelse(tax_analysis_16s_dnak1$family_genome_1== 'f__Lachnospiraceae', "Lachnospiraceae","Other")))
tax_analysis_16s_dnak1$interest_fam_2 <- as.factor(ifelse(tax_analysis_16s_dnak1$family_genome_2 =='f__Ruminococcaceae' ,"Ruminococcaceae", 
                                                          ifelse(tax_analysis_16s_dnak1$family_genome_2== 'f__Lachnospiraceae', "Lachnospiraceae","Other")))
tax_analysis_16s_dnak1$interest_fam <- as.factor(ifelse(tax_analysis_16s_dnak1$interest_fam_1 == tax_analysis_16s_dnak1$interest_fam_2, as.character(tax_analysis_16s_dnak1$interest_fam_1), "Other")) 


#Preliminar graph example: similarity plot 16s vs dnak1 color taxonomy analysis 
plot_16s_dnak1_similiarity_tax_species <- ggplot(tax_analysis_16s_dnak1, aes(x=identity_16s, y=identity_dnak1)) + 
  geom_point(aes(identity_16s, identity_dnak1, color = same_specie, shape= interest_fam), size= 0.5) + 
  theme(legend.title = element_text(size = "10")) +
  theme(legend.position = "right", legend.text = element_text(size = 10, colour = "black")) +
  labs(title= "Pairwise identity 16s rRNA vs. dnaK1") +
  labs(colour = "Belong to same species") +
  xlab ("16s rRNA identity (%)") + ylab ("dnaK1 identity (%)") +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.8) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_dnak1_similiarity_tax_species


#Desviacion estandar 
library(dplyr)
sd16<-sd(tax_analysis_16s_dnak1$identity_16s)
sdkin<-sd(tax_analysis_16s_dnak1$identity_dnak1)

tax_analysis_16s_dnak1_same_species <-filter(tax_analysis_16s_dnak1, same_specie == "true")
sd16<-sd(tax_analysis_16s_dnak1_same_species$identity_16s)
sdkin<-sd(tax_analysis_16s_dnak1_same_species$identity_dnak1)
av16<- mean(tax_analysis_16s_dnak1_same_species$identity_16s)
avkin<- mean(tax_analysis_16s_dnak1_same_species$identity_dnak1)
range16<- range(tax_analysis_16s_dnak1_same_species$identity_16s)
rangekin<-range(tax_analysis_16s_dnak1_same_species$identity_dnak1)

tax_analysis_16s_dnak1_same_genus <-filter(tax_analysis_16s_dnak1, same_genus == "true")
sd16<-sd(tax_analysis_16s_dnak1_same_genus$identity_16s)
sdkin<-sd(tax_analysis_16s_dnak1_same_genus$identity_dnak1)

tax_analysis_16s_dnak1_same_family <-filter(tax_analysis_16s_dnak1, same_family == "true")
sd16<-sd(tax_analysis_16s_dnak1_same_family$identity_16s)
sdkin<-sd(tax_analysis_16s_dnak1_same_family$identity_dnak1)

tax_analysis_16s_dnak1_same_order <-filter(tax_analysis_16s_dnak1, same_order == "true")
sd16<-sd(tax_analysis_16s_dnak1_same_order$identity_16s)
sdkin<-sd(tax_analysis_16s_dnak1_same_order$identity_dnak1)

#Desviacion estandar families Ruminococcaceae and Lachnospiraceae
tax_analysis_16s_dnak1_rumino <-filter(tax_analysis_16s_dnak1, interest_fam == "Ruminococcaceae")
sd16<-sd(tax_analysis_16s_dnak1_rumino$identity_16s)
sdkin<-sd(tax_analysis_16s_dnak1_rumino$identity_dnak1)

tax_analysis_16s_dnak1_lachno <-filter(tax_analysis_16s_dnak1, interest_fam == "Lachnospiraceae")
sd16<-sd(tax_analysis_16s_dnak1_lachno$identity_16s)
sdkin<-sd(tax_analysis_16s_dnak1_lachno$identity_dnak1)


#Contruction of Panel with 4 graphs (4 tax levels)
library(gridExtra)
grupo.lab <- c("False", "True") 
cbp1 <- c( "#C5C5C5" , "#22CEC6")
# "#540EC8" MORADO
# "#C5C5C5" GRIS

#Pairwise identity 16s rRNA vs. dnaK1"
plot_16s_dnak1_similiarity_tax_species <- ggplot(tax_analysis_16s_dnak1, aes(x=identity_16s, y=identity_dnak1)) + 
  geom_point(aes(identity_16s, identity_dnak1, color = same_specie), size= 1) +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))+
  labs(colour = "Belong to same species") +
  theme(axis.title = element_text(size = "15")) +
  scale_colour_manual(values=cbp1, labels = grupo.lab) +
  xlab ("16S rRNA identity (%)") + ylab ("dnaK1 identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_dnak1_similiarity_tax_species

plot_16s_dnak1_similiarity_tax_genus <- ggplot(tax_analysis_16s_dnak1, aes(x=identity_16s, y=identity_dnak1)) + 
  geom_point(aes(identity_16s, identity_dnak1, color = same_genus), size= 1) +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))+
  labs(colour = "Belong to same genus") +
  theme(axis.title = element_text(size = "15")) +
  scale_colour_manual(values=cbp1, labels = grupo.lab) +
  xlab ("16S rRNA identity (%)") + ylab ("dnaK1 identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_dnak1_similiarity_tax_genus

plot_16s_dnak1_similiarity_tax_family <- ggplot(tax_analysis_16s_dnak1, aes(x=identity_16s, y=identity_dnak1)) + 
  geom_point(aes(identity_16s, identity_dnak1, color = same_family), size= 1) +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))+
  labs(colour = "Belong to same family") +
  theme(axis.title = element_text(size = "15")) +
  scale_colour_manual(values=cbp1, labels = grupo.lab) +
  xlab ("16S rRNA identity (%)") + ylab ("dnaK1 identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_dnak1_similiarity_tax_family

plot_16s_dnak1_similiarity_tax_order <- ggplot(tax_analysis_16s_dnak1, aes(x=identity_16s, y=identity_dnak1)) + 
  geom_point(aes(identity_16s, identity_dnak1, color = same_order), size= 1) +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))+
  labs(colour = "Belong to same order") +
  theme(axis.title = element_text(size = "15")) +
  scale_colour_manual(values=cbp1, labels = grupo.lab) +
  xlab ("16S rRNA identity (%)") + ylab ("dnaK1 identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_dnak1_similiarity_tax_order

plot_grid(
  nrow= 2, ncol=2,
  plot_16s_dnak1_similiarity_tax_species, plot_16s_dnak1_similiarity_tax_genus, plot_16s_dnak1_similiarity_tax_family, plot_16s_dnak1_similiarity_tax_order, 
  labels = c('A', 'B', 'C', 'D'),
  align="hv"
)
#export as TIFF (2000x1000)

#Color interest families
tax_analysis_16s_dnak1$interest_fam <- factor(tax_analysis_16s_dnak1$interest_fam, levels = c("Ruminococcaceae", "Lachnospiraceae", "Other"))
cbp2 <- c( "#540EC8", "#B2E65F" , "#C5C5C5" )
plot_16s_dnak1_similiarity_interest_fam <- ggplot(tax_analysis_16s_dnak1, aes(x=identity_16s, y=identity_dnak1)) + 
  geom_point(aes(identity_16s, identity_dnak1, color = interest_fam), size= 1) +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) + 
  guides(color= guide_legend(override.aes = list(size=3)))+
  labs(colour= "Families of interest") +
  theme(axis.title = element_text(size = "15")) +
  scale_colour_manual(values=cbp2) +
  xlab ("16S rRNA identity (%)") + ylab ("dnaK1 identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_dnak1_similiarity_interest_fam
#export as TIFF (1000x500)



#####PLOTS 16s vs. GYRB#####################################################
#ORGANIZE DATA FOR PLOTING
##Open tables of common pairwise comparisons with taxonomy analysis 
tax_analysis_16s_gyrb <- as.data.frame(read.csv(file="/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_gyrb_common_taxonomy_analysis.tsv", sep = "\t", header = TRUE))


#filtrar "same genomes" para no incluirlas en la grafica
tax_analysis_16s_gyrb_true <-tax_analysis_16s_gyrb[tax_analysis_16s_gyrb$same_specie == 'true',]
tax_analysis_16s_gyrb_false <-tax_analysis_16s_gyrb[tax_analysis_16s_gyrb$same_specie == 'false',]
tax_analysis_16s_gyrb_no_diag <- rbind(tax_analysis_16s_gyrb_true, tax_analysis_16s_gyrb_false)


#Eliminate duplicate records 
#First separate contrast
tax_analysis_16s_gyrb_no_diag$contrast_char <- as.character(tax_analysis_16s_gyrb_no_diag$contrast) 
tax_analysis_16s_gyrb_no_diag$assembly1 <- sapply(strsplit(tax_analysis_16s_gyrb_no_diag$contrast_char, "_compared_to_"),`[`, 1) 
tax_analysis_16s_gyrb_no_diag$assembly2 <- sapply(strsplit(tax_analysis_16s_gyrb_no_diag$contrast_char, "_compared_to_"),`[`, 2)
#Export file with separated cols for assemblies 
write.table(tax_analysis_16s_gyrb_no_diag, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_gyrb_common_taxonomy_analysis_2.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)


#Bring back file with separated cols for assemblies to eliminate duplicates
conn_input <- file("/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_gyrb_common_taxonomy_analysis_2.tsv");
conn_output <- file("/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_gyrb_common_taxonomy_analysis_3.tsv", open="wt");
open(conn_input);
open(conn_output);
## escribir el encabezado al archivo de salida
line <- readLines(conn_input, n=1);
writeLines(line,conn_output, sep='\n')  ;
data <- list()
line <- readLines(conn_input, n=1);
while(length(line) > 0) {
  fields = strsplit(line,"\t");
  key = paste(fields[[1]][17],'-',fields[[1]][18], sep='');
  invkey = paste(fields[[1]][18],'-',fields[[1]][17], sep='');
  if (invkey %in% names(data)) {
    line <- readLines(conn_input, n=1);
    next;
  }
  data[[key]] <- line;
  writeLines(line,conn_output, sep='\n');
  line <- readLines(conn_input, n=1);
}
close(conn_input);
close(conn_output)




#PLOTS 1 (16s vs. gyrb):plot color by density: 
#Alternative option check here: https://slowkow.com/notes/ggplot2-color-by-density/

#Import the table with eliminated duplicates
tax_analysis_16s_gyrb <- read.csv(file="/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_gyrb_common_taxonomy_analysis_3.tsv", sep = "\t", header = TRUE)

#Calculate density 
x <- densCols(tax_analysis_16s_gyrb$identity_16s, tax_analysis_16s_gyrb$identity_gyrb, colramp=colorRampPalette(c("black", "white")))
tax_analysis_16s_gyrb$dens <- col2rgb(x)[1,] + 1L

#similarity plot 16s vs gyrb
plot_16s_gyrb_similiarity <- ggplot(tax_analysis_16s_gyrb, aes(x=identity_16s, y=identity_gyrb, colour=dens)) +
  geom_point(aes(identity_16s, identity_gyrb), size= 0.5) +
  scale_color_viridis(option = "D") +
  labs(colour = "Number of pairwise \ncomparisons") +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) +
  xlab ("16S rRNA identity (%)") + ylab ("gyrB identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_gyrb_similiarity 
#Export as TIFF (1000 x 500) 


#PLOTS 2 (16s vs. gyrb): next we are going to plot coloring by the TAXONOMY ANALYSIS
##Organize tables to add col ("interest_fam") of f__Ruminococcaceae and f__Lachnospiraceae (to add shapes in graphs)
tax_analysis_16s_gyrb$interest_fam_1 <- as.factor(ifelse(tax_analysis_16s_gyrb$family_genome_1 =='f__Ruminococcaceae' ,"Ruminococcaceae", 
                                                          ifelse(tax_analysis_16s_gyrb$family_genome_1== 'f__Lachnospiraceae', "Lachnospiraceae","Other")))
tax_analysis_16s_gyrb$interest_fam_2 <- as.factor(ifelse(tax_analysis_16s_gyrb$family_genome_2 =='f__Ruminococcaceae' ,"Ruminococcaceae", 
                                                          ifelse(tax_analysis_16s_gyrb$family_genome_2== 'f__Lachnospiraceae', "Lachnospiraceae","Other")))
tax_analysis_16s_gyrb$interest_fam <- as.factor(ifelse(tax_analysis_16s_gyrb$interest_fam_1 == tax_analysis_16s_gyrb$interest_fam_2, as.character(tax_analysis_16s_gyrb$interest_fam_1), "Other")) 


#Preliminar graph example: similarity plot 16s vs dnak1 color taxonomy analysis 
plot_16s_gyrb_similiarity_tax_species <- ggplot(tax_analysis_16s_gyrb, aes(x=identity_16s, y=identity_gyrb)) + 
  geom_point(aes(identity_16s, identity_gyrb, color = same_specie, shape= interest_fam), size= 0.5) + 
  theme(legend.title = element_text(size = "10")) +
  theme(legend.position = "right", legend.text = element_text(size = 10, colour = "black")) +
  labs(title= "Pairwise identity 16s rRNA vs. gyrB") +
  labs(colour = "Belong to same species") +
  xlab ("16S rRNA identity (%)") + ylab ("gyrB identity (%)") +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.8) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_gyrb_similiarity_tax_species


#Desviacion estandar 
library(dplyr)
sd16<-sd(tax_analysis_16s_gyrb$identity_16s)
sdgyr<-sd(tax_analysis_16s_gyrb$identity_gyrb)

tax_analysis_16s_gyrb_same_species <-filter(tax_analysis_16s_gyrb, same_specie == "true")
sd16<-sd(tax_analysis_16s_gyrb_same_species$identity_16s)
sdgyr<-sd(tax_analysis_16s_gyrb_same_species$identity_gyrb)
av16<- mean(tax_analysis_16s_gyrb_same_species$identity_16s)
avgyr<- mean(tax_analysis_16s_gyrb_same_species$identity_gyrb)
range16<- range(tax_analysis_16s_gyrb_same_species$identity_16s)
rangegyr<-range(tax_analysis_16s_gyrb_same_species$identity_gyrb)

tax_analysis_16s_gyrb_same_genus <-filter(tax_analysis_16s_gyrb, same_genus == "true")
sd16<-sd(tax_analysis_16s_gyrb_same_genus$identity_16s)
sdgyr<-sd(tax_analysis_16s_gyrb_same_genus$identity_gyrb)

tax_analysis_16s_gyrb_same_family <-filter(tax_analysis_16s_gyrb, same_family == "true")
sd16<-sd(tax_analysis_16s_gyrb_same_family$identity_16s)
sdgyr<-sd(tax_analysis_16s_gyrb_same_family$identity_gyrb)

tax_analysis_16s_gyrb_same_order <-filter(tax_analysis_16s_gyrb, same_order == "true")
sd16<-sd(tax_analysis_16s_gyrb_same_order$identity_16s)
sdgyr<-sd(tax_analysis_16s_gyrb_same_order$identity_gyrb)

#Desviacion estandar families Ruminococcaceae and Lachnospiraceae
tax_analysis_16s_gyrb_rumino <-filter(tax_analysis_16s_gyrb, interest_fam == "Ruminococcaceae")
sd16<-sd(tax_analysis_16s_gyrb_rumino$identity_16s)
sdgyr<-sd(tax_analysis_16s_gyrb_rumino$identity_gyrb)

tax_analysis_16s_gyrb_lachno <-filter(tax_analysis_16s_gyrb, interest_fam == "Lachnospiraceae")
sd16<-sd(tax_analysis_16s_gyrb_lachno$identity_16s)
sdgyr<-sd(tax_analysis_16s_gyrb_lachno$identity_gyrb)


#Contruction of Panel with 4 graphs (4 tax levels)
library(gridExtra)
grupo.lab <- c("False", "True") 
cbp1 <- c( "#C5C5C5" , "#22CEC6")
# "#540EC8" MORADO
# "#C5C5C5" GRIS

#Pairwise identity 16s rRNA vs. gyrb"
plot_16s_gyrb_similiarity_tax_species <- ggplot(tax_analysis_16s_gyrb, aes(x=identity_16s, y=identity_gyrb)) + 
  geom_point(aes(identity_16s, identity_gyrb, color = same_specie), size= 1) +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))+
  labs(colour = "Belong to same species") +
  theme(axis.title = element_text(size = "15")) +
  scale_colour_manual(values=cbp1, labels = grupo.lab) +
  xlab ("16S rRNA identity (%)") + ylab ("gyrB identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_gyrb_similiarity_tax_species

plot_16s_gyrb_similiarity_tax_genus <- ggplot(tax_analysis_16s_gyrb, aes(x=identity_16s, y=identity_gyrb)) + 
  geom_point(aes(identity_16s, identity_gyrb, color= same_genus), size= 1) +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))+
  labs(colour = "Belong to same genus") +
  theme(axis.title = element_text(size = "15")) +
  scale_colour_manual(values=cbp1, labels = grupo.lab) +
  xlab ("16S rRNA identity (%)") + ylab ("gyrB identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_gyrb_similiarity_tax_genus

plot_16s_gyrb_similiarity_tax_family <- ggplot(tax_analysis_16s_gyrb, aes(x=identity_16s, y=identity_gyrb)) + 
  geom_point(aes(identity_16s, identity_gyrb, color = same_family), size= 1) +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))+
  labs(colour = "Belong to same family") +
  theme(axis.title = element_text(size = "15")) +
  scale_colour_manual(values=cbp1, labels = grupo.lab) +
  xlab ("16S rRNA identity (%)") + ylab ("gyrB identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_gyrb_similiarity_tax_family

plot_16s_gyrb_similiarity_tax_order <- ggplot(tax_analysis_16s_gyrb, aes(x=identity_16s, y=identity_gyrb)) + 
  geom_point(aes(identity_16s, identity_gyrb, color = same_order), size= 1) +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))+
  labs(colour = "Belong to same order") +
  theme(axis.title = element_text(size = "15")) +
  scale_colour_manual(values=cbp1, labels = grupo.lab) +
  xlab ("16S rRNA identity (%)") + ylab ("gyrB identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_gyrb_similiarity_tax_order


plot_grid(
  nrow= 2, ncol=2,
  plot_16s_gyrb_similiarity_tax_species, plot_16s_gyrb_similiarity_tax_genus, plot_16s_gyrb_similiarity_tax_family, plot_16s_gyrb_similiarity_tax_order, 
  labels = c('A', 'B', 'C', 'D'),
  align="hv"
)
#export as PNG (2000x1000)
#export as TIFF (2000x1000)

#Color interest families
tax_analysis_16s_gyrb$interest_fam <- factor(tax_analysis_16s_gyrb$interest_fam, levels = c("Ruminococcaceae", "Lachnospiraceae", "Other"))
cbp2 <- c( "#540EC8", "#B2E65F" , "#C5C5C5" )
plot_16s_gyrb_similiarity_interest_fam <- ggplot(tax_analysis_16s_gyrb, aes(x=identity_16s, y=identity_gyrb)) + 
  geom_point(aes(identity_16s, identity_gyrb, color = interest_fam), size= 1) +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) + 
  guides(color= guide_legend(override.aes = list(size=3)))+
  labs(colour= "Families of interest") +
  theme(axis.title = element_text(size = "15")) +
  scale_colour_manual(values=cbp2) +
  xlab ("16S rRNA identity (%)") + ylab ("gyrB identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_gyrb_similiarity_interest_fam
#export as TIFF (1000x500)





#####PLOTS DNAK1 vs. GYRB#####################################################
#ORGANIZE DATA FOR PLOTING
##Open tables of common pairwise comparisons with taxonomy analysis 
tax_analysis_dnak1_gyrb <- as.data.frame(read.csv(file="/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_dnak1_gyrb_common_taxonomy_analysis.tsv", sep = "\t", header = TRUE))


#filtrar "same genomes" para no incluirlas en la grafica
tax_analysis_dnak1_gyrb_true <-tax_analysis_dnak1_gyrb[tax_analysis_dnak1_gyrb$same_specie == 'true',]
tax_analysis_dnak1_gyrb_false <-tax_analysis_dnak1_gyrb[tax_analysis_dnak1_gyrb$same_specie == 'false',]
tax_analysis_dnak1_gyrb_no_diag <- rbind(tax_analysis_dnak1_gyrb_true, tax_analysis_dnak1_gyrb_false)


#Eliminate duplicate records 
#First separate contrast
tax_analysis_dnak1_gyrb_no_diag$contrast_char <- as.character(tax_analysis_dnak1_gyrb_no_diag$contrast) 
tax_analysis_dnak1_gyrb_no_diag$assembly1 <- sapply(strsplit(tax_analysis_dnak1_gyrb_no_diag$contrast_char, "_compared_to_"),`[`, 1) 
tax_analysis_dnak1_gyrb_no_diag$assembly2 <- sapply(strsplit(tax_analysis_dnak1_gyrb_no_diag$contrast_char, "_compared_to_"),`[`, 2)
#Export file with separated cols for assemblies 
write.table(tax_analysis_dnak1_gyrb_no_diag, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_dnak1_gyrb_common_taxonomy_analysis_2.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)


#Bring back file with separated cols for assemblies to eliminate duplicates
conn_input <- file("/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_dnak1_gyrb_common_taxonomy_analysis_2.tsv");
conn_output <- file("/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_dnak1_gyrb_common_taxonomy_analysis_3.tsv", open="wt");
open(conn_input);
open(conn_output);
## escribir el encabezado al archivo de salida
line <- readLines(conn_input, n=1);
writeLines(line,conn_output, sep='\n')  ;
data <- list()
line <- readLines(conn_input, n=1);
while(length(line) > 0) {
  fields = strsplit(line,"\t");
  key = paste(fields[[1]][17],'-',fields[[1]][18], sep='');
  invkey = paste(fields[[1]][18],'-',fields[[1]][17], sep='');
  if (invkey %in% names(data)) {
    line <- readLines(conn_input, n=1);
    next;
  }
  data[[key]] <- line;
  writeLines(line,conn_output, sep='\n');
  line <- readLines(conn_input, n=1);
}
close(conn_input);
close(conn_output)




#PLOTS 1 (dnak1 vs. gyrb):plot color by density: 
#Alternative option check here: https://slowkow.com/notes/ggplot2-color-by-density/

#Import the table with eliminated duplicates
tax_analysis_dnak1_gyrb <- read.csv(file="/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_dnak1_gyrb_common_taxonomy_analysis_3.tsv", sep = "\t", header = TRUE)

#Calculate density 
x <- densCols(tax_analysis_dnak1_gyrb$identity_dnak1, tax_analysis_dnak1_gyrb$identity_gyrb, colramp=colorRampPalette(c("black", "white")))
tax_analysis_dnak1_gyrb$dens <- col2rgb(x)[1,] + 1L

#similarity plot 16s vs gyrb
plot_dnak1_gyrb_similiarity <- ggplot(tax_analysis_dnak1_gyrb, aes(x=identity_dnak1, y=identity_gyrb, colour=dens)) +
  geom_point(aes(identity_dnak1, identity_gyrb), size= 0.5) +
  scale_color_viridis(option = "D") +
  labs(colour = "Number of pairwise \ncomparisons") +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) +
  xlab ("dnaK1 identity (%)") + ylab ("gyrB identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(0,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_dnak1_gyrb_similiarity 

#Export as TIFF (1000 x 500) 


#PLOTS 2 (dnak1 vs. gyrb): next we are going to plot coloring by the TAXONOMY ANALYSIS
##Organize tables to add col ("interest_fam") of f__Ruminococcaceae and f__Lachnospiraceae (to add shapes in graphs)
tax_analysis_dnak1_gyrb$interest_fam_1 <- as.factor(ifelse(tax_analysis_dnak1_gyrb$family_genome_1 =='f__Ruminococcaceae' ,"Ruminococcaceae", 
                                                         ifelse(tax_analysis_dnak1_gyrb$family_genome_1== 'f__Lachnospiraceae', "Lachnospiraceae","Other")))
tax_analysis_dnak1_gyrb$interest_fam_2 <- as.factor(ifelse(tax_analysis_dnak1_gyrb$family_genome_2 =='f__Ruminococcaceae' ,"Ruminococcaceae", 
                                                         ifelse(tax_analysis_dnak1_gyrb$family_genome_2== 'f__Lachnospiraceae', "Lachnospiraceae","Other")))
tax_analysis_dnak1_gyrb$interest_fam <- as.factor(ifelse(tax_analysis_dnak1_gyrb$interest_fam_1 == tax_analysis_dnak1_gyrb$interest_fam_2, as.character(tax_analysis_dnak1_gyrb$interest_fam_1), "Other")) 


#Preliminar graph example: similarity plot dnak1 vs dnak1 color taxonomy analysis 
plot_dnak1_gyrb_similiarity_tax_species <- ggplot(tax_analysis_dnak1_gyrb, aes(x=identity_dnak1, y=identity_gyrb)) + 
  geom_point(aes(identity_dnak1, identity_gyrb, color = same_specie, shape= interest_fam), size= 0.5) + 
  theme(legend.title = element_text(size = "10")) +
  theme(legend.position = "right", legend.text = element_text(size = 10, colour = "black")) +
  labs(title= "Pairwise identity dnaK1 vs. gyrB") +
  labs(colour = "Belong to same species") +
  xlab ("dnaK1 identity (%)") + ylab ("gyrB identity (%)") +
  xlim(0,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.8) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_dnak1_gyrb_similiarity_tax_species


#Desviacion estandar 
library(dplyr)
sdkin<-sd(tax_analysis_dnak1_gyrb$identity_dnak1)
sdgyr<-sd(tax_analysis_dnak1_gyrb$identity_gyrb)

tax_analysis_dnak1_gyrb_same_species <-filter(tax_analysis_dnak1_gyrb, same_specie == "true")
sdkin<-sd(tax_analysis_dnak1_gyrb_same_species$identity_dnak1)
sdgyr<-sd(tax_analysis_dnak1_gyrb_same_species$identity_gyrb)
avkin<- mean(tax_analysis_dnak1_gyrb_same_species$identity_dnak1)
avgyr<- mean(tax_analysis_dnak1_gyrb_same_species$identity_gyrb)
rangekin<- range(tax_analysis_dnak1_gyrb_same_species$identity_dnak1)
rangegyr<-range(tax_analysis_dnak1_gyrb_same_species$identity_gyrb)

tax_analysis_dnak1_gyrb_same_genus <-filter(tax_analysis_dnak1_gyrb, same_genus == "true")
sdkin<-sd(tax_analysis_dnak1_gyrb_same_genus$identity_dnak1)
sdgyr<-sd(tax_analysis_dnak1_gyrb_same_genus$identity_gyrb)

tax_analysis_dnak1_gyrb_same_family <-filter(tax_analysis_dnak1_gyrb, same_family == "true")
sdkin<-sd(tax_analysis_dnak1_gyrb_same_family$identity_dnak1)
sdgyr<-sd(tax_analysis_dnak1_gyrb_same_family$identity_gyrb)

tax_analysis_dnak1_gyrb_same_order <-filter(tax_analysis_dnak1_gyrb, same_order == "true")
sdkin<-sd(tax_analysis_dnak1_gyrb_same_order$identity_dnak1)
sdgyr<-sd(tax_analysis_dnak1_gyrb_same_order$identity_gyrb)

#Desviacion estandar families Ruminococcaceae and Lachnospiraceae
tax_analysis_dnak1_gyrb_rumino <-filter(tax_analysis_dnak1_gyrb, interest_fam == "Ruminococcaceae")
sdkin<-sd(tax_analysis_dnak1_gyrb_rumino$identity_dnak1)
sdgyr<-sd(tax_analysis_dnak1_gyrb_rumino$identity_gyrb)

tax_analysis_dnak1_gyrb_lachno <-filter(tax_analysis_dnak1_gyrb, interest_fam == "Lachnospiraceae")
sdkin<-sd(tax_analysis_dnak1_gyrb_lachno$identity_dnak1)
sdgyr<-sd(tax_analysis_dnak1_gyrb_lachno$identity_gyrb)


#Contruction of Panel with 4 graphs (4 tax levels)
library(gridExtra)
grupo.lab <- c("False", "True") 
cbp1 <- c( "#C5C5C5" , "#22CEC6")
# "#540EC8" MORADO
# "#C5C5C5" GRIS

#Pairwise identity dnaK1 vs. gyrb"
plot_dnak1_gyrb_similiarity_tax_species <- ggplot(tax_analysis_dnak1_gyrb, aes(x=identity_dnak1, y=identity_gyrb)) + 
  geom_point(aes(identity_dnak1, identity_gyrb, color = same_specie), size= 1) +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))+
  labs(colour = "Belong to same species") +
  theme(axis.title = element_text(size = "15")) +
  scale_colour_manual(values=cbp1, labels = grupo.lab) +
  xlab ("dnaK1 identity (%)") + ylab ("gyrB identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(0,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_dnak1_gyrb_similiarity_tax_species

plot_dnak1_gyrb_similiarity_tax_genus <- ggplot(tax_analysis_dnak1_gyrb, aes(x=identity_dnak1, y=identity_gyrb)) + 
  geom_point(aes(identity_dnak1, identity_gyrb, color= same_genus), size= 1) +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))+
  labs(colour = "Belong to same genus") +
  theme(axis.title = element_text(size = "15")) +
  scale_colour_manual(values=cbp1, labels = grupo.lab) +
  xlab ("dnaK1 identity (%)") + ylab ("gyrB identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(0,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_dnak1_gyrb_similiarity_tax_genus

plot_dnak1_gyrb_similiarity_tax_family <- ggplot(tax_analysis_dnak1_gyrb, aes(x=identity_dnak1, y=identity_gyrb)) + 
  geom_point(aes(identity_dnak1, identity_gyrb, color = same_family), size= 1) +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))+
  labs(colour = "Belong to same family") +
  theme(axis.title = element_text(size = "15")) +
  scale_colour_manual(values=cbp1, labels = grupo.lab) +
  xlab ("dnaK1 identity (%)") + ylab ("gyrB identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(0,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_dnak1_gyrb_similiarity_tax_family

plot_dnak1_gyrb_similiarity_tax_order <- ggplot(tax_analysis_dnak1_gyrb, aes(x=identity_dnak1, y=identity_gyrb)) + 
  geom_point(aes(identity_dnak1, identity_gyrb, color = same_order), size= 1) +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) +
  guides(color = guide_legend(override.aes = list(size=3)))+
  labs(colour = "Belong to same order") +
  theme(axis.title = element_text(size = "15")) +
  scale_colour_manual(values=cbp1, labels = grupo.lab) +
  xlab ("dnaK1 identity (%)") + ylab ("gyrB identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(0,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_dnak1_gyrb_similiarity_tax_order


plot_grid(
  nrow= 2, ncol=2,
  plot_dnak1_gyrb_similiarity_tax_species, plot_dnak1_gyrb_similiarity_tax_genus, plot_dnak1_gyrb_similiarity_tax_family, plot_dnak1_gyrb_similiarity_tax_order, 
  labels = c('A', 'B', 'C', 'D'),
  align="hv"
)
#export as PNG (2000x1000)
#export as TIFF (2000x1000)


#Color interest families
tax_analysis_dnak1_gyrb$interest_fam <- factor(tax_analysis_dnak1_gyrb$interest_fam, levels = c("Ruminococcaceae", "Lachnospiraceae", "Other"))
cbp2 <- c( "#540EC8", "#B2E65F" , "#C5C5C5" )
plot_dnak1_gyrb_similiarity_interest_fam <- ggplot(tax_analysis_dnak1_gyrb, aes(x=identity_dnak1, y=identity_gyrb)) + 
  geom_point(aes(identity_dnak1, identity_gyrb, color = interest_fam), size= 1) +
  theme(legend.title = element_text(size = "15")) +
  theme(legend.position = "right", legend.text = element_text(size = 15, colour = "black")) + 
  guides(color= guide_legend(override.aes = list(size=3)))+
  labs(colour= "Families of interest") +
  theme(axis.title = element_text(size = "15")) +
  scale_colour_manual(values=cbp2) +
  xlab ("dnaK1 identity (%)") + ylab ("gyrB identity (%)") +
  theme(axis.text = element_text(size = "15")) +
  xlim(0,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.5) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..), size = 5) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..), size = 5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_dnak1_gyrb_similiarity_interest_fam
#export as TIFF (1000x500)

####Panel: 3 density graphs####
plot_grid(
  nrow= 3,
  plot_16s_dnak1_similiarity, plot_16s_gyrb_similiarity, plot_dnak1_gyrb_similiarity,
  labels = c('A', 'B', 'C'),
  align="hv"
)
#export as PNG and TIFF (2000x600) HORIZONTAL PANEL
#export as PNG and TIFF (1000x1000) VERTICAL PANEL

####Panel: 3 interes families graphs####
plot_grid(
  nrow= 1,
  plot_16s_dnak1_similiarity_interest_fam, plot_16s_gyrb_similiarity_interest_fam, plot_dnak1_gyrb_similiarity_interest_fam,
  labels = c('A', 'B', 'C'),
  align="hv"
)

#export as PNG and TIFF (2000x600) HORIZAONTAL PANEL
#export as PNG (1000x1000) VERTICAL PANEL 

