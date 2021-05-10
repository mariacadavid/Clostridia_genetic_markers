#Obtain pairwise distance matrices and frequency histograms (dnak1, gyrb, 16s)
#Alignment needed 
#By: Maria Cadavid Jan2021
#Bassed on: https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter5.html

#libraries
library(seqinr) #read.alignment
library(ape) #dist.dna
library (spaa) #dist2list


####Distance matrix for dnak1####
gene_aln  <- read.alignment(file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Alineamientos/Alineamientos_genes_GTDB/genes_dnak1_GTDB_tblastn_aligned_MAFFT_rename_duplicates.fasta", format = "fasta") #Input alignment
#dist.dna() function requires the input alignment to be in a special format known as “DNAbin” format, so we must use the as.DNAbin() function
aln_bin <- as.DNAbin(gene_aln) # Convert the alignment to "DNAbin" format
dist_matrix <- dist.dna(aln_bin, model= "K80", as.matrix = TRUE, pairwise.deletion = TRUE) # Calculate the genetic distance matrix
write.table(dist_matrix, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/K80/indentity_matrix_ape_dnak1.tsv", sep = "\t",
            row.names = TRUE, col.names = TRUE)

####Frequency histogram dnak1####
dist_object <- dist.dna(aln_bin, model= "K80", as.matrix = FALSE, pairwise.deletion = TRUE) # Calculate the genetic distance in dist object
dist_list_dnak1 <- dist2list(dist_object)
dist_frame <- as.data.frame(dist_list_dnak1)
colnames(dist_frame)[which(names(dist_frame) == "col")] <- "org1"
colnames(dist_frame)[which(names(dist_frame) == "row")] <- "org2"
colnames(dist_frame)[which(names(dist_frame) == "value")] <- "distances"
dist_frame_dnak1 <- dist_frame
freq_dnak1=hist(dist_frame_dnak1$distances, breaks= 100, col="#d8ffc5", xlab="Pairwise distance (dnaK1)", main="", xlim=c(0,1))
tiff(file="/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/K80/frequency_distances_ape_dnak1.tiff", units="cm" ,width=22, height=13, res=600)
hist(dist_frame_dnak1$distances, breaks= 100, col="#d8ffc5", xlab="Pairwise distance (dnaK1)", main="", xlim=c(0,1))
dev.off()



####Distance matrix for gyrb cut 1281_3955####
gene_aln  <- read.alignment(file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Alineamientos/Alineamientos_genes_GTDB/genes_gyrb_GTDB_tblastn_aligned_MAFFT_cut1281_3955.fasta", format = "fasta") #Input alignment
#dist.dna() function requires the input alignment to be in a special format known as “DNAbin” format, so we must use the as.DNAbin() function
aln_bin <- as.DNAbin(gene_aln) # Convert the alignment to "DNAbin" format
dist_matrix <- dist.dna(aln_bin, model= "K80", as.matrix = TRUE, pairwise.deletion = TRUE) # Calculate the genetic distance matrix
write.table(dist_matrix, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/K80/indentity_matrix_ape_gyrb_cut.tsv", sep = "\t",
            row.names = TRUE, col.names = TRUE)

####Frequency histogram gyrb####
dist_object <- dist.dna(aln_bin, model= "K80", as.matrix = FALSE, pairwise.deletion = TRUE) # Calculate the genetic distance in dist object
dist_list_gyrb <- dist2list(dist_object)
dist_frame <- as.data.frame(dist_list_gyrb)
colnames(dist_frame)[which(names(dist_frame) == "col")] <- "org1"
colnames(dist_frame)[which(names(dist_frame) == "row")] <- "org2"
colnames(dist_frame)[which(names(dist_frame) == "value")] <- "distances"
dist_frame_gyrb <- dist_frame
freq_gyrb=hist(dist_frame_gyrb$distances, breaks= 100, col="#d8ffc5", xlab="Pairwise distance (gyrB)", main="", xlim=c(0,1))
tiff(file="/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/K80/frequency_distances_ape_gyrb.tiff", units="cm" ,width=22, height=13, res=600)
hist(dist_frame_gyrb$distances, breaks= 100, col="#d8ffc5", xlab="Pairwise distance (gyrB)", main="", xlim=c(0,1))
dev.off()



####Distance matrix for 16s####
gene_aln  <- read.alignment(file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Alineamientos/Alineamientos_genes_GTDB/genes_16s_GTDB_blastn_aligned_MAFFT_sin_comas.fasta", format = "fasta") #Input alignment
#dist.dna() function requires the input alignment to be in a special format known as “DNAbin” format, so we must use the as.DNAbin() function
aln_bin <- as.DNAbin(gene_aln) # Convert the alignment to "DNAbin" format
dist_matrix<- dist.dna(aln_bin, model= "K80", as.matrix = TRUE, pairwise.deletion = TRUE) # Calculate the genetic distance matrix
write.table(dist_matrix, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/K80/indentity_matrix_ape_16s.tsv", sep = "\t",
            row.names = TRUE, col.names = TRUE)

####Frequency histogram 16s####
dist_object <- dist.dna(aln_bin, model= "K80", as.matrix = FALSE, pairwise.deletion = TRUE) # Calculate the genetic distance in dist object
dist_list_16s <- dist2list(dist_object)
dist_frame <- as.data.frame(dist_list_16s)
colnames(dist_frame)[which(names(dist_frame) == "col")] <- "org1"
colnames(dist_frame)[which(names(dist_frame) == "row")] <- "org2"
colnames(dist_frame)[which(names(dist_frame) == "value")] <- "distances"
dist_frame_16s <- dist_frame
freq_16s=hist(dist_frame_16s$distances, breaks= 100, col="#d8ffc5", xlab="Pairwise distance (16S rRNA)", main="", xlim=c(0,1))
tiff(file="/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/K80/frequency_distances_ape_16s.tiff", units="cm" ,width=22, height=13, res=600)
hist(dist_frame_16s$distances, breaks= 100, col="#d8ffc5", xlab="Pairwise distance (16S rRNA)", main="", xlim=c(0,1))
dev.off()


####Panel all graphs####

par(mfrow=c(3,1))
hist(dist_frame_16s$distances, breaks=   50, col="#d8ffc5", xlab="Pairwise distance (16S rRNA)", main="", xlim=c(0,1), cex.lab=1.5, cex.axis=1.5)
hist(dist_frame_dnak1$distances, breaks= 50, col="#d8ffc5", xlab="Pairwise distance (dnaK1)", main="", xlim=c(0,1), cex.lab=1.5, cex.axis=1.5)
hist(dist_frame_gyrb$distances, breaks= 100, col="#d8ffc5", xlab="Pairwise distance (gyrB)", main="", xlim=c(0,1), cex.lab=1.5, cex.axis=1.5)
#export as PNG (500x900)
