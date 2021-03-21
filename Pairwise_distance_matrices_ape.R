#Obtain pairwise distance matrices and frequency histograms (dnak1, gyrb, 16s)
#Alignment needed 
#By: Maria Cadavid Jan2021
#Bassed on: https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter5.html

#libraries
library(seqinr) #read.alignment
library(ape) #dist.dna
library (spaa) #dist2list


#Distance matrix for dnak1
gene_aln  <- read.alignment(file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Alineamientos/Alineamientos_genes_GTDB/genes_dnak1_GTDB_tblastn_aligned_MAFFT", format = "fasta") #Input alignment
#dist.dna() function requires the input alignment to be in a special format known as “DNAbin” format, so we must use the as.DNAbin() function
aln_bin <- as.DNAbin(gene_aln) # Convert the alignment to "DNAbin" format
dist_matrix <- dist.dna(aln_bin, model= "K80", as.matrix = TRUE, pairwise.deletion = TRUE) # Calculate the genetic distance matrix
write.table(dist_matrix, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/indentity_matrix_ape_dnak1.tsv", sep = "\t",
            row.names = TRUE, col.names = TRUE)

#Frequency histogram dnak1
dist_object <- dist.dna(aln_bin, model= "K80", as.matrix = FALSE, pairwise.deletion = TRUE) # Calculate the genetic distance in dist object
dist_list_dnak1 <- dist2list(dist_object)
dist_frame <- as.data.frame(dist_list_dnak1)
colnames(dist_frame)[which(names(dist_frame) == "col")] <- "org1"
colnames(dist_frame)[which(names(dist_frame) == "row")] <- "org2"
colnames(dist_frame)[which(names(dist_frame) == "value")] <- "distances"
freq_hist = hist(dist_frame$distances, breaks= 100, col="pink", xlab="Distances between reference genes dnak1", main="Frequency of distances calculated by Ape (Clostridia dnak1)",xlim=c(0,1))
tiff(file="/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/frequency_distances_ape_dnak1.tiff", units="cm" ,width=22, height=13, res=600)
hist(dist_frame$distances, breaks= 100, col="pink", xlab="Distances between reference genes dnak1", main="Frequency of distances calculated by Ape (Clostridia dnak1)",xlim=c(0,1))
dev.off()


#Distance matrix for gyrb (cut 1281_3955)
gene_aln  <- read.alignment(file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Alineamientos/Alineamientos_genes_GTDB/genes_gyrb_GTDB_tblastn_aligned_MAFFT.fa", format = "fasta") #Input alignment
aln_bin <- as.DNAbin(gene_aln) # Convert the alignment to "DNAbin" format
dist_matrix <- dist.dna(aln_bin, model= "K80", as.matrix = TRUE, pairwise.deletion = TRUE) # Calculate the genetic distance matrix
write.table(dist_matrix, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/indentity_matrix_ape_gyrb_cut.tsv", sep = "\t",
            row.names = TRUE, col.names = TRUE)

#Frequency histogram gyrb
dist_object <- dist.dna(aln_bin, model= "K80", as.matrix = FALSE, pairwise.deletion = TRUE) # Calculate the genetic distance in dist object
dist_list_gyrb <- dist2list(dist_object)
dist_frame <- as.data.frame(dist_list_gyrb)
colnames(dist_frame)[which(names(dist_frame) == "col")] <- "org1"
colnames(dist_frame)[which(names(dist_frame) == "row")] <- "org2"
colnames(dist_frame)[which(names(dist_frame) == "value")] <- "distances"
freq_hist = hist(dist_frame$distances, breaks= 100, col="pink", xlab="Distances between reference genes gyrb", main="Frequency of distances calculated by Ape (Clostridia gyrb)",xlim=c(0,1))
tiff(file="/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/frequency_distances_ape_gyrb.tiff", units="cm" ,width=22, height=13, res=600)
hist(dist_frame$distances, breaks= 100, col="pink", xlab="Distances between reference genes gyrb", main="Frequency of distances calculated by Ape (Clostridia gyrb)",xlim=c(0,1))
dev.off()


#Distance matrix for 16s
gene_aln  <- read.alignment(file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Alineamientos/Alineamientos_genes_GTDB/genes_16s_GTDB_blastn_aligned_MAFFT.fasta", format = "fasta") #Input alignment
aln_bin <- as.DNAbin(gene_aln) # Convert the alignment to "DNAbin" format
dist_matrix<- dist.dna(aln_bin, model= "K80", as.matrix = TRUE, pairwise.deletion = TRUE) # Calculate the genetic distance matrix
write.table(dist_matrix, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/indentity_matrix_ape_16s.tsv", sep = "\t",
            row.names = TRUE, col.names = TRUE)


#Frequency histogram 16s
dist_object <- dist.dna(aln_bin, model= "K80", as.matrix = FALSE, pairwise.deletion = TRUE) # Calculate the genetic distance in dist object
dist_list_16s <- dist2list(dist_object)
dist_frame <- as.data.frame(dist_list_16s)
colnames(dist_frame)[which(names(dist_frame) == "col")] <- "org1"
colnames(dist_frame)[which(names(dist_frame) == "row")] <- "org2"
colnames(dist_frame)[which(names(dist_frame) == "value")] <- "distances"
freq_hist = hist(dist_frame$distances, breaks= 100, col="pink", xlab="Distances between reference genes 16s", main="Frequency of distances calculated by Ape (Clostridia 16s)", xlim=c(0,1))
tiff(file="/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/frequency_distances_ape_16s.tiff", units="cm" ,width=22, height=13, res=600)
hist(dist_frame$distances, breaks= 100, col="pink", xlab="Distances between reference genes 16s", main="Frequency of distances calculated by Ape (Clostridia 16s)",xlim=c(0,1))
dev.off()


