## Evaluation of genetic markers for class Clostridia 

**Goals:**
Evaluate the suitability of dna Kinase1 and gyrase subunit B genes as genetic markers for the study of community diversity of the class Clostridia.
-	Extract dnaK1 and gyrB gene sequences by identity with a reference from complete genomes available in public databases for microorganisms of the class Clostridia.
-	Estimate the genetic variability of dnaK1 and gyrB genes in the Clostridia class compared to the universal marker 16s rRNA.
-	Compare taxonomic resolution found with dnaK1 and gyrB genes against the universal marker 16s rRNA in the class Clostridia.


**Methods:**

**1. Download reference genomes from database: GTDB**  
A search of reference genomes was performed through the recent Genome Taxonomy Database (GTDB). This database is an initiative of the Australian Research Council Laurate Fellowship to establish a standardized microbial taxonomy based on genome phylogeny (Parks et al., 2020). GTDB uses genomes published to Genbank and RefSeq to construct a revised phylogeny. To build our dataset of reference genomes, we selected accession numbers of genomes that belong to the Clostridia class according to the GTDB taxonomy, which  have a high completeness (>= 95%), low contamination (<1%) and have fewer than 10 SSU copies (often a sign of mis-assembly). Using the selected accession numbers, the corresponding genome sequences were subsequently downloaded from the National Center for Biotechnology Information (NCBI). A total of 4069 genomes were downloaded for the dataset on 28/10/2020.   
(R Script for retrieval of accession numbers can be found here as “Reference_genomes_retrieval.R”, Download was done following ncbi-genome-download pipeline: https://github.com/kblin/ncbi-genome-download)


**2.	Search of dnaK1, gyrB and 16s rRNA genes in reference genomes**  
Once the dataset with reference genomes was built, we searched within each genome the genes of interest: dna Kinase1 (dnaK1) and gyrase subunit B (gyrB) and also the widely used universal marker gene 16s rRNA.   
To retrieve  the  gene sequences within the reference genomes we used two Basic Local Alignment Search Tools (BLAST) provided by NCBI. For the 16s rRNA gene we used the traditional BLASTn program, which carries out a search of a nucleotide query sequence within a nucleotide reference database. The 16s rRNA query sequence used was downloaded from NCBI, reported as ‘16s ribosomal RNA partial sequence’ (NR_074399.1), size 1500 base pairs, from the source organism Ruminococcus albus 7 = DSM 20455. The BLASTn was executed with ‘max_target_seqs’ parameter set to 5000 to indicate the number of aligned sequences to keep and ‘max_hsps’ parameter set to 1 to keep a single hit/alignment for any query-subject pair. For the protein coding genes dnaK1 and gyrB we used tBLASTn, a translated version of BLAST which finds regions of local similarity between a query protein sequence compared to the six-frame translations of sequences in a nucleotide reference database (Wheeler & Bhagwat, 2007).   
Query sequences for tBLASTn were downloaded from UniProt database. For gene dnaK1 we selected a protein sequence classified as Chaperone protein DnaK (A0A143WYF3), gene dnaK1, size: 594 amino acids, from the organism Clostridiales bacterium CHKCI006. For gene gyrB we selected a protein sequence classified as DNA gyrase subunit B (R6D2S9), gene gyrB, size: 680 amino acids, from the organism Ruminococcus sp. CAG:579. The tBLASTn was also executed with the ‘max_target_seqs’ parameter set to 5000.   
(Bash Script for BLAST commands can be found here as “BLAST_commands.txt”)  


**3.	Extraction of genes’ sequence**  
We accepted BLAST hits with expected values (e-value) equal to 0 as a significant threshold. The significant hits contain information of the position where the gene of interest is located in each  reference genome and the frame on which it is read. We used this information to extract the nucleotide sequences corresponding to our genes of interest and load them into 3 independent files, one per gene (dnaK1, gyrB and 16S rRNA) containing all its hit sequences.   
(Python Script for the creation of new file with extracted genes can be found here as “Extract_blasted_genes.py”)  
Note: in the final file we had to change the spaces for “_”, the “:” for “-“ and erase all parenthesis() and square brackets[] to avoid downstream errors.   


**4.	Alignment and pairwise distances**  
For each gene we carried out an alignment. Each gene’s hit sequences were aligned using MAFFT  (Katoh & Standley, 2013) with default parameters (mafft/7.402-with-extensions_intel-17.0.1) as it proved to be computationally efficient and accurate.   
Example command: >> mafft --quiet --auto --thread $SLURM_NTASKS genes_dnak1_NCBI_tblastn.txt > genes_dnak1_NCBI_tblastn_aligned_MAFFT.fasta  
Additionally, we did translated alignments using MEGA X to check that the extracted sequences were working fine and matching. Alignments were reviewed manually for a quality check, in this process, the gyrB alignment was trimmed at the ends keeping the section the nucleotide positions 1281 and 3955.    
Once we had the alignments, we built pairwise distance matrices to compare how different one organism is from another in a given gene. To calculate the distances between every pair of sequences we used the R  ‘Analyses of Phylogenetics and Evolution’ (APE) package with Kimura evolution model (k80). Distances were expressed as a number between 0 and 1; the greatest the distance, the more different the sequences and the closest to cero, the more similar the sequences. Finally, we evaluated the frequency distribution of the pairwise distances for each gene.    
(R Script for the pairwise distances calculations and frequency graphs can be found here as “Pairwise_distance_matrices_ape.r”)  


**5.	Comparison of marker genes**   
With the pairwise distances we calculated identity between a pair of reference genomes in a given gen. Having the pairwise identities for each gene, we contrasted these pairs across the three different genes.  We calculated a linear regression for the pairwise identities of every two genes, hence we obtained a linear regression for the comparison of dnaK1 vs. 16s rRNA, a second regression for the comparison of gyrB vs. 16s rRNA and a third regression for the comparison of dnaK1 vs. gyrB. These regressions allow us to compare variability between two genes within the same genomes.    
(R Script for the comparison of genes and graphs can be found here as “Compare_gene_distances.r”)  
To understand how the variability of the genes behaves at different taxonomic levels, we analyzed how data belonging to the same order, family, genus and specie are distributed throughout the linear regression. To assign the taxonomic classification of the genome of origin to each gene sequence we matched the nucleotide ID of the sequence where the gene was found and its assembly ID which has an associated GTDB taxonomic classification.   


**6.	Phylogenic trees**  
For the phylogenetic analysis we used Randomized Axelerated Maximum Likelihood (RAxML) version 8 (Stamatakis, 2014). We built a phylogenetic tree per gene with the aligned sequences using GRT+Gamma as substitution model and algorithms of rapid Bootstrap analysis and searched for best-scoring ML tree in one program run. Finally, we compared the phylogenetic trees obtained using the different genes and this way, reaffirm the genes’ suitability as genetic markers for the class Clostridia, the most variable gene is expected to be a better marker for these closely related taxa.   


**References**  
Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: Improvements in performance and usability. Molecular Biology and Evolution, 30(4), 772–780. https://doi.org/10.1093/molbev/mst010  
Parks, D. H., Chuvochina, M., Chaumeil, P. A., Rinke, C., Mussig, A. J., & Hugenholtz, P. (2020). A complete domain-to-species taxonomy for Bacteria and Archaea. Nature Biotechnology, 38(9), 1079–1086. https://doi.org/10.1038/s41587-020-0501-8  
Stamatakis, A. (2014). RAxML version 8: A tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics, 30(9), 1312–1313. https://doi.org/10.1093/bioinformatics/btu033  
Wheeler, D., & Bhagwat, M. (2007). BLAST QuickStart: example-driven web-based BLAST tutorial. In Methods in molecular biology (Clifton, N.J.) (Vol. 395). https://doi.org/10.1007/978-1-59745-514-5_9  

