## Evaluation of genetic markers for the class Clostridia
Repository for project: Comparison of the variation of gut bacterial diversity marker genes in public genomes: a study case of the class Clostridia.
Maria Cadavid-Velez, 2021

**Goals**  
Evaluate the suitability of dna Kinase1 and gyrase subunit B genes as genetic markers for the study of community diversity of the class Clostridia.
-	Extract dnaK1 and gyrB gene sequences by identity with a reference from genomes available in public databases for microorganisms of the class Clostridia.
-	Estimate the genetic variability of dnaK1 and gyrB genes in the Clostridia class compared to the universal marker 16s rRNA.



**Methods**  
 
![image](https://user-images.githubusercontent.com/37601806/117590497-5a190f80-b0f5-11eb-83e1-cbf0e4e8277b.png)  
Figure 1. General methods workflow: steps (white boxes) and detail of written scripts and used programs for each purpose (gray boxes).  

**Download reference genomes from public database: GTDB**  
In order to build a local database of reference genomes of the class Clostridia, a search of genome assemblies was performed through the GTDB (https://gtdb.ecogenomic.org). This public database is an initiative of the Australian Research Council Laurate Fellowship to establish a standardized microbial taxonomy based on genome phylogeny (Parks et al., 2020). GTDB uses genome sequences published in Genbank and RefSeq to construct a revised phylogeny. To determine our set of reference genome assemblies, we selected accession numbers of assemblies that belong to the Clostridia class according to the GTDB taxonomy, which have a high completeness (>= 95%), low contamination (<1%) and fewer than ten 16S rRNA copies to avoid mis-assemblies (Větrovský & Baldrian, 2013). Using the selected accession numbers, the corresponding genome assembly sequences were subsequently downloaded from the National Center for Biotechnology Information (NCBI). A total of 4073 assemblies were downloaded on 28/10/2020. Within the downloaded assemblies, we had both types of data: culture independent metagenome assembled genomes (MAGs) and genome assemblies from cultured and isolated strains; both with different assembly/fragmentation levels. With this set of selected assemblies, a local database was built with the application Makeblastdb version 2.6.0 (Fig. 1). (R Script for retrieval of accession numbers can be found here as “Reference_genomes_retrieval.R”, Download was done following ncbi-genome-download pipeline: https://github.com/kblin/ncbi-genome-download).

**Search of dnaK1, gyrB and 16S rRNA genes in the local database**  
Once the local database was built, we proceeded to identify and retrieve, within each assembly, the genes of interest: dnaK1 and gyrB, and the widely used universal marker 16S rRNA. The process of identifying such genes in the assemblies’ sequences was performed through similarity with a query sequence. Query sequences for dnaK1 and gyrB were downloaded from the UniProt database. For dnaK1, we selected a protein sequence classified as Chaperone protein DnaK (A0A143WYF3), gene ‘dnaK1’, size: 594 amino acids, from the organism Clostridiales bacterium CHKCI006. For gyrB, we selected a protein sequence classified as DNA gyrase subunit B (R6D2S9), gene ‘gyrB’, size: 680 amino acids, from the organism Ruminococcus sp. CAG:579. The 16S rRNA query sequence was downloaded from NCBI, reported as ‘16S ribosomal RNA partial sequence’ (NR_074399.1), size 1500 base pairs, from the source organism Ruminococcus albus 7 = DSM 20455.  

Each query sequence was used to find the coordinates of the queried gene on each whole assembly sequence in the local database. For this, we used two Basic Local Alignment Search Tools (BLAST) provided by the NCBI (Fig. 1). For the 16S rRNA gene, we used the BLASTn program, which carries out a search of a nucleotide query sequence within a nucleotide reference database. The BLASTn was executed with ‘max_target_seqs’ parameter set to 5000 to indicate the number of aligned sequences to keep and ‘max_hsps’ parameter set to 1 to keep a single alignment for any query-subject pair. For the protein coding genes dnaK1 and gyrB, we used tBLASTn, a translated version of BLAST which finds regions of local similarity between a query protein sequence compared to the six-frame translations of sequences in a nucleotide reference database (Wheeler & Bhagwat, 2007). tBLASTn was also executed with the ‘max_target_seqs’ parameter set to 5000. (Bash Script for BLAST commands can be found here as “BLAST_commands.txt”).

**Extraction of genes’ sequences**  
The BLAST results contained a list of records that aligned with the query. The records, here called hits, contained information of their coordinates in each genome assembly and the frame on which it is read. We accepted BLAST hits with expected values (e-value) equal to 0 as a significant threshold. We took the significant hits’ coordinates to extract the nucleotide sequences corresponding to our genes of interest from the reference assemblies. The sequences corresponding to the significant hits for each gene were saved as three independent FASTA files (National Center of Biotechnology Information, 2021), one per gene (dnaK1, gyrB and 16S rRNA) (Fig. 1). (Python Script for the creation of new file with extracted genes can be found here as “Extract_blasted_genes.py”).   
Note: in the final file we had to change the spaces for “_”, the “:” for “-“ and erase all parenthesis() and square brackets[] to avoid downstream errors.     
Although we purposely chose house-keeping genes (dnaK1 and gyrB) with one expected copy per genome, when extracting the genes from the reference genome assemblies we checked for possible occurrences of multi-copies. (Python Script to check repeated hits can be found here as "Check_repeated_orgs_gtdb.py").

**Alignment and pairwise distances**  
For each gene, we carried out an alignment with the retrieved sequences. We tested four different aligning programs: PRANK, MUSCLE, MAFFT and Clustal (Edgar, 2004; Katoh & Standley, 2013; Löytynoja, 2014; Sievers et al., 2011). The alignments were run in the Apolo Scientific Computing Center at Universidad EAFIT and were visually inspected for quality check (Fig. 1). During this process, the gyrB alignment was trimmed at the ends keeping the section between 1281 and 3955 nucleotide positions.  

In order to quantify the variation for each marker, we built pairwise distance matrices from the alignments (Fig. 2). To compute a distance matrix, each pair of sequences was considered separately to calculate a respective distance value. In brief, two aligned sequences were assumed to be homologous (i.e., to share a common ancestor) and each nucleotide was considered as a character. Similarities mean that the character was conserved, whereas differences were considered as derived traits that originated when substitutions, deletions or insertions occurred in the nucleotide position (Yang & Bielawski, 2000). The pairwise distances give an estimation of the degree of similarity (0) or dissimilarity (1) of two sequences from two different genome assemblies. We calculated the distances between every pair of sequences using the R ‘Analyses of Phylogenetics and Evolution’ (APE 5.0) package (Paradis & Schliep, 2019) with Kimura 1980 evolution model (K80). This two-parameter model assumes equal base frequencies, one transition probability and one transversion probability (Kimura, 1980). 

A frequency distribution of all the pairwise distances for each gene allowed us to contrast how similar or dissimilar were the sequences of the different markers. Additionally, with the pairwise distances, we calculated the reciprocal identities, used later on, with values closer to one (1) indicating similarity between the sequences. (R Script for the pairwise distances calculations and frequency graphs can be found here as “Pairwise_distance_matrices_ape.R").

![image](https://user-images.githubusercontent.com/37601806/117590507-6ef5a300-b0f5-11eb-87ad-d8ee86de30ff.png)  
Figure 2. Workflow diagram for alignment and pairwise distances calculation.  

**Comparison of marker genes**  
To understand in more detail how does the genetic variation correlate between marker genes, we compared matrices for every two genes (16S vs. dnaK1, 16S vs. gyrB and dnaK1 vs. gyrB) (Fig. 1). We plotted the pairwise identities of one gene against the pairwise identities of the second gene and calculated linear regressions, which resulted in one regression for the comparison of dnaK1 vs. 16S rRNA, a second regression for the comparison of gyrB vs. 16S rRNA and a third one for the comparison of dnaK1 vs. gyrB. These regressions allowed us to compare variability between two genes within the same genome assembly. (R and Python Scripts for the comparison of genes and graphs can be found here as “Compare_gene_distances.R” and "Merge_df_distances.py").  

Next, we visualized how the variability of the genes behaved at different taxonomic levels, by highlighting data belonging to the same order, family, genus and species throughout the linear regression. To do this, it was necessary to associate each gene sequence to the respective taxonomic classification of the assembly of origin. To assign the corresponding taxonomies, a match was made between the nucleotide ID of the single sequence fragment where the gene was found with its corresponding assembly ID, which has an associated GTDB taxonomic classification. We calculated standard deviations of data at all taxonomic levels to better understand the dispersion of identities of both marker genes (dnaK1 and gyrB) compared to that of the 16S rRNA gene.

Finally, we also highlighted the comparisons and calculated standard deviations for data belonging to two important families that were relevant in the characterization of the Colombian microbial gut population: Ruminococcaceae and Lachnospiraceae.  


**References**  
Edgar, R. C. (2004). MUSCLE: Multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research, 32(5), 1792–1797. https://doi.org/10.1093/nar/gkh340  
Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: Improvements in performance and usability. Molecular Biology and Evolution, 30(4), 772–780. https://doi.org/10.1093/molbev/mst010  
Kimura, M. (1980). A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. Journal of Molecular Evolution, 16(2), 111–120. https://doi.org/10.1007/BF01731581  
Löytynoja, A. (2014). Phylogeny-aware alignment with PRANK. In Methods in Molecular Biology (Vol. 1079, pp. 155–170). https://doi.org/10.1007/978-1-62703-646-7_10  
National Center of Biotechnology Information. (2021). BLAST Topics. Retrieved April 24, 2021, from https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp  
Paradis, E., & Schliep, K. (2019). Ape 5.0: An environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics, 35(3), 526–528. https://doi.org/10.1093/bioinformatics/bty633  
Parks, D. H., Chuvochina, M., Chaumeil, P. A., Rinke, C., Mussig, A. J., & Hugenholtz, P. (2020). A complete domain-to-species taxonomy for Bacteria and Archaea. Nature Biotechnology, 38(9), 1079–1086. https://doi.org/10.1038/s41587-020-0501-8  
Sievers, F., Wilm, A., Dineen, D., Gibson, T. J., Karplus, K., Li, W., … Higgins, D. G. (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology, 7(539). https://doi.org/10.1038/msb.2011.75  
Větrovský, T., & Baldrian, P. (2013). The Variability of the 16S rRNA Gene in Bacterial Genomes and Its Consequences for Bacterial Community Analyses. PLoS ONE, 8(2), 1–10. https://doi.org/10.1371/journal.pone.0057923  
Wheeler, D., & Bhagwat, M. (2007). BLAST QuickStart: example-driven web-based BLAST tutorial. In Methods in molecular biology (Clifton, N.J.) (Vol. 395). https://doi.org/10.1007/978-1-59745-514-5_9  
Yang, Z., & Bielawski, J. R. (2000). Statistical methods for detecting molecular adaptation. Trends in Ecology and Evolution, 15(12), 496–503. https://doi.org/10.1016/S0169-5347(00)01994-7  



