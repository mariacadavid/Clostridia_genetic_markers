#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 15:32:42 2020
@author: mariacadavid
"""



def obtain_sequence(blast_output_file, complete_fasta_file):
    
   #Open files
   blastdata = open(blast_output_file,"r")
   line = blastdata.readline()[:-1]
   
   referencedata = open(complete_fasta_file,"r")
   refs= referencedata.readlines()
   referencedata.close() 
   print ("ARCHIVO LISTO")
   genes= []
   
   while line != '':
        if line.split(",")[1] == "0.0":
            header=">"+line.split(",")[0] 
            start_position= int( line.split(",")[2])
            end_position=int (line.split(",")[3])
            line= blastdata.readline()
            idxline= 0                    #como tengo en memoria el archivo no lo tengo que volver a abrir y cerrar
            line_fasta = refs[idxline][:-1]
            while idxline < len(refs): #controlando que no se termino la lista de lineas en memoria
                if header == line_fasta:
                    seqnuc=""
                    idxline += 1
                    line_fasta = refs[idxline][:-1]
                    while idxline < len(refs):
                        if line_fasta[0] != ">":
                            seqnuc += str(line_fasta)
                            idxline += 1
                            line_fasta = refs[idxline][:-1]
                        else:
                            break
                
                    gene= seqnuc[start_position -1 :end_position]
                    genes.append(header)
                    genes.append(gene)
        

                else:
                    if idxline +1 != len(refs):
                        idxline += 1
                        line_fasta = refs[idxline][:-1]
                    else:
                        break
            
        else:
            line= blastdata.readline()
   
            
   blastdata.close()
   return (genes)
  

   
def write_genes_to_fasta(genes:list, output_filename:str): 
    genes_fasta= open(output_filename,"w")
    genes_fasta.writelines("%s\n" % s for s in genes)
    genes_fasta.close()

genes= obtain_sequence("/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Reference_genomes_ncbi/DNA_clostridia_GTDB/blastout_dnak1_tblastn_ref1uniprot_GTDBgenomes_csv.txt", "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Reference_genomes_ncbi/DNA_clostridia_GTDB/GTDB_reference_genomes/gtdb_all_genomes.fna")
write_genes_to_fasta(genes, "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Reference_genomes_ncbi/DNA_clostridia_GTDB/genes_dnak1_GTDB_tblastn.txt")