#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 16:14:32 2020

@author: mariacadavid
"""

import sqlite3 as sql



    

def obtain_sequence(blast_output_file, complete_fasta_file, db):
    
   #Open files
   blastdata = open(blast_output_file,"r")
   line = blastdata.readline()[:-1]
   
   referencedata = open(complete_fasta_file,"r")
   refs= referencedata.readlines()
   referencedata.close() 
   print ("ARCHIVO LISTO")
   genes= []
   
   
   con = sql.connect(db)
   cur = con.cursor()
   
   while line != '':
        if line.split(",")[1] == "0.0":
            
            line_splited=line.split("l|")[1]
            line_splited_1= line_splited.split("_cds_")[0]
            partial_locus= line_splited_1.split(".")[0]
            
            #Add taxonomy form databse
            print(partial_locus)
            query= 'SELECT organism_name from Organism WHERE locus like (?)'
            cur.execute(query,("%"+partial_locus+"%",))
            tax = cur.fetchall()
            
            #Add locus from database
            #query= 'SELECT locus from Organism WHERE locus like (?)'
            #cur.execute(query,("%"+partial_locus+"%",))
            #loc = cur.fetchall()
            
            #taxonomy= (tax[0])
            name=">"+tax[0][0] + "_" + line.split(",")[0]  #name in fasta be: tax + locus
            start_position= int(line.split(",")[2])
            end_position=int (line.split(",")[3])
            line= blastdata.readline()
            idxline= 0                    #como tengo en memoria el archivo no lo tengo que volver a abrir y cerrar
            line_fasta = refs[idxline][:-1]
            header=">"+line.split(",")[0]
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
                    print (name)
                    genes.append(name)
                    genes.append(gene)
        

                else:
                    if idxline +1 != len(refs):
                        idxline += 1
                        line_fasta = refs[idxline][:-1]
                    else:
                        break
            
        else:
            line= blastdata.readline()
   
   con.close()        
   blastdata.close()
   return (genes)
  

   
def write_genes_to_fasta(genes:list, output_filename:str): 
    genes_fasta= open(output_filename,"w")
    genes_fasta.writelines("%s\n" % s for s in genes)
    genes_fasta.close()

genes= obtain_sequence("/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Reference_genomes_ncbi/CDS_clostridia_NCBI/genome_assemblies_cds_fasta/blastout_dnak1_tblastn_ref1uniprot_ncbiCDS_csv.txt", "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Reference_genomes_ncbi/CDS_clostridia_NCBI/genome_assemblies_cds_fasta/ncbi-genomes-2020-09-01/ncbi_all_genomes.fna","/Users/mariacadavid/Google Drive/Universidad /BIOLOGÍA/Enfasis_computacional/Programación_Enfasis/Gene_Bank_SQL_db/Gene_Bank_Data_Clostridia2.db")
write_genes_to_fasta(genes, "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Reference_genomes_ncbi/CDS_clostridia_NCBI/genes_taxonomy_dnak1_NCBI_tblastn.txt")