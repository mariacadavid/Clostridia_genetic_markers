#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 16:06:06 2021

@author: mariacadavid
"""

#It needs to run from in the directory where the .fna files of interest are

import os
import pandas as pd
from Bio import SeqIO


list_assemblies = []
assembly_files = []

#test 
#dir_assemblies = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Reference_genomes_ncbi/DNA_clostridia_GTDB/test_folder"

#real
dir_assemblies = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Reference_genomes_ncbi/DNA_clostridia_GTDB/GTDB_reference_genomes"

for f in os.listdir(dir_assemblies):
    if os.path.splitext(f)[1] == ".fna":
       list_assemblies.append(os.path.splitext(f)[0])
       assembly_files.append(f)

assembly_id=[]
nucelotide_id=[]


for assembly in assembly_files:
    for seq_record in SeqIO.parse(assembly,"fasta"):
        part_1= assembly.split("_")[0]
        part_2=assembly.split("_")[1]
        a_id= part_1 + "_" + part_2
        assembly_id.append(a_id)
        nucelotide_id.append(seq_record.id)
        
        
id_full_table = pd.DataFrame(list(zip(assembly_id, nucelotide_id)),
               columns =['assembly', 'nucleotide'])  

id_full_table.to_csv("/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Reference_genomes_ncbi/DNA_clostridia_GTDB/ids_of_selected_clostridia_genomes_full.tsv", sep='\t', header=True, index=False)
