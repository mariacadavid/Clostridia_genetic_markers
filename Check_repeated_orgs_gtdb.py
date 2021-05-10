#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 16:27:59 2020

@author: mariacadavid
"""
#CHECK REPEATED FOR GTDB DATA

import os 
import io #tipo para un archivo "io.TextIOWrapper"
import collections

#Funcion para identificar id repetidos en NCBI data
def check_repeated_orgs_GTDB(file: io.TextIOWrapper):
    list_org = []
    set_org= ()
    #Validaciones 
    assert os.path.isfile(file) , "archivo no existe"
    
    #Abrir archivo 
    data = open(file,"r") 
    
    #logica
    for line in data:
        if line[0] == '>':   #Buscar en los encabezados del fasta
            header= line[0:-1]
            list_org.append(header)
    data.close()
    
    set_org = set(list_org)
    
    if (len(list_org))!= (len(set_org)):
        print("contains replicates")
        print([item for item, count in collections.Counter(list_org).items() if count > 1])
    print (len(list_org), "total genes in the list")
    print (len (set_org), "total organisms")
    print ((len(list_org)-len (set_org)), "repeated hits = means two or more hits for the same organisms")


#GTDB dnak1
check_repeated_orgs_GTDB("/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Reference_genomes_ncbi/DNA_clostridia_GTDB/OUT_genes_dnak1_GTDB_tblastn.txt")

#GTDB gyrb
check_repeated_orgs_GTDB("/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Reference_genomes_ncbi/DNA_clostridia_GTDB/OUT_genes_gyrb_GTDB_tblastn.txt")

#GTDB gyrb
check_repeated_orgs_GTDB("/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Reference_genomes_ncbi/DNA_clostridia_GTDB/OUT_genes_16s_GTDB_blastn.txt")

