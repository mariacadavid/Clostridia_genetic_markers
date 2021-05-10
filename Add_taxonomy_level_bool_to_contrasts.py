#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 09:54:16 2021

@author: mariacadavid
"""

#I have a table with contrast between genomes and the distances in two marker genes
#Here i aim to retrieve the taxonomy of both genomes and check if the correspond to the same fam/gen/specie. 
#Add a bool cols where i can tell if they correspond (true) or not (false) to the same fam/gen/specie. 


#Function to get taxonomy and compare taxonomy between contrasts
def obtain_taxonomy_bool(contrasts_distance_table, taxonomy_metadata):
   
   #Open file contrasts
   contrasts_table = open(contrasts_distance_table,"r")
   line = contrasts_table.readline()[:-1]
   line = contrasts_table.readline()[:-1]
   
   #create empty list for contrasts and identities
   contrast_col= []
   identity_gen1 = []
   identity_gen2=  []
   
   #Create empty list for real and boolean data taxonomy
   order_genome_1 = []
   order_genome_2 = []
   order_bool= []
   family_genome_1 =[]
   family_genome_2 =[]
   family_bool= []
   genus_genome_1= []
   genus_genome_2= []
   genus_bool= []
   specie_genome_1= []
   specie_genome_2= []
   specie_bool= []
   
   
   #Check for similarity in taxonomy
   tracking_files = 0
   while line != '':
       
        identity_1 = line.split("\t")[3] 
        identity_gen1.append (identity_1) #list for one maker distance between contrast
        identity_2= line.split("\t")[4]
        identity_gen2.append (identity_2) #list for the other marker distance between contrast
       
        contrast = line.split("\t")[0]
        contrast_col.append (contrast [1:-1])
        header_genome_1 = contrast.split("_compared_to_")[0]
        nucleotide_id_1 = header_genome_1.split("_")[0]
        nucleotide_id_1= nucleotide_id_1[1:]
        header_genome_2 = contrast.split("_compared_to_")[1]
        nucleotide_id_2 = header_genome_2.split("_")[0]
        
        
        #recorrer las filas del tax table
        taxonomy_table = open(taxonomy_metadata,"r")
        line_tax = taxonomy_table.readline()[:-1]
        line_tax = taxonomy_table.readline()[:-1]
    
        
        while line_tax != '':
            nucleotide_id = line_tax.split(",")[0]
            nucleotide_id= nucleotide_id [1:-1]
    
            
            if nucleotide_id == nucleotide_id_1:
                order_1 = line_tax.split(",")[5]
                family_1 = line_tax.split(",")[6]
                genus_1 = line_tax.split(",")[7] 
                specie_1 = line_tax.split(",")[8]
    
            if nucleotide_id == nucleotide_id_2:
                order_2 = line_tax.split(",")[5]
                family_2 = line_tax.split(",")[6]
                genus_2 = line_tax.split(",")[7] 
                specie_2 = line_tax.split(",")[8]
            
            line_tax = taxonomy_table.readline()[:-1]
      
        order_genome_1.append (order_1)
        family_genome_1.append (family_1)
        genus_genome_1.append (genus_1)
        specie_genome_1.append (specie_1)
        order_genome_2.append (order_2)
        family_genome_2.append (family_2)
        genus_genome_2.append (genus_2)
        specie_genome_2.append (specie_2)
            


        if nucleotide_id_1 == nucleotide_id_2:
            order_bool.append("same genome")
            family_bool.append("same genome")
            genus_bool.append("same genome")
            specie_bool.append("same genome")
                
        else:
            if order_1 == order_2: #bool order
                order_bool.append("true")
            else: 
                order_bool.append("false")
                
            if family_1 == family_2: #bool family
                family_bool.append("true")
            else: 
                family_bool.append("false")
            
            if genus_1 == genus_2: #bool genera
                genus_bool.append("true")
            else: 
                genus_bool.append("false")
            
            if specie_1 == specie_2: #bool specie
                specie_bool.append("true")
            else: 
                specie_bool.append("false")

        line = contrasts_table.readline()[:-1]
        tracking_files= tracking_files + 1
        print (tracking_files)
        
   return (contrast_col, identity_gen1, identity_gen2,
           order_genome_1, order_genome_2, order_bool, 
           family_genome_1,family_genome_2, family_bool,
           genus_genome_1, genus_genome_2, genus_bool,
           specie_genome_1, specie_genome_2, specie_bool)
        


#function to build dataframe from lists
def write_lists_to_df_file(lists: list , output_filename:str): 
     contrast_col= lists[0]
     identity_gen1= lists[1]
     identity_gen2= lists[2]
     order_genome_1= lists[3]
     order_genome_2= lists[4]
     order_bool= lists[5]
     family_genome_1= lists[6]
     family_genome_2= lists[7]
     family_bool= lists[8]
     genus_genome_1= lists[9]
     genus_genome_2= lists[10]
     genus_bool= lists[11]
     specie_genome_1= lists[12]
     specie_genome_2= lists[13]
     specie_bool= lists[14]
     
     import pandas as pd 
     df = pd.DataFrame(list(zip(contrast_col, identity_gen1, identity_gen2,
           order_genome_1, order_genome_2, order_bool, 
           family_genome_1,family_genome_2, family_bool,
           genus_genome_1, genus_genome_2, genus_bool,
           specie_genome_1, specie_genome_2, specie_bool)), 
               columns =['contrast', 'identity_16s', 'identity_dnak1',   #Aca es donde tengo que cambiar el nombre de las columnas segun el gen que este usando 
                         'order_genome_1', 'order_genome_2', 'same_order',
                         'family_genome_1','family_genome_2', 'same_family',
                         'genus_genome_1', 'genus_genome_2', 'same_genus',
                         'specie_genome_1', 'specie_genome_2', 'same_specie'])  
    
     df.to_csv(output_filename, sep='\t', header=True, index=False)
###


if __name__ == "__main__":
    import sys
    args = (sys.argv)
    contrasts_distance_table = args[1]
    taxonomy_metadata = args [2]
    output_filename= args [3]
    lists = obtain_taxonomy_bool(contrasts_distance_table, taxonomy_metadata)
    write_lists_to_df_file(lists , output_filename)


#Command line: 
    #general= python add_taxonomy_level_bool_to_contrasts.py contrasts_distance_table taxonomy_metadata output_filename
        #command for dnak1 vs gyrb= python add_taxonomy_level_bool_to_contrasts.py "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_dnak1_gyrb_common.tsv" "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Reference_genomes_ncbi/DNA_clostridia_GTDB/summary_genomes_nucleotides_blast_separate_taxonomy.csv" "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_dnak1_gyrb_common_taxonomy_analysis.tsv"
        #command for 16s vs dnak1= python add_taxonomy_level_bool_to_contrasts.py "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_dnak1_common.tsv" "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Reference_genomes_ncbi/DNA_clostridia_GTDB/summary_genomes_nucleotides_blast_separate_taxonomy.csv" "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_dnak1_common_taxonomy_analysis.tsv"
        #command for 16s vs gyrb= python add_taxonomy_level_bool_to_contrasts.py "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_gyrb_common.tsv" "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Reference_genomes_ncbi/DNA_clostridia_GTDB/summary_genomes_nucleotides_blast_separate_taxonomy.csv" "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_gyrb_common_taxonomy_analysis.tsv"