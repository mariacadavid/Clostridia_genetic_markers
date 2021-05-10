#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 16:15:56 2021

@author: mariacadavid
"""

import pandas as pd 

#####open dataframe per gene
df_dnak1= pd.read_csv('/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/compare_dataframe_dnak1.tsv', sep='\t')
df_dnak1.head()

df_16s= pd.read_csv('/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/compare_dataframe_16s.tsv', sep='\t')
df_16s.head()

df_gyrb= pd.read_csv('/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Indentity_matrix_R/compare_dataframe_gyrb.tsv', sep='\t')
df_gyrb.head()

#####create dataframe of comparisions
#16s vs dnak1
df_compare_16s_dnak1= df_16s.set_index('comparison_16s').join(df_dnak1.set_index('comparison_dnak1'))
df_compare_16s_dnak1.head()
df_compare_16s_dnak1.to_csv("/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/distances_16s_dnak1.csv", sep='\t')

#16s vs gyrb 
df_compare_16s_gyrb= df_16s.set_index('comparison_16s').join(df_gyrb.set_index('comparison_gyrb'))
df_compare_16s_gyrb.head()
df_compare_16s_gyrb.to_csv("/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/distances_16s_gyrb.csv", sep='\t')

#dnak1 vs gyrb 
df_compare_dnak1_gyrb= df_dnak1.set_index('comparison_dnak1').join(df_gyrb.set_index('comparison_gyrb'))
df_compare_dnak1_gyrb.head()
df_compare_dnak1_gyrb.to_csv("/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/distances_dnak1_gyrb.csv", sep='\t')


#####create dataframe for the 3 markers

df_dnak1_new = df_dnak1
df_dnak1_new.head()
df_dnak1_new.columns = ['distance_dnak1', 'pairs']  

df_16s_new = df_16s
df_16s_new.head()
df_16s_new.columns = ['distance_16s', 'pairs']  

df_gyrb_new= df_gyrb
df_gyrb_new.head()
df_gyrb_new.columns = ['distance_gyrb', 'pairs'] 


merge1= pd.merge(df_dnak1_new,df_16s_new,on= "pairs",how='outer')
merge_final= pd.merge(merge1,df_gyrb_new,on= "pairs",how='outer')
merge_final.to_csv("/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/distances_complete.csv", sep='\t')




