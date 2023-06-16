# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 14:41:29 2023

@author: lawashburn
"""

import pandas as pd
import os
import scipy
from scipy.spatial import distance
import numpy as np
import pysptools
from pysptools import distance as dist_2
from distance_metrics_mcda import distance_metrics
import csv
import re


final_rep = r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\KD_Training_Spectra\Perfect_spectra_compile\perfect_spectra_list_.csv"
dir_to_search = r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230214\KD_search_results_v5\SG_Unlabeled_DDA_TR3\xcorr_data"
output_directoy = r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230215\iteration_test"
sample = 'KD_' + 'SG_Unlabeled_DDA_TR3'
interation_numbers = 25

number_frag_charges = 4
peaks_per_peak = 5
neutral_loss = 3
number_ions = 2

final_rep = pd.read_csv(final_rep)
final_rep["Identifier"] = final_rep['Sequence'].astype(str) +"_"+ final_rep["Scan"].astype(str)

final_rep_seq = final_rep['Identifier'].values.tolist()
final_rep_seq = set(final_rep_seq)

def get_file_names_with_strings(str_list):
    full_list = os.listdir(dir_to_search)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

seqs = []
iteration_storage = []
dotttt = []


denominator = number_frag_charges*peaks_per_peak*neutral_loss*number_ions
fraction = 1/denominator

for seq in final_rep_seq:
    file_query = seq
    fragment_list = (get_file_names_with_strings([file_query]))
    print(file_query)
    print(fragment_list)
    
    mod_seq = re.sub("[\(\[].*?[\)\]]", "", seq)
    mod_seq = re.sub("_", "", mod_seq)
    mod_seq = re.sub(r'[0-9]+', "", mod_seq)

    sequence_length = len(mod_seq)

    if len(fragment_list) == 2:
        
        exp_df = pd.DataFrame()
        theo_df = pd.DataFrame()
        
        exp_query = 'exp'
        theo_query = 'theo'
        
        bin_storage = []
        bin_storage_no_dups = []
        for a in fragment_list:
            if exp_query in a:
                path = dir_to_search + '\\' + a
                
                read_file = pd.read_csv(path)
                read_file = read_file.round({'experimental m/z': 0})
                read_file = read_file.sort_values(by=['experimental intensity'],ascending=False)
                read_file = read_file.drop_duplicates(subset='experimental m/z')
                exp_df = pd.concat([exp_df,read_file])
                exp_df_bins = exp_df['experimental m/z'].values.tolist()
                bin_storage.extend(exp_df_bins)
            if theo_query in a:
                path = dir_to_search + '\\' + a
                
                read_file = pd.read_csv(path)
                read_file = read_file.round({'theoretical m/z': 0})
                read_file = read_file.sort_values(by=['theoretical intensity'],ascending=False)
                read_file = read_file.drop_duplicates(subset='theoretical m/z')
                theo_df = pd.concat([theo_df,read_file])
                theo_df_bins = theo_df['theoretical m/z'].values.tolist()
                bin_storage.extend(theo_df_bins)
            else:
                pass
        else:
            pass
            
        for bins in bin_storage:
            if bins not in bin_storage_no_dups:
                bin_storage_no_dups.append(bins)
            else:
                pass
        
        for iteration in range(0,interation_numbers):
            theo_exp_array_mod = theo_df.sample(n=sequence_length)
            theo_exp_array_mod = theo_exp_array_mod.rename(columns={'theoretical m/z': 'pseudo-experimental m/z', 'theoretical intensity': 'pseudo-experimental intensity'})
            theo_exp_array_mod['pseudo-experimental intensity'] = 50
            
            output_path = output_directoy + '\\'+ sample +'_' + seq + '_' + str(iteration) + '_iteration_test.csv'
            with open(output_path,'w',newline='') as filec:
                    writerc = csv.writer(filec)
                    theo_exp_array_mod.to_csv(filec,index=False)   

            merge_array = pd.DataFrame()
            merge_array['Bin'] = bin_storage_no_dups
            merge_array = merge_array.merge(theo_exp_array_mod,left_on='Bin', right_on='pseudo-experimental m/z',how='left')
            merge_array = merge_array.merge(theo_df,left_on='Bin', right_on='theoretical m/z',how='left')
            
            merge_array = merge_array.fillna(0)
            
            array1 = merge_array['pseudo-experimental intensity'].values.tolist()
            array2 = merge_array['theoretical intensity'].values.tolist()
            
            dot_prod = np.dot(array1, array2) #dot product calculation
            
            dotttt.append(dot_prod)
            seqs.append(seq)
            iteration_storage.append(iteration)
    else:
        pass
results_report = pd.DataFrame()
results_report['Sequence'] =  seqs      
results_report['Iteration'] =  iteration_storage
results_report['Dot product'] =  dotttt
    
output_path = output_directoy + '\\'+ sample + '_iteration_test.csv'
with open(output_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        results_report.to_csv(filec,index=False)   