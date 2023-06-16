# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 15:41:25 2023

@author: lawashburn
"""
import csv
import pandas as pd
import re
import os
from itertools import permutations
from Bio.SeqIO.FastaIO import SimpleFastaParser
import random
import collections
import time
import numpy as np
from scipy import signal
from datetime import datetime
import scipy
from scipy.spatial import distance
import numpy as np

predefined_db_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale_weighting_hyperscore\weighting_input_data\database_search_input\target_decoy_df_full_validated.csv" #predefined DB path
base_file_path = r"C:\Users\lawashburn\Documents\DBpep_v2\validation\Nhu_RAW_files"
output_parent_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\PEAKS_oursoftware_compare_brain_only\DB_search_optimize\03_precursor_AMM_optimize_v2"

raw_conv_mzml_storage = [[r"C:\Users\lawashburn\Documents\DBpep_v2\validation\Nhu_RAW_files\2021_0817_Brain_1.mzML",
                          r"C:\Users\lawashburn\Documents\DBpep_v2\validation\Nhu_RAW_files\2021_0817_Brain_1_formatted.txt"]]

h_mass = 1.00784
precursor_error_cutoff = 20 #ppm
spectra_segments = 50

###Definitions###
def raw_file_detail_extraction(raw_file_path):
    raw_file_sample_name1 = raw_converter_path.replace(base_file_path,'')
    raw_file_sample_name2 = raw_file_sample_name1.replace('_formatted','')
    raw_file_sample_name3 = raw_file_sample_name2.replace('\\','')
    sample_name = raw_file_sample_name3.replace('.txt','')
    output_folder = output_parent_directory+'\\'+sample_name
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    ### generate output folder ###
    return sample_name, output_folder

def raw_file_data_extraction(raw_file_path):
    raw_converter = pd.read_csv(raw_converter_path, sep=",",skiprows=[0], names= ["m/z","resolution","charge","intensity", "MS2", "scan_number","precursor_charge",'null','Sample','Identifier','Iteration'])
    raw_converter = raw_converter.rename(columns={'m/z':'Fragment actual m/z',
                                                  'charge': 'Fragment actual charge',
                                                  'intensity':'Fragment actual intensity',
                                                  'MS2':'Precursor actual m/z',
                                                  "precursor_charge":'Precursor actual charge',
                                                  'scan_number':'Scan'})

    raw_converter['Fragment actual charge'] = raw_converter['Fragment actual charge'].replace(to_replace=0,value=1) #assume z=0 is z=1   
    exp_precursor = raw_converter.drop_duplicates() 
    exp_precursor = exp_precursor.copy()
    exp_precursor['Precursor actual monoisotopic mass'] =  ((exp_precursor['Precursor actual m/z']) * 
                                                            (exp_precursor['Precursor actual charge']))-(h_mass*(exp_precursor['Precursor actual charge']))
    return exp_precursor

def precursor_amm(db,exp_precursor):
    merge_match=pd.DataFrame()
    precursor_temp_cutoff = precursor_error_cutoff*3 #rough estimate of ppm to Da to minimize extra search space
    
    if len(db)<1: #throws an error if database file is empty
        raise ValueError('Database file is empty')
    
    db_sorted = db.sort_values(by = 'Precursor theoretical monoisotopic mass') #sort database monoisotopic mass and experimental for mergeasof
    exp_precursor_sorted = exp_precursor.sort_values(by = 'Precursor actual monoisotopic mass')
        
    for a in range(0,(len(db_sorted)-1)):
        print(a)
        db_filtered = db_sorted.iloc[[a]]
        
        merge_match2 = pd.merge_asof(exp_precursor_sorted,db_filtered, left_on='Precursor actual monoisotopic mass', 
                                    right_on='Precursor theoretical monoisotopic mass',
                                    tolerance = precursor_temp_cutoff, allow_exact_matches=True,direction='forward') 

        merge_match3 = pd.merge_asof(exp_precursor_sorted,db_filtered, left_on='Precursor actual monoisotopic mass', 
                                    right_on='Precursor theoretical monoisotopic mass',
                                    tolerance = precursor_temp_cutoff, allow_exact_matches=True,direction='backward') 
        merge_match = pd.concat([merge_match,merge_match2, merge_match3])

    merge_match = merge_match.drop_duplicates()
    merge_match = merge_match.sort_values(by='scan')
    merge_match = merge_match[merge_match['sequence'].notna()]
    
    
    scans_missing_check = []
    
    # for xxxx in scans_check:
    #     if xxxx not in scans_accounted:
    #         scans_missing_check.append(xxxx)
    #     else:
    #         pass
    
    # if len(scans_missing_check)>0:
    #     exp_precursor_segment = exp_precursor_sorted[exp_precursor_sorted['Scan'].isin(scans_missing_check)]

    #     merge_match2 = pd.merge_asof(exp_precursor_segment,db_sorted, left_on='Precursor actual monoisotopic mass', 
    #                                 right_on='Precursor theoretical monoisotopic mass',
    #                                 tolerance = precursor_error_cutoff, allow_exact_matches=True,direction='forward') 
    #     merge_match3 = pd.merge_asof(exp_precursor_segment,db_sorted, left_on='Precursor actual monoisotopic mass', 
    #                                 right_on='Precursor theoretical monoisotopic mass',
    #                                 tolerance = precursor_error_cutoff, allow_exact_matches=True,direction='backward') 
        
    #     merge_match = pd.concat([merge_match, merge_match2, merge_match3])
    # else:
    #     pass
    
    # merge_match= merge_match.drop_duplicates()

    # merge_match_filtered = merge_match.dropna(subset=['Sequence','Precursor theoretical monoisotopic mass'])
    # merge_match_filtered = merge_match_filtered.copy() #silences pesky warning
    # merge_match_filtered['Precursor error (ppm)'] = ((abs((merge_match_filtered['Precursor theoretical monoisotopic mass'])-
    #                                                       (merge_match_filtered['Precursor actual monoisotopic mass'])))/
    #                                                   (merge_match_filtered['Precursor theoretical monoisotopic mass'])) * 1E6

    # precursor_amm_results = merge_match_filtered[merge_match_filtered['Precursor error (ppm)'] <= precursor_error_cutoff]

    file_path = peptide_report_output + '\\precursor_amm_results.csv'
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            merge_match.to_csv(filec,index=False)   

    return merge_match
###

fasta_w_mass = pd.read_csv(predefined_db_path)

for file_category in raw_conv_mzml_storage:
    start = time.time()
    raw_converter_path = file_category[1]
    mzml_path_input = file_category[0]
    rounds_number = [1]
    unique_IDS_number = []
    details = raw_file_detail_extraction(raw_converter_path)
    sample_name = details[0]

    peptide_report_output = details[1]
    exp_precursor = raw_file_data_extraction(raw_converter_path)
  #%%  
    fasta_w_mass = fasta_w_mass[fasta_w_mass["Sequence"].str.contains("\(") == False]
    precursor_amm_results = precursor_amm(fasta_w_mass,exp_precursor)
    print(len(precursor_amm_results))
    finish = time.time()
    elapsed_time = start-finish
    elapsed_time_report = 'Time elapsed : ' + str(elapsed_time) 