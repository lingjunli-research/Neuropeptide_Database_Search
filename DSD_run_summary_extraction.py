# -*- coding: utf-8 -*-
"""
Created on Sat May 27 16:43:31 2023

@author: lawashburn
"""

import pandas as pd

results_dir = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\DSD_eval\motif_1_NFDEID"

reps = [1,2,3]
samples = ['Brain','CoG','PO','TG','SG']
runs = 25+1

run_log = []
sample_log = []
rep_log = []
num_unique_log = []

for run in range(1,runs):
    for sample in samples:
        for rep in reps:
            results_path = results_dir + '\\' + 'final_results2021_0817_' + sample + '_' + str(rep) + '_DSD_run_' + str(run) + '_target.csv'
            results = pd.read_csv(results_path)
            unique = results.drop_duplicates(subset='Peptide')
            
            run_log.append(run)
            sample_log.append(sample)
            rep_log.append(rep)
            num_unique_log.append(len(unique))
            
dsd_results = pd.DataFrame()
dsd_results['Sample'] = sample_log
dsd_results['Replicate'] = rep_log
dsd_results['Run #'] = run_log
dsd_results['# Unique IDs'] = num_unique_log

dsd_results.to_csv((results_dir + '\\summary_of_runs.csv'), index=False)  