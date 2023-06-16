# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 17:29:43 2023

@author: lawashburn
"""

import pandas as pd
import numpy as np

correlation_results_path = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\PEAKS_oursoftware_compare_brain_only\DB_search_optimize\03_precursor_AMM_optimize_v2\2021_0817_Brain_1\correlation_results.csv"
correlation_results = pd.read_csv(correlation_results_path)

standard_err_percent = 0.1

max_seq_cov = correlation_results['Sequence coverage'].max()
max_thresh = max_seq_cov * (1-standard_err_percent)

standard_err_subset = correlation_results[correlation_results['Sequence coverage'] >= max_thresh]
standard_err = standard_err_subset['Correlation value'].std()

correlation_results['count'] = correlation_results.groupby('Scan')['Scan'].transform(len)

final_psm = correlation_results[correlation_results['count'] == 1]

psm_candidate = correlation_results[correlation_results['count'] > 1]
psm_candidate = correlation_results[correlation_results['count'] > 1]

scan_candidate = psm_candidate.drop_duplicates(subset='Scan')
scan_candidate_list = scan_candidate['Scan'].values.tolist()

for scan in scan_candidate_list:
    scan_filtered_psm_candidate = psm_candidate[psm_candidate['Scan'] == scan]
    top_corr_score = scan_filtered_psm_candidate['Correlation value'].max()
    score_scan_filtered_psm_candidate = scan_filtered_psm_candidate[scan_filtered_psm_candidate['Correlation value'] >= (top_corr_score-standard_err)]
    if len(score_scan_filtered_psm_candidate) == 1:
        final_psm = pd.concat([final_psm,score_scan_filtered_psm_candidate])
    else:
        peptide_candidates = score_scan_filtered_psm_candidate['Sequence'].values.tolist()
        peptide_counts = []
        
        for value in peptide_candidates:
            instance_assess = final_psm[final_psm['Sequence'] == (value)]
            peptide_counts.append(len(instance_assess))
        
        peptide_min_instance = min(peptide_counts)
        count_min_instance = peptide_counts.count(peptide_min_instance)
        
        peptide_instance_assess = pd.DataFrame()
        peptide_instance_assess['Sequence'] = peptide_candidates
        peptide_instance_assess['PSM Count'] = peptide_counts
        
        if count_min_instance == 0:
            sequence_to_report = peptide_instance_assess[peptide_instance_assess['PSM Count'] == (top_corr_score-count_min_instance)]
            sequence_psm = sequence_to_report['Sequence'].values[0]
            
            psm_entry = score_scan_filtered_psm_candidate[score_scan_filtered_psm_candidate['Sequence'] == (sequence_psm)]
            final_psm = pd.concat([final_psm,psm_entry])
            
        if count_min_instance > 0:
            score_scan_filtered_psm_candidate = score_scan_filtered_psm_candidate.sort_values(by='Correlation value',ascending=False)
            psm_entry = score_scan_filtered_psm_candidate.head(1)
            final_psm = pd.concat([final_psm,psm_entry])
#%%
final_psm_no_dups = final_psm.drop_duplicates(subset='Sequence')