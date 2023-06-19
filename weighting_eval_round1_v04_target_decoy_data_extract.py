# -*- coding: utf-8 -*-
"""
Created on Tue May 16 14:52:42 2023

@author: lawashburn
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
plt.style.use('seaborn-deep')

weighted_metrics_results_path = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\PEAKS_oursoftware_compare_brain_only\DB_search_optimize\03_precursor_AMM_optimize_v2\2021_0817_Brain_1\metrics_extracted.csv"
target_list_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\target_list.csv"
round1_dsd_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale_weighting_hyperscore\weighting_dsd_eval\round1_DSD_0_10_w_hyperscore.csv"
output_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\PEAKS_oursoftware_compare_brain_only\DB_search_optimize\03_precursor_AMM_optimize_v2\2021_0817_Brain_1"

round1_fdr = 0.01
power = 2

target_list = pd.read_csv(target_list_path)
target_list['Sequence'] = target_list['Sequence'].str.replace(r"\([^)]*\)","")
target_list_str = target_list['Sequence'].values.tolist()

weighted_metrics_results = pd.read_csv(weighted_metrics_results_path)
weighted_metrics_results['Unmodified sequence'] = weighted_metrics_results['Peptide']
weighted_metrics_results['Unmodified sequence'] = weighted_metrics_results['Unmodified sequence'].str.replace(r"\([^)]*\)","")
weighted_metrics_results['Status'] = weighted_metrics_results['Unmodified sequence'].apply(lambda x: any([k in x for k in target_list_str]))

sample_list = weighted_metrics_results['Sample'].values.tolist()
sample_nodups = []
for sample in sample_list:
    if sample not in sample_nodups:
        sample_nodups.append(sample)


#DSD Values (from run 31)
avg_frag_err = 0
prec_err = 0
consec_b = 10
consec_y = 0
percent_seq_cov = 0
no_annot_peak = 10
avg_matched_ions_per_AA = 10
avg_matched_ions_per_AA_non_neut = 10
hyperscore = 10

for sample_ID in sample_nodups:
    weighted_metrics_results_filtered = weighted_metrics_results[weighted_metrics_results['Sample'] == sample_ID]
    avg_frag_err_norm_score = 1/(max(weighted_metrics_results_filtered['Average Fragment Error']))
    prec_err_norm_score = 1/(max(weighted_metrics_results_filtered['Precursor Error']))
    consec_b_ions_norm_score = 1/(max(weighted_metrics_results_filtered['# Consecutive b-ions']))
    consec_y_ions_norm_score = 1/(max(weighted_metrics_results_filtered['# Consecutive y-ions']))
    avg_annotation_norm_score = 1/(max(weighted_metrics_results_filtered['Average annotations/fragment ']))
    avg_frag_ions_norm_score = 1/(max(weighted_metrics_results_filtered['Average number of fragment ions per AA']))
    avg_non_neut_frag_ions_norm_score = 1/(max(weighted_metrics_results_filtered['Number of non-neutral-loss fragment ions per AA']))
    hyperscore_norm_score = 1/(max(weighted_metrics_results_filtered['Hyperscore']))
    precent_sec_cov_norm_score = 1/(max(weighted_metrics_results_filtered['Motif_Score']))
    
    weighted_metrics_results_filtered['peptide_score'] = (
        (((weighted_metrics_results_filtered['Average Fragment Error']**power)*avg_frag_err_norm_score)*avg_frag_err)+
        (((weighted_metrics_results_filtered['Precursor Error']**power)*prec_err_norm_score)*prec_err)+
        (((weighted_metrics_results_filtered['# Consecutive b-ions']**power)*consec_b_ions_norm_score)*consec_b)+
        (((weighted_metrics_results_filtered['# Consecutive y-ions']**power)*consec_y_ions_norm_score)*consec_y)+
        (((weighted_metrics_results_filtered['Motif_Score']**power)*precent_sec_cov_norm_score)*percent_seq_cov)+
        (((weighted_metrics_results_filtered['Average annotations/fragment ']**power)*avg_annotation_norm_score)*no_annot_peak)+
        (((weighted_metrics_results_filtered['Average number of fragment ions per AA']**power)*avg_frag_ions_norm_score)*avg_matched_ions_per_AA)+
        (((weighted_metrics_results_filtered['Number of non-neutral-loss fragment ions per AA']**power)*avg_non_neut_frag_ions_norm_score)*avg_matched_ions_per_AA_non_neut)+
        (((weighted_metrics_results_filtered['Hyperscore']**power)*hyperscore_norm_score)*hyperscore)
        )
    
    target_results = weighted_metrics_results_filtered[weighted_metrics_results_filtered['Status'] == True]
    decoy_results = weighted_metrics_results_filtered[weighted_metrics_results_filtered['Status'] == False]
    
    target_scores = target_results['peptide_score'].values.tolist()
    decoy_scores = decoy_results['peptide_score'].values.tolist()
    
    all_scores = weighted_metrics_results_filtered['peptide_score'].values.tolist()
    min_scores = 0
    max_scores = (max(all_scores)) +10
    
    bins = np.linspace(min_scores, max_scores, 50)
            
    int_list_no_dups = []
    
    for aa in target_scores:
        aaa = int(aa)
        if aaa not in int_list_no_dups:
            int_list_no_dups.append(aaa)
    
    
    score_list_sorted = sorted(int_list_no_dups)
    
    score_cutoff = []
    fdr = []
    num_target_IDs = []
    num_decoy_IDs = []
    
    for score in score_list_sorted:
        target_results_filtered = target_results[target_results['peptide_score'] >= score]
        decoy_results_filtered = decoy_results[decoy_results['peptide_score'] >= score]
        
        fdr_report = len(decoy_results_filtered)/len(target_results_filtered)
        
        score_cutoff.append(score)
        fdr.append(fdr_report)
        num_target_IDs.append(len(target_results_filtered))
        num_decoy_IDs.append(len(decoy_results_filtered))
    
    fdr_table = pd.DataFrame()
    fdr_table['Score Threshold'] = score_cutoff
    fdr_table['FDR'] = fdr
    fdr_table['# Target IDs'] = num_target_IDs
    fdr_table['# Decoy IDs'] = num_decoy_IDs
    
    file_path = output_directory + '\\FDR_eval_table' + sample_ID + '.csv'
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            fdr_table.to_csv(filec,index=False)
    
    fdr_results_filtered = fdr_table[fdr_table['FDR'] <= round1_fdr]
    
    best_num_targets = fdr_results_filtered['# Target IDs'].max()
    
    fdr_table_filtered = fdr_table[fdr_table['# Target IDs'] == best_num_targets]
    score_threshold_to_apply = fdr_table_filtered['Score Threshold'].min()
    
    target_results_final = target_results[target_results['peptide_score'] >= score_threshold_to_apply]
    decoy_results_final = decoy_results[decoy_results['peptide_score'] >= score_threshold_to_apply]
    
    file_path = output_directory + '\\final_results' + sample_ID + '_decoy.csv'
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            decoy_results_final.to_csv(filec,index=False)
    
    file_path = output_directory + '\\final_results' + sample_ID + '_target.csv'
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            target_results_final.to_csv(filec,index=False)