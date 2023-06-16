# -*- coding: utf-8 -*-
"""
Created on Mon May  8 17:09:44 2023

@author: lawashburn
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
plt.style.use('seaborn-deep')

weighted_metrics_results_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale_weighting_hyperscore\weighting_metric_extraction\metrics_out_v7_20230515.csv"
target_list_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\target_list.csv"
round1_dsd_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale_weighting_hyperscore\weighting_dsd_eval\round1_DSD_0_100_w_hyperscore.csv"
output_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\finale_weighting_hyperscore\dsd_output\1_100_1_FDR_hyperscore_v03_20230515"

round1_fdr = 0.01

weighted_metrics_results = pd.read_csv(weighted_metrics_results_path)
round1_dsd = pd.read_csv(round1_dsd_path)

target_list = pd.read_csv(target_list_path)
target_list['Sequence'] = target_list['Sequence'].str.replace(r"\([^)]*\)","")
target_list_str = target_list['Sequence'].values.tolist()

weighted_metrics_results['Unmodified sequence'] = weighted_metrics_results['Peptide']
weighted_metrics_results['Unmodified sequence'] = weighted_metrics_results['Unmodified sequence'].str.replace(r"\([^)]*\)","")
weighted_metrics_results['Status'] = weighted_metrics_results['Unmodified sequence'].apply(lambda x: any([k in x for k in target_list_str]))

sample_list = weighted_metrics_results['Sample'].values.tolist()

sample_nodups = []

for sample in sample_list:
    if sample not in sample_nodups:
        sample_nodups.append(sample)



run_numbers = round1_dsd['Run'].values.tolist()
dsd_eval = pd.DataFrame()
dsd_eval['Run #'] = run_numbers

for sample_ID in sample_nodups:
    run_log = []
    num_target_IDs_log = []
    for ind in round1_dsd.index:

        weighted_metrics_results_filtered = weighted_metrics_results[weighted_metrics_results['Sample'] == sample_ID]
        
        avg_frag_err = round1_dsd['average fragment error'][ind]
        prec_err = round1_dsd['precursor error'][ind]
        consec_b = round1_dsd['# consecutive b-ions'][ind]
        consec_y = round1_dsd['# consecutive y-ions'][ind]
        c_term_amid = round1_dsd['C-terminal amidation'][ind]
        percent_seq_cov = round1_dsd['% sequence coverage'][ind]
        no_annot_peak = round1_dsd['average # annotations per peak'][ind]
        percent_b_ions = round1_dsd['% ions are b'][ind]
        percent_y_ions = round1_dsd['% ions are y'][ind]
        avg_matched_ions_per_AA = round1_dsd['average # matched ions per AA'][ind]
        avg_matched_ions_per_AA_non_neut = round1_dsd['# AAs annotated w/ non-neutral-loss ion'][ind]
        hyperscore = round1_dsd['hyperscore'][ind]
        
        weighted_metrics_results_filtered['peptide_score'] = (
            (weighted_metrics_results_filtered['Average Fragment Error']*avg_frag_err)+
            (weighted_metrics_results_filtered['Precursor Error']*prec_err)+
            (weighted_metrics_results_filtered['# Consecutive b-ions']*consec_b)+
            (weighted_metrics_results_filtered['# Consecutive y-ions']*consec_y)+
            (weighted_metrics_results_filtered['C-termini amidation']*c_term_amid)+
            (weighted_metrics_results_filtered['% Sequence coverage']*percent_seq_cov)+
            (weighted_metrics_results_filtered['Average annotations/fragment ']*no_annot_peak)+
            (weighted_metrics_results_filtered['% Fragment ions are b']*percent_b_ions)+
            (weighted_metrics_results_filtered['% Fragment ions are y']*percent_y_ions)+
            (weighted_metrics_results_filtered['Average number of fragment ions per AA']*avg_matched_ions_per_AA)+
            (weighted_metrics_results_filtered['Number of non-neutral-loss fragment ions per AA']*avg_matched_ions_per_AA_non_neut)+
            (weighted_metrics_results_filtered['Hyperscore']*hyperscore)
            )
        
        target_results = weighted_metrics_results_filtered[weighted_metrics_results_filtered['Status'] == True]
        decoy_results = weighted_metrics_results_filtered[weighted_metrics_results_filtered['Status'] == False]
        
        target_scores = target_results['peptide_score'].values.tolist()
        decoy_scores = decoy_results['peptide_score'].values.tolist()
        
        all_scores = weighted_metrics_results_filtered['peptide_score'].values.tolist()
        min_scores = 0
        max_scores = (max(all_scores)) +10
        
        bins = np.linspace(min_scores, max_scores, 50)
        
        plt.hist([target_scores, decoy_scores], bins, label=['Target', 'Decoy'])
        plt.legend(loc='upper right')
        plt.show()
        
        int_list_no_dups = []
        
        for aa in target_scores:
            aaa = int(aa)
            if aaa not in int_list_no_dups:
                int_list_no_dups.append(aaa)
        
        
        score_list_sorted = sorted(int_list_no_dups)
        
        score_cutoff = []
        fdr = []
        num_target_IDs = []
        
        for score in score_list_sorted:
            target_results_filtered = target_results[target_results['peptide_score'] >= score]
            decoy_results_filtered = decoy_results[decoy_results['peptide_score'] >= score]
            
            fdr_report = len(decoy_results_filtered)/len(target_results_filtered)
            
            score_cutoff.append(score)
            fdr.append(fdr_report)
            num_target_IDs.append(len(target_results_filtered))
        
        fdr_table = pd.DataFrame()
        fdr_table['Score Threshold'] = score_cutoff
        fdr_table['FDR'] = fdr
        fdr_table['# Target IDs'] = num_target_IDs
        
        fdr_results_filtered = fdr_table[fdr_table['FDR'] <= round1_fdr]
        
        best_num_targets = fdr_results_filtered['# Target IDs'].max()
        
        run_log.append(ind+1)
        num_target_IDs_log.append(best_num_targets)
        
        dsd_merge_table = pd.DataFrame()
        dsd_merge_table['Run #'] = run_log
        dsd_merge_table['# Target IDs: ' + sample_ID] = num_target_IDs_log
        
    file_path = output_directory + '\\dsd_eval_round1_' + sample_ID + '.csv'
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            dsd_merge_table.to_csv(filec,index=False)
