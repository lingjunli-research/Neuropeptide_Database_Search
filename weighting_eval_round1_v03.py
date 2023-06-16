# -*- coding: utf-8 -*-
"""
Created on Tue May 16 14:47:00 2023

@author: lawashburn
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
plt.style.use('seaborn-deep')

weighted_metrics_results_path = r"C:\Users\lawashburn\Documents\DBpep_v2\validation\weighting_extraction_results_20230518\metrics_out_v9_20230516.csv"
target_list_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\target_list.csv"
round1_dsd_path = r"C:\Users\lawashburn\Documents\DBpep_v2\Weighting_Assessment\Definitive Screening Design_v3_20230524_short.csv"
output_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\finale_weighting_hyperscore\PSM_unique_ID_validation"

round1_fdr = 0.01
power = 2

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
    num_unique_IDs = []
    for ind in round1_dsd.index:

        weighted_metrics_results_filtered = weighted_metrics_results[weighted_metrics_results['Sample'] == sample_ID]

        avg_frag_err = round1_dsd['average fragment error'][ind]
        prec_err = round1_dsd['precursor error'][ind]
        consec_b = round1_dsd['# consecutive b-ions'][ind]
        consec_y = round1_dsd['# consecutive y-ions'][ind]
        # c_term_amid = round1_dsd['C-terminal amidation'][ind]
        percent_seq_cov = round1_dsd['% sequence coverage'][ind]
        no_annot_peak = round1_dsd['average # annotations per peak'][ind]
        # percent_b_ions = round1_dsd['% ions are b'][ind]
        # percent_y_ions = round1_dsd['% ions are y'][ind]
        avg_matched_ions_per_AA = round1_dsd['average # matched ions per AA'][ind]
        avg_matched_ions_per_AA_non_neut = round1_dsd['# AAs annotated w/ non-neutral-loss ion'][ind]
        hyperscore = round1_dsd['hyperscore'][ind]
        
        avg_frag_err_norm_score = 1/(max(weighted_metrics_results_filtered['Average Fragment Error']))
        prec_err_norm_score = 1/(max(weighted_metrics_results_filtered['Precursor Error']))
        consec_b_ions_norm_score = 1/(max(weighted_metrics_results_filtered['# Consecutive b-ions']))
        consec_y_ions_norm_score = 1/(max(weighted_metrics_results_filtered['# Consecutive y-ions']))
        # precent_sec_cov_norm_score = 1/(max(weighted_metrics_results_filtered['% Sequence coverage']))
        avg_annotation_norm_score = 1/(max(weighted_metrics_results_filtered['Average annotations/fragment ']))
        avg_frag_ions_norm_score = 1/(max(weighted_metrics_results_filtered['Average number of fragment ions per AA']))
        avg_non_neut_frag_ions_norm_score = 1/(max(weighted_metrics_results_filtered['Number of non-neutral-loss fragment ions per AA']))
        hyperscore_norm_score = 1/(max(weighted_metrics_results_filtered['Hyperscore']))

        weighted_metrics_results_filtered['peptide_score'] = (
            (((weighted_metrics_results_filtered['Average Fragment Error']**power)*avg_frag_err_norm_score)*avg_frag_err)+
            (((weighted_metrics_results_filtered['Precursor Error']**power)*prec_err_norm_score)*prec_err)+
            (((weighted_metrics_results_filtered['# Consecutive b-ions']**power)*consec_b_ions_norm_score)*consec_b)+
            (((weighted_metrics_results_filtered['# Consecutive y-ions']**power)*consec_y_ions_norm_score)*consec_y)+
            # (((weighted_metrics_results_filtered['% Sequence coverage']**power)*precent_sec_cov_norm_score)*percent_seq_cov)+
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
        
        unique_peptide_score_thresh_filter = fdr_table[fdr_table['# Target IDs'] == best_num_targets]
        unique_peptide_score_thresh = unique_peptide_score_thresh_filter['Score Threshold'].min()
        
        filter_target_score_thresh = target_results[target_results['peptide_score'] >= unique_peptide_score_thresh]
        filter_target_score_thresh_unique = filter_target_score_thresh.drop_duplicates(subset=['Sample','Peptide'])
        
        run_log.append(ind+1)
        num_target_IDs_log.append(best_num_targets)
        num_unique_IDs.append(len(filter_target_score_thresh_unique))
        
        dsd_merge_table = pd.DataFrame()
        dsd_merge_table['Run #'] = run_log
        dsd_merge_table['# Target IDs: ' + sample_ID] = num_target_IDs_log
        dsd_merge_table['# Unique Target IDs:' + sample_ID] = num_unique_IDs
        
    file_path = output_directory + '\\dsd_eval_round1_' + sample_ID + '.csv'
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            dsd_merge_table.to_csv(filec,index=False)