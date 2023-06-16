# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 17:29:43 2023

@author: lawashburn
"""

import pandas as pd
import numpy as np

correlation_results_path = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\PEAKS_oursoftware_compare_brain_only\DB_search_optimize\03_precursor_AMM_optimize_v2\2021_0817_Brain_1\correlation_results.csv"
all_correlation_results = pd.read_csv(correlation_results_path)

standard_err_percent = 0.1

max_seq_cov = all_correlation_results['Sequence coverage'].max()
max_thresh = max_seq_cov * (1-standard_err_percent)

standard_err_subset = all_correlation_results[all_correlation_results['Sequence coverage'] >= max_thresh]
standard_err = standard_err_subset['Correlation value'].std()

all_correlation_results['count'] = all_correlation_results.groupby('Scan')['Scan'].transform(len)

final_psm = all_correlation_results[all_correlation_results['count'] == 1]

psm_candidate = all_correlation_results[all_correlation_results['count'] > 1]


scan_candidate = psm_candidate.drop_duplicates(subset='Scan')
scan_candidate_list = scan_candidate['Scan'].values.tolist()

future_storage = []

for scan in scan_candidate_list:
    final_psm = final_psm.drop_duplicates()
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

        
        if count_min_instance == 1:
            sequence_to_report = peptide_instance_assess[peptide_instance_assess['PSM Count'] == peptide_min_instance]
            sequence_psm = sequence_to_report['Sequence'].values[0]
            
            psm_entry = score_scan_filtered_psm_candidate[score_scan_filtered_psm_candidate['Sequence'] == (sequence_psm)]
            final_psm = pd.concat([final_psm,psm_entry])


        elif count_min_instance > 1:

            sequence_to_report = peptide_instance_assess[peptide_instance_assess['PSM Count'] == peptide_min_instance]
            sequence_to_report_list = sequence_to_report['Sequence'].values.tolist()
            filtered_scans = score_scan_filtered_psm_candidate[score_scan_filtered_psm_candidate['Sequence'].isin(sequence_to_report_list)]
            score_score_scan_filtered_psm_candidate = filtered_scans[filtered_scans['Correlation value'] > 0]
            if len(score_score_scan_filtered_psm_candidate) == 1:
                final_psm = pd.concat([final_psm,score_score_scan_filtered_psm_candidate])
            if len(score_score_scan_filtered_psm_candidate) == 0:
                max_seq_cov = filtered_scans['Sequence coverage'].max()
                cov_score_scan_filtered_psm_candidate = filtered_scans[filtered_scans['Sequence coverage'] == max_seq_cov]
                if len(cov_score_scan_filtered_psm_candidate) == 1:
                    final_psm = pd.concat([final_psm,cov_score_scan_filtered_psm_candidate])
                if len(cov_score_scan_filtered_psm_candidate) > 1:
                    future_storage.append(cov_score_scan_filtered_psm_candidate)
            if len(score_score_scan_filtered_psm_candidate) > 1:
                max_seq_cov = score_score_scan_filtered_psm_candidate['Sequence coverage'].max()
                cov_score_scan_filtered_psm_candidate = score_score_scan_filtered_psm_candidate[score_score_scan_filtered_psm_candidate['Sequence coverage'] == max_seq_cov]
                if len(cov_score_scan_filtered_psm_candidate) == 1:
                    final_psm = pd.concat([final_psm,cov_score_scan_filtered_psm_candidate])
                if len(cov_score_scan_filtered_psm_candidate) > 1:
                    future_storage.append(cov_score_scan_filtered_psm_candidate)

first_round_psm_no_dups = final_psm.drop_duplicates(subset='Sequence')
print('Round: 1')
print('# PSMs: ',len(final_psm))
print('# Unique IDs',len(first_round_psm_no_dups))

max_round = 5
for rounds in range(0,max_round):
    if len(future_storage) > 0:
        future_storage_df = pd.concat(future_storage,ignore_index=True)   
        future_storage.clear()
        
        repeat_scan_candidate = future_storage_df.drop_duplicates(subset='Scan')
        repeat_scan_candidate_list = repeat_scan_candidate['Scan'].values.tolist()
        
        for repeat_scan in repeat_scan_candidate_list:
            final_psm = final_psm.drop_duplicates()
            repeat_scan_filtered_psm_candidate = repeat_scan_candidate[repeat_scan_candidate['Scan'] == repeat_scan]
            repeat_top_corr_score = repeat_scan_filtered_psm_candidate['Correlation value'].max()
            repeat_score_scan_filtered_psm_candidate = repeat_scan_filtered_psm_candidate[repeat_scan_filtered_psm_candidate['Correlation value'] >= (repeat_top_corr_score-standard_err)]
            if len(repeat_score_scan_filtered_psm_candidate) == 1:
                final_psm = pd.concat([final_psm,repeat_score_scan_filtered_psm_candidate])
            else:
                repeat_peptide_candidates = repeat_score_scan_filtered_psm_candidate['Sequence'].values.tolist()
                repeat_peptide_counts = []
                
                for repeat_value in repeat_peptide_candidates:
                    repeat_instance_assess = final_psm[final_psm['Sequence'] == (repeat_value)]
                    repeat_peptide_counts.append(len(repeat_instance_assess))
                
                repeat_peptide_min_instance = min(repeat_peptide_counts)
                repeat_count_min_instance = repeat_peptide_counts.count(repeat_peptide_min_instance)
                
                repeat_peptide_instance_assess = pd.DataFrame()
                repeat_peptide_instance_assess['Sequence'] = repeat_peptide_candidates
                repeat_peptide_instance_assess['PSM Count'] = repeat_peptide_counts

                
                if repeat_count_min_instance == 1:
                    repeat_sequence_to_report = repeat_peptide_instance_assess[repeat_peptide_instance_assess['PSM Count'] == repeat_peptide_min_instance]
                    repeat_sequence_psm = repeat_sequence_to_report['Sequence'].values[0]
                    
                    repeat_psm_entry = repeat_score_scan_filtered_psm_candidate[repeat_score_scan_filtered_psm_candidate['Sequence'] == (repeat_sequence_psm)]
                    final_psm = pd.concat([final_psm,repeat_psm_entry])


                elif repeat_count_min_instance > 1:

                    repeat_sequence_to_report = repeat_peptide_instance_assess[repeat_peptide_instance_assess['PSM Count'] == repeat_peptide_min_instance]
                    repeat_sequence_to_report_list = repeat_sequence_to_report['Sequence'].values.tolist()
                    repeat_filtered_scans = repeat_score_scan_filtered_psm_candidate[repeat_score_scan_filtered_psm_candidate['Sequence'].isin(repeat_sequence_to_report_list)]
                    repeat_score_score_scan_filtered_psm_candidate = repeat_filtered_scans[repeat_filtered_scans['Correlation value'] > 0]
                    if len(repeat_score_score_scan_filtered_psm_candidate) == 1:
                        final_psm = pd.concat([final_psm,repeat_score_score_scan_filtered_psm_candidate])
                    if len(repeat_score_score_scan_filtered_psm_candidate) == 0:
                        repeat_max_seq_cov = repeat_filtered_scans['Sequence coverage'].max()
                        repeat_cov_score_scan_filtered_psm_candidate = repeat_filtered_scans[repeat_filtered_scans['Sequence coverage'] == repeat_max_seq_cov]
                        if len(repeat_cov_score_scan_filtered_psm_candidate) == 1:
                            final_psm = pd.concat([final_psm,repeat_cov_score_scan_filtered_psm_candidate])
                        if len(repeat_cov_score_scan_filtered_psm_candidate) > 1:
                            future_storage.append(repeat_cov_score_scan_filtered_psm_candidate)
                    if len(repeat_score_score_scan_filtered_psm_candidate) > 1:
                        repeat_max_seq_cov = repeat_score_score_scan_filtered_psm_candidate['Sequence coverage'].max()
                        repeat_cov_score_scan_filtered_psm_candidate = repeat_score_score_scan_filtered_psm_candidate[repeat_score_score_scan_filtered_psm_candidate['Sequence coverage'] == repeat_max_seq_cov]
                        if len(repeat_cov_score_scan_filtered_psm_candidate) == 1:
                            final_psm = pd.concat([final_psm,repeat_cov_score_scan_filtered_psm_candidate])
                        if len(repeat_cov_score_scan_filtered_psm_candidate) > 1:
                            if rounds == max_round:
                                repeat_cov_score_scan_filtered_psm_candidate_select = repeat_cov_score_scan_filtered_psm_candidate.sample(n=1)
                                final_psm = pd.concat([final_psm,repeat_cov_score_scan_filtered_psm_candidate_select])
                            else:
                                future_storage.append(repeat_cov_score_scan_filtered_psm_candidate)
        prelim_final_psm_no_dups = final_psm.drop_duplicates(subset='Sequence')
        print('Round: ',(rounds+2))
        print('# PSMs: ',len(final_psm))
        print('# Unique IDs',len(prelim_final_psm_no_dups))
    
    else:
        pass
    
final_psm_no_dups = final_psm.drop_duplicates(subset='Sequence')
print('Total:')
print('# PSMs: ',len(final_psm))
print('# Unique IDs',len(final_psm_no_dups))

