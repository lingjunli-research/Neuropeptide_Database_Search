# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 17:48:07 2023

@author: lawashburn
"""
import pandas as pd
import re
import numpy as np
import csv

psm_assignments_path = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\PEAKS_oursoftware_compare_brain_only\DB_search_optimize\03_precursor_AMM_optimize_v2\2021_0817_Brain_1\psm_assignments.csv"
motif_database_path = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\PEAKS_oursoftware_compare_brain_only\DB_search_optimize\03_precursor_AMM_optimize_v2\2021_0817_Brain_1\motif_db_20230616.csv"
fragment_match_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\PEAKS_oursoftware_compare_brain_only\DB_search_optimize\03_precursor_AMM_optimize_v2\2021_0817_Brain_1\fragment_matches"
output_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\PEAKS_oursoftware_compare_brain_only\DB_search_optimize\03_precursor_AMM_optimize_v2\2021_0817_Brain_1"

psm_assignments = pd.read_csv(psm_assignments_path)
psm_assignments['Sequence format'] = psm_assignments['Sequence'].str.replace(r"\(.*\)","")

motif_database = pd.read_csv(motif_database_path)

motif_sequences = motif_database['Sequence'].values.tolist()

motif_match_motif_storage = []
motif_match_sequence_storage = []
motif_match_sequence_modded_storage = []
motif_match_scan_storage = []
motif_match_coverage_storage = []

for psm in range(0,len(psm_assignments)):
    sequence_modded = psm_assignments['Sequence'].values[psm]
    sequence = psm_assignments['Sequence format'].values[psm]
    scan = psm_assignments['Scan'].values[psm]
    for motif in motif_sequences:
        if motif in sequence:
            motif_len = len(motif)
            seq_len = len(sequence)
            motif_index = [match.span() for match in re.finditer(motif, sequence)]
            for index in motif_index:
                b_ion_list = []
                y_ion_list = []
                
                start, end = index
                
                for val in range(start,end):
                    b_ion_list.append('b' + str(val+1))
                    y_ion_list.append('y' + str(seq_len - val))

                exp_ions_path = fragment_match_directory + '\\' + sequence + '_' + str(scan) + '_fragment_report.csv'
                exp_ions = pd.read_csv(exp_ions_path)
                
                theo_ion_list = b_ion_list + y_ion_list

                exp_b_y = exp_ions['ion'].values.tolist() 
                exp_b_y_noPTM = [] 
                for c in exp_b_y:
                    c = c.replace('-H2O','')
                    c = c.replace('-NH3','')
                    if c not in exp_b_y_noPTM:
                        if c in theo_ion_list:
                            exp_b_y_noPTM.append(c)
                        else:
                            pass
                    else:
                        pass

                motif_seq_cov_calc = len(exp_b_y_noPTM) / len(theo_ion_list)

                motif_match_motif_storage.append(motif)
                motif_match_sequence_storage.append(sequence)
                motif_match_sequence_modded_storage.append(sequence_modded)
                motif_match_scan_storage.append(scan)
                motif_match_coverage_storage.append(motif_seq_cov_calc)

            
motif_matches = pd.DataFrame()
motif_matches['Sequence format'] = motif_match_sequence_storage
motif_matches['Sequence'] = motif_match_sequence_modded_storage
motif_matches['Matched motif'] = motif_match_motif_storage
motif_matches['Scan'] = motif_match_scan_storage
motif_matches['Motif Sequence Coverage'] = motif_match_coverage_storage
motif_matches['Motif Score (1)'] = motif_matches.apply(lambda row: len(row['Sequence']) / len(motif_matches['Matched motif']), axis=1)
motif_matches['Sequence Length'] = motif_matches['Sequence'].str.len()
motif_matches['Norm. Motif Score (2)'] = motif_matches['Motif Score (1)'] * np.sqrt(motif_matches['Sequence Length'])
motif_matches['Final motif score'] = motif_matches['Norm. Motif Score (2)']*(motif_matches['Motif Sequence Coverage'])

psm_w_motif_report = pd.merge(psm_assignments, motif_matches, on=['Sequence','Scan','Sequence format'])
psm_w_motif_report = psm_w_motif_report.sort_values(by='Final motif score', ascending=False)
psm_w_motif_report = psm_w_motif_report.drop_duplicates(subset='Scan')

file_path = output_directory + '\\psm_results_w_motif.csv'
with open(file_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        psm_w_motif_report.to_csv(filec,index=False)

