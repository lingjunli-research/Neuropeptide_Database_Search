# -*- coding: utf-8 -*-
"""
Created on Wed May 31 10:12:17 2023

@author: lawashburn
"""

import pandas as pd
import csv
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt

peaks_results_directory = r"H:\PEAKS_exports\individual_tissue_rep_results"
db_search_inital_output = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\PEAKS_oursoftware_compare_brain_only\DB_search_optimize\01_baseline"
FDR_filtered_db_search_results = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\DSD_eval\motif_2"
np_database_path = r"C:\Users\lawashburn\Desktop\ALC50_Mass_Search_Files\duplicate_removed_crustacean_database_validated_formatted20220725.csv"
output_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\PEAKS_oursoftware_compare_brain_only\DB_search_optimize\01_baseline\software_compare"
np_database = pd.read_csv(np_database_path)
samples = ['Brain_1']
#samples = ['Brain_1','Brain_2','Brain_3','CoG_1','CoG_2','CoG_3','PO_1','PO_2','PO_3','SG_1','SG_2','SG_3','TG_1','TG_2','TG_3']
#samples = ['Brain_1']

unique_to_peaks_all = []
unique_to_peaks_all_sample = []
unique_to_peaks_from_FDR_filtered = []
unique_to_peaks_from_FDR_filtered_sample = []

for sample in samples:
    sample_format = sample.replace('_',"")
    peaks_results_file_path = peaks_results_directory + '\\' + sample_format + '\\protein-peptides.csv'
    peaks_results = pd.read_csv(peaks_results_file_path)
    peaks_results_filtered = peaks_results[peaks_results['Unique'] == 'Y']
    peaks_unique_peptide_filtered = peaks_results_filtered.drop_duplicates(subset='Protein Accession')
    peaks_unique_peptides = peaks_unique_peptide_filtered['Protein Accession'].values.tolist() #PEAKS results

    initial_db_search_results_path = db_search_inital_output + '\\2021_0817_' + sample + '\\results_with_correlation.csv'
    initial_db_search_results = pd.read_csv(initial_db_search_results_path)
    initial_db_search_results['Sequence (no mods)'] = initial_db_search_results['Sequence'].str.replace(r"\([^)]*\)","")
    initial_db_search_results_filtered = initial_db_search_results.drop_duplicates(subset='Sequence (no mods)')
    unfiltered_unqiue_peptide_results = initial_db_search_results_filtered['Sequence (no mods)'].values.tolist() #Results not FDR filtered

    sample_type = sample[:(len(sample)-2)]
    sample_rep = sample[(len(sample)-1):] 
    
    all_FDR_filtered_results_path = FDR_filtered_db_search_results + '\\summary_of_runs.csv'
    all_FDR_filtered_results = pd.read_csv(all_FDR_filtered_results_path)
    all_FDR_filtered_results_sample_filtered = all_FDR_filtered_results[all_FDR_filtered_results['Sample'] == sample_type]
    all_FDR_filtered_results_sample_rep_filtered = all_FDR_filtered_results_sample_filtered[all_FDR_filtered_results_sample_filtered['Replicate'] == int(sample_rep)]
    best_ID_count = max(all_FDR_filtered_results_sample_rep_filtered['# Unique IDs'])
    best_FDR_filtered_results_sample_rep_filtered = all_FDR_filtered_results_sample_rep_filtered[all_FDR_filtered_results_sample_rep_filtered['# Unique IDs'] == best_ID_count]
    best_FDR_filtered_results_sample_rep_filtered = best_FDR_filtered_results_sample_rep_filtered.reset_index(drop=True)
    target_DSD = best_FDR_filtered_results_sample_rep_filtered['Run #'][0]

    fdr_filtered_peptide_results_path  = FDR_filtered_db_search_results + '\\final_results2021_0817_' + sample + "_DSD_run_" + str(target_DSD) + "_target.csv"
    fdr_filtered_peptide_results = pd.read_csv(fdr_filtered_peptide_results_path)
    fdr_filtered_peptide_results = fdr_filtered_peptide_results.drop_duplicates(subset='Unmodified sequence')
    fdr_filtered_unqiue_peptide_results = fdr_filtered_peptide_results['Unmodified sequence'].values.tolist() #Results at 1% FDR
    
    
    db_seqs = np_database['Sequence'].values.tolist()
    
    results_unfiltered_accession = [] #Accession of unfiltered results
    for a in unfiltered_unqiue_peptide_results:
        if a in db_seqs: # need to filter because the unfiltered results will include decoys not in DB
            np_database_filtered = np_database[np_database['Sequence'] == a]
            np_database_filtered = np_database_filtered.reset_index(drop=True)
            accession_to_report = np_database_filtered['Accession'][0]
            results_unfiltered_accession.append(int(accession_to_report))
        else:
            pass

    results_FDR_filtered_accession = [] #Accession of results @ 1% FDR
    for b in fdr_filtered_unqiue_peptide_results:
        if b in db_seqs:
            np_database_filtered = np_database[np_database['Sequence'] == b]
            np_database_filtered = np_database_filtered.reset_index(drop=True)
            accession_to_report = np_database_filtered['Accession'][0]
            results_FDR_filtered_accession.append(int(accession_to_report))
        else:
            pass
    

    
    #make a Venn Diagram
    #sample A = PEAKS results; sample B = unfiltered results; sample C = results @ 1% FDR
    
    AB_overlap = set(peaks_unique_peptides) & set(results_unfiltered_accession)  #compute intersection of set A & set B
    AC_overlap = set(peaks_unique_peptides) & set(results_FDR_filtered_accession)
    BC_overlap = set(results_unfiltered_accession) & set(results_FDR_filtered_accession)
    ABC_overlap = set(peaks_unique_peptides) & set(results_unfiltered_accession) & set(results_FDR_filtered_accession)
    A_rest = set(peaks_unique_peptides) - AB_overlap - AC_overlap #see left graphic
    B_rest = set(results_unfiltered_accession) - AB_overlap - BC_overlap
    C_rest = set(results_FDR_filtered_accession) - AC_overlap - BC_overlap
    AB_only = AB_overlap - ABC_overlap   #see right graphic
    AC_only = AC_overlap - ABC_overlap
    BC_only = BC_overlap - ABC_overlap
    
    sets = Counter()               #set order A, B, C   
    sets['100'] = len(A_rest)      #100 denotes A on, B off, C off 
    sets['010'] = len(B_rest)      #010 denotes A off, B on, C off
    sets['001'] = len(C_rest)      #001 denotes A off, B off, C on 
    sets['110'] = len(AB_only)     #110 denotes A on, B on, C off
    sets['101'] = len(AC_only)     #101 denotes A on, B off, C on 
    sets['011'] = len(BC_only)     #011 denotes A off, B on, C on 
    sets['111'] = len(ABC_overlap) #011 denotes A on, B on, C on
    labels = ('PEAKS', 'Ours: unfiltered', 'Ours: 1% FDR')  
    plt.figure(figsize=(7,7)) 
    ax = plt.gca() 
    v = venn3(subsets=sets, set_labels=labels, ax=ax,set_colors=    
          ('#001244','#318fb5','#f7d6bf'),alpha=0.7)    
    venn3_circles(subsets=sets, linewidth=1)
    
    for text in v.set_labels:
          text.set_font('Arial')
          text.set_fontsize(16)
          text.set_fontweight('bold')
 
    for text in v.subset_labels:
        if text is None:
            pass
        else:
             text.set_font('Arial')
             text.set_fontsize(16)
             text.set_fontweight('bold')

    plt_title = "ID comparison for sample: " + sample
    plt.title(plt_title, fontsize=24, fontweight='bold',font='Arial')
    plt.tight_layout()
    #plt.show()

    plt_out_path = output_directory + '\\PEAKS_our_software_compare_3_way_' + sample + '.png'
    plt.savefig(plt_out_path)
    plt.clf()
    
    #2 way venn only

    AC_overlap_2 = set(peaks_unique_peptides) & set(results_FDR_filtered_accession)
    A_rest_2 = set(peaks_unique_peptides) - AC_overlap_2 #see left graphic
    C_rest_2 = set(results_FDR_filtered_accession) - AC_overlap_2

    
    sets2 = Counter()               #set order A, B, C   
    sets2['10'] = len(A_rest_2)      #100 denotes A on, B off, C off 
    sets2['01'] = len(C_rest_2)      #001 denotes A off, B off, C on 
    sets2['11'] = len(AC_overlap_2)     #101 denotes A on, B off, C on 
    labels2 = ('PEAKS', 'Ours: 1% FDR')  
    plt.figure(figsize=(7,7)) 
    ax = plt.gca() 
    v2 = venn2(subsets=sets2, set_labels=labels2, ax=ax,set_colors=    
          ('#001244','#318fb5'),alpha=0.7)    
    venn2_circles(subsets=sets2, linewidth=1)
    
    for text in v2.set_labels:
          text.set_font('Arial')
          text.set_fontsize(16)
          text.set_fontweight('bold')
 
    for text in v2.subset_labels:
        if text is None:
            pass
        else:
              text.set_font('Arial')
              text.set_fontsize(16)
              text.set_fontweight('bold')

    plt_title = "ID comparison for sample: " + sample
    plt.title(plt_title, fontsize=24, fontweight='bold',font='Arial')
    plt.tight_layout()
    #plt.show()
    
    plt_out_path = output_directory + '\\PEAKS_our_software_compare_2_way_' + sample + '.png'
    plt.savefig(plt_out_path)
    plt.clf()

    
    for c in peaks_unique_peptides:
        if c in results_FDR_filtered_accession:
            pass
        else:
            unique_to_peaks_from_FDR_filtered.append(c)
            unique_to_peaks_from_FDR_filtered_sample.append(sample)
            
    
    for d in peaks_unique_peptides:
        if d in results_unfiltered_accession:
            pass
        else:
            unique_to_peaks_all.append(d)
            unique_to_peaks_all_sample.append(sample)


unique_to_peaks_all_no_dups = []
for j in unique_to_peaks_all:
    if j not in unique_to_peaks_all_no_dups:
        unique_to_peaks_all_no_dups.append(j)
    else:
        pass



sample_retrieve_df_all = pd.DataFrame()
sample_retrieve_df_all['Accession'] = unique_to_peaks_all
sample_retrieve_df_all['Sample'] = unique_to_peaks_all_sample


all_results_storage_accession = []
all_results_storage_peptide = []
all_results_storage_instances = []
all_results_storage_samples = []

for k in unique_to_peaks_all_no_dups:
    instances = unique_to_peaks_all.count(k)
    np_database_filtered = np_database[np_database['Accession'] == k]
    np_database_filtered = np_database_filtered.reset_index(drop=True)
    sequence_to_report = np_database_filtered['Sequence'][0]
    
    samples_filter_df = sample_retrieve_df_all[sample_retrieve_df_all['Accession'] == k]
    samples_impacted = samples_filter_df['Sample'].values.tolist()
    
    all_results_storage_peptide.append((sequence_to_report))
    all_results_storage_accession.append(k)
    all_results_storage_instances.append(instances)
    all_results_storage_samples.append(samples_impacted)

##

sample_retrieve_df_filtered = pd.DataFrame()
sample_retrieve_df_filtered['Accession'] = unique_to_peaks_from_FDR_filtered
sample_retrieve_df_filtered['Sample'] = unique_to_peaks_from_FDR_filtered_sample

unique_to_peaks_filtered_no_dups = []
for l in unique_to_peaks_from_FDR_filtered:
    if l not in unique_to_peaks_all_no_dups:
        unique_to_peaks_filtered_no_dups.append(l)
    else:
        pass

filtered_results_storage_accession = []
filtered_results_storage_peptide = []
filtered_results_storage_instances = []
filtered_results_storage_samples = []

for m in unique_to_peaks_filtered_no_dups:
    instances = unique_to_peaks_from_FDR_filtered.count(m)
    np_database_filtered = np_database[np_database['Accession'] == m]
    np_database_filtered = np_database_filtered.reset_index(drop=True)
    sequence_to_report = np_database_filtered['Sequence'][0]
    
    samples_filter_df = sample_retrieve_df_filtered[sample_retrieve_df_filtered['Accession'] == m]
    samples_impacted = samples_filter_df['Sample'].values.tolist()
    
    filtered_results_storage_peptide.append((sequence_to_report))
    filtered_results_storage_accession.append(m)
    filtered_results_storage_instances.append(instances)
    filtered_results_storage_samples.append(samples_impacted)

unfiltered_unique_to_PEAKS = pd.DataFrame()
unfiltered_unique_to_PEAKS['Sequence'] = all_results_storage_peptide
unfiltered_unique_to_PEAKS['Accession'] = all_results_storage_accession
unfiltered_unique_to_PEAKS['Instances'] = all_results_storage_instances
unfiltered_unique_to_PEAKS['Samples'] = all_results_storage_samples


file_path = output_directory + '\\unique_to_PEAKS_from_unfiltered_db_search_results.csv'
with open(file_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        unfiltered_unique_to_PEAKS.to_csv(filec,index=False)
##

filtered_unique_to_PEAKS = pd.DataFrame()
filtered_unique_to_PEAKS['Sequence'] = filtered_results_storage_peptide
filtered_unique_to_PEAKS['Accession'] = filtered_results_storage_accession
filtered_unique_to_PEAKS['Instances'] = filtered_results_storage_instances
filtered_unique_to_PEAKS['Samples'] = filtered_results_storage_samples


file_path = output_directory + '\\unique_to_PEAKS_from_filtered_db_search_results_1_FDR.csv'
with open(file_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        filtered_unique_to_PEAKS.to_csv(filec,index=False)
