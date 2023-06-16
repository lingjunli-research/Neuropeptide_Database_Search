# -*- coding: utf-8 -*-
"""
Created on Tue May 16 14:47:00 2023

@author: lawashburn
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib_venn import venn3, venn3_circles
from collections import Counter
plt.style.use('seaborn-deep')

weighted_metrics_results_path = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\weighting_metric_extraction_motif_full\metrics_out_motif_score_1_NFDEID.csv"
target_list_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\target_list.csv"
round1_dsd_path = r"C:\Users\lawashburn\Documents\DBpep_v2\Weighting_Assessment\Definitive Screening Design_v3_20230527_w_motif.csv"
output_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\DSD_eval\motif_1_NFDEID"
validation_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\NV_output_DB_search_results"
sample_types = ['Brain','CoG','PO','SG','TG']
#sample_types = ['Brain']
round1_fdr = 0.01
power = 2
weighting_integer = 1
num_rounds = 10
reps = [1,2,3]

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
    num_decoy_IDs = []
    for ind in round1_dsd.index:

        weighted_metrics_results_filtered = weighted_metrics_results[weighted_metrics_results['Sample'] == sample_ID]

        avg_frag_err = round1_dsd['Avg Fragment Error'][ind]
        prec_err = round1_dsd['Precursor Error'][ind]
        consec_b = round1_dsd['# Consecutive b-ions'][ind]
        consec_y = round1_dsd['# Consecutive y-ions'][ind]
        # c_term_amid = round1_dsd['C-terminal amidation'][ind]
        percent_seq_cov = round1_dsd['Motif Score'][ind]
        no_annot_peak = round1_dsd['Avg annotations/fragment'][ind]
        # percent_b_ions = round1_dsd['% ions are b'][ind]
        # percent_y_ions = round1_dsd['% ions are y'][ind]
        avg_matched_ions_per_AA = round1_dsd['Avg frag ions/AA'][ind]
        avg_matched_ions_per_AA_non_neut = round1_dsd['Avg non-neut frag/AA'][ind]
        hyperscore = round1_dsd['Hyperscore'][ind]
        
        avg_frag_err_norm_score = weighting_integer/(max(weighted_metrics_results_filtered['Average Fragment Error']))
        prec_err_norm_score = weighting_integer/(max(weighted_metrics_results_filtered['Precursor Error']))
        consec_b_ions_norm_score = weighting_integer/(max(weighted_metrics_results_filtered['# Consecutive b-ions']))
        consec_y_ions_norm_score = weighting_integer/(max(weighted_metrics_results_filtered['# Consecutive y-ions']))
        precent_sec_cov_norm_score = weighting_integer/(max(weighted_metrics_results_filtered['Motif_Score']))
        avg_annotation_norm_score = weighting_integer/(max(weighted_metrics_results_filtered['Average annotations/fragment ']))
        avg_frag_ions_norm_score = weighting_integer/(max(weighted_metrics_results_filtered['Average number of fragment ions per AA']))
        avg_non_neut_frag_ions_norm_score = weighting_integer/(max(weighted_metrics_results_filtered['Number of non-neutral-loss fragment ions per AA']))
        hyperscore_norm_score = weighting_integer/(max(weighted_metrics_results_filtered['Hyperscore']))

        weighted_metrics_results_filtered['peptide_score'] = (
            # (((weighted_metrics_results_filtered['Average Fragment Error']**power)*avg_frag_err_norm_score)*avg_frag_err)+
            # (((weighted_metrics_results_filtered['Precursor Error']**power)*prec_err_norm_score)*prec_err)+
            # # (((weighted_metrics_results_filtered['# Consecutive b-ions']**power)*consec_b_ions_norm_score)*consec_b)+
            # # (((weighted_metrics_results_filtered['# Consecutive y-ions']**power)*consec_y_ions_norm_score)*consec_y)+
            (((weighted_metrics_results_filtered['Motif_Score']**power)*precent_sec_cov_norm_score)*percent_seq_cov)+
            # (((weighted_metrics_results_filtered['Average annotations/fragment ']**power)*avg_annotation_norm_score)*no_annot_peak)+
            # (((weighted_metrics_results_filtered['Average number of fragment ions per AA']**power)*avg_frag_ions_norm_score)*avg_matched_ions_per_AA)+
            # (((weighted_metrics_results_filtered['Number of non-neutral-loss fragment ions per AA']**power)*avg_non_neut_frag_ions_norm_score)*avg_matched_ions_per_AA_non_neut)+
            (((weighted_metrics_results_filtered['Hyperscore']**power)*hyperscore_norm_score)*hyperscore)
            )
        
        # weighted_metrics_results_filtered = weighted_metrics_results_filtered.rename(columns={"peptide_score": "inital peptide score"})
        # weighted_metrics_results_filtered['Pep no mod'] = weighted_metrics_results_filtered['Peptide'].str.replace(r"\(.*\)","")
        # weighted_metrics_results_filtered['Pep len'] = weighted_metrics_results_filtered['Pep no mod'].astype(str).map(len)
        # weighted_metrics_results_filtered['peptide_score'] = weighted_metrics_results_filtered['inital peptide score']/weighted_metrics_results_filtered['Pep len']
        
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
            
            if len(target_results_filtered) >0:
                fdr_report = len(decoy_results_filtered)/len(target_results_filtered)
            else:
                fdr_report = 1
            
            score_cutoff.append(score)
            fdr.append(fdr_report)
            if len(target_results_filtered)>0:
                num_target_IDs.append(len(target_results_filtered))
            else:
                num_target_IDs.append(0)
           
            if len(decoy_results_filtered)>0:
                num_decoy_IDs.append(len(decoy_results_filtered))
            else:
                num_decoy_IDs.append(0)
        
        fdr_table = pd.DataFrame()
        fdr_table['Score Threshold'] = score_cutoff
        fdr_table['FDR'] = fdr
        fdr_table['# Target IDs'] = num_target_IDs
        #fdr_table['# Decoy IDs'] = num_decoy_IDs
        
        file_path = output_directory + '\\FDR_eval_table' + sample_ID + '_DSD_run_'+ str(ind+1) + '.csv'
        with open(file_path,'w',newline='') as filec:
                writerc = csv.writer(filec)
                fdr_table.to_csv(filec,index=False)
        
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
        
        fdr_table_filtered = fdr_table[fdr_table['# Target IDs'] == best_num_targets]
        score_threshold_to_apply = fdr_table_filtered['Score Threshold'].min()
        
        target_results_final = target_results[target_results['peptide_score'] >= score_threshold_to_apply]
        decoy_results_final = decoy_results[decoy_results['peptide_score'] >= score_threshold_to_apply]
        
        file_path = output_directory + '\\final_results' + sample_ID +  '_DSD_run_'+ str(ind+1) + '_decoy.csv'
        with open(file_path,'w',newline='') as filec:
                writerc = csv.writer(filec)
                decoy_results_final.to_csv(filec,index=False)
        
        file_path = output_directory + '\\final_results' + sample_ID +  '_DSD_run_'+ str(ind+1) + '_target.csv'
        with open(file_path,'w',newline='') as filec:
                writerc = csv.writer(filec)
                target_results_final.to_csv(filec,index=False)
        
    file_path = output_directory + '\\dsd_eval_round1_' + sample_ID + '.csv'
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            dsd_merge_table.to_csv(filec,index=False)


sample_type_storage = []
run_number_storage = []
number_unique_IDs_storage = []
number_quantifiable_unique_IDs_storage = []
 #%%    
for l in sample_types:
    run_number_storage_inner = []
    number_unique_IDs_storage_inner = []
    number_quantifiable_unique_IDs_storage_inner = []
    #%%
    for j in run_numbers:
        
        sample_1_results_path = output_directory + '\\final_results2021_0817_' + l + '_1_DSD_run_' + str(j) + '_target.csv'
        sample_2_results_path = output_directory + '\\final_results2021_0817_' + l + '_2_DSD_run_' + str(j) + '_target.csv'
        sample_3_results_path = output_directory + '\\final_results2021_0817_' + l + '_3_DSD_run_' + str(j) + '_target.csv'
        
        sample1_results = pd.read_csv(sample_1_results_path)
        sample2_results = pd.read_csv(sample_2_results_path)
        sample3_results = pd.read_csv(sample_3_results_path)
        
        if len(sample1_results)>0:
            if len(sample2_results)>0:
                if len(sample3_results)>0:
                    sample1_results_no_dups = sample1_results.drop_duplicates(subset='Peptide')
                    sample1_unique_peps = set(sample1_results_no_dups['Peptide'].values.tolist())
                    
                    
                    sample2_results_no_dups = sample2_results.drop_duplicates(subset='Peptide')
                    sample2_unique_peps = set(sample2_results_no_dups['Peptide'].values.tolist())
                    
                    
                    sample3_results_no_dups = sample3_results.drop_duplicates(subset='Peptide')
                    sample3_unique_peps = set(sample3_results_no_dups['Peptide'].values.tolist())
                    
                    AB_overlap = sample1_unique_peps & sample2_unique_peps  #compute intersection of set A & set B
                    AC_overlap = sample1_unique_peps & sample3_unique_peps
                    BC_overlap = sample2_unique_peps & sample3_unique_peps
                    ABC_overlap = sample1_unique_peps & sample2_unique_peps & sample3_unique_peps
                    A_rest = sample1_unique_peps - AB_overlap - AC_overlap #see left graphic
                    B_rest = sample2_unique_peps - AB_overlap - BC_overlap
                    C_rest = sample3_unique_peps - AC_overlap - BC_overlap
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
                    labels = ('Rep 1', 'Rep 2', 'Rep 3')  
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
            
                    plt_title = l + ' unique IDs' + ' DSD run ' + str(j)
                    plt.title(plt_title, fontsize=24, fontweight='bold',font='Arial')
                    plt.tight_layout()
                    plt.show()
            #%%
                    plt_out_path = output_directory + '\\' + plt_title + '.png'
                    plt.savefig(plt_out_path)
                    plt.clf()
  #%%          
                    all_unique_peps = []
                    quantifiable_peps = []
                    
                    for a in sample1_unique_peps:
                        if a not in all_unique_peps:
                            all_unique_peps.append(a)
                            
                    
                    for b in sample2_unique_peps:
                        if b not in all_unique_peps:
                            all_unique_peps.append(b)
                    
                    for c in sample3_unique_peps:
                        if c not in all_unique_peps:
                            all_unique_peps.append(c)
                            
                    for d in sample1_unique_peps:
                        if d in sample2_unique_peps:
                            if d in sample3_unique_peps:
                                if d not in quantifiable_peps:
                                    quantifiable_peps.append(d)
                                else:
                                    pass
                            else:
                                pass
                        else:
                            pass
                    
                    sample_type_storage.append(l)
                    run_number_storage.append(str(j))
                    number_unique_IDs_storage.append(len(all_unique_peps))
                    number_quantifiable_unique_IDs_storage.append(len(quantifiable_peps))
                    
                    run_number_storage_inner.append(str(j))
                    number_unique_IDs_storage_inner.append(len(all_unique_peps))
                    number_quantifiable_unique_IDs_storage_inner.append(len(quantifiable_peps))
            
                    for x in reps:
                        unique_id_count = []
                        unique_id_round_log = []
                        unique_psm_count = []
                        unique_psm_round_log = []
                        
                        
                        results_directory_rounds = validation_directory + '\\2021_0817_' + l + '_' + str(x)
                        for i in range(0,num_rounds):
                            unique_ID_storage = []
                            unique_psm_storage = []
                            inital_round_storage = []
                            psm_round_storage = []
                            
                            psm_report_path = results_directory_rounds + '\\final_psm_report' + str(i) + '.csv'
                            
                            psm_report = pd.read_csv(psm_report_path)
                            psm_report_sequence = psm_report['Sequence'].values.tolist()
                            
                            seq_scan_concat = []
                            psm_report_scan = psm_report['Scan'].values.tolist()
                            
                            for m in range(0,(len(psm_report_scan))):
                                seq_ind = psm_report_sequence[m]
                                scan_ind = psm_report_scan[m]
                                concat_values = seq_ind + '_' + str(scan_ind)
                                seq_scan_concat.append(concat_values)
                            
                            
                            sample_results_path = output_directory + '\\final_results2021_0817_' + l + '_' + str(x) + '_DSD_run_' + str(j) + '_target.csv'
                            sample_results = pd.read_csv(sample_results_path)
                            sample_results_no_dups = sample_results.drop_duplicates(subset='Peptide')
                            sample_unique_peps = set(sample_results_no_dups['Peptide'].values.tolist())
                            
                            round_seq_report = sample_results['Peptide'].values.tolist()
                            round_scan_report = sample_results['Scan'].values.tolist()
                            
                            seq_scan_round_concat = []
                            
                            for p in range(0,(len(round_seq_report))):
                                seq_ind_round = round_seq_report[p]
                                scan_ind_round = round_scan_report[p]
                                concat_values_round = seq_ind_round + '_' + str(scan_ind_round)
                                seq_scan_round_concat.append(concat_values_round)
                            
                            for y in sample_unique_peps:
                                if y not in unique_ID_storage:
                                    if y in psm_report_sequence:
                                        unique_ID_storage.append(y)
                                        inital_round_storage.append(i)
                            unique_id_count.append(len(unique_ID_storage))
                            unique_id_round_log.append(i)
                            
                            for w in seq_scan_concat:
                                if w not in unique_psm_storage:
                                    if w in seq_scan_round_concat:
                                        unique_psm_storage.append(w)
                                        psm_round_storage.append(i)
                            unique_psm_count.append(len(unique_psm_storage))
                            unique_psm_round_log.append(i)
            
                        colors = ['#b0cac7']
                        rounds_evaluation = pd.DataFrame()
                        rounds_evaluation['# Unique IDs'] = unique_id_count
                        rounds_evaluation['Round #'] = unique_id_round_log
                        df_groups = rounds_evaluation.groupby(['Round #'])['# Unique IDs'].sum()
                        plt_name2 = l + ' Rep ' + str(x) + ' DSD ' + str(j) + ' Unique Peptides'
                        v = df_groups.plot(kind='bar', title=plt_name2,
                           ylabel='Count', xlabel='Round', figsize=(10, 6),color=colors)
                        
                        plt.title(plt_name2, fontsize=24, fontweight='bold',font='Arial')
                        v.set_xlabel('Round #', rotation=0, fontsize=18,font='Arial',fontweight='bold')
                        v.set_ylabel('Count', rotation=90, fontsize=18,font='Arial',fontweight='bold')
                        plt.xticks(rotation=360,fontsize=16,font='Arial')
                        plt.yticks(rotation=360,fontsize=16,font='Arial')
            
                        plt.tight_layout()
                        plt_out_path = output_directory + '\\' + plt_name2 + '.png'
                        plt.savefig(plt_out_path)
                        plt.clf()
            
                        colors = ['#001244']
                        rounds_evaluation_psm = pd.DataFrame()
                        rounds_evaluation_psm['# PSMs'] = unique_psm_count
                        rounds_evaluation_psm['Round #'] = unique_psm_round_log
                        df_groups_rounds = rounds_evaluation_psm.groupby(['Round #'])['# PSMs'].sum()
                        plt_name3 = l + ' Rep ' + str(x) + ' DSD ' + str(j) + ' PSMs'
                        v = df_groups_rounds.plot(kind='bar', title=plt_name3,
                            ylabel='Count', xlabel='Round', figsize=(10, 6),color=colors)
                        plt.title(plt_name3, fontsize=24, fontweight='bold',font='Arial')
                        v.set_xlabel('Round #', rotation=0, fontsize=18,font='Arial',fontweight='bold')
                        v.set_ylabel('Count', rotation=90, fontsize=18,font='Arial',fontweight='bold')
                        plt.xticks(rotation=360,fontsize=16,font='Arial')
                        plt.yticks(rotation=360,fontsize=16,font='Arial')
            
                        plt.tight_layout()
                        plt_out_path = output_directory + '\\' + plt_name3 + '.png'
                        plt.savefig(plt_out_path)
                        plt.clf()
            
                        rounds_eval_grouped = pd.DataFrame()
                        rounds_eval_grouped['Round #'] = rounds_evaluation_psm['Round #']
                        rounds_eval_grouped['# Unique IDs'] = rounds_evaluation['# Unique IDs']
                        rounds_eval_grouped['# PSMs'] = rounds_evaluation_psm['# PSMs']
            
                        plt.style.use('default')
                        colors = ['#b0cac7', '#001244']
                        v = rounds_eval_grouped[['Round #', '# Unique IDs', '# PSMs']].plot(x='Round #', kind='bar',stacked=True, color=colors,
                                                                                        ylabel='Count', xlabel='Round')
                        #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                        plt.legend(loc='lower right',fontsize=10,prop={"family":"Arial"})
                        plt.title((l + ' Rep ' + str(x) + ' DSD ' + str(j) + ' PSMs & Unique IDs'),fontsize=24, fontweight='bold',font='Arial')
                        plt_name4 = l + ' Rep ' + str(x) + ' DSD ' + str(j) + ' PSMs_Peptides'
                        v.set_xlabel('Round #', rotation=0, fontsize=18,font='Arial',fontweight='bold')
                        v.set_ylabel('Count', rotation=90, fontsize=18,font='Arial',fontweight='bold')
                        plt.xticks(rotation=360,fontsize=16,font='Arial')
                        plt.yticks(rotation=360,fontsize=16,font='Arial')
                        plt_out_path = output_directory + '\\' + plt_name4 + '.png'
                        plt.savefig(plt_out_path)
                        plt.clf()

    inner_bar_plot_generate = pd.DataFrame()
    inner_bar_plot_generate['Run #'] = run_number_storage_inner
    
    inner_bar_plot_generate['Quantifiable'] = number_quantifiable_unique_IDs_storage_inner
    inner_bar_plot_generate['Unique IDs'] = number_unique_IDs_storage_inner
    if len(inner_bar_plot_generate)>0:
        ax = inner_bar_plot_generate.plot.bar(x='Run #', stacked=True, color=['#b0cac7','#001244'], figsize=(8,6))
        ax.set_title(l, fontsize=24,fontweight='bold',font='Arial')
        ax.set_xticklabels(run_number_storage_inner, rotation=0, fontsize=16,font='Arial')
        plt.xticks(rotation=360,fontsize=16,font='Arial')
        plt.yticks(rotation=360,fontsize=16,font='Arial')
        ax.set_xlabel('Run #', rotation=0, fontsize=20,font='Arial',fontweight='bold')
        ax.set_ylabel('Count', rotation=90, fontsize=20,font='Arial',fontweight='bold')
        plt.legend(fontsize=20,prop={"family":"Arial"})
        plt.tight_layout()
    
        plt_out_path = output_directory + '\\' + l + '_quantifiable_IDs.png'
        plt.savefig(plt_out_path)
        plt.clf()
    else:
        pass

final_results_df = pd.DataFrame()
final_results_df['Sample Type'] = sample_type_storage
final_results_df['Run #'] = run_number_storage
final_results_df['# Unique IDs'] = number_unique_IDs_storage
final_results_df['# Quantifiable IDs'] = number_quantifiable_unique_IDs_storage   
file_path = output_directory + '\\summary_all_DSDs.csv'
with open(file_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        final_results_df.to_csv(filec,index=False)