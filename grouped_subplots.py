# -*- coding: utf-8 -*-
"""
Created on Thu May 25 17:24:11 2023

@author: lawashburn
"""

# Data and imports

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
import matplotlib

df_path = r"C:\Users\lawashburn\Documents\DBpep_v2\validation\dbsearch_results_20230520\pos_neg_PSM_unique_ID_motif_weightinginteger_10\summary_all_DSDs.csv"
output_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\validation\dbsearch_results_20230520\pos_neg_PSM_unique_ID_motif_weightinginteger_10"
rounds_eval_grouped = pd.read_csv(df_path)
sample_names = ['Brain', 'CoG', 'PO', 'SG', 'TG']

runs = rounds_eval_grouped['Run #'].values.tolist()

runs_no_dups = []
for a in runs:
    if a not in runs_no_dups:
        runs_no_dups.append(a)
print(len(runs_no_dups))

# fig = plt.figure(1,figsize=(16,18))
# y_max = (max(rounds_eval_grouped['# PSMs']) + max(rounds_eval_grouped['# Unique IDs']) + max(rounds_eval_grouped['# Quantifiable IDs']))*1.1

# for i,k in enumerate(runs_no_dups,1):
#     fig.add_subplot(6,5,i,)

#     rounds_eval_grouped_filter = rounds_eval_grouped[rounds_eval_grouped['Run #'] == k]
#     if len(rounds_eval_grouped_filter) != len(sample_names):
#         for name in sample_names:
#             rounds_eval_grouped_filter2 = rounds_eval_grouped_filter[rounds_eval_grouped_filter['Sample Type'] == name]
#             if len(rounds_eval_grouped_filter2) == 0:
#                 df2 = {'Sample Type': name, 'Run #': k, '# Unique IDs': 0, '# Quantifiable IDs':0,'# PSMs':0}
#                 rounds_eval_grouped_filter = rounds_eval_grouped_filter.append(df2, ignore_index = True)
#             else:
#                 pass
#     else:
#         pass
                
    
#     plt.style.use('default')
#     colors = ['#b0cac7', '#001244','#318fb5']
#     v = rounds_eval_grouped_filter[['Sample Type', '# Quantifiable IDs', '# Unique IDs','# PSMs']].plot(x='Sample Type', kind='bar',stacked=True, color=colors,
#                                                                     ylabel='Count',ax=plt.gca(),legend=False)
#     #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#     v.set_ylim([0, y_max])
#     plt.title((str(k)),fontsize=24, fontweight='bold',font='Arial')
#     # plt_name4 = l + ' Rep ' + str(x) + ' DSD ' + str(j) + ' PSMs_Peptides'
#     list_of_x_axes = [17,18,20,22,24]
#     if k in list_of_x_axes:
#         v.set_xlabel('Run #', rotation=0, fontsize=18,font='Arial',fontweight='bold')
#         plt.xticks(rotation=45,fontsize=16,font='Arial')
#     else:
#         v.set_xlabel('', rotation=0, fontsize=18,font='Arial',fontweight='bold')
#         plt.xticks(rotation=360,fontsize=0,font='Arial')
    
#     list_of_y_axes = [1,12,18]
#     if k in list_of_y_axes:
#         v.set_ylabel('Count', rotation=90, fontsize=18,font='Arial',fontweight='bold')
#         plt.yticks(rotation=360,fontsize=16,font='Arial')
#     else:
#         v.set_ylabel('', rotation=0, fontsize=18,font='Arial',fontweight='bold')
#         plt.yticks(rotation=360,fontsize=0,font='Arial')

# plt.legend(loc = 'lower right',bbox_to_anchor = (2,.7,.6,0))
# plt.tight_layout(pad=20.0)
# plt.subplots_adjust(left=0.1,
#                     bottom=0.3,
#                     right=0.9,
#                     top=1.9,
#                     wspace=0.6,
#                     hspace=0.4)
# plt.show()
# plt.clf()
#%%
fig2 = plt.figure(1,figsize=(16,20))
y_max2 = (max(rounds_eval_grouped['# Unique IDs']))*1.1

for i,k in enumerate(runs_no_dups,1):
    fig2.add_subplot(6,5,i,)

    rounds_eval_grouped_filter = rounds_eval_grouped[rounds_eval_grouped['Run #'] == k]
    if len(rounds_eval_grouped_filter) != len(sample_names):
        for name in sample_names:
            rounds_eval_grouped_filter2 = rounds_eval_grouped_filter[rounds_eval_grouped_filter['Sample Type'] == name]
            if len(rounds_eval_grouped_filter2) == 0:
                df2 = {'Sample Type': name, 'Run #': k, '# Unique IDs': 0, '# Quantifiable IDs':0,'# PSMs':0}
                rounds_eval_grouped_filter = rounds_eval_grouped_filter.append(df2, ignore_index = True)
            else:
                pass
    else:
        pass
                
    
    plt.style.use('default')
    colors = ['#318fb5']
    v = rounds_eval_grouped_filter[['Sample Type', '# Unique IDs']].plot(x='Sample Type', kind='bar',stacked=True, color=colors,
                                                                    ylabel='Count',ax=plt.gca(),legend=False)
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    v.set_ylim([0, y_max2])
    plt.title((str(k)),fontsize=24, fontweight='bold',font='Arial')
    # plt_name4 = l + ' Rep ' + str(x) + ' DSD ' + str(j) + ' PSMs_Peptides'
    list_of_x_axes = [26,17,23,24,25]
    if k in list_of_x_axes:
        v.set_xlabel('Run #', rotation=0, fontsize=18,font='Arial',fontweight='bold')
        plt.xticks(rotation=45,fontsize=16,font='Arial')
    else:
        v.set_xlabel('', rotation=0, fontsize=18,font='Arial',fontweight='bold')
        plt.xticks(rotation=360,fontsize=0,font='Arial')
    
    list_of_y_axes = [5,16,26]
    if k in list_of_y_axes:
        v.set_ylabel('Count', rotation=90, fontsize=18,font='Arial',fontweight='bold')
        plt.yticks(rotation=360,fontsize=16,font='Arial')
    else:
        v.set_ylabel('', rotation=0, fontsize=18,font='Arial',fontweight='bold')
        plt.yticks(rotation=360,fontsize=0,font='Arial')

plt.legend(loc = 'lower right',bbox_to_anchor = (2,.7,.6,0))
plt.tight_layout(pad=20.0)
plt.subplots_adjust(left=0.1,
                    bottom=0.3,
                    right=0.9,
                    top=1.9,
                    wspace=0.6,
                    hspace=0.4)
plt.show()
