# -*- coding: utf-8 -*-
"""
Created on Wed May 10 15:53:33 2023

@author: lawashburn
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-deep')

weighted_metrics_results_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale_weighting_hyperscore\weighting_metric_extraction\metrics_out_v3_TD.csv"
target_list_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\target_list.csv"
round1_dsd_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale_weighting_hyperscore\weighting_dsd_eval\round1_DSD_test.csv"

weighted_metrics_results = pd.read_csv(weighted_metrics_results_path)
round1_dsd = pd.read_csv(round1_dsd_path)

target_list = pd.read_csv(target_list_path)
target_list['Sequence'] = target_list['Sequence'].str.replace(r"\([^)]*\)","")
target_list_str = target_list['Sequence'].values.tolist()

weighted_metrics_results['Unmodified sequence'] = weighted_metrics_results['Peptide']
weighted_metrics_results['Unmodified sequence'] = weighted_metrics_results['Unmodified sequence'].str.replace(r"\([^)]*\)","")
weighted_metrics_results['Status'] = weighted_metrics_results['Unmodified sequence'].apply(lambda x: any([k in x for k in target_list_str]))

target_results = weighted_metrics_results[weighted_metrics_results['Status'] == True]
decoy_results = weighted_metrics_results[weighted_metrics_results['Status'] == False]

#%%

plt.hist([x, y], bins, label=['x', 'y'])
plt.legend(loc='upper right')
plt.show()