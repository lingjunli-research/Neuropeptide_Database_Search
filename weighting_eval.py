# -*- coding: utf-8 -*-
"""
Created on Mon May  8 17:09:44 2023

@author: lawashburn
"""

import pandas as pd

weighted_metrics_results_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale_weighting_hyperscore\weighting_metric_extraction\metrics_out_v2.csv"
round1_dsd_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale_weighting_hyperscore\weighting_dsd_eval\round1_DSD_test.csv"

weighted_metrics_results = pd.read_csv(weighted_metrics_results_path)
round1_dsd = pd.read_csv(round1_dsd_path)

for ind in round1_dsd.index:
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
    
    weighted_metrics_results['peptide_score'] = (
        (weighted_metrics_results['Average Fragment Error']*avg_frag_err)+
        (weighted_metrics_results['Precursor Error']*prec_err)+
        (weighted_metrics_results['# Consecutive b-ions']*consec_b)+
        (weighted_metrics_results['# Consecutive y-ions']*consec_y)+
        (weighted_metrics_results['C-termini amidation']*c_term_amid)+
        (weighted_metrics_results['% Sequence coverage']*percent_seq_cov)+
        (weighted_metrics_results['Average annotations/fragment ']*no_annot_peak)+
        (weighted_metrics_results['% Fragment ions are b']*percent_b_ions)+
        (weighted_metrics_results['% Fragment ions are y']*percent_y_ions)+
        (weighted_metrics_results['Average number of fragment ions per AA']*avg_matched_ions_per_AA)+
        (weighted_metrics_results['Number of non-neutral-loss fragment ions per AA']*avg_matched_ions_per_AA_non_neut)
        )
    
