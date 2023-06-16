# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 13:00:08 2023

@author: lawashburn
"""

import pandas as pd
import os
import scipy
from scipy.spatial import distance
import numpy as np
import pysptools
from pysptools import distance as dist_2
from distance_metrics_mcda import distance_metrics
import csv


final_rep = r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\KD_Training_Spectra\Perfect_spectra_compile\perfect_spectra_list_.csv"
dir_to_search = r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230214\KD_search_results_v5\SG_Unlabeled_DDA_TR3\xcorr_data"
output_directoy = r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230214\KD_perterb_results_v2"
sample = 'KD_' + 'SG_Unlabeled_DDA_TR3'

final_rep = pd.read_csv(final_rep)
final_rep["Identifier"] = final_rep['Sequence'].astype(str) +"_"+ final_rep["Scan"].astype(str)

final_rep_seq = final_rep['Identifier'].values.tolist()
final_rep_seq = set(final_rep_seq)

def get_file_names_with_strings(str_list):
    full_list = os.listdir(dir_to_search)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

seqs = []
perterb_level_rep = []
# euc = []
# bray = []
# sam = []
# spec_div = []
# cheb = []
# can = []
# manhattan = []
# cosine = []
# jen_shan = []
# minko = []
# sq_euclid = []
# diced = []
# hamming2 = []
# jaccard2 = []
# kulsin = []
# roger = []
# rus = []
# sokal2 = []
# sokal3 = []
# yule2 = []
# haus = []
# sq_cho = []
# ent = []
dotttt = []



for seq in final_rep_seq:
    file_query = seq
    fragment_list = (get_file_names_with_strings([file_query]))
    print(file_query)
    print(fragment_list)

    if len(fragment_list) == 2:
        
        exp_df = pd.DataFrame()
        theo_df = pd.DataFrame()
        
        exp_query = 'exp'
        theo_query = 'theo'
        
        bin_storage = []
        bin_storage_no_dups = []
        for a in fragment_list:
            if exp_query in a:
                path = dir_to_search + '\\' + a
                
                read_file = pd.read_csv(path)
                read_file = read_file.round({'experimental m/z': 0})
                read_file = read_file.sort_values(by=['experimental intensity'],ascending=False)
                read_file = read_file.drop_duplicates(subset='experimental m/z')
                exp_df = pd.concat([exp_df,read_file])
                exp_df_bins = exp_df['experimental m/z'].values.tolist()
                bin_storage.extend(exp_df_bins)
            if theo_query in a:
                path = dir_to_search + '\\' + a
                
                read_file = pd.read_csv(path)
                read_file = read_file.round({'theoretical m/z': 0})
                read_file = read_file.sort_values(by=['theoretical intensity'],ascending=False)
                read_file = read_file.drop_duplicates(subset='theoretical m/z')
                theo_df = pd.concat([theo_df,read_file])
                theo_df_bins = theo_df['theoretical m/z'].values.tolist()
                bin_storage.extend(theo_df_bins)
            else:
                pass
       
        
        for bins in bin_storage:
            if bins not in bin_storage_no_dups:
                bin_storage_no_dups.append(bins)
            else:
                pass 
        perterb_list = [0,10,20,30,40,50,60,70,80,90,100]
        for perterb in perterb_list:
                perterbation_level = perterb/100
                perterb_level_rep.append(perterbation_level)
                exp_array_mod = exp_df.sample(frac=perterbation_level)

                merge_array = pd.DataFrame()
                merge_array['Bin'] = bin_storage_no_dups
                merge_array = merge_array.merge(exp_array_mod,left_on='Bin', right_on='experimental m/z',how='left')
                merge_array = merge_array.merge(theo_df,left_on='Bin', right_on='theoretical m/z',how='left')
                
                merge_array = merge_array.fillna(0)
                
                array1 = merge_array['experimental intensity'].values.tolist()
                array2 = merge_array['theoretical intensity'].values.tolist()
                
                # euc_distance = distance.euclidean(array1, array2) #euclidean distance calculation
                # bray_curt_distance = distance.braycurtis(array1, array2) #bray-curtis distance calculation
                dot_prod = np.dot(array1, array2) #dot product calculation
                # #normXcorr = dist_2.NormXCorr(array1, array2) #normalized cross correlation
                # sam_map = pysptools.distance.SAM(array1, array2) #Computes the spectral angle mapper between two vectors (in radians)
                # spectral_divergence = pysptools.distance.SID(array1, array2) #Computes the spectral information divergence between two vectors
                # chebychev_distance = distance.chebyshev(array1, array2) #Chebychev distance
                # canberra = distance.canberra(array1, array2) #canberra distance
                # nyc = distance.cityblock(array1, array2) #manhattan/city block distance
                # cos = distance.cosine(array1, array2) #cosine distance
                # jen_sha = distance.jensenshannon(array1, array2) #jensen-shannon distance
                # #mahala = distance.mahalanobis(array1, array2) #Mahalanobis distance
                # minkowski = distance.minkowski(array1, array2,p=1) #Minkowski distance
                # sq_euc = distance.sqeuclidean(array1, array2) #squared euclidean distance
                # dice = distance.dice(array1, array2) #Dice dissimilarity distance
                # hamming = distance.hamming(array1, array2) #Hamming distance
                # jaccard = distance.jaccard(array1, array2) #Jaccard-Needham dissimilary
                # kulsinski = distance.kulsinski(array1, array2) #Kulsinski dissimilarity
                # #kulczynski1 = distance.kulczynski1(array1, array2) #Kulczynski 1 dissimilarity
                # rogerstanimoto = distance.rogerstanimoto(array1, array2) #Rogers-Tanimoto dissimilarity
                # russelrao = distance.russellrao(array1, array2) #Russell-Rao dissimilarity
                # sokal = distance.sokalmichener(array1, array2) #Sokal-Michener dissimilarity
                # sokal_sneath = distance.sokalsneath(array1, array2) #Sokal-Sneath dissimilarity
                # yule = distance.yule(array1, array2) #Yule dissimilarity
                # hausdorff = distance_metrics.hausdorff(array1, array2) #Hausdorff distance
                # #lorentzian = distance_metrics.lorentzian(array1, array2) #Lorentzian distance
                # #bhattacharyya = distance_metrics.bhattacharyya(array1, array2) #bhattacharyya distance
                # #hellinger = distance_metrics.hellinger(array1, array2) #hellinger distance
                # #matusita = distance_metrics.matusita(array1, array2) #matusita distance
                # sq_chord = distance_metrics.squared_chord(array1, array2) #Squared-chrod distance
                # #pearson_chi_square = distance_metrics.pearson_chi_square(array1, array2) #pearson_chi_square distance
                # #squared_chi_square = distance_metrics.pearson_chi_square(array1, array2) #squared_chi_square distance
                # #entropy = scipy.stats.entropy(array1, array2) #Shannon entropy calculation
                
                seqs.append(seq)
                # euc.append(euc_distance)
                # bray.append(bray_curt_distance)
                dotttt.append(dot_prod)
                # sam.append(sam_map)
                # spec_div.append(spectral_divergence)
                # cheb.append(chebychev_distance)
                # can.append(canberra)
                # manhattan.append(nyc)
                # cosine.append(cos)
                # jen_shan.append(jen_sha)
                # minko.append(minkowski)
                # sq_euclid.append(sq_euc)
                # diced.append(dice)
                # hamming2.append(hamming)
                # jaccard2.append(jaccard)
                # kulsin.append(kulsinski)
                # roger.append(rogerstanimoto)
                # rus.append(russelrao)
                # sokal2.append(sokal)
                # sokal3.append(sokal_sneath)
                # yule2.append(yule)
                # haus.append(hausdorff)
                # sq_cho.append(sq_chord)
                #ent.append(entropy)

                
results_report = pd.DataFrame()
results_report['Sequence'] =  seqs      
results_report['Pertermation_level'] =  perterb_level_rep
# results_report['Euclidean distance'] =  euc  
# results_report['Bray-Curtis'] =  bray  
# results_report['Spectral angle'] =  sam  
# results_report['Spectral divergence'] =  spec_div  
# results_report['Chebychev'] =  cheb  
# results_report['Canberra'] =  can  
# results_report['Manhattan'] =  manhattan  
# results_report['Cosine distance'] =  cosine  
# results_report['Jensen-shannon'] =  jen_shan  
# results_report['Minkowski'] =  minko  
# results_report['Squared euclidean'] =  sq_euclid 
# results_report['Dice'] =  diced
# results_report['Hamming'] =  hamming2
# results_report['Jaccord'] =  jaccard2
# results_report['Kulsinski'] =  kulsin
# results_report['Rogers tanimoto'] =  roger
# results_report['Russel Rao'] =  rus
# results_report['Sokal'] =  sokal2
# results_report['Sokal Sneath'] =  sokal3
# results_report['Yule'] =  yule2
# results_report['Hausdorff'] =  haus
# results_report['Square Chord'] =  sq_cho
#results_report['Entropy'] =  ent
results_report['Dot product'] =  dotttt
    
output_path = output_directoy + '\\'+ sample + '_results_test_peterbs.csv'
with open(output_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        results_report.to_csv(filec,index=False)    
