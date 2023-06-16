# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 15:46:12 2023

@author: lawashburn
"""


import scipy
from scipy.spatial import distance
import numpy as np
import pysptools
from pysptools import distance as dist_2
from distance_metrics_mcda import distance_metrics
import pandas as pd

csv_path = r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\array_test.csv"
array_input = pd.read_csv(csv_path)
array_input = array_input.fillna(0)

array_input = array_input.sample(frac=.5)

array1 = array_input['theoretical intensity'].values.tolist()
array2 = array_input['experimental intensity'].values.tolist()



euc_distance = distance.euclidean(array1, array2) #euclidean distance calculation
bray_curt_distance = distance.braycurtis(array1, array2) #bray-curtis distance calculation
dot_prod = np.dot(array1, array2) #dot product calculation
#normXcorr = dist_2.NormXCorr(array1, array2) #normalized cross correlation
sam_map = pysptools.distance.SAM(array1, array2) #Computes the spectral angle mapper between two vectors (in radians)
spectral_divergence = pysptools.distance.SID(array1, array2) #Computes the spectral information divergence between two vectors
chebychev_distance = distance.chebyshev(array1, array2) #Chebychev distance
canberra = distance.canberra(array1, array2) #canberra distance
nyc = distance.cityblock(array1, array2) #manhattan/city block distance
cos = distance.cosine(array1, array2) #cosine distance
jen_sha = distance.jensenshannon(array1, array2) #jensen-shannon distance
#mahala = distance.mahalanobis(array1, array2) #Mahalanobis distance
minkowski = distance.minkowski(array1, array2,p=1) #Minkowski distance
sq_euc = distance.sqeuclidean(array1, array2) #squared euclidean distance
dice = distance.dice(array1, array2) #Dice dissimilarity distance
hamming = distance.hamming(array1, array2) #Hamming distance
jaccard = distance.jaccard(array1, array2) #Jaccard-Needham dissimilary
kulsinski = distance.kulsinski(array1, array2) #Kulsinski dissimilarity
#kulczynski1 = distance.kulczynski1(array1, array2) #Kulczynski 1 dissimilarity
rogerstanimoto = distance.rogerstanimoto(array1, array2) #Rogers-Tanimoto dissimilarity
russelrao = distance.russellrao(array1, array2) #Russell-Rao dissimilarity
sokal = distance.sokalmichener(array1, array2) #Sokal-Michener dissimilarity
sokal_sneath = distance.sokalsneath(array1, array2) #Sokal-Sneath dissimilarity
yule = distance.yule(array1, array2) #Yule dissimilarity
hausdorff = distance_metrics.hausdorff(array1, array2) #Hausdorff distance
#lorentzian = distance_metrics.lorentzian(array1, array2) #Lorentzian distance
#bhattacharyya = distance_metrics.bhattacharyya(array1, array2) #bhattacharyya distance
#hellinger = distance_metrics.hellinger(array1, array2) #hellinger distance
#matusita = distance_metrics.matusita(array1, array2) #matusita distance
sq_chord = distance_metrics.squared_chord(array1, array2) #Squared-chrod distance
#pearson_chi_square = distance_metrics.pearson_chi_square(array1, array2) #pearson_chi_square distance
#squared_chi_square = distance_metrics.pearson_chi_square(array1, array2) #squared_chi_square distance
entropy = scipy.stats.entropy(array1, array2) #Shannon entropy calculation