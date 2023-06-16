# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 14:55:14 2023

@author: lawashburn
"""

from urllib.request import urlretrieve
from pyopenms import *
import pandas as pd
import csv
import time
import pickle

start = time.time()

spectra_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\perfect_spectra_mzML\20180524_PO_DDA_top20_TR3.mzML"
# sequences = ['YKIFEPLR(Amidated)','IPLRYEKF(Amidated)','LEFIKPYR(Amidated)',
#              'FKELPYIR(Amidated)','ELKRFPIY(Amidated)','RIFPKELY(Amidated)',
#              'ELKYIRFP(Amidated)','EIFPLYRK(Amidated)','FKILRPEY(Amidated)']
sequences = ['YKIFEPLR(Amidated)']
output_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\HyperScore_Evaluation"
target_sequence = 'YKIFEPLR(Amidated)'
exp_name = 'PO3_top20'
scan = 8785
#gh = "https://raw.githubusercontent.com/OpenMS/pyopenms-docs/master"
#urlretrieve(gh + "/src/data/SimpleSearchEngine_1.mzML", "searchfile.mzML")

# seq_list = []
# score_list = []
# tsg = TheoreticalSpectrumGenerator()
# thspec = MSSpectrum()
# p = Param()
# p.setValue("add_metainfo", "true")
# tsg.setParameters(p)

# for seq in sequences:
    
#     peptide = AASequence.fromString(seq)
#     tsg.getSpectrum(thspec, peptide, 1, 1)
#     e = MSExperiment()
#     MzMLFile().load(spectra_path, e)
#     spectrum_of_interest = e[scan-1]
#     print("Spectrum native id", spectrum_of_interest.getNativeID())
#     mz, i = spectrum_of_interest.get_peaks()
#     peaks = [(mz, i) for mz, i in zip(mz, i) if i > 1000 and mz > 50]
#     print(type(peaks))
#     with open('peaks.pkl', 'wb') as f:
#         pickle.dump(peaks, f)

with open('peaks.pkl', 'rb') as f:
    peaks = pickle.load(f)
#     #    print(peak[0], "mz", peak[1], "int")

    hscore = HyperScore()
    fragment_mass_tolerance = 20.0
    is_tol_in_ppm = True
    result = hscore.compute(
        fragment_mass_tolerance, is_tol_in_ppm, spectrum_of_interest, thspec
    )
    print(result)
#     seq_list.append(seq)
#     score_list.append(result)

# results_table = pd.DataFrame()
# results_table['Sequence'] = seq_list
# results_table['Score'] = score_list
# results_table

# file_path = output_directory + '\\' + exp_name + '_target_' + target_sequence + '_'+ str(scan) + '_hyperscore_results.csv'
# with open(file_path,'w',newline='') as filec:
#             writerc = csv.writer(filec)
#             results_table.to_csv(filec,index=False) 

# end = time.time()
# print(end - start)