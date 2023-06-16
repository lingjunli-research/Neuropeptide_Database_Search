# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 17:01:14 2023

@author: lawashburn
"""
import csv
import pandas as pd
import re
import os
from itertools import permutations
from Bio.SeqIO.FastaIO import SimpleFastaParser
import random
import collections
import time
import numpy as np
from scipy import signal
from datetime import datetime
import scipy
from scipy.spatial import distance
import numpy as np
from pyopenms import *

output_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\PEAKS_oursoftware_compare_brain_only\DB_search_optimize\03_precursor_AMM_optimize_v2\2021_0817_Brain_1"
unfiltered_psm_path = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\PEAKS_oursoftware_compare_brain_only\sample_input_data\unfiltered_psms.csv"
mzml_path = r"C:\Users\lawashburn\Documents\DBpep_v2\validation\Nhu_RAW_files\2021_0817_Brain_1.mzML"



def hyperscore_calc(ss,peptide):        
    ss = int(ss)
    tsg = TheoreticalSpectrumGenerator()
    thspec = MSSpectrum()
    p = Param()
    p.setValue("add_metainfo", "true")
    tsg.setParameters(p)
    peptide = AASequence.fromString(peptide)
    tsg.getSpectrum(thspec, peptide, 1, 1)

    spectrum_of_interest = e[ss-1]
    print("Spectrum native id", spectrum_of_interest.getNativeID())
    mz, i = spectrum_of_interest.get_peaks()
    peaks = [(mz, i) for mz, i in zip(mz, i) if i > 1000 and mz > 50]

    hscore = HyperScore()
    fragment_mass_tolerance = 20.0
    is_tol_in_ppm = True
    result = hscore.compute(
        fragment_mass_tolerance, is_tol_in_ppm, spectrum_of_interest, thspec)

    return result

e = MSExperiment()
MzMLFile().load(mzml_path, e)

unfiltered_psm = pd.read_csv(unfiltered_psm_path)
unfiltered_psm_scan_only = unfiltered_psm.drop_duplicates(subset='Scan')

correlation_results = []

for scan in range(0,len(unfiltered_psm_scan_only)):

    scan_isolate = unfiltered_psm_scan_only.iloc[[scan]]
    scan_value = scan_isolate['Scan'].values[0]
    scan_filtered_full_results = unfiltered_psm[unfiltered_psm['Scan'] == scan_value]

    for peptide in range(0,len(scan_filtered_full_results)):
        peptide_isolate = scan_filtered_full_results.iloc[[peptide]]
        peptide_value = peptide_isolate['Sequence'].values[0]
        correlation_value = hyperscore_calc(scan_value,peptide_value)
        peptide_isolate['Correlation value'] = correlation_value
        correlation_results.append(peptide_isolate)
        
all_correlation_results = pd.concat(correlation_results,ignore_index=True)   

output_path_rep = output_directory + '\\correlation_results.csv'

with open(output_path_rep,'w',newline='') as filec:
        writerc = csv.writer(filec)
        all_correlation_results.to_csv(filec,index=False)
