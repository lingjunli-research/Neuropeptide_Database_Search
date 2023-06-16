# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 14:30:42 2023

@author: lawashburn
"""
import csv
import pandas as pd

spectra_path = r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\training_spectra.txt"

spectra = pd.read_csv(spectra_path, sep=",",skiprows=[0], names= ["m/z","resolution","charge","intensity", "MS2", "scan_number","precursor_charge",'null','Sample'])

spectra["Identifier"] = spectra['Sample'].astype(str) +"_"+ spectra["scan_number"].astype(str)

identifiers = set(spectra['Identifier'].values.tolist())

all_peaks = pd.DataFrame()

for identifier in identifiers:
    single_scan = spectra[spectra['Identifier'] == identifier]
    number_peaks = len(single_scan)
    for x in range(100,-1,-1):
        percentage = x / 100
        peaks_to_keep = round(percentage*number_peaks)
        single_scan = single_scan.sample(n=peaks_to_keep)
        single_scan['Iteration'] = x
        all_peaks = pd.concat([all_peaks,single_scan])
output_path =  r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\iterated_spectra.txt"
spectra_export = all_peaks.to_csv(output_path, header=["m/z","resolution","charge","intensity", "MS2", "scan_number","precursor_charge",'null','Sample','Identifier','Iteration'], index=None, sep=',', mode='a')        