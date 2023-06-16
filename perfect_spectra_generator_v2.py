# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 13:22:44 2023

@author: lawashburn
"""

import csv
import pandas as pd

spectra_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\Ref_MS2_formatted\20180524_PO_DDA_top10_TR1_180525095121.txt"
output_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\perfect_spectra_formatted_w_RT"

old_scan_number = [6461,9444,9577]



dict_scans = {6461:9,9444:10,9577:11}

spectra = pd.read_csv(spectra_path, sep=",",skiprows=[0], names= ["m/z","resolution","charge","intensity", "precursor m/z","Scan #","Precursor Charge","RT"])

spectra_filtered = spectra[spectra['Scan #'].isin(old_scan_number)]


spectra_filtered2=spectra_filtered.replace({"Scan #": dict_scans})

csv_out_path2 = output_directory  + '\\perfect_spectra.txt'
with open(csv_out_path2,'a',newline='') as filec:
        writerc = csv.writer(filec)
        spectra_filtered2.to_csv(filec,index=False,header=None)
