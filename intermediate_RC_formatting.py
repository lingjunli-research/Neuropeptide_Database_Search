# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 11:09:18 2023

@author: lawashburn
"""

import csv
import pandas as pd
import os

working_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\perfect_spectra_MS2_w_RT\v3"
output_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\perfect_spectra_MS2_w_RT\formatted\round2"

def get_file_names_with_strings_list(str_list): #definition for finding a file containing a string in filename in specified directory
    full_list = os.listdir(working_directory)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

query = '.txt' #search for ion list pertaining to the sequence
ion_file = (get_file_names_with_strings_list([query])) #search for the file based on query

# select_ion = ion_file[0]
# ion_file = []
# ion_file.append(select_ion)

for a in ion_file:
    raw_converter_path = working_directory+'\\'+a
    #raw_converter_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\perfect_spectra_MS2_w_RT\20180507_Brain_DDA_top20_TR2.txt"
    raw_converter = pd.read_csv(raw_converter_path, sep=",",skiprows=[0], names= ["m/z","resolution","charge","intensity", "Scan #","RT"])
    raw_converter_filtered = raw_converter[pd.to_numeric(raw_converter['m/z'], errors='coerce').notnull()]
    
    csv_out_path2 = output_directory  + '\\' + a
    with open(csv_out_path2,'w',newline='') as filec:
            writerc = csv.writer(filec)
            raw_converter_filtered.to_csv(filec,index=False)