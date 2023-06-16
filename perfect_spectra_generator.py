# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 13:15:10 2023

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
start = time.time()

##User input##
output_parent_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs" #folder in which all output directories will be generated
db_path = r"C:\Users\lawashburn\Desktop\ALC50_Mass_Search_Files\duplicate_removed_crustacean_database_validated_formatted20220725.fasta" #database fasta path
base_file_path = r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\formatted_MS2"
key = r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\perfect_spectra_list.csv"

raw_converter_path_input = [r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\formatted_MS2\Brain1_formatted.txt",
                            r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\formatted_MS2\Brain2_formatted.txt",
                            r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\formatted_MS2\Brain3_formatted.txt",
                            r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\formatted_MS2\CoG1_formatted.txt",
                            r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\formatted_MS2\CoG2_formatted.txt",
                            r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\formatted_MS2\CoG3_formatted.txt",
                            r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\formatted_MS2\PO1_formatted.txt",
                            r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\formatted_MS2\PO2_formatted.txt",
                            r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\formatted_MS2\PO3_formatted.txt",
                            r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\formatted_MS2\SG1_formatted.txt",
                            r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\formatted_MS2\SG2_formatted.txt",
                            r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\formatted_MS2\SG3_formatted.txt",
                            r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\formatted_MS2\TG1_formatted.txt",
                            r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\formatted_MS2\TG2_formatted.txt",
                            r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230207\NV_Training_Spectra_w_PTMs\formatted_MS2\TG3_formatted.txt"]

source_key = pd.read_csv(key)

spectra_out = pd.DataFrame()

for file in raw_converter_path_input:
    raw_converter = pd.read_csv(file, sep=",",skiprows=[0], names= ["m/z","resolution","charge","intensity", "MS2", "scan_number","precursor_charge",'null'])

    
    file_name = file.replace(base_file_path,'')
    file_name = file_name.replace('_formatted.txt','')
    file_name = file_name.replace('\\','')
    
    raw_converter['Sample'] = file_name
    
    key_filter = source_key[source_key['Sample'] == file_name]
    scan_list = key_filter['Scan'].values.tolist()
    
    spectra_filtered = raw_converter[raw_converter['scan_number'].isin(scan_list)]
    
    spectra_out = pd.concat([spectra_out,spectra_filtered])

output_path = output_parent_directory + '\\training_spectra.txt'
spectra_export = spectra_out.to_csv(output_path, header=["m/z","resolution","charge","intensity", "MS2", "scan_number","precursor_charge",'null','Sample'], index=None, sep=',', mode='a')
    