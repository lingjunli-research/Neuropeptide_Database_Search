# -*- coding: utf-8 -*-
"""
Created on Thu May 18 11:09:46 2023

@author: lawashburn
"""

import csv
import pandas as pd
import os


working_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\validation\Nhu_RAW_files"
output_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\validation\Nhu_RAW_files" #output folder

def get_file_names_with_strings_list(str_list): #definition for finding a file containing a string in filename in specified directory
    full_list = os.listdir(working_directory)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

query = '.ms2' #search for ion list pertaining to the sequence
ion_file = (get_file_names_with_strings_list([query])) #search for the file based on query

for a in ion_file:
    print(a)
    MS2_path = working_directory + '\\' + a
    tissue_type = a.replace('.ms2','')
    with open(MS2_path) as input:
        lst = [line.strip() for line in input]
    
    new_list= []
    final_lst = []
    final_lst.append(['m/z', 'resolution', 'charge', 'intensity', 'MS2','scan_number','precursor_charge'])
    ms2_list = []
    
    new = lst
    
    for i in new:
        new_list.append(i.split())
        if '@' in i:
            x = i.split()
            for y in x:
                if '@' in y:
                    ms2 = y[0:y.index('@')]
                    ms2_list.append(str(ms2))
    
    header_list = new_list[0:26]
    new_list = new_list[26:] # starts from line 26 to remove the first few header lines so that program could proceed
    seperation_list = []
    scan_number_list = []    
    precursor_charge_list = []
    
    for i in header_list:
        if 'S' in i:
            scan_number_list.append(i[1])
        if 'Z' in i:
            precursor_charge_list.append(i[1])
    
    for i in range(len(new_list)):
        #print(i, new_list[i])
        if 'RetTime' in new_list[i]:
            seperation_list.append(i-1)
        if 'PrecursorInt' in new_list[i]:
            seperation_list.append(i+2)
        if 'S' in new_list[i]:
            scan_number_list.append(new_list[i][1])
        if 'Z' in new_list[i]:
            precursor_charge_list.append(new_list[i][1])
    
    seperation_pairs = []
    start = 0
    for i in range(int(len(seperation_list)/2)):
        seperation_pairs.append((seperation_list[i+start],seperation_list[i+start+1]))
        start +=1 
     
    update_index = 0
    for start,end in seperation_pairs:
        start += update_index
        end += update_index
        new_list[start:end] = '-'
        update_index -= (end-start-1)
    
    ms2_list_index = 0
    scan_number_index = 0
    precursor_charge_index = 0
    
    for element in new_list:
        if element == '-':
            ms2_list_index+=1
            scan_number_index+=1
            precursor_charge_index+=1
            continue   
        element.append(ms2_list[ms2_list_index])
        element.append(scan_number_list[scan_number_index])
        element.append(precursor_charge_list[precursor_charge_index])
        final_lst.append(element)
    
    out_name = working_directory + '\\'+tissue_type+'_formatted.txt'
    with open(out_name,'w') as output:
        for i in final_lst:
            for j in i:
                output.write(str(j + ','))
            output.write('\n')