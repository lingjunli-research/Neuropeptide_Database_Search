# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 13:11:34 2023

@author: 15156
"""
#packages 
import pandas as pd
import math
import csv
import pdb 
import numpy as np
from random import sample
import statistics 
import matplotlib.pyplot as plt
working_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\opt_predict_test\20230301"
TBY_path = r"D:\Tina\random_fragment_ion\Theoretical_b_y_ion_APQRNFLRF(Amidated).txt"
predicted = r"D:\Tina\random_fragment_ion\APQRNFLRF(Amidated)_8233_theo_rep.csv"
experimental_path = r"D:\Tina\random_fragment_ion\APQRNFLRF(Amidated)_8233_exp_rep.csv"

hmass = 1.00784
bin_tolerance = 0.02
rounds = 5
exp_intensity = 27
fragment_charge = [1,2,3,4] #list of z 
max_scores = [93000,86750]

##main code##
theoreticalBY = pd.read_csv(TBY_path) #Reading txt file 
predicted_intensity = pd.read_csv(predicted) #reading theoretical DF 
experimental_spectra = pd.read_csv(experimental_path)
experimental_spectra['experimental intensity'] = exp_intensity
experimental_spectra = experimental_spectra.sort_values(by='experimental m/z')

theoreticalBY.drop([16,41,50], axis=0, inplace=True)#removing MH from theoretical 

theoreticalBY_charges = pd.DataFrame()
for a in fragment_charge:
    theoreticalBY['m/z'] = (theoreticalBY['mass'] + (a * hmass)) / a
    theoreticalBY['m/z'] = theoreticalBY['m/z'].round(4)
    theoreticalBY['charge'] = a
    theoreticalBY_charges = pd.concat([theoreticalBY_charges,theoreticalBY])

#merging columns i'm comparing 
mass_charge_merged = predicted_intensity.merge(theoreticalBY_charges, left_on = 'theoretical m/z', right_on = 'm/z',how='outer')

#create empty df that will have ions identities 
merge_match = pd.DataFrame()
ions_in_order = mass_charge_merged['ion'].values.tolist()
ions_in_order_no_nan = [x for x in ions_in_order if str(x) != 'nan'] #remove nan from list
updated_list = []
#filling in where nan used to be with ion identity 
for aa in  ions_in_order_no_nan:
    for zz in range(0,4):
        updated_list.append(aa) #append each entry 5 times (0,1,2,3,4)
merge_match['ion'] = updated_list #use list to add values in DF for ion identity 

# replace NaN values in column 'ion' with the list, without overwriting non-NaN values
mass_charge_merged['ion'] = mass_charge_merged['ion'].apply(lambda x: x if pd.notna(x) else updated_list.pop(0) if len(updated_list) > 0 else x)

#dropping last 3 columns
mass_charge_merged.drop(['mass', 'm/z', 'charge'], axis=1, inplace=True)

#getting rid of ion types with h20 or nh3
mass_charge_merged[['ion type', 'delete']] = mass_charge_merged['ion'].str.split('-', 1, expand=True)
#drop column that isnt ion type without adducts 
mass_charge_merged.drop(['delete', 'ion'], axis=1, inplace=True)
#getting ion to a list
ion_list = mass_charge_merged['ion type'].tolist()
#remove duplicates in ion_list 
ion_nodup = list(set(ion_list)) #got 16 entries yay


run_storage = []
list_of_scores_storage = []
average_storage = []
std_storage = []

for xx in range(1,rounds+1):
        cor_score = []
        for jj in range(0,xx):
            #creating pseudo spectra 
            pseudo_spectra = pd.DataFrame() #empty df to add stuff to later 
            
            for x in ion_nodup:
                mass_charge_merged2 = mass_charge_merged[mass_charge_merged['ion type'] == x]
                pseudo = mass_charge_merged2.sample(1)
                pseudo_spectra = pd.concat([pseudo_spectra,pseudo])
            
            #filter by column substring
            pseudo_b = pseudo_spectra[pseudo_spectra['ion type'].str.contains("b")]
            pseudo_y = pseudo_spectra[pseudo_spectra['ion type'].str.contains("y")]    
            
            #getting rid of b and y before number on ion type column
            pseudo_y = pseudo_y.copy()
            pseudo_y['num'] = pseudo_y['ion type'].str[1:]
            pseudo_b = pseudo_b.copy()
            pseudo_b['num'] = pseudo_b['ion type'].str[1:]
            pseudo_y = pseudo_y.sort_values(by=['num'])
            #converting low to high y ions into list
            low_high_y = pseudo_y['num'].values.tolist()
            high_low_y = list(reversed(low_high_y))
            #copy list over to pseudo y df 
            pseudo_y['num'] = high_low_y
            #sort b by ascending order
            pseudo_b = pseudo_b.sort_values(by=['num'])
            
            #concat y and b into new df 
            pseudo_by_concat = pd.concat([pseudo_b,pseudo_y])
            #export aa to list 
            aa_num = pseudo_by_concat['num'].values.tolist()
            #remove duplicates from aa list 
            aa_nodupe = list(set(aa_num))
            
            
            pseudo_spec1 = pd.DataFrame() #empty df to add stuff to later 
            
            for x in aa_nodupe:
                by_con = pseudo_by_concat[pseudo_by_concat['num'] == x]
                pseudo2 = by_con.sample(1)
                pseudo_spec1 = pd.concat([pseudo_spec1,pseudo2])
            
            #write over pseudo_spec1 intensity to = 50 
            pseudo_spec1['theoretical intensity'] = 50 
            pseudo_spec1 = pseudo_spec1.rename(columns={"theoretical intensity": "pseudo intensity"})
            
            predicted_intensity = predicted_intensity.sort_values(by=['theoretical m/z'])
            pseudo_spec1 = pseudo_spec1.sort_values(by=['theoretical m/z'])
            
            xcorr_array = pd.merge_asof(predicted_intensity,pseudo_spec1,left_on='theoretical m/z',
                                                        right_on='theoretical m/z', tolerance=bin_tolerance,allow_exact_matches=True, direction='nearest')
                # replace nan with 0 for placeholder for new list 
            xcorr_array2 = pd.merge_asof(xcorr_array, experimental_spectra,left_on='theoretical m/z',
                                                        right_on='experimental m/z', tolerance=bin_tolerance,allow_exact_matches=True, direction='nearest')

            xcorr_array2 = xcorr_array2.replace(np.nan, 0)
            
            experimental_array = []
            
            pseudo_intensity_report = xcorr_array2['pseudo intensity'].values.tolist()
            experimental_intensity_report = xcorr_array2['experimental intensity'].values.tolist()
            
            for i in range(0,len(experimental_intensity_report)):
                pseudo_rep = pseudo_intensity_report[i]
                exp_rep = experimental_intensity_report[i]
                
                if pseudo_rep > exp_rep:
                    experimental_array.append(pseudo_rep)
                elif exp_rep > pseudo_rep:
                    experimental_array.append(exp_rep)
                else:
                    experimental_array.append(pseudo_rep)

            theoretical_array = xcorr_array['theoretical intensity'].values.tolist()
            #experimental_array = xcorr_array['pseudo intensity'].values.tolist()
            
            correlation_score = np.dot(theoretical_array,experimental_array)
            cor_score.append(correlation_score)
        avg = np.mean(cor_score)
        std = np.std(cor_score)
        
        run_storage.append(xx)
        list_of_scores_storage.append(cor_score)
        average_storage.append(avg)
        std_storage.append(std)
            

results_df = pd.DataFrame()
results_df['Run #'] = run_storage
results_df['Score list'] = list_of_scores_storage
results_df['Average Score'] = average_storage
results_df['Std Score'] = std_storage
## make a plot to quickly illustrate results ##
mean_1 = np.array(average_storage)
std_1 = np.array(std_storage)
#%%
x = np.arange(len(mean_1))
plt.plot(x, mean_1, '#420D09', label='calculated mean')
plt.fill_between(x, mean_1 - std_1, mean_1 + std_1, color='#A45A52', alpha=0.2)

for a in max_scores:
    plt.axhline(y=a,linewidth=4, color='#949494',label=('observed mean ' + str(max_scores.index(a)+1)))

plt.xticks(np.arange(0, len(x)+1, 10))
plt.legend()
plt.show()