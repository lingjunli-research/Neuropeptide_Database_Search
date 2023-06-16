import csv
import pandas as pd
import os

working_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\Ref_MS2"
output_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\Ref_MS2_formatted"

def get_file_names_with_strings_list(str_list): #definition for finding a file containing a string in filename in specified directory
    full_list = os.listdir(working_directory)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

query = '.ms2' #search for ion list pertaining to the sequence
ion_file = (get_file_names_with_strings_list([query])) #search for the file based on query

# select_ion = ion_file[0]
# ion_file = []
# ion_file.append(select_ion)

for a in ion_file:
    out_its = []
    out_exps = []
    print(a)
    sample_name = a[:-4]
    print(sample_name)
    path = working_directory + '\\' + a
    print(path)
    with open(path) as input:
        lst = [line.strip() for line in input]

    new_list= []
    final_lst = []
    final_lst.append(['m/z', 'Resolution', 'Charge', 'Intensity', 'MS2','Scan_Number','Precursor_Charge', 'RetTime'])
    #final_lst.append(['m/z', 'Resolution', 'Charge', 'Intensity','Scan_Number', 'Ion_Injection_Time'])
    ms2_list = []

    new = lst
    #print(new)
       
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
    ion_inject_time_list = []

    for i in header_list:
        if 'S' in i:
            scan_number_list.append(i[1])
        if 'Z' in i:
            precursor_charge_list.append(i[1])
        if 'RetTime' in i:
            ion_inject_time_list.append(i[2])
        
    for i in range(len(new_list)):
        #print(i, new_list[i])
        if 'ActivationType' in new_list[i]:
            seperation_list.append(i-1)
        if 'PrecursorInt' in new_list[i]:
            seperation_list.append(i+2)
        if 'S' in new_list[i]:
            scan_number_list.append(new_list[i][1])
        if 'Z' in new_list[i]:
            precursor_charge_list.append(new_list[i][1])
        if 'RetTime' in new_list[i]:
            ion_inject_time_list.append(new_list[i][2])

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
    ion_injection_time_index = 0

    for element in new_list:
        if element == '-':
           # ms2_list_index+=1
            scan_number_index+=1
            #precursor_charge_index+=1
            ion_injection_time_index +=1
            continue   
        element.append(ms2_list[ms2_list_index])
        element.append(scan_number_list[scan_number_index])
        element.append(precursor_charge_list[precursor_charge_index])
        element.append(ion_inject_time_list[ion_injection_time_index])
    
        final_lst.append(element)

    out_path = output_directory + '\\' + sample_name + '.txt'

    with open(out_path,'w') as output:
        for i in final_lst:
            for j in i:
                output.write(str(j + ' '))
                output.write('\n')

    out_suf = '_ms2_output_list.txt'
    out_name = working_directory + '\\'  + out_suf
    with open(out_name,'w') as output:
        for i in final_lst:
            for j in i:
                output.write(str(j + ' '))
            output.write('\n')
    csv_suf = '_ms2_output_list.csv'
    csv_name =  working_directory + '\\'  + out_suf
    read_txt = pd.read_csv(out_name)
    read_txt.to_csv(csv_name, index=None)
    read_csv_path = csv_name
    csv_out_path = output_directory  + '\\' + sample_name + '.csv'
    SIM_result_imp = pd.read_csv(read_csv_path, sep = ' ')
    #with open(csv_out_path,'w',newline='') as file:
    #    writer = csv.writer(file)
    #    SIM_result_imp.to_csv(file,index=False)
    SIM_result = pd.DataFrame()
    SIM_result['m/z'] = SIM_result_imp['m/z']
    SIM_result['resolution'] = SIM_result_imp['Resolution']
    SIM_result['charge'] = SIM_result_imp['Charge']
    SIM_result['intensity'] = SIM_result_imp['Intensity']
    SIM_result['precursor m/z'] = SIM_result_imp['MS2']
    SIM_result['Scan #'] = SIM_result_imp['Scan_Number']
    SIM_result['Precursor Charge'] = SIM_result_imp['Precursor_Charge']
    SIM_result['RT'] = SIM_result_imp['RetTime']
    SIM_result = SIM_result[pd.to_numeric(SIM_result['m/z'], errors='coerce').notnull()]
    csv_path = working_directory + '\\'  + csv_suf
    csv_out_path2 = output_directory  + '\\' + sample_name + '.csv'
    #with open(csv_out_path2,'w',newline='') as filec:
    #        writerc = csv.writer(filec)
    #        SIM_result.to_csv(filec,index=False)
    csv_path = working_directory + '\\'  + csv_suf
    csv_out_path2 = output_directory  + '\\' + sample_name + '.txt'
    with open(csv_out_path2,'w',newline='') as filec:
            writerc = csv.writer(filec)
            SIM_result.to_csv(filec,index=False)
    print('Analysis Complete')
