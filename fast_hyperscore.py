# -*- coding: utf-8 -*-
"""
Created on Fri May 12 15:54:03 2023

@author: lawashburn
"""
from scipy.spatial import distance
import numpy as np
from pyopenms import *
import time


mzml_path_input = r"C:\Users\lawashburn\Documents\DBpep_v2\finale_weighting_hyperscore\weighting_input_data\database_search_input\20190326_SG_Unlabeled_DDA_TR3.mzML"
# scan = 3819
# peptide = 'AVLLPKKTEKK'

scan_peptide = [[28749,'E(Glu->pyro-Glu)STNWMSSLRSAW']]

e = MSExperiment()
MzMLFile().load(mzml_path_input, e)

for aa in scan_peptide:
    start = time.time()
    scan = aa[0]
    peptide = aa[1]

    
    
    ss = int(scan)
    tsg = TheoreticalSpectrumGenerator()
    
    
    p = Param()
    p.setValue("add_metainfo", "true")
    tsg.setParameters(p)
    thspec = MSSpectrum()
    
    peptide = AASequence.fromString(peptide)
    
    tsg.getSpectrum(thspec, peptide, 1, 1)
    
    
    
    
    spectrum_of_interest = e[ss-1]
    
    
    mz, i = spectrum_of_interest.get_peaks()
    
    peaks = [(mz, i) for mz, i in zip(mz, i) if i > 1000 and mz > 50]
    
    
    hscore = HyperScore()
    fragment_mass_tolerance = 5.0
    is_tol_in_ppm = True
    
    result = hscore.compute(
        fragment_mass_tolerance, is_tol_in_ppm, spectrum_of_interest, thspec)
    
    end = time.time()
    print('Elapsed time: ' + str(end-start))
    print(result)