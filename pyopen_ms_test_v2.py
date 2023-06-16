# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 16:14:07 2023

@author: lawashburn
"""

from pyopenms import *
import time
start = time.time()

spectrum_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\perfect_spectra_mzML\20180524_PO_DDA_top20_TR3.mzML"
peptide_sequence = 'YKIFEPLR(Amidated)'
scan = 8785
fragment_mass_tolerance = 5.0
is_tol_in_ppm = True

tsg = TheoreticalSpectrumGenerator()
thspec = MSSpectrum()
print(thspec)
p = Param()
p.setValue("add_metainfo", "true")
tsg.setParameters(p)
peptide = AASequence.fromString(peptide_sequence)
tsg.getSpectrum(thspec, peptide, 1, 1)
# Iterate over annotated ions and their masses
# for ion, peak in zip(thspec.getStringDataArrays()[0], thspec):
#     print(ion, peak.getMZ())

e = MSExperiment()
MzMLFile().load(spectrum_path, e)
spectrum_of_interest = e[scan-1]
print(spectrum_of_interest)
print("Spectrum native id", spectrum_of_interest.getNativeID())
# mz, i = spectrum_of_interest.get_peaks()
# peaks = [(mz, i) for mz, i in zip(mz, i) if i > 1000 and mz > 50]
# for peak in peaks:
#     print(peak[0], "mz", peak[1], "int")

hscore = HyperScore()

result = hscore.compute(
    fragment_mass_tolerance, is_tol_in_ppm, spectrum_of_interest, thspec
)
print(result)

end = time.time()
print(end - start)