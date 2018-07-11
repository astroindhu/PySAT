import numpy as np
import sys
from bokeh.plotting import output_file,show
from bokeh.layouts import column
from libpysat.plotting import plot_psl
import copy
sys.path.append("/Users/vara_id/Mercury/PlanetSpec/")
# from PlanetSpec.deconvolution import endmember
# from PlanetSpec.deconvolution import mixture
from libpysat.spectral.unmixing import Deconvolution_1D as decon

telescope={'IP_22E - Sprague et al (1994)': np.genfromtxt("/Users/vara_id/PhD/NASA IRTF MIRSI/MIR Mercury/telescope_spectra_screenshots_spectra/Spargue et al 1998/sprague1994_8Dec1990_22deg_scaled.csv", delimiter=",", dtype=float)}

endmember = {'Diopside':np.loadtxt("/Users/vara_id/PSL/emiss_minerals_mertis_analogues/emiss_00000756_diopside_0_125_padova_450c.txt")}
endmember['Labradorite']=np.loadtxt("/Users/vara_id/PSL/emiss_minerals_mertis_analogues/emiss_00001571_labradorite0_125_450C.txt")
endmember['Orthopyroxenite']=np.loadtxt("/Users/vara_id/PSL/emiss_minerals_mertis_analogues/emiss_00001734_orthopyroxenite_125_250_700K.txt")

for k in list(endmember.keys()):
    endmember[k]=np.delete(endmember[k],[0],axis=1)

for k in list(endmember.keys()):
    endmember[k]=endmember[k][(endmember[k][:,0]>=7)&(endmember[k][:,0]<=14)]

for k in list(telescope.keys()):
    telescope[k]=telescope[k][np.argsort(telescope[k][:,0])]

for k in list(telescope.keys()): #Normalising spectra for CF = 1 for deconvolution
    cf= telescope[k][(telescope[k][:,0]>=7) & (telescope[k][:,0]<=9)]
    cf_max = np.nanmax(cf[:,1])
    if cf_max >= 1:
        telescope[k][:,1]=telescope[k][:,1]-(cf_max-1)
    else:
        telescope[k][:,1]=telescope[k][:,1]+(1-cf_max)


result_path = "/Users/vara_id/PhD/NASA IRTF MIRSI/MIR Mercury/telescope_spectra_screenshots_spectra/"
result_fname = "TelSpec_decon_combination_of_three_endmembers"
unmixing_class = deconvolution.Deconvolution_1D(telescope,endmember,result_path, result_fname)
unmixing = unmixing_class.linear_unmixing_emi()

