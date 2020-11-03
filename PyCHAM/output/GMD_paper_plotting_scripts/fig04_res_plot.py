'''code to make plot that illustrates effect of temporal resolution on gas-phase photochemistry using the alpha-pinene ozonolysis in presence of NOx simulation using for photo_chem_res_plot'''
# introduction
# use the AtChem2_apinene_scheme.txt in the Results folder of the GMD_paper for chemical 
# scheme input
# use the xml file from the PyCHAM inputs folder
# use the Photo_chem_inputs_hiNOx.txt as Model Variables inputs

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from retrieve_PyCHAM_outputs import retrieve_outputs as retr


# ----------------------------------------------------------------------------------------

# file name
# empty dictionary of results from each simulation
num_sb_dict = {}
num_speci_dict = {}
Cfac_dict = {}
y_dict = {}
N_dict = {}
sbb_dict = {}
x_dict = {}
thr_dict = {}

# get current working directory
cwd = os.getcwd()
# 60 s intervals
(num_sb_dict['num_sb0'], num_speci_dict['num_speci0'], Cfac_dict['Cfac0'], y_dict['y0'], N_dict['N0'], sbb_dict['sbb0'], x_dict['x0'], thr_dict['thr0'], _, _, _, _, _, _) = retr(str(cwd + '/photo_chem_data/PyCHAM_time_res/PyCHAM_time_res60s'))
# 600 s intervals
(num_sb_dict['num_sb1'], num_speci_dict['num_speci1'], Cfac_dict['Cfac1'], y_dict['y1'], N_dict['N1'], sbb_dict['sbb1'], x_dict['x1'], thr_dict['thr1'], _, _, _, _, _, _) = retr(str(cwd + '/photo_chem_data/PyCHAM_time_res/PyCHAM_time_res600s'))
# 6000 s intervals
(num_sb_dict['num_sb2'], num_speci_dict['num_speci2'], Cfac_dict['Cfac2'], y_dict['y2'], N_dict['N2'], sbb_dict['sbb2'], x_dict['x2'], thr_dict['thr2'], _, _, _, _, _, _) = retr(str(cwd + '/photo_chem_data/PyCHAM_time_res/PyCHAM_time_res6000s'))

# make plot with gas-phase concentration deviations shown --------------------------------
compnum = 0 # count on components
fig, (ax0) = plt.subplots(1, 1, figsize=(8,6))


# alpha-pinene
ax0.plot(thr_dict['thr1'], (y_dict['y1'][:,312]-y_dict['y0'][0::10,312])/np.max(np.abs(y_dict['y0'][0::10,312]))*100.0, '-k', label=r'$\mathrm{\alpha}$-pinene - $\mathrm{6x10^{2}\, s}$')
ax0.plot(thr_dict['thr2'][0:-1], (y_dict['y2'][0:-1,312]-y_dict['y0'][0::100,312])/np.max(np.abs(y_dict['y0'][0::100,312]))*100.0, '--k' , label=r'$\mathrm{\alpha}$-pinene - $\mathrm{6x10^{3}\, s}$')
# O3
ax0.plot(thr_dict['thr1'], (y_dict['y1'][:,1]-y_dict['y0'][0::10,1])/np.max(np.abs(y_dict['y0'][0::10,1]))*100.0, '-r', label=r'$\mathrm{O3}$ - $\mathrm{6x10^{2}\, s}$')
ax0.plot(thr_dict['thr2'][0:-1], (y_dict['y2'][0:-1,1]-y_dict['y0'][0::100,1])/np.max(np.abs(y_dict['y0'][0::100,1]))*100.0, '--r', label=r'$\mathrm{O3}$ - $\mathrm{6x10^{3}\, s}$')
# OH
ax0.plot(thr_dict['thr1'], (y_dict['y1'][:,7]-y_dict['y0'][0::10,7])/np.max(np.abs(y_dict['y0'][0::10,7]))*100.0, '-b', label=r'$\mathrm{OH}$ - $\mathrm{6x10^{2}\, s}$')
ax0.plot(thr_dict['thr2'][0:-1], (y_dict['y2'][0:-1,7]-y_dict['y0'][0::100,7])/np.max(np.abs(y_dict['y0'][0::100,7]))*100.0, '--b', label=r'$\mathrm{OH}$ - $\mathrm{6x10^{3}\, s}$')

ax0.set_ylabel(r'Deviation from $\mathrm{6x10^{1}\, s}$ (%)', fontsize=14)
ax0.set_xlabel(r'Time of day (hours)', fontsize=14)
ax0.yaxis.set_tick_params(size=14)
ax0.xaxis.set_tick_params(size=14)
ax0.legend(fontsize=12)
fig.savefig('fig05.png')
plt.show()