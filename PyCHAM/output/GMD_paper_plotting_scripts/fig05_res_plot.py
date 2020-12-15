'''code to make plot that illustrates effect of temporal resolution on gas-phase photochemistry using the alpha-pinene ozonolysis in presence of NOx simulation using for photo_chem_res_plot'''
# introduction
# use the fig03_scheme.txt in the Results folder of the GMD_paper for chemical 
# scheme input
# use the xml file from the PyCHAM inputs folder
# use the fig03_mod_var_hiNOx.txt as Model Variables inputs

# can call from the GMD paper Results folder of the PyCHAM home directory (if relevant data saved there)

import numpy as np
import sys
import os
import matplotlib.pyplot as plt



# ----------------------------------------------------------------------------------------

# get current working directory
cwd = os.getcwd()
import retr_out

try: # if calling from the GMD paper results folder

	# 60 s intervals
	# file name
	Pyfname = str(cwd + '/fig03_data/PyCHAM_time_res/PyCHAM_time_res60s')
	(num_sb, num_comp, Cfac, y0, Ndry, rbou_rec, xfm, thr0, PyCHAM_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(Pyfname)

	# 600 s intervals
	Pyfname = str(cwd + '/fig03_data/PyCHAM_time_res/PyCHAM_time_res600s')
	(num_sb, num_comp, Cfac, y1, Ndry, rbou_rec, xfm, thr1, PyCHAM_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(Pyfname)

	# 6000 s intervals
	Pyfname = str(cwd + '/fig03_data/PyCHAM_time_res/PyCHAM_time_res6000s')
	(num_sb, num_comp, Cfac, y2, Ndry, rbou_rec, xfm, thr2, PyCHAM_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(Pyfname)

except: # if calling from the PyCHAM home folder
	
	# 60 s intervals
	Pyfname =str(cwd+'/PyCHAM/output/GMD_paper_plotting_scripts/fig03_data/PyCHAM_time_res/PyCHAM_time_res60s')
	(num_sb, num_comp, Cfac, y0, Ndry, rbou_rec, xfm, thr0, PyCHAM_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(Pyfname)
	
	# 600 s intervals
	Pyfname =str(cwd+'/PyCHAM/output/GMD_paper_plotting_scripts/fig03_data/PyCHAM_time_res/PyCHAM_time_res600s')
	(num_sb, num_comp, Cfac, y1, Ndry, rbou_rec, xfm, thr1, PyCHAM_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(Pyfname)
		
	# 6000 s intervals
	Pyfname =str(cwd+'/PyCHAM/output/GMD_paper_plotting_scripts/fig03_data/PyCHAM_time_res/PyCHAM_time_res6000s')
	(num_sb, num_comp, Cfac, y2, Ndry, rbou_rec, xfm, thr2, PyCHAM_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(Pyfname)

# make plot with gas-phase concentration deviations shown --------------------------------
compnum = 0 # count on components
fig, (ax0) = plt.subplots(1, 1, figsize=(8, 6))


# alpha-pinene
ax0.plot(thr1, (y1[:, 312]-y0[0::10, 312])/np.max(np.abs(y0[0::10, 312]))*100.0, '-k', label=r'$\mathrm{\alpha}$-pinene - $\mathrm{6x10^{2}\, s}$')
ax0.plot(thr2[0:-1], (y2[0:-1, 312]-y0[0::100, 312])/np.max(np.abs(y0[0::100, 312]))*100.0, '--k' , label=r'$\mathrm{\alpha}$-pinene - $\mathrm{6x10^{3}\, s}$')
# O3
ax0.plot(thr1, (y1[:,1]-y0[0::10,1])/np.max(np.abs(y0[0::10,1]))*100.0, '-r', label=r'$\mathrm{O3}$ - $\mathrm{6x10^{2}\, s}$')
ax0.plot(thr2[0:-1], (y2[0:-1, 1]-y0[0::100, 1])/np.max(np.abs(y0[0::100, 1]))*100.0, '--r', label=r'$\mathrm{O3}$ - $\mathrm{6x10^{3}\, s}$')
# OH
ax0.plot(thr1, (y1[:, 7]-y0[0::10, 7])/np.max(np.abs(y0[0::10, 7]))*100.0, '-b', label=r'$\mathrm{OH}$ - $\mathrm{6x10^{2}\, s}$')
ax0.plot(thr2[0:-1], (y2[0:-1, 7]-y0[0::100, 7])/np.max(np.abs(y0[0::100,7]))*100.0, '--b', label=r'$\mathrm{OH}$ - $\mathrm{6x10^{3}\, s}$')

ax0.set_ylabel(r'Deviation from $\mathrm{6x10^{1}\, s}$ (%)', fontsize=14)
ax0.set_xlabel(r'Time of day (hours)', fontsize=14)
ax0.yaxis.set_tick_params(size=14)
ax0.xaxis.set_tick_params(size=14)
ax0.legend(fontsize=12)
fig.savefig('fig05.png')
plt.show()