'''plots a replication of mass spectrum as reported by a chemical ionisation mass spectrometer (CIMS)'''
# simulation results are represented graphically

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import retr_out
import numpy as np
import scipy.constants as si
import importlib

def plotter_CIMS(dir_path, mcr_res, tn, iont, sens_func):

	# inputs: -----------------
	# dir_path - path to folder containing results files to plot
	# mcr_res - resolution of molar mass to charge ratio (g/mol/charge)
	# tn - time through experiment to plot at (s)
	# iont - type of ionisation
	# sens_func - sensitivity to molar mass function
	# ---------------------------

	# retrieve results, note that num_sb (number of size bins)
	# includes wall if wall turned on
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_MW, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, H2Oi, _, _, _, RO2i) = retr_out.retr_out(dir_path)
	
	# convert to 2D numpy array
	y_MW = np.array((y_MW)).reshape(-1, 1)
	
	# get index of time wanted
	ti = (np.where(np.abs(timehr-tn/3600.) == np.min(np.abs(timehr-tn/3600.))))[0][0]
	
	# convert yrec from 1D to 2D with times in rows, then select time wanted
	yrec = (yrec.reshape(len(timehr), num_comp*(num_sb+1)))[ti, :]
	
	# get gas-phase concentrations (ppt, note starting concentration is ppb)
	gp = yrec[0:num_comp]*1.e3  # conversion to ug/m3 (if wanted): /Cfac[ti]/si.N_A*y_MW[:, 0]*1.e12
	
	# get particle-phase concentrations (molecules/cm3)
	pp = yrec[num_comp:num_comp*(num_sb+1-wall_on)]
	
	# sum each component over size bins (molecules/cm3)
	pp = np.sum(pp.reshape(num_sb-wall_on, num_comp), axis=0)
	
	# convert to ppt
	pp = (pp/si.N_A)*Cfac[ti]*1.e6 # or convert to ug/m3: *y_MW[:, 0]*1.e12
	
	# correct for sensitivity to molar mass
	fac_per_comp = write_sens2mm(0, sens_func, y_MW)

	gp = gp*fac_per_comp[:, 0]
	pp = pp*fac_per_comp[:, 0]
	
	# remove water
	gp = np.append(gp[0:H2Oi], gp[H2Oi+1::])
	pp = np.append(pp[0:H2Oi], pp[H2Oi+1::])
	y_MW = np.append(y_MW[0:H2Oi, 0], y_MW[H2Oi+1::, 0])

	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
	ax0.semilogy(y_MW, gp, '+m', markersize = 10,  label = str('gas-phase'))
	ax0.semilogy(y_MW, pp, 'xb', markersize = 10, label = str('particle-phase'))
	
	ax0.set_title(str('Mass spectrum at ' + str(timehr[ti]) + ' hours'), fontsize = 14)
	ax0.set_xlabel(r'Mass/charge (Th)', fontsize = 14)
	ax0.set_ylabel(r'Concentration (ppt)', fontsize = 14)
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	ax0.legend(fontsize = 14)
	plt.show()

	return()

# function for plotting weighting of particles by age due to instrument response time
def write_sens2mm(caller, sens_func, y_MW):

	import datetime

	# inputs: --------------------------
	# caller - flag for the calling function
	# sens_func - the sensitivity (Hz/ppt) function to molar mass (g/mol)
	# y_MW - molar mass of components (g/mol)
	# ------------------------------------

	# create new  file - will contain module response time weighting function
	f = open('PyCHAM/sens2mm.py', mode='w')

	f.write('\'\'\'solving the sensitivity (Hz/ppt) of instrument to molar mass (g/mol)\'\'\'\n')
	f.write('# module to estimate the sensitivity of an instrument to the molar mass of components, for example a Chemical Ionisiation Mass Spectrometer\n')
	f.write('# File Created at %s\n' %(datetime.datetime.now()))	
	f.write('\n')
	f.write('import numpy as np\n')
	f.write('\n')
	f.write('# function for sensitivity\n')
	f.write('def sens2mm(caller, y_MW):\n')
	f.write('	\n')
	f.write('	# inputs: -----------------\n')
	f.write('	# caller - flag for the calling function\n')
	f.write('	# y_MW - molar mass (g/mol) of components in question\n')
	f.write('	# ---------------------------\n')
	f.write('	\n')
	f.write('	fac_per_comp = %s # sensitivity (Hz/ppt) per molar mass (g/mol) \n' %(sens_func))
	f.write('	fac_per_comp = np.array((fac_per_comp)).reshape(-1, 1) # reshape \n')
	f.write('	\n')
	f.write('	if (caller == 3): # called on to plot sensitivity to molar mass\n')
	f.write('		import matplotlib.pyplot as plt \n')
	f.write('		plt.ion()\n')
	f.write('		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))\n')
	f.write('		ax0.plot(y_MW, fac_per_comp)\n')
	f.write('		ax0.set_title(\'Sensitivity of instrument to molar mass\')\n')
	f.write('		ax0.set_ylabel(\'Sensitivity (fraction (0-1))\', size = 14)\n')
	f.write('		ax0.yaxis.set_tick_params(labelsize = 14, direction = \'in\', which = \'both\')\n')
	f.write('		ax0.set_xlabel(\'Molar Mass ($\mathrm{g\,mol^{-1}}$)\', fontsize=14)\n')
	f.write('		ax0.xaxis.set_tick_params(labelsize = 14, direction = \'in\', which = \'both\')\n')
	f.write('	\n')
	f.write('	return(fac_per_comp)')
	f.close() # close file
	
	# get sensitivity for each component
	import sens2mm
	importlib.reload(sens2mm)
	fac_per_mass = sens2mm.sens2mm(caller, y_MW)

	return(fac_per_mass)