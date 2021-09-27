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

def plotter_CIMS(dir_path, res_in, tn, iont, sens_func):
	
	# inputs: -----------------
	# dir_path - path to folder containing results files to plot
	# res_in - inputs for resolution of molar mass to charge ratio (g/mol/charge)
	# tn - time through experiment to plot at (s)
	# iont - type of ionisation
	# sens_func - sensitivity to molar mass function
	# ---------------------------

	# retrieve results, note that num_sb (number of size bins)
	# includes wall if wall turned on
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_MW, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, H2Oi, _, _, _, RO2i, _) = retr_out.retr_out(dir_path)
	
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
	
	gp = gp*fac_per_comp[:]
	pp = pp*fac_per_comp[:]

	# if ionisation source molar mass to be added 
	# (e.g. because not corrected for in measurment software), then add
	if (int(iont[1]) == 1):
		if (iont[0] == 'I'): # (https://pubchem.ncbi.nlm.nih.gov/compound/Iodide-ion)
			y_MW += 126.9045
		if (iont[0] == 'N'): # (https://pubchem.ncbi.nlm.nih.gov/compound/nitrate)
			y_MW += 62.005
	
	# remove water
	gp = np.append(gp[0:H2Oi], gp[H2Oi+1::])
	pp = np.append(pp[0:H2Oi], pp[H2Oi+1::])
	y_MW = np.append(y_MW[0:H2Oi, 0], y_MW[H2Oi+1::, 0])

	# account for mass to charge resolution
	[pdf, comp_indx, comp_prob, mm_all] = write_mzres(1, res_in, y_MW)
	gpres = np.zeros((len(comp_indx)))
	ppres = np.zeros((len(comp_indx)))
	
	for pdfi in range(len(comp_indx)): # loop through resolution intervals
		gpres[pdfi] = np.sum(gp[comp_indx[pdfi]]*comp_prob[pdfi])
		ppres[pdfi] = np.sum(pp[comp_indx[pdfi]]*comp_prob[pdfi])
		
	plt.ion() # disply plot in interactive mode

	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
	ax0.semilogy(mm_all, gpres, '+m', markersize = 14, markeredgewidth = 5,  label = str('gas-phase'))
	ax0.semilogy(mm_all, ppres, 'xb', markersize = 14, markeredgewidth = 5, label = str('particle-phase'))
	
	ax0.set_title(str('Mass spectrum at ' + str(timehr[ti]) + ' hours'), fontsize = 14)
	ax0.set_xlabel(r'Mass/charge (Th)', fontsize = 14)
	ax0.set_ylabel(r'Concentration (ppt)', fontsize = 14)
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	ax0.legend(fontsize = 14)
	
	return()

# function for plotting sensitivity to components
def write_sens2mm(caller, sens_func, y_MW):

	import datetime

	# inputs: --------------------------
	# caller - flag for the calling function
	# sens_func - the sensitivity (Hz/ppt) function to molar mass (g/mol)
	# y_MW - molar mass of components (g/mol)
	# ------------------------------------

	# create new  file
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
	f.write('	fac_per_comp = np.array((fac_per_comp)).reshape(-1) # reshape \n')
	f.write('	if (len(fac_per_comp) == 1): # if just a single value then tile across components \n')
	f.write('		fac_per_comp = np.tile(fac_per_comp, len(y_MW)) # if just a single value then tile across components \n')
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

# function for plotting mass:charge resolution
def write_mzres(caller, res_in, y_mw):

	import datetime

	# inputs: --------------------------
	# caller - flag for the calling function
	# res_in - inputs for the mass:charge resolution (start point and width descriptors)
	# y_mw - the molar mass range to describe (g/mol)
	# ------------------------------------

	# create new  file
	f = open('PyCHAM/mzres.py', mode='w')

	f.write('\'\'\'solving probability density function of mass:charge resolution\'\'\'\n')
	f.write('# module to estimate the probability density function that is demonstrative of an instrument\'s mass:charge resolution, for example a Chemical Ionisiation Mass Spectrometer\n')
	f.write('# File Created at %s\n' %(datetime.datetime.now()))	
	f.write('\n')
	f.write('import numpy as np\n')
	f.write('import scipy.stats as st\n')
	f.write('\n')
	f.write('# function for sensitivity\n')
	f.write('def mzres(caller, res_in, y_mw):\n')
	f.write('	\n')
	f.write('	# inputs: -----------------\n')
	f.write('	# caller - flag for the calling function\n')
	f.write('	# res_in - inputs for the mass:charge resolution (start point and width descriptors)\n')
	f.write('	# y_mw - molar mass (g/mol) of components in question\n')
	f.write('	# ---------------------------\n')
	f.write('	\n')
	f.write('	if (caller == 3): # called on to plot\n')
	f.write('		import matplotlib.pyplot as plt \n')
	f.write('		plt.ion()\n')
	f.write('		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))\n')
	f.write('	\n')
	f.write('	y_mw = np.array((y_mw)) # ensure numpy array rather than list\n')
	f.write('	comp_indx = [] # empty list to hold results\n')
	f.write('	comp_prob = [] # empty list to hold results\n')
	f.write('	maxmm = np.max(y_mw) + res_in[0]\n')
	f.write('	mm_acc = 0. + %s # count on accumulated molar mass (g/mol) \n' %(res_in[0]))
	f.write('	# get maximum probability possible\n')
	f.write('	pdfm = st.norm.pdf(mm_acc, mm_acc, %s)\n' %(res_in[1]))
	f.write('	# loop through until upper end of molar mass reached \n')
	f.write('	while (mm_acc < maxmm):\n')
	f.write('		pdf = st.norm.pdf(y_mw, mm_acc, %s)\n' %(res_in[1]))
	f.write('		try: # in case a maximum can be identified\n')
	f.write('			pdf = pdf/pdfm # ensure that probability at distribution peak is 1\n')
	f.write('			# minimum and maximum molar masses covered significantly by this resolution interval\n')
	f.write('			mm = [np.min(y_mw[pdf>1.e-2]), np.max(y_mw[pdf>1.e-2])]\n')
	
	f.write('			if (caller == 1): # if called from non-test module\n')
	f.write('				if (len(y_mw[pdf > 1.e-2]) > 0): # if components do contribute to this interval\n')
	f.write('					# store indices of components contributing to this mass:charge interval\n')
	f.write('					ci = (np.where((y_mw >= mm[0])*(y_mw <= mm[1]) == 1))[0][:]\n')
	f.write('					comp_indx.append(ci)\n')
	f.write('					# store probability of contribution to this resolution interval\n')
	f.write('					comp_prob.append(pdf[ci])\n')
	f.write('				else: # if components do not contribute to this interval\n')
	f.write('					comp_indx.append([])\n')
	f.write('					comp_prob.append([])\n')
	f.write('		except: # no maximum, so no components contributing to this interval\n')
	f.write('			if (caller == 1):\n')
	f.write('				comp_indx.append([])\n')
	f.write('				comp_prob.append([])\n')
	f.write('		if (caller == 3): # called on to plot\n')
	f.write('			ax0.plot(y_mw[pdf>1.e-2], pdf[pdf>1.e-2])\n')
	f.write('		mm_acc += res_in[0] # keep count on accumulated molar mass \n')
	f.write('	\n')
	f.write('	if (caller == 3): # called on to plot\n')
	f.write('		ax0.set_title(\'Sensitivity of instrument due to mass:charge resolution\')\n')
	f.write('		ax0.set_ylabel(\'Probability of inclusion in resolution interval (fraction (0-1))\', size = 14)\n')
	f.write('		ax0.yaxis.set_tick_params(labelsize = 14, direction = \'in\', which = \'both\')\n')
	f.write('		ax0.set_xlabel(\'Molar Mass ($\mathrm{g\,mol^{-1}}$)\', fontsize=14)\n')
	f.write('		ax0.xaxis.set_tick_params(labelsize = 14, direction = \'in\', which = \'both\')\n')
	f.write('	\n')
	f.write('	# remember the range of molar masses representing mass:charge resolution\n')
	f.write('	mm_all = np.arange((0. + res_in[0]), (mm_acc), res_in[0]) \n')
	f.write('	return(pdf, comp_indx, comp_prob, mm_all)')
	f.close() # close file
	
	# get resolution result for mass spectrum
	import mzres
	importlib.reload(mzres)
	[pdf, comp_indx, comp_prob, mm_all] = mzres.mzres(caller, res_in, y_mw)

	return(pdf, comp_indx, comp_prob, mm_all)