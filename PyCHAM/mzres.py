'''solving probability density function of mass:charge resolution'''
# module to estimate the probability density function that is demonstrative of an instrument's mass:charge resolution, for example a Chemical Ionisiation Mass Spectrometer
# File Created at 2021-09-20 12:00:11.313553

import numpy as np
import scipy.stats as st

# function for sensitivity
def mzres(caller, res_in, y_mw):
	
	# inputs: -----------------
	# caller - flag for the calling function
	# res_in - inputs for the mass:charge resolution (start point and width descriptors)
	# y_mw - molar mass (g/mol) of components in question
	# ---------------------------
	
	if (caller == 3): # called on to plot
		import matplotlib.pyplot as plt 
		plt.ion()
		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
	y_mw = np.array((y_mw)) # ensure numpy array rather than list
	comp_indx = [] # empty list to hold results
	comp_prob = [] # empty list to hold results
	maxmm = np.max(y_mw) + res_in[0]
	mm_acc = 0. + 1.0 # count on accumulated molar mass (g/mol) 
	# get maximum probability possible
	pdfm = st.norm.pdf(mm_acc, mm_acc, 0.3)
	# loop through until upper end of molar mass reached 
	while (mm_acc < maxmm):
		pdf = st.norm.pdf(y_mw, mm_acc, 0.3)
		try: # in case a maximum can be identified
			pdf = pdf/pdfm # ensure that probability at distribution peak is 1
			# minimum and maximum molar masses covered significantly by this resolution interval
			mm = [np.min(y_mw[pdf>1.e-2]), np.max(y_mw[pdf>1.e-2])]
			if (caller == 1): # if called from non-test module
				if (len(y_mw[pdf > 1.e-2]) > 0): # if components do contribute to this interval
					# store indices of components contributing to this mass:charge interval
					ci = (np.where((y_mw >= mm[0])*(y_mw <= mm[1]) == 1))[0][:]
					comp_indx.append(ci)
					# store probability of contribution to this resolution interval
					comp_prob.append(pdf[ci])
				else: # if components do not contribute to this interval
					comp_indx.append([])
					comp_prob.append([])
		except: # no maximum, so no components contributing to this interval
			if (caller == 1):
				comp_indx.append([])
				comp_prob.append([])
		if (caller == 3): # called on to plot
			ax0.plot(y_mw[pdf>1.e-2], pdf[pdf>1.e-2])
		mm_acc += res_in[0] # keep count on accumulated molar mass 
	
	if (caller == 3): # called on to plot
		ax0.set_title('Sensitivity of instrument due to mass:charge resolution')
		ax0.set_ylabel('Probability of inclusion in resolution interval (fraction (0-1))', size = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
		ax0.set_xlabel('Molar Mass ($\mathrm{g\,mol^{-1}}$)', fontsize=14)
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	
	# remember the range of molar masses representing mass:charge resolution
	mm_all = np.arange((0. + res_in[0]), (mm_acc), res_in[0]) 
	return(pdf, comp_indx, comp_prob, mm_all)