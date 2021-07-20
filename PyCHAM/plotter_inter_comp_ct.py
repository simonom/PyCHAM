'''function for comparing simulated concentrations of individual components from multiple models'''
# for example, PyCHAM, FACSIMILE and EASY

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import os
import retr_out
import numpy as np
import scipy.constants as si

# define function
def plotter_inter_comp():

	
	
	# PyCHAM --------------------------------------

	# list containing components of interest
	comp_of_int = ['HO2', 'NO2']

	# list containing pythonic reaction number of interest
	reac_of_int = [24]
	
	# path to results
	dir_path = '/Users/Psymo/Documents/PyCHAM_vW/PyCHAM/PyCHAM/output/ic_chem_scheme/Flow_Reactor_gas_phase_Intercomparison_APINENE_20NOx_HO2'
	
	# get required information from PyCHAM
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_MW, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, H2Oi, _, _, cham_env, RO2i) = retr_out.retr_out(dir_path)
	

	# reshape so that time in rows and components per size bin in columns
	PCrec = yrec.reshape(len(timehr), num_comp*(num_sb+1))
	# isolate just gas-phase concentrations (ppb)
	PCrec = PCrec[:, 0:num_comp]
	
	Cfac = (np.array(Cfac)).reshape(-1, 1)# convert to numpy array from list

	# get change tendency (molecules/cm3/s)
	ct_of_int = ['HO2'] # component with change tendency of interest
	fname = str(dir_path+ '/' + ct_of_int[0] +'_rate_of_change')
	dydt = np.loadtxt(fname, delimiter = ',', skiprows = 1) # skiprows = 1 omits header	

	# note that penultimate column in dydt is gas-particle 
	# partitioning and final column is gas-wall partitioning, whilst
	# the first row contains chemical reaction numbers
		
	# prepare to store results of change tendency due to chemical reactions
	res = np.zeros((dydt.shape[0], dydt.shape[1]-2))
	res[:, :] = dydt[:, 0:-2] # get chemical reaction numbers and change tendencies (molecules/cm3/s)	

	# convert from molecules/cm3/s to ppb/cc/s
	#res[1::, :] = res[1::, :]/Cfac.reshape(-1, 1)	

	# get index of reaction interested in
	reac_indx = (res[0, :] == reac_of_int)
	
	# get NO2 gas-phase concentration (molecules/cm3)
	for i in range(len(comp_names)):
		if comp_names[i] == comp_of_int[0]:
			comp0 = ((yrec[:, i]).reshape(-1, 1))*Cfac
		if comp_names[i] == comp_of_int[1]:
			comp1 = ((yrec[:, i]).reshape(-1, 1))*Cfac

	# EASY ----------------------------------------------	
	# list containing components of interest
	comp_of_int = ['HO2', 'KMT09']
	dir_path = '/Users/Psymo/Documents/PyCHAM_vW/PyCHAM/PyCHAM/output/ic_chem_scheme/EASY/data.APINENE.CS.nc'
	
	# get required information from EASY
	[Etime_s, Ecomp_names, ECrec, [], []] = retr_out.retr_out_noncsv(dir_path, comp_of_int)
	
	# remove repition of final time from EASY
	ECrec = ECrec[0:-1, :]
	Etime_s = Etime_s[0:-1]

	# get reaction rate (molecules/cm3/s)
	#EASY_rr = ECrec[:, 0]#*ECrec[:, 1]

	# plot ---------------------------------------------- 
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7)) # prepare plot

	# reaction rate coefficients
	ax0.plot(timehr*3600., -2*res[1::, reac_indx]/(comp0*comp1), label = 'P_rrc')
	ax0.plot(Etime_s, ECrec[:, 1], '--', label = 'E_rrc')	

	# components
	#ax0.plot(timehr*3600., comp0, label = 'P_')
	#ax0.plot(timehr*3600., comp1, label = 'P_O')
	#ax0.plot(Etime_s, ECrec[:, 0], label = 'E_')
	#ax0.plot(Etime_s, ECrec[:, 1], label = 'E_O')

			
	# details of plot
	ax0.set_ylabel(r'Flux (ppb/cc/s)', fontsize = 14)
	ax0.set_xlabel(r'Time through simulation (seconds)', fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.set_title(r'EASY and PyCHAM when alpha-pinene 20 ppb, NOx 20 ppb, lit', fontsize = 14)
	ax0.legend(fontsize = 14, loc = 'lower right')
	plt.show()
	
plotter_inter_comp() # run function if module called