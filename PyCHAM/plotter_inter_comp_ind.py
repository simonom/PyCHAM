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

	# list containing components of interest
	comp_of_int = ['APINENE', 'H2O2', 'NO', 'NO2', 'OH', 'RO2']
	
	# PyCHAM --------------------------------------
	#dir_path = '/Users/Psymo/Documents/PyCHAM_vW/PyCHAM/PyCHAM/output/ic_chem_scheme/PyCHAM/Flow_Reactor_gas_phase_Intercomparison_APINENE_20N2O5_dark_400s'
	# get required information from PyCHAM
	#(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
	#	y_MW, _, comp_names, y_MV, _, wall_on, space_mode, 
	#	_, _, _, PsatPa, OC, H2Oi, _, _, _, group_indx) = retr_out.retr_out(dir_path)
	
	#RO2i = group_indx['RO2i']

	# reshape so that time in rows and components per size bin in columns
	#PCrec = yrec.reshape(len(timehr), num_comp*(num_sb+1))
	# isolate just gas-phase concentrations (ppb)
	#PCrec = PCrec[:, 0:num_comp]
	
	#Cfac = (np.array(Cfac)).reshape(-1, 1)# convert to numpy array from list

	# PyCHAM --------------------------------------
	dir_path = '/Users/Psymo/Documents/PyCHAM_vW/PyCHAM/PyCHAM/output/ic_chem_scheme/PyCHAM/Flow_Reactor_gas_phase_Intercomparison_APINENE_20N2O5_dark'
	# get required information from PyCHAM
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_MW, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, H2Oi, _, _, _, group_indx) = retr_out.retr_out(dir_path)
	
	RO2i = group_indx['RO2i']

	# reshape so that time in rows and components per size bin in columns
	PCrec = yrec.reshape(len(timehr), num_comp*(num_sb+1))
	# isolate just gas-phase concentrations (ppb)
	PCrec = PCrec[:, 0:num_comp]
	
	Cfac = (np.array(Cfac)).reshape(-1, 1)# convert to numpy array from list

	# FACSIMILE ----------------------------------
	dir_path = '/Users/Psymo/Documents/PyCHAM_vW/PyCHAM/PyCHAM/output/ic_chem_scheme/FACSIMILE/APINENEdarkNewresults.dat'

	# get required information from facsimile
	[Ftime_s, Fcomp_names, FCrec, [], []] = retr_out.retr_out_noncsv(dir_path, comp_of_int)
	
	# convert FACSIMILE to 60 s intervals
	FCrec = np.append(FCrec[0::60, :], FCrec[-1, :].reshape(1, -1), axis = 0)
	Ftime_s = np.append(Ftime_s[0::60], 3600.)
	
	# EASY --------------------------------------------
	dir_path = '/Users/Psymo/Documents/PyCHAM_vW/PyCHAM/PyCHAM/output/ic_chem_scheme/EASY/data.APINENEdark.CS.nc'
	
	# get required information from EASY
	[Etime_s, Ecomp_names, ECrec, [], []] = retr_out.retr_out_noncsv(dir_path, comp_of_int)

	# convert EASY to 60 s intervals
	#ECrec = np.append(ECrec[0::60, :], ECrec[-1, :].reshape(1, -1), axis = 0)
	#Etime_s = np.append(Etime_s[0::60], 3600.)

	# remove repition of final time from EASY
	ECrec = ECrec[0:-1, :]
	Etime_s = Etime_s[0:-1]
	# -----------------------------------------------------------------
	
	
	# convert PyCHAM concentrations from ppb to molecules/cm3
	PCrec = PCrec*Cfac[0]
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7)) # prepare plot
	
	for i in comp_of_int: # loop through components of interest
		
		Fi = Fcomp_names.index(i) # FACSIMILE index
		Ei = Ecomp_names.index(i) # EASY index
		
		
		if (i != 'RO2'): # individual components
			
			Pi = comp_names.index(i) # PyCHAM index
			Ti = FCrec[:, Fi] > 0. # allowed values

			# plot absolute values			
			#ax0.plot(Ftime_s, (FCrec[:, Fi]), '-x', linewidth = 2., label = str('F_'+i))
			#ax0.plot(timehr*3600., (PCrec[:, Pi]), '-x', linewidth = 2., label = str('P_hires_'+i))
			#ax0.plot(timehr2*3600., (PCrec2[:, Pi2]), '-x', linewidth = 2., label = str('P_lores_'+i))
			#ax0.plot(Etime_s, (ECrec[:, Ei]), '--x', linewidth = 2., label = str('E_'+i))
			
			# plot deviation
			ax0.plot(Ftime_s[Ti]/3600., ((PCrec[Ti, Pi]-FCrec[Ti, Fi])/FCrec[Ti, Fi])*100., '-x', linewidth = 2., label = str('P_'+i))
			ax0.plot(Etime_s[Ti]/3600., ((ECrec[Ti, Ei]-FCrec[Ti, Fi])/FCrec[Ti, Fi])*100., '-x', linewidth = 2., label = str('E_'+i))
			
		if (i == 'RO2'): # sum of organic peroxy radical components
			Ti = FCrec[:, Fi] > 0. # allowed values
			PCrecn =  np.sum(PCrec[:, RO2i], axis=1)# PyCHAM index

			# plot absolute values			
			#ax0.plot(Ftime_s, (FCrec[:, Fi]), '-x', linewidth = 2., label = str('F_'+i))
			#ax0.plot(timehr2*3600., (PCrecn), '-x', linewidth = 2., label = str('P_'+i))
			#ax0.plot(Etime_s, (ECrec[:, Ei]), '-x', linewidth = 2., label = str('E_'+i))

			# plot deviation
			ax0.plot(Etime_s[Ti]/3600., ((ECrec[Ti, Ei]-FCrec[Ti, Fi])/FCrec[Ti, Fi])*100., '--x', linewidth = 2., label = str('E_'+i))
			ax0.plot(Ftime_s[Ti]/3600., ((PCrecn[Ti]-FCrec[Ti, Fi])/FCrec[Ti, Fi])*100., '--^', linewidth = 2., label = str('P_'+i))
			
			

	# details of plot
	#ax0.set_ylabel(r'Concentration (molecules/cm3)', fontsize = 10)
	ax0.set_ylabel(r'Deviation ((((PyCHAM or EASY)-FACSIMILE)/FACSIMILE)*100) (%)', fontsize = 10)
	ax0.set_xlabel(r'Time through simulation (s)', fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	#ax0.set_title(r'Deviation from FACSIMILE for EASY (E_) and PyCHAM (P_) when alpha-pinene 20 ppb, 20 ppb N2O5, 0 ppb H2O2', fontsize = 13)
	ax0.legend(fontsize = 14, loc = 'lower right')
	plt.show()
	
plotter_inter_comp() # run function if module called