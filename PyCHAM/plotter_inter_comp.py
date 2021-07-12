'''function for comparing simulated concentrations from multiple models'''
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


	# inputs : ----------------------------------------------
	# 
	# ---------------------------------------------------------

	# list containing components of interest
	comp_of_int = ['CRESOL', 'H2O2', 'NO', 'NO2', 'OH', 'RO2']

	# PyCHAM --------------------------------------
	dir_path = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v301/PyCHAM/PyCHAM/output/ic_chem_scheme/Flow_Reactor_gas_phase_Intercomparison_Ocres_10NOx'
	# get required information from PyCHAM
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_MW, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, H2Oi, _, _, _, RO2i) = retr_out.retr_out(dir_path)
	
	# reshape so that time in rows and components per size bin in columns
	PCrec = yrec.reshape(len(timehr), num_comp*(num_sb+1))
	# isolate just gas-phase concentrations (ppb)
	PCrec = PCrec[:, 0:num_comp]
	
	Cfac = (np.array(Cfac)).reshape(-1, 1)# convert to numpy array from list

	# FACSIMILE ----------------------------------
	dir_path = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v301/PyCHAM/PyCHAM/output/ic_chem_scheme/CRESOLwithNOx.dat'
	
	# get required information from facsimile
	[Ftime_s, Fcomp_names, FCrec, [], []] = retr_out.retr_out_noncsv(dir_path, comp_of_int)
	
	# convert FACSIMILE to 60 s intervals
	FCrec = np.append(FCrec[0::60, :], FCrec[-1, :].reshape(1, -1), axis = 0)
	Ftime_s = np.append(Ftime_s[0::60], 3600.)
	
	# EASY --------------------------------------------
	dir_path = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v301/PyCHAM/PyCHAM/output/ic_chem_scheme/data.CRESOL.CS.nc'
	
	# get required information from EASY
	[Etime_s, Ecomp_names, ECrec, [], []] = retr_out.retr_out_noncsv(dir_path, comp_of_int)
	
	# remove repetition of final time in EASY
	Etime_s = Etime_s[0:-1]
	ECrec = ECrec[0:-1, :]
	# -----------------------------------------------------------------
	
	
	# convert PyCHAM concentrations from ppb to molecules/cm3
	PCrec = PCrec*Cfac[0]

	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7)) # prepare plot

	for i in comp_of_int: # loop through components of interest
		
		Fi = Fcomp_names.index(i) # FACSIMILE index
		Ei = Ecomp_names.index(i) # EASY index
		
		
		if (i != 'RO2'): # individual components
			ax0.plot(Etime_s[1::]/3600., (ECrec[1::, Ei]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '-x', linewidth = 2., label = str('E_'+i))
			Pi = comp_names.index(i) # PyCHAM index
			ax0.plot(timehr[1::], (PCrec[1::, Pi]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '--+', linewidth = 2., label = str('P_'+i))
		if (i == 'RO2'): # sum of organic peroxy radical components
			ax0.plot(Etime_s[1::]/3600., (ECrec[1::, Ei]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '-^', linewidth = 2., label = str('E_'+i))
			PCrecn =  np.sum(PCrec[:, RO2i], axis=1)# PyCHAM index
			ax0.plot(timehr[1::], (PCrecn[1::]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '--^', linewidth = 2., label = str('P_'+i))
			
	# details of plot
	ax0.set_ylabel(r'% Deviation from FACSIMILE', fontsize = 14)
	ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
	ax0.set_title(r'Deviation from FACSIMILE for EASY and PyCHAM when O-cresol 20 ppb, NOx 10 ppb', fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.legend(fontsize = 14, loc = 'lower right')
	plt.show()
	
	
	# O-cresol no NOx ----------------------------------------------------
	
	# list containing components of interest
	comp_of_int = ['CRESOL', 'H2O2', 'OH', 'RO2']

	# PyCHAM --------------------------------------
	dir_path = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v301/PyCHAM/PyCHAM/output/ic_chem_scheme/Flow_Reactor_gas_phase_Intercomparison_Ocres_0NOx'
	# get required information from PyCHAM
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_MW, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, H2Oi, _, _, _, RO2i) = retr_out.retr_out(dir_path)
	
	# reshape so that time in rows and components per size bin in columns
	PCrec = yrec.reshape(len(timehr), num_comp*(num_sb+1))
	# isolate just gas-phase concentrations (ppb)
	PCrec = PCrec[:, 0:num_comp]
	
	Cfac = (np.array(Cfac)).reshape(-1, 1)# convert to numpy array from list

	# FACSIMILE ----------------------------------
	dir_path = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v301/PyCHAM/PyCHAM/output/ic_chem_scheme/CRESOLwithoutNOx_orig.dat'
	
	# get required information from facsimile
	[Ftime_s, Fcomp_names, FCrec, [], []] = retr_out.retr_out_noncsv(dir_path, comp_of_int)
	
	# convert FACSIMILE to 60 s intervals
	FCrec = np.append(FCrec[0::60, :], FCrec[-1, :].reshape(1, -1), axis = 0)
	Ftime_s = np.append(Ftime_s[0::60], 3600.)
	
	# EASY --------------------------------------------
	dir_path = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v301/PyCHAM/PyCHAM/output/ic_chem_scheme/data.CRESOLnoNOx.CS.nc'
	
	# get required information from EASY
	[Etime_s, Ecomp_names, ECrec, [], []] = retr_out.retr_out_noncsv(dir_path, comp_of_int)
	
	# convert EASY to 60 s intervals
	#ECrec = np.append(ECrec[0::60, :], ECrec[-1, :].reshape(1, -1), axis = 0)
	#Etime_s = np.append(Etime_s[0::60], 3600.)
	# remove repetition of final time in EASY
	Etime_s = Etime_s[0:-1]
	ECrec = ECrec[0:-1, :]
	
	# -----------------------------------------------------------------
	
	
	# convert PyCHAM concentrations from ppb to molecules/cm3
	PCrec = PCrec*Cfac[0]

	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7)) # prepare plot

	for i in comp_of_int: # loop through components of interest
		
		Fi = Fcomp_names.index(i) # FACSIMILE index
		Ei = Ecomp_names.index(i) # EASY index
		
		
		if (i != 'RO2'): # individual components
			ax0.plot(Etime_s[1::]/3600., (ECrec[1::, Ei]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '-x', linewidth = 2., label = str('E_'+i))
			Pi = comp_names.index(i) # PyCHAM index
			ax0.plot(timehr[1::], (PCrec[1::, Pi]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '--+', linewidth = 2., label = str('P_'+i))
		if (i == 'RO2'): # sum of organic peroxy radical components
			ax0.plot(Etime_s[1::]/3600., (ECrec[1::, Ei]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '-^', linewidth = 2., label = str('E_'+i))
			PCrecn =  np.sum(PCrec[:, RO2i], axis=1) # PyCHAM index
			ax0.plot(timehr[1::], (PCrecn[1::]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '--^', linewidth = 2., label = str('P_'+i))
			
	# details of plot
	ax0.set_ylabel(r'% Deviation from FACSIMILE', fontsize = 14)
	ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
	ax0.set_title(r'Deviation from FACSIMILE for EASY and PyCHAM when O-cresol 20 ppb, NOx 0 ppb', fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.legend(fontsize = 14, loc = 'lower right')
	plt.show()
	
	# alpha-pinene with NOx ----------------------------------------------------
	
	# list containing components of interest
	comp_of_int = ['APINENE', 'H2O2', 'NO', 'NO2', 'OH', 'RO2']

	# PyCHAM --------------------------------------
	dir_path = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v301/PyCHAM/PyCHAM/output/ic_chem_scheme/Flow_Reactor_gas_phase_Intercomparison_APINENE_20NOx'
	# get required information from PyCHAM
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_MW, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, H2Oi, _, _, _, RO2i) = retr_out.retr_out(dir_path)
	
	# reshape so that time in rows and components per size bin in columns
	PCrec = yrec.reshape(len(timehr), num_comp*(num_sb+1))
	# isolate just gas-phase concentrations (ppb)
	PCrec = PCrec[:, 0:num_comp]
	
	Cfac = (np.array(Cfac)).reshape(-1, 1)# convert to numpy array from list

	# FACSIMILE ----------------------------------
	dir_path = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v301/PyCHAM/PyCHAM/output/ic_chem_scheme/APINENEwithNOx.dat'
	
	# get required information from facsimile
	[Ftime_s, Fcomp_names, FCrec, [], []] = retr_out.retr_out_noncsv(dir_path, comp_of_int)
	
	# convert FACSIMILE to 60 s intervals
	FCrec = np.append(FCrec[0::60, :], FCrec[-1, :].reshape(1, -1), axis = 0)
	Ftime_s = np.append(Ftime_s[0::60], 3600.)
	
	# EASY --------------------------------------------
	dir_path = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v301/PyCHAM/PyCHAM/output/ic_chem_scheme/data.APINENE.CS.nc'
	
	# get required information from EASY
	[Etime_s, Ecomp_names, ECrec, [], []] = retr_out.retr_out_noncsv(dir_path, comp_of_int)
	
	# convert EASY to 60 s intervals
	#ECrec = np.append(ECrec[0::60, :], ECrec[-1, :].reshape(1, -1), axis = 0)
	#Etime_s = np.append(Etime_s[0::60], 3600.)
	# remove repetition of final time in EASY
	Etime_s = Etime_s[0:-1]
	ECrec = ECrec[0:-1, :]
	# -----------------------------------------------------------------
	
	
	# convert PyCHAM concentrations from ppb to molecules/cm3
	PCrec = PCrec*Cfac[0]

	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7)) # prepare plot

	for i in comp_of_int: # loop through components of interest
		
		Fi = Fcomp_names.index(i) # FACSIMILE index
		Ei = Ecomp_names.index(i) # EASY index
		
		
		if (i != 'RO2'): # individual components
			ax0.plot(Etime_s[1::]/3600., (ECrec[1::, Ei]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '-x', linewidth = 2., label = str('E_'+i))
			Pi = comp_names.index(i) # PyCHAM index
			ax0.plot(timehr[1::], (PCrec[1::, Pi]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '--+', linewidth = 2., label = str('P_'+i))
		if (i == 'RO2'): # sum of organic peroxy radical components
			ax0.plot(Etime_s[1::]/3600., (ECrec[1::, Ei]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '-^', linewidth = 2., label = str('E_'+i))
			PCrecn =  np.sum(PCrec[:, RO2i], axis=1)# PyCHAM index
			ax0.plot(timehr[1::], (PCrecn[1::]-ECrec[1::, Ei])/ECrec[1::, Ei]*100., '--^', linewidth = 2., label = str('P_'+i))
			
	# details of plot
	ax0.set_ylabel(r'% Deviation from FACSIMILE', fontsize = 14)
	ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
	ax0.set_title(r'Deviation from FACSIMILE for EASY and PyCHAM when alpha-pinene 20 ppb, NOx 20 ppb', fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.legend(fontsize = 14, loc = 'lower right')
	plt.show()
	
	# alpha-pinene without NOx ----------------------------------------------------
	
	# list containing components of interest
	comp_of_int = ['APINENE', 'H2O2', 'OH', 'RO2']

	# PyCHAM --------------------------------------
	dir_path = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v301/PyCHAM/PyCHAM/output/ic_chem_scheme/Flow_Reactor_gas_phase_Intercomparison_APINENE_0NOx'
	# get required information from PyCHAM
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_MW, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, H2Oi, _, _, _, RO2i) = retr_out.retr_out(dir_path)
	
	# reshape so that time in rows and components per size bin in columns
	PCrec = yrec.reshape(len(timehr), num_comp*(num_sb+1))
	# isolate just gas-phase concentrations (ppb)
	PCrec = PCrec[:, 0:num_comp]
	
	Cfac = (np.array(Cfac)).reshape(-1, 1)# convert to numpy array from list

	# FACSIMILE ----------------------------------
	dir_path = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v301/PyCHAM/PyCHAM/output/ic_chem_scheme/APINENEwithoutNOx.dat'
	
	# get required information from facsimile
	[Ftime_s, Fcomp_names, FCrec, [], []] = retr_out.retr_out_noncsv(dir_path, comp_of_int)
	
	# convert FACSIMILE to 60 s intervals
	FCrec = np.append(FCrec[0::60, :], FCrec[-1, :].reshape(1, -1), axis = 0)
	Ftime_s = np.append(Ftime_s[0::60], 3600.)
	
	# EASY --------------------------------------------
	dir_path = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v301/PyCHAM/PyCHAM/output/ic_chem_scheme/data.APINENEnoNOx.CS.nc'
	
	# get required information from EASY
	[Etime_s, Ecomp_names, ECrec, [], []] = retr_out.retr_out_noncsv(dir_path, comp_of_int)
	
	# convert EASY to 60 s intervals
	#ECrec = np.append(ECrec[0::60, :], ECrec[-1, :].reshape(1, -1), axis = 0)
	#Etime_s = np.append(Etime_s[0::60], 3600.)
	# remove repetition of final time in EASY
	Etime_s = Etime_s[0:-1]
	ECrec = ECrec[0:-1, :]
	# -----------------------------------------------------------------
	
	
	# convert PyCHAM concentrations from ppb to molecules/cm3
	PCrec = PCrec*Cfac[0]

	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7)) # prepare plot

	for i in comp_of_int: # loop through components of interest
		
		Fi = Fcomp_names.index(i) # FACSIMILE index
		Ei = Ecomp_names.index(i) # EASY index
		
		
		if (i != 'RO2'): # individual components
			ax0.plot(Etime_s[1::]/3600., (ECrec[1::, Ei]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '-x', linewidth = 2., label = str('E_'+i))
			Pi = comp_names.index(i) # PyCHAM index
			ax0.plot(timehr[1::], (PCrec[1::, Pi]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '--+', linewidth = 2., label = str('P_'+i))
		if (i == 'RO2'): # sum of organic peroxy radical components
			ax0.plot(Etime_s[1::]/3600., (ECrec[1::, Ei]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '-^', linewidth = 2., label = str('E_'+i))
			PCrecn =  np.sum(PCrec[:, RO2i], axis=1)# PyCHAM index
			ax0.plot(timehr[1::], (PCrecn[1::]-ECrec[1::, Ei])/ECrec[1::, Ei]*100., '--^', linewidth = 2., label = str('P_'+i))
			
	# details of plot
	ax0.set_ylabel(r'% Deviation from FACSIMILE', fontsize = 14)
	ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
	ax0.set_title(r'Deviation from FACSIMILE for EASY and PyCHAM when alpha-pinene 20 ppb, NOx 0 ppb', fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.legend(fontsize = 14, loc = 'lower right')
	plt.show()

	# alpha-pinene dark with N2O5 ----------------------------------------------------
	
	# list containing components of interest
	comp_of_int = ['APINENE', 'H2O2', 'OH', 'RO2', 'N2O5']

	# PyCHAM --------------------------------------
	dir_path = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v301/PyCHAM/PyCHAM/output/ic_chem_scheme/Flow_Reactor_gas_phase_Intercomparison_APINENE_20N2O5_dark'
	# get required information from PyCHAM
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_MW, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, H2Oi, _, _, _, RO2i) = retr_out.retr_out(dir_path)
	
	# reshape so that time in rows and components per size bin in columns
	PCrec = yrec.reshape(len(timehr), num_comp*(num_sb+1))
	# isolate just gas-phase concentrations (ppb)
	PCrec = PCrec[:, 0:num_comp]
	
	Cfac = (np.array(Cfac)).reshape(-1, 1)# convert to numpy array from list

	# FACSIMILE ----------------------------------
	dir_path = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v301/PyCHAM/PyCHAM/output/ic_chem_scheme/APINENElightoff.dat'
	
	# get required information from facsimile
	[Ftime_s, Fcomp_names, FCrec, [], []] = retr_out.retr_out_noncsv(dir_path, comp_of_int)
	
	# convert FACSIMILE to 60 s intervals
	FCrec = np.append(FCrec[0::60, :], FCrec[-1, :].reshape(1, -1), axis = 0)
	Ftime_s = np.append(Ftime_s[0::60], 3600.)
	
	# EASY --------------------------------------------
	dir_path = '/Users/Simon_OMeara/Documents/Manchester/postdoc/box/PyCHAM_v301/PyCHAM/PyCHAM/output/ic_chem_scheme/data.APINENEdark.CS.nc'
	
	# get required information from EASY
	[Etime_s, Ecomp_names, ECrec, [], []] = retr_out.retr_out_noncsv(dir_path, comp_of_int)
	
	# convert EASY to 60 s intervals
	#ECrec = np.append(ECrec[0::60, :], ECrec[-1, :].reshape(1, -1), axis = 0)
	#Etime_s = np.append(Etime_s[0::60], 3600.)
	# remove repetition of final time in EASY
	Etime_s = Etime_s[0:-1]
	ECrec = ECrec[0:-1, :]
	# -----------------------------------------------------------------
	
	
	# convert PyCHAM concentrations from ppb to molecules/cm3
	PCrec = PCrec*Cfac[0]

	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7)) # prepare plot

	for i in comp_of_int: # loop through components of interest
		
		Fi = Fcomp_names.index(i) # FACSIMILE index
		Ei = Ecomp_names.index(i) # EASY index
		
		
		if (i != 'RO2'): # individual components
			ax0.plot(Etime_s[1::]/3600., (ECrec[1::, Ei]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '-x', linewidth = 2., label = str('E_'+i))
			#ax0.plot(Etime_s[1::]/3600., (ECrec[1::, Ei]), '-x', linewidth = 2., label = str('E_'+i))
			#ax0.plot(Ftime_s[1::]/3600., (FCrec[1::, Fi]), '-x', linewidth = 2., label = str('F_'+i))

			Pi = comp_names.index(i) # PyCHAM index
			ax0.plot(timehr[1::], (PCrec[1::, Pi]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '--+', linewidth = 2., label = str('P_'+i))
			#ax0.plot(timehr[1::], (PCrec[1::, Pi]), '--+', linewidth = 2., label = str('P_'+i))
		if (i == 'RO2'): # sum of organic peroxy radical components
			ax0.plot(Etime_s[1::]/3600., (ECrec[1::, Ei]-FCrec[1::, Fi])/FCrec[1::, Fi]*100., '-^', linewidth = 2., label = str('E_'+i))
			PCrecn =  np.sum(PCrec[:, RO2i], axis=1)# PyCHAM index
			ax0.plot(timehr[1::], (PCrecn[1::]-ECrec[1::, Ei])/ECrec[1::, Ei]*100., '--^', linewidth = 2., label = str('P_'+i))
			
	# details of plot
	ax0.set_ylabel(r'% Deviation from FACSIMILE', fontsize = 14)
	ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
	ax0.set_title(r'Deviation from FACSIMILE for EASY and PyCHAM when alpha-pinene 20 ppb, N2O5 20 ppb', fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.legend(fontsize = 14, loc = 'lower right')
	plt.show()

		
plotter_inter_comp() # run function if module called