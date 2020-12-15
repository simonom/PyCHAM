'''function to initiate concentrations of components'''
# based on inputs, initial concentrations and their holding arrays 
# are set

import numpy as np
import scipy.constants as si
import math
from water_calc import water_calc
import write_dydt_rec
import pybel


def init_conc(num_comp, Comp0, init_conc, TEMP, RH, PInit, Pybel_objects,
	testf, pconc, dydt_trak, end_sim_time, save_step, 
	rindx, pindx, num_eqn, nreac, nprod, 
	comp_namelist, Compt, seed_name, seed_mw,
	core_diss, nuc_comp, comp_xmlname, comp_smil):
		
	# inputs:------------------------------------------------------
	
	# num_comp - number of unique components
	# Comp0 - chemical scheme names of components present at start of experiment
	# init_conc - initial concentrations of components (molecules/cc)	
	# TEMP - temperature in chamber at start of experiment (K)
	# RH - relative humidity in chamber (dimensionless fraction 0-1)
	# PInit - initial pressure (Pa)
	# init_SMIL - SMILES of components present at start of experiment (whose 
	# concentrations are given in init_conc)
	# testf - flag for whether in normal mode (0) or testing mode (1/2)
	# pconc - initial concentration of particles (# particles/cc (air))
	# dydt_trak - chemical scheme name of components for which user wants the tendency to  
	#			change tracked
	# end_sim_time - total simulation time (s)
	# save_step - recording frequency (s)
	# num_eqn - number of equations
	# comp_namelist - list of components' names in chemical equation file
	# Compt - name of component injected after start of experiment
	# seed_name - name of core component (input by user)
	# seed_mw - molecular weight of seed material (g/mol)
	# core_diss - dissociation constant of seed material
	# nuc_comp - name of nucleating component (input by user, or defaults to 'core')
	# comp_xmlname - component names in xml file
	# comp_smil - SMILE strings in xml file
	# -----------------------------------------------------------

	if testf==1: # testing mode
		# return dummies
		return(0,0,0,0,0,0,0,0)

	NA = si.Avogadro # Avogadro's number (molecules/mol)
	# empty array for storing species' concentrations, must be an array
	y = np.zeros((num_comp))
	y_mw = np.zeros((num_comp, 1)) # species' molecular weight (g/mol)
	# empty array for storing index of interesting gas-phase components
	y_indx_plot = []
	
	# convert concentrations
	# total number of molecules in 1 cc air using ideal gas law.  R has units cc.Pa/K.mol
	ntot = PInit*(NA/(8.3144598e6*TEMP))
	# one billionth of number of molecules in chamber unit volume
	Cfactor = ntot*1.0e-9 # ppb-to-molecules/cc
	

	# prepare dictionary for tracking tendency to change of user-specified components
	dydt_vst = {}

	# insert initial concentrations where appropriate
	for i in range(len(Comp0)):
    		# index of where initial components occur in list of components
		try: # in case components already listed via interpretation of the chemical scheme
			y_indx = comp_namelist.index(Comp0[i])
			
		# if component not already listed via interpretation of the chemical scheme
		# then add to list
		except:
			print('Note: component specified in initial gas-phase concentrations of model variables has not been registered in the chemical scheme, but will try to be added to the component list inside the module for initialising gas-phase concentrations')
			comp_namelist.append(Comp0[i]) # add to name list
			
			# convert MCM chemical names to SMILES
			if (Comp0[i] in comp_xmlname):
				# index where xml file name matches reaction component name
				name_indx = comp_xmlname.index(Comp0[i])
				name_SMILE = comp_smil[name_indx] # SMILES of component
				# Generate pybel
				Pybel_object = pybel.readstring('smi', name_SMILE)
				# append to Pybel object list
				Pybel_objects.append(Pybel_object)
				
				comp_num += 1 # number of unique species
				
				y_indx = len(comp_namelist)-1
				y = y.append(np.zeros((1)))
			else:
				print(str('Error: inside eqn_parser, chemical scheme name '+str(name_only)+' not found in xml file'))
				sys.exit()
			
		y[y_indx] = init_conc[i]*Cfactor # convert from ppb to molecules/cc (air)
		# remember index for plotting gas-phase concentrations later
		y_indx_plot.append(y_indx)

	# number of recording steps
	nrec_steps = int(math.ceil(end_sim_time/save_step)+1)

	# get index of user-specified components for tracking their dydt due to model 
	# mechanisms
	if len(dydt_trak)>0:
		
		dydt_traki = [] # empty list for indices of these components
		
		for i in range (len(dydt_trak)):
			reac_index = [] # indices of reactions involving this species
			# index of where initial species occurs in SMILE string
			y_indx = comp_namelist.index(dydt_trak[i])

			# remember index for plotting gas-phase concentrations later
			dydt_traki.append(int(y_indx))
			# search through reactions to see where this component is reactant or product
			for ri in range(num_eqn):
				if sum(rindx[ri,0:nreac[ri]]==y_indx)>0:
					reac_index.append(int(ri)) # append reaction index
			for ri in range(num_eqn): # repeat for products
				if sum(pindx[ri,0:nprod[ri]]==y_indx)>0:
					reac_index.append(int(ri)) # append reaction index
					
	
			# save reaction indices in dictionary value for this component,
			# when creating empty rec_array, add two rows onto the end for particle- and 
			# wall-partitioning, respectively.  Note the extra row to hold the reacion 
			# indices
			rec_array = np.zeros((nrec_steps+1, len(reac_index)+2))
			rec_array[0, 0:-2] = (reac_index)
			dydt_vst[y_indx] = rec_array # dictionary entry to hold results

		dydt_vst['comp_index'] = dydt_traki # dictionary entry to hold component names
		
		# call on write_dydt_rec to generate the module that will process
		# the tendency to change during the simulation
		write_dydt_rec.write_dydt_rec()
	
	for i in range(num_comp): # loop through all components
		y_mw[i] = Pybel_objects[i].molwt # molecular weight (g/mol)
	
	# ------------------------------------------------------------------------------------
	# account for water's properties
	H2Oi = num_comp # index for water
	num_comp += 1 # update number of species to account for water
	
	# update gas-phase concentration (molecules/cc (air)) and vapour pressure
	# of water (log10(atm))
	[C_H2O, Psat_water, H2O_mw] = water_calc(TEMP, RH, si.N_A)
	
	# append empty element to y and y_mw to hold water values
	y = np.append(y, C_H2O)
	y_mw = (np.append(y_mw, H2O_mw)).reshape(-1, 1)
	comp_namelist.append('H2O') # append water's name to component name list

	# ------------------------------------------------------------------------------------
	# account for seed properties - note that even if no seed particle, this code ensures
	# that an index is provided for core material

	# empty array for index of core component
	seedi = (np.zeros((len(seed_name)))).astype(int)
	
	comp_namelist.append('core') # append name of core to component name list
	corei = [num_comp] # index for core component
	# increase number of components to account for 'core' component
	num_comp += 1

	# append core gas-phase concentration (molecules/cc (air)) and molecular 
	# weight (g/mol) (needs to have a 1 length in second dimension for the kimt 
	# calculations)
	y = np.append(y, 0.)
	y_mw = (np.append(y_mw, seed_mw)).reshape(-1, 1)
		
	# if nucleating component formed of core component
	if (nuc_comp[0] == 'core'):
		nuci = num_comp-1 # index of core component
	else:
		nuci = -1 # filler
	
	# get indicices of seed particle component(s)
	indx = 0 # count on seed component(s)
	for sname in seed_name:
		# index of core component
		seedi[indx] = int(comp_namelist.index(sname))
		indx += 1 # count on seed component(s)
	
	# get index of component with latter injections
	if len(Compt)>0:
		inj_indx = np.zeros((len(Compt)))
		for i in range(len(Compt)):
			# index of where initial species occurs in SMILE string
			inj_indx[i] = comp_namelist.index(Compt[i])
	else:
		inj_indx = np.zeros((1)) # dummy

	# ensure index arrays are integer type
	inj_indx = inj_indx.astype('int')
	corei = np.array((corei)).astype('int')
	
	return (y, H2Oi, y_mw, num_comp, Cfactor, y_indx_plot, corei, dydt_vst, 
				comp_namelist, inj_indx, core_diss,
				Psat_water, nuci, nrec_steps, seedi)
