########################################################################
#								       #
# Copyright (C) 2018-2024					       #
# Simon O'Meara : simon.omeara@manchester.ac.uk			       #
#								       #
# All Rights Reserved.                                                 #
# This file is part of PyCHAM                                          #
#                                                                      #
# PyCHAM is free software: you can redistribute it and/or modify it    #
# under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or    #
# (at  your option) any later version.                                 #
#                                                                      #
# PyCHAM is distributed in the hope that it will be useful, but        #
# WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# General Public License for more details.                             #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with PyCHAM.  If not, see <http://www.gnu.org/licenses/>.      #
#                                                                      #
########################################################################
'''function to initiate concentrations of components'''
# based on inputs, initial concentrations and their holding arrays 
# are set

import numpy as np
import scipy.constants as si
import math
from water_calc import water_calc
import write_dydt_rec
import openbabel.pybel as pybel


def init_conc(num_comp, init_conc, PInit,
	testf, num_eqn, Compt, seed_mw,
	core_diss, self):
		
	# inputs:------------------------------------------------------
	
	# num_comp - number of unique components
	# self.comp0 - chemical scheme names of components present at 
	# 	start of experiment
	# init_conc - initial concentrations of components (ppb)	
	# self.TEMP[0] - temperature in chamber at start of 
	# 	experiment (K)
	# self.RH - relative humidity in chamber (dimensionless 
	#	fraction 0-1)
	# PInit - initial pressure (Pa)
	# init_SMIL - SMILES of components present at start of 
	#	experiment (whose concentrations are given in init_conc)
	# testf - flag for whether in normal mode (0) or testing 
	# mode (1/2)
	# self.dydt_trak - chemical scheme name of components for 
	#	which user wants the tendency to  change tracked
	# self.tot_time - total simulation time (s)
	# self.save_step - recording frequency (s)
	# self.rindx_g - indices of reactants per equation
	# self.pindx_g - indices of products per equation
	# num_eqn - number of equations
	# self.comp_namelist - list of names of components as 
	#	presented in the chemical scheme file
	# Compt - name of components injected instantaneously after 
	#	start of experiment
	# self.seed_name - name of core component (input by user)
	# seed_mw - molecular weight of seed material (g/mol)
	# core_diss - dissociation constant of seed material
	# self.nuc_comp - name of nucleating component (input by user, 
	#	or defaults to 'core')
	# self.comp_xmlname - component names in xml file
	# self.comp_smil - all SMILES strings in xml file
	# self.self.rel_SMILES - only the SMILES strings of 
	# components present in the chemical scheme file
	# self.RO_indx - RO chemical scheme indices of alkoxy radicals
	# self.RO2_indices - RO2 list indices and chemical 
	# scheme indices of non-HOM-RO2 molecules
	# self.HOMRO2_indx - chemical scheme indices of HOM-RO2 molecules
	# self.rstoi_g - stoichiometry of reactants per equation
	# self.pstoi_g - stoichiometry of products per equation
	# self - reference to program
	# -----------------------------------------------------------
	
	# start by assuming no error
	erf = 0
	err_mess = ''
	
	if (testf == 1): # testing mode
		# return dummies
		return(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

	NA = si.Avogadro # Avogadro's number (# molecules/mol)
	# empty array for storing component gas-phase concentrations, must be an array
	y = np.zeros((num_comp))
	# empty array for storing component surface-phase concentrations, must be an array,
	# note that 2 components added on to prepare for water and core and that if water
	# included in chemical scheme then we adjust for that below
	y_w = np.zeros(((num_comp+2)*self.wall_on))
	y_mw = np.zeros((num_comp, 1)) # empty array for component molar mass (g/mol)
	# empty array for storing index of interesting gas-phase components
	y_indx_plot = []
	
	# convert concentrations
	# total number of molecules in 1 cm3 air using ideal gas law.  
	# R has units m3.Pa/K.mol (changing to cm3.Pa/K.mol when 
	# multiplied by 1e6)
	ntot = PInit*(NA/((si.R*1.e6)*self.TEMP[0]))
	# one billionth of number of # molecules in chamber unit volume
	Cfactor = ntot*1.e-9 # ppb to # molecules/cm3 conversion factor
	self.Cfactor = Cfactor
	Cfac_flag = 1 # flag to use Cfactor on initial concentrations
	# prepare dictionary for tracking tendency to change of user-specified components
	self.dydt_vst = {}

	# if no initial concentrations given in model variables file,
	# then check if the observed concentrations file has 
	# concentrations at the starting time
	if hasattr(self, 'obs'):
		if (len(self.comp0) == 0 and sum(self.obs[:, 0] == 0)>0):
			self.comp0 = (np.array((self.comp_namelist))[self.obs_comp_i]).tolist()
			# molecules/cm3
			init_conc = np.squeeze(self.obs[self.obs[:, 0] == 0, 1::])
			# flag to not use Cfactor on initial concentrations
			Cfac_flag = 0

	

	# insert initial concentrations where appropriate
	for i in range(len(self.comp0)):

		if ('_wall' in self.comp0[i]): # in case component concentration relates to wall
			wall_flag = 1
			# get index of string where _wall stated
			str_cnt = 0 # count way through string
			for ii in self.comp0[i]:
				
				if ii == '_':
					if self.comp0[i][str_cnt:str_cnt+5] == '_wall':
						# component name
						comp_now = self.comp0[i][0:str_cnt]
						# wall numbers provided by user start from 1
						wall_number = int(self.comp0[i][str_cnt+5])-1
					break
				str_cnt += 1 # count way through string	
			
		else: # in case component concentration does not relate to wall
			wall_flag = 0

		if (wall_flag == 0):
    			# index of where initial components occur in list of components
			# in case components already listed via 
			# interpretation of the chemical scheme
			try:
				y_indx = self.comp_namelist.index(self.comp0[i])
			
			# if component not already listed via 
			# interpretation of the chemical scheme
			# then send error message
			except:
				erf = 1
				err_mess = str('Error: component called ' + str(self.comp0[i]) + ', which has an initial concentration specified in the model variables input file has not been found in the chemical scheme.  Please check the scheme and associated chemical scheme markers, which are stated in the model variables input file.')
				return (0, 0, 0, 0, 0, 0, 0, 
					0, 0, 0,
					0, 0, erf, err_mess, 0, 0, 0, 0)
			
			# set initial concentration
			if (Cfac_flag == 1):
				# convert from ppb to # molecules/cm3
				y[y_indx] = init_conc[i]*Cfactor
			if (Cfac_flag == 0):
				# # molecules/cm3
				y[y_indx] = init_conc[i]

			# remember index for plotting gas-phase concentrations later
			y_indx_plot.append(y_indx)
	
		if (wall_flag == 1):

			# index of where initial components occur in list of components
			# in case components already listed via 
			# interpretation of the chemical scheme
			try:
				y_indx = self.comp_namelist.index(comp_now)
				
			# if component not already listed via interpretation 
			# of the chemical scheme
			# then send error message
			except:
				erf = 1
				err_mess = str('Error: component called ' + str(comp_now) + ', which has an initial concentration specified in the model variables input file has not been found in the chemical scheme.  Please check the scheme and associated chemical scheme markers, which are stated in the model variables input file.')
				return (0, 0, 0, 0, 0, 0, 0, 
					0, 0, 0,
					0, 0, erf, err_mess, 0, 0, 0, 0, 0)

			# insert surface concentration into surface array 
			# convert from ppb to # molecules/cm3
			y_w[(num_comp+2)*wall_number+y_indx] = init_conc[i]*Cfactor
			y_indx_plot.append(y_indx) # remember for plotting
		
	# check on whether O3 isopleth due to be made
	# if isopleth due to be made overide any 
	# originally provided initial conditions
	if (self.testf == 5):
		y[self.VOCi] = self.VOCequil # VOC concentration
		y[self.NOi] = self.NOxequil/2. # NO concentration
		y[self.NO2i] = self.NOxequil/2. # NO2 concentration
		# if a previous equilibrium ozone concentration has been found, 
		# then use this as a first guess to help speed-up iteration
		#if hasattr(self, 'O3equil'):
		#	y[self.O3i] = self.O3equil
		
		#else: # if no existing guess of [O3], then estimate of [NOx] and [VOC]
		y[self.O3i] = self.VOCequil+self.NOxequil

		self.O3equil = y[self.O3i] # register first guess at O3
		# note that for NOx, the NO:NO2 ratio is allowed to change during 
		# integration steps,
		# however, the total NO+NO2 concentration is held constant through commands
		# contained in ode_updater, likewise inside ode solver the VOC concentration is 
		# not allowed to change

	# number of recording steps
	nrec_steps = int(math.ceil(self.tot_time/self.save_step)+1)
	
	for i in range(num_comp): # loop through all components to get molar weights
		y_mw[i] = self.Pybel_objects[i].molwt # molecular mass (g/mol)
	
	# --------------------------------------------------------------
	# account for water's properties
	
	# get initial gas-phase concentration (# molecules/cm3 (air)) 
	# and vapour pressure
	# of water (log10(atm))
	[C_H2O, Psat_water, H2O_mw] = water_calc(self.TEMP[0], 
		self.RH[0], si.N_A)
	
	# if skipping parsing and property estimation
	if (self.pars_skip == 1):
		H2Oi = self.H2Oi # index for water
		y_mw = self.y_mw
		num_comp = self.num_comp
		if (y.shape[0] == self.y.shape[0]-2):
			y = np.append(y, np.zeros((2)))
			y[self.H2Oi] = self.y[self.H2Oi]
			y[self.seedi] = self.y[self.seedi]

	if (self.pars_skip == 0 or self.pars_skip == 2):
		# holder for water index (will be used if not identified in chemical scheme)
		H2Oi = num_comp # index for water
		self.H2Oi = H2Oi # index for water
	
		# check for water presence in chemical scheme via its SMILE string
		# count on components
		indx = -1
		for single_chem in self.rel_SMILES:
			indx += 1
			# ensure this is water rather than single 
			# oxygen (e.g. due to ozone photolysis 
			# (O is the MCM chemical scheme name for single oxygen))
			
			if (single_chem == 'O' and self.comp_namelist[indx] != 'O'):
				
				H2Oi = indx
				self.H2Oi = H2Oi # index for water
				# include initial concentration of water
				# (# molecules/cm3)
				y[H2Oi] = C_H2O
				
				# include molar weight of water (g/mol)
				y_mw[H2Oi] = H2O_mw
			
				# remove the addition of water in the 
				# surface concentrations
				# first rearrange matrix so that 
				# components in rows, surface number in 
				#columns
				y_w = y_w.reshape(self.wall_on, 
					num_comp+2)
				# remove the excess water column
				y_w = np.concatenate((y_w[:, 0:-2], 
				y_w[:, -1].reshape(-1, 1)), axis=1)
				# then flatten back to 1D array
				y_w = y_w.flatten()

		
		# if not included in chemical scheme file, then add 
		# water to end of component list
		if (H2Oi == num_comp):
	
			# update number of components to account for 
			# water
			num_comp += 1
			# append empty element to y and y_mw to hold 
			# water values
			y = np.append(y, C_H2O)
			# append molar weight of water (g/mol)
			y_mw = (np.append(y_mw, H2O_mw)).reshape(-1, 1)
			# append water's name to component name list
			self.comp_namelist.append('H2O')
			# add to SMILES list
			H2O_SMILES = 'HOH'
			self.rel_SMILES.append(H2O_SMILES)
			# generate pybel object
			Pybel_object = pybel.readstring('smi', 'O')
				
			# append to Pybel object list
			self.Pybel_objects.append(Pybel_object)

	# ------------------------------------------------------------------------------------
	# seed components

	# if this information already gained in previous run then skip
	if (self.pars_skip == 0 or self.pars_skip == 2): 

		# empty array for index of core component
		self.seedi = (np.zeros((len(self.seed_name)))).astype(int)

		# append any seed names that are not seen in the 
		# chemical scheme file to the list of components
		for seed_cnt in range(len(self.seed_name)):
			if self.seed_name[seed_cnt] == 'core':
				self.seedi[seed_cnt] = num_comp
			else:
				self.seedi[seed_cnt] = self.comp_namelist.index(self.seed_name[seed_cnt])

		# append name of core to component name list
		self.comp_namelist.append('core')
		
		# increase number of components to account for 'core' component
		num_comp += 1
		# prepare for skipping of parsing in following simulations
		self.num_comp = num_comp
		# add to SMILES list
		self.rel_SMILES.append('[NH4+].[NH4+].[O-]S(=O)(=O)[O-]')

		# append core gas-phase concentration (molecules/cm3 (air)) and molar 
		# mass (g/mol) (needs to have a 1 length in second dimension for the kimt 
		# calculations)
		y = np.append(y, 0.)
		y_mw = (np.append(y_mw, seed_mw)).reshape(-1, 1)
		# prepare for skipping of parsing in future simulations
		self.y_mw = y_mw
		self.y = y # prepare for skipping of parsing in following simulations 
	
	# account for seed properties - note that even if no seed particle, this code ensures
	# that an index is provided for core material
	corei = [num_comp-1] # index for core component
	
	# finally append surface concentration array to gas-phase array
	y = np.append(y, y_w)
	
	# tracked components (dydt)--------------------------------------
	# check for tracking of all alkyl peroxy radicals
	if ('RO2_ind' in self.dydt_trak):
		# append existing list and list of RO2 radicals 
		self.dydt_trak = self.dydt_trak + ((np.array((
			self.comp_namelist)))[self.RO2_indices[:, 
			1]]).tolist()
	
		# remove the RO2_ind item from the tracked component 
		# list
		self.dydt_trak.remove('RO2_ind')

	# check for tracking of all alkoxy radicals
	if ('RO_ind' in self.dydt_trak):
		# append existing list and list of RO radicals 
		self.dydt_trak = self.dydt_trak + ((np.array((
		self.comp_namelist)))[self.RO_indx]).tolist()
	
		# remove the RO_ind item from the tracked component list
		self.dydt_trak.remove('RO_ind')

	# get index of user-specified components for tracking their 
	# change tendencies (dydt) due to modelled mechanisms
	if (len(self.dydt_trak) > 0):
		
		# empty list for indices of these components
		dydt_traki = [] 	

		for i in range (len(self.dydt_trak)):

			# indices of reactions involving this component
			reac_index = [] 
			# value for whether this component is reactant 
			# or product in a reaction and stoichiometry
			reac_sign = []
			
			if (self.dydt_trak[i] != 'RO2' and 
				self.dydt_trak[i] != 'HOMRO2'):
			
				# index of components in component list
				try:
					y_indx = self.comp_namelist.\
					index(self.dydt_trak[i])
				# if component not already listed via 
				# interpretation of the chemical scheme
				# then send error message
				except:
					erf = 1
					err_mess = str('''Error: 
					component called ''' + 
					str(self.dydt_trak[i]
					) + ''', which is specified to 
					be tracked in the model 
					variables input file has not 
					been found in the chemical 
					scheme.  Please check the scheme
					 and associated chemical scheme 
					markers, which are stated in 
					the model variables input 
					file.''')
					return (0, 0, 0, 0, 0, 0, 0, 
					0, 0, 0,
					0, 0, erf, err_mess, 0, 0, 0, 0)
				# remember index for plotting gas-phase
				# concentrations later
				dydt_traki.append([int(y_indx)])
			
				# search through reactions to see 
				# where this component is reactant or 
				# product
				for ri in range(num_eqn):
					if sum(self.rindx_g[ri, 
						0:self.nreac_g[ri]] == 
						y_indx) > 0:
						# append reaction index
						reac_index.append(int(
						ri)) 
						reac_place = np.where(self.rindx_g[ri, 0:self.nreac_g[ri]] == y_indx)[0]
						reac_sign.append(-1*self.rstoi_g[int(ri), reac_place])
					if sum(self.pindx_g[ri, 0:self.nprod_g[ri]] == y_indx) > 0:
						reac_index.append(int(ri)) # append reaction index
						reac_place = np.where(self.pindx_g[ri, 0:self.nprod_g[ri]] == y_indx)[0]
						reac_sign.append(1*self.pstoi_g[int(ri), reac_place])
			
			if (self.dydt_trak[i] == 'RO2'): # for non-HOM-RO2 (alkyl peroxy radicals)
			
				# remember index for plotting gas-phase concentrations later
				dydt_traki.append(list(self.RO2_indices[:, 1]))
				
				# loop through non-HOM-RO2 components to get their reaction indices
				for y_indx in self.RO2_indices[:, 1]:
					
					# search through reactions to see where this component is reactant or product
					for ri in range(num_eqn):
						if sum(self.rindx_g[ri, 0:self.nreac_g[ri]] == y_indx) > 0:
							reac_index.append(int(ri)) # append reaction index
							reac_place = np.where(self.rindx_g[ri, 0:self.nreac_g[ri]] == y_indx)[0]
							reac_sign.append(-1*self.rstoi_g[int(ri), reac_place])
							
						if sum(self.pindx_g[ri, 0:self.nprod_g[ri]] == y_indx) > 0:
							reac_index.append(int(ri)) # append reaction index
							reac_place = np.where(self.pindx_g[ri, 0:self.nprod_g[ri]] == y_indx)[0]
							reac_sign.append(1*self.pstoi_g[int(ri), reac_place])
					
					y_indx = self.RO2_indices[:, 1] # ready for storing below
			
			if (self.dydt_trak[i] == 'HOMRO2'): # for HOM-RO2 (highly oxidised molecule alkyl peroxy radicals)
			
				# remember index for plotting gas-phase concentrations later
				dydt_traki.append(list(np.squeeze(self.HOMRO2_indx)))
				
				# loop through non-HOM-RO2 components to get their reaction indices
				for y_indx in self.HOMRO2_indx[:]:
					
					# search through reactions to see where this component is reactant or product
					for ri in range(num_eqn):
						if sum(self.rindx_g[ri, 0:self.nreac_g[ri]] == y_indx) > 0:
							reac_index.append(int(ri)) # append reaction index
							reac_place = np.where(self.rindx_g[ri, 0:self.nreac_g[ri]] == y_indx)[0]
							reac_sign.append(-1*self.rstoi_g[int(ri), reac_place])
						if (sum(self.pindx_g[ri, 0:self.nprod_g[ri]] 
							== y_indx) > 0):
							# append reaction index
							reac_index.append(int(ri))
							reac_place = np.where(self.pindx_g[
							ri, 0:self.nprod_g[ri]] == y_indx)[0]
							reac_sign.append(1*self.pstoi_g[
							int(ri), reac_place])
				y_indx = self.HOMRO2_indx[:] # ready for storing below
				
			# save reaction indices in dictionary value 
			# for this component,
			# when creating empty rec_array, add three 
			# columns onto the end for 
			# gas-particle-partitioning, gas-wall-partitioning and 
			# dilution, and include a top row for reaction indices
			rec_array = np.zeros((nrec_steps+1, 
				len(reac_index)+3))
			rec_array[0, 0:-3] = reac_index

			comp_indx_str = str(self.dydt_trak[i] + 
				'_comp_indx')
			res_string = str(self.dydt_trak[i] + '_res')
			reac_string = str(self.dydt_trak[i] + 
				'_reac_sign')
	
			if hasattr(self, 'sim_ci_file'):
				ci_string = str(self.dydt_trak[i] + 
					'_ci')
				ci_arrayi = np.zeros((nrec_steps-1, 1))
			# dictionary entry to hold index of tracked 
			# component	
			self.dydt_vst[comp_indx_str] = y_indx
			# dictionary entry to hold reaction indices 
			# and results 
			self.dydt_vst[res_string] = rec_array
			# dictionary entry to hold sign (source or sink
			# per reaction) 
			self.dydt_vst[reac_string] = reac_sign 
			if hasattr(self, 'sim_ci_file'):
				# dictionary entry to hold continuous 
				#influx rate
				self.dydt_vst[ci_string] = ci_arrayi

		# dictionary entry to hold component names of 
		# components to track
		self.dydt_vst['comp_names'] = self.dydt_trak
		
		# call on write_dydt_rec to generate the module that 
		# will process
		# the tendency to change during the simulation
		write_dydt_rec.write_dydt_rec(self)
	
	# --------------------------------------
	
	# if nucleating component formed of core component
	if (self.nuc_comp[0] == 'core'):
		nuci = num_comp-1 # index of core component
	else:
		nuci = -1 # filler
	
	# get indices of seed particle component(s)
	indx = 0 # count on seed component(s)
	for sname in self.seed_name:
		# index of core component
		self.seedi[indx] = int(self.comp_namelist.index(sname))
		indx += 1 # count on seed component(s)
	
	# get index of component with latter injections
	if len(Compt)>0:
		inj_indx = np.zeros((len(Compt)))
		for i in range(len(Compt)):
			# index of where instantaneously injected components 
			# occur in SMILES string
			inj_indx[i] = self.comp_namelist.index(Compt[i])
	else:
		inj_indx = np.zeros((1)) # dummy

	# ensure index arrays are integer type
	inj_indx = inj_indx.astype('int')
	
	corei = np.array((corei)).astype('int')
	
	# get indices of NO, HO2 and NO3 (for reaction rate calculations)
	
	try:
		NOi = self.comp_namelist.index(NO)
	except:
		NOi = 0 # filler
	try:
		HO2i = self.comp_namelist.index(HO2)
	except:
		HO2i = 0 # filler
	try:
		NO3i = self.comp_namelist.index(NO3)
	except:
		NO3i = 0 # filler

	# if user wants to see molar masses of all components
	if (self.testf == 3.1):

		import matplotlib.pyplot as plt
		from matplotlib.colors import BoundaryNorm
		from matplotlib.ticker import MaxNLocator
		# for customised colormap
		from matplotlib.colors import LinearSegmentedColormap
		# set colormap tick labels to standard notation
		import matplotlib.ticker as ticker 

		plt.ion() # show results to screen and turn on interactive mode
		
		# prepare plot
		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
		# get indices of molar mass in ascending order
		asc_ind = np.argsort(y_mw, axis = 0)
		array_names = (np.array(self.comp_namelist)).reshape(-1, 1)

		# plot molar masses against component names in ascending order
		ax0.plot(np.arange(len(y_mw)), y_mw[asc_ind][:, 0, 0], '+')

		ax0.set_ylabel(r'Molar Mass (g mol$\rm{^{-1}}$)', fontsize = 14)
		ax0.set_xlabel(r'Component name', fontsize = 14)
		# set location of x ticks
		ax0.set_xticks(np.arange(len(self.comp_namelist)))
		ax0.set_xticklabels(array_names[asc_ind], rotation = 45)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
		ax0.set_title(str('Molar masses of all components'), fontsize = 14)
		
		# ensure that software stops on return to middle module
		erf = 1
		err_mess = 'Stop'

	try: # in case called from autorun, in which a finisher simulation is setup
		if (self.param_const['sim_type'] == 'finisher'):
			# gas-phase concentrations (# molecules/cm3)
			y[0:num_comp] = self.param_const['ynow'][0:num_comp]
			
	except: # not called from finisher simulation
		y[:] = y[:]

	return (y, H2Oi, y_mw, num_comp, Cfactor, y_indx_plot, corei, 
			inj_indx, core_diss,
			Psat_water, nuci, nrec_steps, erf, err_mess, NOi, 
			HO2i, NO3i, init_conc, self)
