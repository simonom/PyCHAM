########################################################################
#                                                                      #
# Copyright (C) 2018-2026                                              #
# Simon O'Meara : simon.omeara@manchester.ac.uk                        #
#                                                                      #
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
'''module to save user-defined PyCHAM results to files'''
# module called following simulation in PyCHAM, to store only 
# user-defined information
# for tips on representation of units: 
# https://ncics.org/portfolio/other-resources/udunits2/
# and standard names for CF convention are here:
# https://cfconventions.org/Data/cf-standard-names/
# current/build/cf-standard-name-table.html

# import dependencies
import scipy.constants as si
import numpy as np
import time

def bespoke_saving(savei, y_mat, y_MM, t_out, cham_env, Cfactor_vst, 
				   numsb, Nresult_wet, Nresult_dry, rbou_rec, H2Oi,
				   MV, indx_plot, x, rootgrp, self):
	
	# inputs: -------------------------------------------------------
	# savei - output index to save
	# y_mat - component (columns) concentrations with time (rows) 
	#	(# molecules/cm^3 (air))
	# y_MM - molar masses of components (g/mol)
	# t_out - the times (s) through simulation that outputs
	# 	correspond to
	# cham_env - chamber environmental conditions (temperature (K), 
	# pressure (Pa), relative humidity (0-1), transmission factor of 
	#	light (0-1)
	# rootgrp - the netCDF file object
	# self - the PyCHAM object
	# ---------------------------------------------------------------

	# ensure dimensions with particle number concentration per size bin
	# are a useful shape
	if ((numsb-self.wall_on) == 0): # if zero size bins
		# ensure numpy arrays, rather than float
		Nresult_wet = (np.array((Nresult_wet))).reshape(-1, 1)
		Nresult_dry = (np.array((Nresult_dry))).reshape(-1, 1)
		rbou_rec = (np.array((rbou_rec))).reshape(-1, 1)
		x = (np.array((x))).reshape(-1, 1)
		
	if ((numsb-self.wall_on) == 1): # if one size bin
		# ensure numpy arrays, rather than float
		Nresult_wet = (np.array((Nresult_wet)))
		Nresult_dry = (np.array((Nresult_dry)))
		rbou_rec = (np.array((rbou_rec)))	
		x = (np.array((x)))	

	# if first index, then save the unconditional variables
	if (savei == 0):

		# unconditional saving of: general attributes, time, 
		# physical conditions 
		# (temperature, pressure, RH, transmission factor of light)

		# store general attributes
		# get time now
		now_time = time.time()
		rootgrp.history = str('Saved at (Greenwhich Meridian Time) ' + str(time.gmtime(now_time)))
		rootgrp.source = str('CHemistry with Aerosol Microphysics in Python (PyCHAM), version ' + self.vnum)

		# for CF compliance, time should be a numeric value with
		# a time unit based on a reference starting time. Note that
		# the number immediately after the year is the month, followed
		# by the number for the day, as explained in:
		#https://cfconventions.org/Data/cf-conventions/
		#cf-conventions-1.7/build/ch04s04.html#:~:text=The%
		#20acceptable%20units%20for%20time,Plural%20forms%
		#20are%20also%20acceptable.
		# set the unit for time, note that using a real world 
		# reference time is somewhat arbitrary for simulations
		# that don't necessarily relate to any real world time, therefore,
		# by default the reference time is set to the real world time at
		# which the computer simulation began - time stamped in 
		# ode_updater module
		time_unit = str('seconds since ' + str(self.sim_start_timestamp[0]) + '-' +
		str(self.sim_start_timestamp[1]) + '-' + str(self.sim_start_timestamp[2]) + 
		' ' + str(self.sim_start_timestamp[3]) + ':' + 
		str(self.sim_start_timestamp[4]) + ':' + str(self.sim_start_timestamp[5]) + 
		' +00:00')

		# set dimensions that netCDF variables may use ---------------
		
		# set the time dimension
		rootgrp.createDimension('time', len(t_out))

		# set the particle size bin number dimension
		rootgrp.createDimension('particle_size_bins', self.nasb)

		# set the particle size bin bounds dimension
		rootgrp.createDimension('particle_size_bin_bounds', self.nasb+1)

		# set the surface number dimension
		rootgrp.createDimension('surface_number', self.wall_on)

		# set the component dimension
		rootgrp.createDimension('components', len(self.comp_namelist))

		# set the reaction dimension
		rootgrp.createDimension('reactions', 
						  len(self.reac_coef_g)+len(self.reac_coef_aq)+len(self.reac_coef_su))
		
		# set the constant with time value dimension (i.e. for single numbers)
		rootgrp.createDimension('single_values', 1)

		# set dimension for seed components
		rootgrp.createDimension('seed_components', len(self.seedi[:]))

		# set dimension for organic peroxy radical components
		rootgrp.createDimension('organic_peroxy_radical_components', len(self.RO2_indices[:, -1]))

		# set dimension for component with initial concentration specified
		rootgrp.createDimension('initial_components', len(self.comp0))

		# set dimension for maximum number of reactants in reactions
		rootgrp.createDimension('maximum_reactants', max(self.nreac_g))

		# set dimension for maximum number of products in reactions
		rootgrp.createDimension('maximum_products', max(self.nprod_g))
	
		# end of setting dimensions ---------------------------------

		# create time variable
		timevar = rootgrp.createVariable('time', 'float32', ('time'))
		# set units for time variable
		timevar.setncattr('units', time_unit)
		# set values for time variable (s)
		timevar[:] = t_out

		# create temperature variable
		tempervar = rootgrp.createVariable('temperature', 'float32', ('time'))
		# set units for temperature variable
		tempervar.setncattr('units', 'K')
		# set values for temperature variable (K)
		tempervar[:] = cham_env[:, 0]

		# create pressure variable
		pressvar = rootgrp.createVariable('air_pressure', 'float32', ('time'))
		# set units for pressure variable
		pressvar.setncattr('units', 'pascal')
		# set values for pressure variable (Pa)
		pressvar[:] = cham_env[:, 1]

		# create RH variable
		rhvar = rootgrp.createVariable('relative_humidity', 'float32', ('time'))
		# set units for pressure variable
		rhvar.setncattr('units', 'fraction 0-1')
		# set values for relative humidity variable (0-1)
		rhvar[:] = cham_env[:, 2]

		# create light transmission factor variable
		tfvar = rootgrp.createVariable('transmission factor of light', 
			'float32', ('time'))
		# set units for transmission factor variable
		tfvar.setncattr('units', 'fraction 0-1')
		# set values for tranmission factor of light variable (0-1)
		tfvar[:] = cham_env[:, 3]

		# create component chemical scheme name variable
		cnamesvar = rootgrp.createVariable('component_chemical_scheme_name', 
			'str', ('components'))
		# set units for component chemical scheme name variable
		cnamesvar.setncattr('units', 'chemical scheme name')
		# set values for component chemical scheme name variable
		# note that because the names have variable numbers of 
		# characters, they must be assigned to the netCDF variable 
		# individually
		for stri in range(len(self.comp_namelist)):
			cnamesvar[stri] = self.comp_namelist[stri]

		# create variable for number of components
		ncvar = rootgrp.createVariable('number_of_components', 
			'int', ('single_values'))
		# set units
		ncvar.setncattr('units', 'total number of components')
		# set values
		ncvar[0] = self.nc

		# create wall on flag variable
		wallonvar = rootgrp.createVariable('wall_on_flag', 
			'float32', ('single_values'))
		# set units
		wallonvar.setncattr('units', 'wall_on_flag_0forNO_>0forYES')
		# set values
		wallonvar[0] = self.wall_on

		# create number of particle size bin variable
		nsb_var = rootgrp.createVariable('number_of_particle_size_bins', 
			'float32', ('single_values'))
		# set units
		nsb_var.setncattr('units', 'number of particle size bins')
		# set values
		nsb_var[0] = numsb

		# number concenration of anhydrous particles variable
		Ndry_var = rootgrp.createVariable('particle_number_concentration_anhydrous', 
			'float32', ('time', 'particle_size_bins'))
		# set units
		Ndry_var.setncattr('units', '# particles/cm^3 (air)')
		# set values
		Ndry_var[:, :] = Nresult_dry[:, :]

		# number concenration of hydrous particles variable
		Nwet_var = rootgrp.createVariable('particle_number_concentration_hydrous', 
			'float32', ('time', 'particle_size_bins'))
		# set units
		Nwet_var.setncattr('units', '# particles/cm^3 (air)')
		# set values
		Nwet_var[:, :] = Nresult_wet[:, :]

		# radius bounds of particles variable (um)
		rbou_var = rootgrp.createVariable('particle_radius_bounds', 
			'float32', ('time', 'particle_size_bin_bounds'))
		# set units
		rbou_var.setncattr('units', 'um')
		# set values
		rbou_var[:, :] = rbou_rec[:, :]

		# radius of hydrous particles variable (um)
		r_var = rootgrp.createVariable('hydrous_particle_radius', 
			'float32', ('time', 'particle_size_bins'))
		# set units
		r_var.setncattr('units', 'um')
		# set values
		r_var[:, :] = x[:, :]

		# spacing between particle size bins
		space_mode_var = rootgrp.createVariable('particle_size_bin_spacing_mode', 
			'str', ('single_values'))
		# set units
		space_mode_var.setncattr('units', 'string')
		# set values
		space_mode_var[0] = self.space_mode


		# create SMILES variable for components in chemical scheme
		SMILES_var = rootgrp.createVariable('SMILES_strings', 
			'str', ('components'))
		# set units
		SMILES_var.setncattr('units', 'SMILES strings of components in chemical scheme')
		# set values
		# note that because the strings have variable numbers of 
		# characters, they must be assigned to the netCDF variable 
		# individually
		for stri in range(len(self.rel_SMILES)):
			SMILES_var[stri] = self.rel_SMILES[stri]

		# create molar mass variable for components in chemical scheme
		mm_var = rootgrp.createVariable('molar_masses', 
			'float32', ('components'))
		# set units
		mm_var.setncattr('units', 'g/mol')
		# set values
		mm_var[:] = np.squeeze(y_MM)[:]

		# create molar volume variable for components in chemical scheme
		mm_var = rootgrp.createVariable('molar_volumes', 
			'float32', ('components'))
		# set units
		mm_var.setncattr('units', 'cm^3/mol')
		# set values
		mm_var[:] = np.squeeze(MV)[:]

		# create variable for pure component vapour pressure at 298.15 K
		# for components in the chemical scheme
		PsatPa_var = rootgrp.createVariable('pure_component_saturation_vapour_pressure_of_components_at_298.15_K', 
			'float32', ('components'))
		# set units
		PsatPa_var.setncattr('units', 'Pa')
		# set values
		PsatPa_var[:] = self.Psat_Pa_rec[:]

		# create variable for oxygen to carbon ratios
		# for components in the chemical scheme
		OC_var = rootgrp.createVariable('oxygen_to_carbon_ratios', 
			'float32', ('components'))
		# set units
		OC_var.setncattr('units', 'fraction')
		# set values
		OC_var[:] = self.OC[:]

		# create seed index variable
		h2oi_var = rootgrp.createVariable('water_component_index', 
			'float32', ('single_values'))
		# set units
		h2oi_var.setncattr('units', 'index of water')
		# set values
		h2oi_var[0] = H2Oi

		# create seed index variable
		seedi_var = rootgrp.createVariable('seed_component_indices', 
			'float32', ('seed_components'))
		# set units
		seedi_var.setncattr('units', 'indices of seed')
		# set values
		seedi_var[:] = self.seedi[:]

		# create organic peroxy radical index variable
		ro2i_var = rootgrp.createVariable('organic_peroxy_radical_component_indices', 
			'float32', ('organic_peroxy_radical_components'))
		# set units
		ro2i_var.setncattr('units', 'indices of organic peroxy radicals')
		# set values
		ro2i_var[:] = self.RO2_indices[:, -1]

		# chemical scheme names of components with initial gas-phase concentration variable
		comp0_var = rootgrp.createVariable('names_of_components_with_initial_gas_phase_concentrations_specified', 
			'str', ('initial_components'))
		# set units
		comp0_var.setncattr('units', 'chemical scheme names of components with initial gas-phase concentrations specified')
		# set values
		# note that because the strings have variable numbers of 
		# characters, they must be assigned to the netCDF variable 
		# individually
		for stri in range(len(self.comp0)):
			comp0_var[stri] = self.comp0[stri]
		
		# indices of components with initial gas-phase concentration variable
		indx0_var = rootgrp.createVariable('indices_of_components_with_initial_gas_phase_concentrations_specified', 
			'int', ('initial_components'))
		
		# set units
		indx0_var.setncattr('units', 'indices of components with initial gas-phase concentrations specified')
		# set values
		indx0_var[:] = indx_plot[:]

		return(rootgrp)

	else:
		# get name of user-defined variable to save
		ud_var = str(self.user_output[savei]).strip()

	# move onto conditional saving ---------------------------------

	# the gas-phase concentration(s) (molecules/cm^3)	
	if ('_g' in ud_var):
		
		# save the conversion factor (with time) for converting molecules/cm^3 to ppb
		# create variable inside the concentrations group
		cfac_var = rootgrp.createVariable(str('factor_for_multiplying_ppb_to_get_molec/cm3_with_time'), 
			'float32', ('time'))
		
		# set values for variable
		cfac_var[:] = Cfactor_vst

		# set units
		cfac_var.setncattr('units', 'molecules/cm^3/ppb')

		# get the component to be stored
		comp_name = ud_var[0:ud_var.index('_g')]

		# if all components to be saved
		if (comp_name == 'all'):
			cindx_all = np.arange(self.nc).astype('int')
		
		else: # if individual component to be saved
			# get index of component
			cindx_all = np.array((self.comp_namelist.index(comp_name))).reshape(1)

		# loop through components
		for ci in cindx_all:
			
			# set variable name
			var_name = str('number_concentration_of_' +
				'gas_phase_' + self.comp_namelist[ci] + '_molecules_in_air')

			# if component already done, then continue to next component
			if var_name in rootgrp['concentrations_g'].variables.keys():
				continue

			# create variable inside the concentrations group
			c_gvar = rootgrp.createVariable(str('/concentrations_g/' + var_name), 
				'float32', ('time'))
			# get gas-phase concentration with time (molecules/cm^3)
			c_g = y_mat[:, ci]
			# set values for variable
			c_gvar[:] = c_g

			# set units of concentration
			c_gvar.setncattr('units', 'molecules/cm^3')

	# the particle-phase concentration(s) (molecules/cm^3)	
	if ('_p' in ud_var):
		
		# get the component to be stored
		comp_name = ud_var[0:ud_var.index('_p')]

		# if all components to be saved
		if (comp_name == 'all'):
			cindx_all = np.arange(self.nc).astype('int')
		
		else: # if individual component to be saved
			# get index of component
			cindx_all = np.array((self.comp_namelist.index(comp_name))).reshape(1)

		# loop through components
		for ci in cindx_all:
			
			# set variable name
			var_name = str('number_concentration_of_' +
				'particle_phase_' + self.comp_namelist[ci] + 
				'_molecules_in_air')

			# if component already done, then continue to next component
			if var_name in rootgrp['concentrations_p'].variables.keys():
				continue

			# create variable inside the concentrations group
			c_pvar = rootgrp.createVariable(str('/concentrations_p/' + var_name), 
				'float32', ('time', 'particle_size_bins'))
			
			# prepare to hold paricle-phase concentrations,
			# with times in rows and particle size bins in columns
			c_p = np.zeros((y_mat.shape[0], self.nasb))

			# get particle-phase concentrations with time across 
			# particle size bins (molecules/cm3)
			for sbi in range(0, self.nasb):
				c_p[:, sbi] = y_mat[:, (self.nc*(sbi+1))+ci]
			# set values for variable
			c_pvar[:, :] = c_p[:, :]

			# set units of concentration
			c_pvar.setncattr('units', 'molecules/cm^3')

	# the surface-phase concentration(s) (molecules/cm^3)	
	if ('_s' in ud_var):
		
		# get the component to be stored
		comp_name = ud_var[0:ud_var.index('_s')]

		# if all components to be saved
		if (comp_name == 'all'):
			cindx_all = np.arange(self.nc).astype('int')
		
		else: # if individual component to be saved
			# get index of component
			cindx_all = np.array((self.comp_namelist.index(comp_name))).reshape(1)

		# loop through components
		for ci in cindx_all:
			
			# set variable name
			var_name = str('number_concentration_of_' +
				'surface_phase_' + self.comp_namelist[ci] + 
				'_molecules_in_air')

			# if component already done, then continue to next component
			if var_name in rootgrp['concentrations_s'].variables.keys():
				continue

			# create variable inside the concentrations group
			c_svar = rootgrp.createVariable(str('/concentrations_s/' + var_name), 
				'float32', ('time', 'surface_bins'))
			
			# prepare to hold concentrations,
			# with times in rows and surface bins in columns
			c_s = np.zeros((y_mat.shape[0], self.wall_on))

			# get concentrations with time across 
			# surface bins (molecules/cm3)
			for sbi in range(0, self.wall_on):
				c_s[:, sbi] = y_mat[:, (self.nc*(self.nasb+1)+self.nc*sbi)+ci]
			# set values for variable
			c_svar[:, :] = c_s[:, :]

			# set units of concentration
			c_svar.setncattr('units', 'molecules/cm^3')

	# the mass concentration of secondary organic particulate matter (ug/m^3)	
	if (ud_var == 'SOA'):

		# isolate just particle-phase component concentrations
		# (molecules/cm^3)
		y_pp = np.zeros((y_mat.shape[0], self.nc*self.nasb))
		y_pp[:, :] = y_mat[:, self.nc:self.nc*(self.nasb+1)]

		# sum concentrations of individual components over size bins
		# loop over size bins beyond first
		for sbi in range(1, self.nasb):
			y_pp[:, 0:self.nc] += y_pp[:, self.nc*sbi:self.nc*(sbi+1)]

		# keep just the total concentration in particle phase
		# (molecules/cm^3)
		y_pp = y_pp[:, 0:self.nc]

		# zero the concentrations of components without a carbon 
		# (molecules/cm3), note this will zero water, and so
		# represent dry particulate matter
		y_pp[:, self.HC[0, :] == 0] = 0.

		# zero the seed components (molecules/cm^3)
		y_pp[:, self.seedi] = 0

		# tile molar masses over times
		y_MMtiled = np.tile(y_MM.reshape(1, -1), (y_pp.shape[0], 1))

		# convert concentrations of organics from molecules/cm^3 to ug/m^3
		y_pp[:, :] = ((y_pp[:, :]/si.N_A)*y_MMtiled)*1.e12

		# sum mass concentrations over components to get SOA concentration (ug/m^3)
		y_pp = np.sum(y_pp, axis= 1)

		# create SOA variable
		soavar = rootgrp.createVariable(str('mass_concentration_of_secondary_' + 
		'particulate_organic_matter_dry_aerosol_particles_in_air'), 
		'float32', ('time'))
		# set SOA units
		soavar.setncattr('units', 'ug/m^3')
		# set values for SOA variable
		soavar[:] = y_pp

	# the mass concentration of all particulate matter (ug/m^3)	
	if (ud_var == 'PM'):

		
		# isolate just particle-phase component concentrations
		# (molecules/cm^3)
		y_pp = np.zeros((y_mat.shape[0], self.nc*self.nasb))
		y_pp[:, :] = y_mat[:, self.nc:self.nc*(self.nasb+1)]

		# sum concentrations of individual components over size bins
		# loop over size bins beyond first
		for sbi in range(1, self.nasb):
			y_pp[:, 0:self.nc] += y_pp[self.nc*sbi:self.nc*(sbi+1)]

		# keep just the total concentration in particle phase
		# (molecules/cm^3)
		y_pp = y_pp[:, 0:self.nc]

		# tile molar masses over times
		y_MMtiled = np.tile(y_MM.reshape(1, -1), (y_pp.shape[0], 1))

		# convert concentrations from molecules/cm^3 to ug/m^3
		y_pp[:, :] = ((y_pp[:, :]/si.N_A)*y_MMtiled)*1.e12

		# sum mass concentrations over components to get total
		# particulate matter concentration (ug/m^3)
		y_pp = np.sum(y_pp, axis= 1)

		# create total particulate matter mass concentration variable
		pmvar = rootgrp.createVariable(str('mass_concentration_of_' + 
		'particulate_matter_wet_aerosol_particles_in_air'), 
		'float32', ('time'))
		# set units
		pmvar.setncattr('units', 'ug/m^3')
		# set values for variable
		pmvar[:] = y_pp

	# the oxygen number to carbon number ratio of organics with C>4 
	# in gas-phase
	if (ud_var == 'O:C'):

		# array to hold number concentration weighted O:C 
		# of organics in the gas-phase with C>4 in 
		# gas-phase
		OC = np.zeros((len(t_out), self.nc))

		# index for components with more than four carbons
		cni = (self.Cnum>4).reshape(1, -1)

		# tile over times
		cni = np.tile(cni, (y_mat.shape[0], 1))

		# copy of gas-phase number concentrations (molecules/cm^3)
		y_g = np.zeros((y_mat.shape[0], self.nc))
		y_g[:, :] = y_mat[:, 0:self.nc]

		# zero the components to be ignored (Cn<=4)
		y_g[cni == 0] = 0.

		# sum of considered components (molecules/cm^3)
		y_gsum = np.sum(y_g, axis=1).reshape(-1, 1)

		# times when sum is above zero
		ti = y_gsum[:, 0]>0.

		# normalise the number concentration of components
		# with Cn>4
		y_g[ti, :] = y_g[ti, :]/(y_gsum[ti].reshape(-1, 1))

		# multiply number concentration fractions by O:C ratio,
		# thereby weighting the O:C ratio
		OC = y_g*(self.OC[0, :].reshape(1, -1))

		# sum O:C across components to get 
		# number concentration-weighted O:C ratio
		OC = np.sum(OC, axis = 1)

		# create O:C variable
		OCvar = rootgrp.createVariable(str('O:C_for_C>4_weighted_by_number' + 
				'_concentration_in_air'), 'float32', ('time'))
		# set O:C units
		OCvar.setncattr('units', 'fraction')
		# set values for O:C variable
		OCvar[:] = OC

	# the rate of reactions (molecules/cm^3/s)	
	if ('reactionRate' in ud_var):

		# create variable for reaction rates
		rrvar = rootgrp.createVariable(str('rates_of_reaction_of_all_reactions'), 
								 'float32', ('time', 'reactions'))
		# set units
		rrvar.setncattr('units', str('molecules/cm^3/s'))
		
		# set values
		rrvar[:, :] = self.reactionRates[:, :]
		# ------

		# create variable for indices of components in gas-phase reactions,
		# with reactions in rows and reactants in columns
		rig_var = rootgrp.createVariable(str('reactant_indices'), 
								 'int', ('reactions', 'maximum_reactants'))
		# set units
		rig_var.setncattr('units', str('component indices for reactions (rows) and reactants (columns)'))
		
		# set values
		rig_var[:, :] = self.rindx_g[:, :]
		# ------

		# create variable for stoichiometry of components in gas-phase reactions,
		# with reactions in rows and reactants in columns
		stoi_g_var = rootgrp.createVariable(str('reactant_stoichiometries'), 
								 'float32', ('reactions', 'maximum_reactants'))
		# set units
		stoi_g_var.setncattr('units', str('stoichiometries of reactants for reactions (rows) and reactants (columns)'))
		
		# set values
		stoi_g_var[:, :] = self.rstoi_g[:, :]
		# ------

		# create variable for number of reactants per reaction
		nreac_var = rootgrp.createVariable(str('number_of_reactants_per_reaction'), 
								 'int', ('reactions'))
		# set units
		nreac_var.setncattr('units', str('number of reactants per reaction'))
		
		# set values
		nreac_var[:] = self.nreac_g[:]
		# ------

		# create variable for indices of components as products in gas-phase reactions,
		# with reactions in rows and reactants in columns
		pig_var = rootgrp.createVariable(str('product_indices'), 
								 'int', ('reactions', 'maximum_products'))
		# set units
		pig_var.setncattr('units', str('component indices for reactions (rows) and products (columns)'))
		
		# set values
		pig_var[:, :] = self.pindx_g[:, 0:max(self.nprod_g)]
		# ------

		# create variable for stoichiometry of components in gas-phase reactions,
		# with reactions in rows and products in columns
		stoi_p_var = rootgrp.createVariable(str('product_stoichiometries'), 
								 'float32', ('reactions', 'maximum_products'))
		# set units
		stoi_p_var.setncattr('units', str('stoichiometries of products for reactions (rows) and products (columns)'))
		
		# set values
		stoi_p_var[:, :] = self.pstoi_g[:, 0:max(self.nprod_g)]
		
		# ------

		# create variable for number of products per reaction
		nprod_var = rootgrp.createVariable(str('number_of_products_per_reaction'), 
								 'int', ('reactions'))
		# set units
		nprod_var.setncattr('units', str('number of products per reaction'))
		
		# set values
		nprod_var[:] = self.nprod_g[:]
		# ------

		# create variable for reaction strings
		eq_str_var = rootgrp.createVariable(str('equations_per_reaction'), 
								 'str', ('reactions'))
		# set units
		eq_str_var.setncattr('units', str('chemical equations per reaction'))
		
		# set values
		# note that because of variable numbers of 
		# characters per string, they must be assigned to the netCDF variable 
		# individually
		for stri in range(len(self.eqn_list)):
			eq_str_var[stri] = self.eqn_list[stri]
		# ------

	# end of bespoke saving function
	return(rootgrp)