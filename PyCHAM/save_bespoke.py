########################################################################
#								       #
# Copyright (C) 2018-2025					       #
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

def bespoke_saving(savei, y_mat, y_MM, t_out, cham_env, rootgrp, self):
	
	# inputs: -------------------------------------------------------
	# savei - output index to save
	# y_mat - component (columns) concentrations with time (rows) 
	#	(# molecules/cm3 (air))
	# y_MM - molar masses of components (g/mol)
	# t_out - the times (s) through simulation that outputs
	# 	correspond to
	# cham_env - chamber environmental conditions (temperature (K), 
	# pressure (Pa), relative humidity (0-1), transmission factor of 
	#	light (0-1)
	# rootgrp - the netCDF file object
	# self - the PyCHAM object
	# ---------------------------------------------------------------

	# if first index, then save the unconditional variables
	if (savei == 0):

		# unconditional saving of: time, physical conditions 
		# (temperature, pressure, RH, transmission factor of light)

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

		# set the surface number dimension
		rootgrp.createDimension('surface_number', self.wall_on)

		# set the component dimension
		rootgrp.createDimension('components', len(self.comp_namelist))
	
		# end of setting dimensions ---------------------------------

		# create time variable
		timevar = rootgrp.createVariable('time', 'float32', ('time'))
		# set units for time variable
		timevar.setncattr('units', time_unit)
		# set values for time variable (s)
		timevar[:] = t_out

		# create time variable
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

		return(rootgrp)

	else:
		# get name of user-defined variable to save
		ud_var = str(self.user_output[savei]).strip()

	# move onto conditional saving ---------------------------------

	# the gas-phase concentration(s) (molecules/cm3)	
	if ('_g' in ud_var):
		
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
			# get gas-phase concentration with time (molecules/cm3)
			c_g = y_mat[:, ci]
			# set values for variable
			c_gvar[:] = c_g

			# set units of concentration
			c_gvar.setncattr('units', 'molecules/cm^3')

	# the particle-phase concentration(s) (molecules/cm3)	
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

	# the surface-phase concentration(s) (molecules/cm3)	
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

	# the mass concentration of secondary organic aerosol (ug/m3)	
	if (ud_var == 'SOA'):

		
		# isolate just particle-phase component concentrations
		# (molecules/cm3)
		y_pp = np.zeros((y_mat.shape[0], self.nc*self.nasb))
		y_pp[:, :] = y_mat[:, self.nc:self.nc*(self.nasb+1)]

		# sum concentrations of individual components over size bins
		# loop over size bins beyond first
		for sbi in range(1, self.nasb):
			y_pp[:, 0:self.nc] += y_pp[self.nc*sbi:self.nc*(sbi+1)]

		# keep just the total concentration in particle phase
		# (molecules/cm3)
		y_pp = y_pp[:, 0:self.nc]


		# zero the concentrations of components without a carbon 
		# (molecules/cm3)
		y_pp[:, self.HC[0, :] == 0] = 0.

		# zero the seed components (molecules/cm3)
		y_pp[:, self.seedi] = 0

		# tile molar masses over times
		y_MMtiled = np.tile(y_MM.reshape(1, -1), (y_pp.shape[0], 1))

		# convert concentrations of organics from molecules/cm3 to ug/m3
		y_pp[:, :] = ((y_pp[:, :]/si.N_A)*y_MMtiled)*1.e12

		# sum mass concentrations over components to get SOA concentration (ug/m3)
		y_pp = np.sum(y_pp, axis= 1)

		# create SOA variable
		soavar = rootgrp.createVariable(str('mass_concentration_of_secondary_' + 
		'particulate_organic_matter_dry_aerosol_particles_in_air'), 
		'float32', ('time'))
		# set SOA units
		soavar.setncattr('units', 'ug/m^3')
		# set values for SOA variable
		soavar[:] = y_pp
	
	# end of bespoke saving function
	return(rootgrp)
