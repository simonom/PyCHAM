##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2022 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
#                                                                                       											 #
#    All Rights Reserved.                                                                									 #
#    This file is part of PyCHAM                                                         									 #
#                                                                                        											 #
#    PyCHAM is free software: you can redistribute it and/or modify it under              						 #
#    the terms of the GNU General Public License as published by the Free Software       					 #
#    Foundation, either version 3 of the License, or (at your option) any later          						 #
#    version.                                                                            										 #
#                                                                                        											 #
#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT                						 #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       			 #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              				 #
#    details.                                                                            										 #
#                                                                                        											 #
#    You should have received a copy of the GNU General Public License along with        					 #
#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                 							 #
#                                                                                        											 #
##########################################################################################
'''module to save PyCHAM results to files'''
# module called following simulation in PyCHAM to store generated and input information

import os
import sys
import numpy as np
import csv
from shutil import copyfile

def saving(y_mat, Nresult_dry, Nresult_wet, t_out, savefolder, num_comp, 
	Cfactor_vst, testf, numsb, comp_namelist, y_mw, MV,
	time_taken, seed_name, x2, rbou_rec, wall_on, space_mode, rbou00, upper_bin_rad_amp, 
	indx_plot, comp0, yrec_p2w, rel_SMILES, Psat_Pa_rec, OC, H2Oi,
	siz_str, cham_env, tot_in_res_ft, self):

	# inputs: ----------------------------------------------------------------------------
	
	# y_mat - species (columns) concentrations with time (rows) (# molecules/cm3 (air))
	# Nresult_dry  - number concentration of dry particles per size bin (# particles/cm3 (air))
	# Nresult_wet  - number concentration of dry particles per size bin (# particles/cm3 (air))
	# Cfactor - conversion factor to change gas-phase concentrations from # molecules/cm3 
	# (air) into ppb
	# testf - flag to show whether in normal mode (0) or test mode (1)
	# numsb - number of size bins
	# self.dydt_vst - tendency to change of user-specified components
	# self.dydt_trak - user-input names of components to track
	# comp_namelist - names of components given by the chemical scheme file
	# upper_bin_rad_amp - factor upper bin radius found increased by in 
	#						Size_distributions.py for more than 1 size bin, or in 
	# 						pp_intro.py for 1 size bin
	# Cfactor_vst - one billionth the molecular concentration in a unit volume of chamber
	#				(molecules/cc) per recording time step
	# time_taken - computer time for entire simulation (s)
	# seed_name - name of seed component
	# y_mw - molecular weights (g/mol)
	# MV - molar volumes of all components (cm3/mol)
	# time_taken - simulation computer time (s)
	# seed_name - chemical scheme name of component comprising seed particles
	# x2 - record of size bin radii (um)
	# rbou_rec - radius bounds per size bin (columns) per time step (rows) (um)
	# wall_on - marker for whether wall on or off
	# space_mode - type of spacing between particle size bins (log or lin)
	# rbou00 - initial lower size (radius) bound of particles (um)
	# upper_bin_rad_amp - factor increase in radius of upper bound
	# indx_plot - index of components to plot gas-phase concentration temporal profiles of in 
	# standard results plot
	# comp0 - names of components to plot gas-phase concentration temporal profiles of in 
	# standard results plot
	# yrec_p2w - concentration of components on the wall due to 
	#	particle-wall loss, stacked by component first then by
	#	size bin (molecules/cc)
	# self.sch_name - path to chemical scheme file
	# self.inname - path to model variables file
	# rel_SMILES - SMILES strings for components in chemical scheme
	# Psat_Pa_rec - pure component saturation vapour pressures at 298.15 K
	# OC - oxygen to carbon ratio of components
	# H2Oi - index of water
	# self.seedi - index of seed components
	# siz_str - the size structure
	# cham_env - chamber environmental conditions (temperature (K), 
	# pressure (Pa) and relative humidity
	# self.RO2_indices[:, 1] - alkyl peroxy radical indices
	# self.RO_indx - alkoxy radical indices
	# tot_in_res_ft - record of cumulative inputs of injected components (ug/m3)
	# self - reference to PyCHAM
	# ---------------------------------------------------------------
	
	
	if ((numsb-wall_on) > 0): # correct for changes to size bin radius bounds
		rbou_rec[:, 0] = rbou00
		rbou_rec[:, -1] = rbou_rec[:, -1]/upper_bin_rad_amp

	if (testf == 1):
		return(0) # return dummy

	dir_path = os.getcwd() # current working directory
	output_root = 'PyCHAM/output'
	filename = os.path.basename(self.sch_name)
	filename = os.path.splitext(filename)[0]
	# one folder for one simulation
	output_by_sim = os.path.join(dir_path, output_root, filename, savefolder)
	# create folder to store results
	os.makedirs(output_by_sim)
	
	# create folder to store copies of inputs
	os.makedirs(str(output_by_sim+'/inputs'))
	# making a copy of the chemical scheme and model variables input files
	output_by_sim_sch_ext = str(output_by_sim+'/inputs/'+self.sch_name.split('/')[-1])
	copyfile(self.sch_name, output_by_sim_sch_ext)
	output_by_sim_mv_ext = str(output_by_sim+'/inputs/'+self.inname.split('/')[-1])
	if (self.inname != 'Default'): # if not in default model variables mode 
		copyfile(self.inname, output_by_sim_mv_ext)
	
	# saving dictionary
	# dictionary containing variables for model and components
	const = {}
	const["number_of_size_bins"] = numsb
	const["number_of_components"] = num_comp
	const["molecular_weights_g/mol_corresponding_to_component_names"] = (np.squeeze(y_mw[:, 0]).tolist())
	const["molar_volumes_cm3/mol"] = (MV[:, 0].tolist())
	const["organic_peroxy_radical_index"] = (self.RO2_indices[:, 1].tolist())
	const["organic_alkoxy_radical_index"] = self.RO_indx
	const["chem_scheme_names"] = comp_namelist
	const["SMILES"] = rel_SMILES
	const["factor_for_multiplying_ppb_to_get_molec/cm3_with_time"] = (Cfactor_vst.tolist())
	const["simulation_computer_time(s)"] = time_taken
	const["seed_name"] = seed_name
	const["wall_on_flag_0forNO_1forYES"] = wall_on
	const["space_mode"] = space_mode
	const["pure_component_saturation_vapour_pressures_at_298.15K"] = Psat_Pa_rec.tolist()
	const["oxygen_to_carbon_ratios_of_components"] = OC.tolist()
	const["index_of_water"] = H2Oi
	const["index_of_seed_components"] = self.seedi.tolist()
	const["size_structure_0_for_moving_centre_1_for_full_moving"] = siz_str
	const["output_by_sim_sch_ext"] = output_by_sim_sch_ext
	const["output_by_sim_mv_ext"] = output_by_sim_mv_ext

	with open(os.path.join(output_by_sim,'model_and_component_constants'),'w') as f:
		for key in const.keys():
			f.write("%s,%s\n"%(key, const[key]))
	
	# convert gas-phase concentrations from # molecules/cm3 (air) into ppb
	# leaving any particle-phase concentrations as molecules/cc (air)
	y_mat[:, 0:num_comp] = y_mat[:, 0:num_comp]/(Cfactor_vst.reshape(len(Cfactor_vst), 1))
	
	
	y_header = str('') # prepare header for concentrations with time file 
	x2_header = str('') # prepare header for files relating to size bins
	
	for i in range(numsb+1): # loop through size bins

		if i == 0:
			end = '_g'
		if ((i > 0) and (i < numsb)):
			end = str('_p' + str(i))
			x2_header = str(x2_header + str(i))
		if (i == numsb):
			if (wall_on == 0):
				end = str('_p'  + str(i))
				x2_header = str(x2_header+str(np.repeat(i, num_comp)))
			else:
				end = '_w'
		for ii in range(num_comp):
			if i  == 0 and ii == 0:
				start = ''
			else:
				start = ', '
			y_header = str(y_header + str(start + comp_namelist[ii]) + end)
			
	# saving both gas, particle and wall concentrations of components
	np.savetxt(os.path.join(output_by_sim, 'concentrations_all_components_all_times_gas_particle_wall'), y_mat, delimiter=',', header=str('time changes with rows which correspond to the time output file, components in columns, with _g representing gas phase (ppb), _pi representing particle phase where i is the size bin number (starting at 1) (molecules/cm3 (air)) and _w is the wall phase (molecules/cm3 (air))\n' + y_header)) 		
	 
	# saving time of outputs
	np.savetxt(os.path.join(output_by_sim, 'time'), t_out, delimiter=',', header='time (s), these correspond to the rows in the concentrations_all_components_all_times_gas_particle_wall, particle_number_concentration and size_bin_radius output files')
	
	# saving environmental conditions (temperature, pressure, relative humidity)
	np.savetxt(os.path.join(output_by_sim, 'chamber_environmental_conditions'), cham_env, delimiter=',', header='chamber environmental conditions throughout the simulation, with rows corresponding to the time points in the time output file, first column is temperature (K), second is pressure (Pa) and third is relative humidity (fraction (0-1))')
	
	# saving the index and names of components whose gas-phase temporal profiles can be plotted on the standard results plot
	fname = os.path.join(output_by_sim, 'components_with_initial_gas_phase_concentrations_specified')
	np.savetxt(fname, [indx_plot, comp0], delimiter =', ', header='index (top row) and chemical scheme name (bottom row) of components with initial gas-phase concentrations specified', fmt ='% s') 
	
	# if tracking of tendencies to change requested by user, loop through the components
	# and save the tendency record for each of these (%/hour)
	if (len(self.dydt_vst) > 0):
		compind = 0
		dydtnames = self.dydt_vst['comp_names'] # get names of tracked components
		# loop through components to record the tendency of change \n')
		for comp_name in dydtnames:
			# open relevant dictionary value, to get the 2D numpy array for saving
			dydt_rec = np.array(self.dydt_vst[str(str(comp_name) + '_res')])
			
			# get user-input name of this component
			comp_name = str(self.dydt_trak[compind] +'_rate_of_change')
			# save
			np.savetxt(os.path.join(output_by_sim, comp_name), dydt_rec, delimiter=',', header='tendency to change, top row gives equation number (where number 0 is the first equation), penultimate column is gas-particle partitioning and final column is gas-wall partitioning (molecules/cc.s (air))')
			compind += 1
	
	# saving generation of components
	np.savetxt(os.path.join(output_by_sim, 'component_generation'), self.gen_num, delimiter=',', header='generation number of each component (where the initial unoxidised VOC is generation number 0), with the order corresponding to that of components in the concentrations_all_components_all_times_gas_particle_wall file.')
	
	if ((numsb-wall_on) > 0): # if particles present
	
		# saving the concentration of components on the wall due to particle deposition to wall
		np.savetxt(os.path.join(output_by_sim, 'concentrations_all_components_all_times_on_wall_due_to_particle_deposition_to_wall'), yrec_p2w, delimiter=',', header=str('concentration of components on wall due to particle deposition to wall (molecules/cc (air)) time changes with rows which correspond to the time output file, components in columns and size bin changing with columns with size bin numbers given in the second row of the header\n'+ x2_header)) 
	
		np.savetxt(os.path.join(output_by_sim, 'particle_number_concentration_dry'), Nresult_dry, delimiter=',',
				header=('particle number concentration assuming water removed from particles (#/cc (air)), with time changing with rows (corresponding times given in the time output file) and size bin changing with columns with size bin numbers given in the second row of the header\n'+x2_header))
		
		np.savetxt(os.path.join(output_by_sim, 'particle_number_concentration_wet'), Nresult_wet, delimiter=',',
				header=('particle number concentration assuming water not removed from particles (#/cc (air)), with time changing with rows (corresponding times given in the time output file) and size bin changing with columns with size bin numbers given in the second row of the header\n'+x2_header))	
	
		np.savetxt(os.path.join(output_by_sim, 'size_bin_radius'), x2, delimiter=',',
				header= str('particle radii (um) per size_bin (including water contribution to size), with size bins represented by columns and their number (starting from 1) given in second line of header, per time step which is represented by rows and corresponding times given in the time output file \n'+x2_header))
	
		np.savetxt(os.path.join(output_by_sim, 'size_bin_bounds'), rbou_rec, delimiter=',',
				header=str('particle size bin bounds (um), with size bins represented by columns and their number (starting at 1 and in line with the lower bound) given in second line of header, per time step which is is represented by rows and corresponding times given in the time output file \n'+x2_header))		
	
	
	np.savetxt(os.path.join(output_by_sim, 'total_concentration_of_injected_components'), tot_in_res_ft, delimiter=',',
				header=str('the total concentration (ug/m3) of injected (through initial gas-phase concentration, instantaneous and/or continuous gas-phase influx) components, with component index (relative to all components) in the first row and its cumulative injected concentrations in following rows'))

	# if save name is the default, then remove to ensure no duplication in future
	if (savefolder == 'default_res_name'):
		import shutil	
		shutil.rmtree(output_by_sim)
	return()