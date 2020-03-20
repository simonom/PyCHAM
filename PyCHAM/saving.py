'''module to save PyCHAM results to files'''

import os
import sys
import numpy as np

def saving(filename, y_mat, t_out, Nresult, x2, numsb, y_mw, num_speci, 
			savefolder, rbou, Cfactor, MV, testf, spec_namelist):
			
	#  -------------------------------
	# inputs:
	
	# filename - name of inputs file name
	# y_mat - species (columns) concentrations with time (rows) (molecules/cc (air))
	# x2 - particle radius per size bin (columns), per time step (row) (um) including 
	#		water
	# Cfactor - conversion factor to multiply gas-phase concentrations in ppb to get 
	# molecules/cc (air)
	# MV - molar volume (cm3/mol) (set in pp_intro)
	# testf - flag to show whether in normal mode (0) or test mode (1)
	# spec_namelist - names of species given by the equation file
	# -------------------------------
	
	if testf==1:
		return(0) # return dummy
	dir_path = os.getcwd()
	  
	output_root = 'PyCHAM/output'
	filename = os.path.basename(filename)
	filename = os.path.splitext(filename)[0]
	# one folder for one simulation
	output_by_sim = os.path.join(dir_path, output_root, filename, savefolder)
	# create folder to store results
	os.makedirs(output_by_sim)
	
	# dictionary containing constants for model and components
	const = {}
	const["number_of_size_bins"] = numsb
	const["number_of_components"] = num_speci
	const["molecular_weights_g/mol_corresponding_to_component_names"] = (np.squeeze(y_mw[:, 0]).tolist())
	const["molecular_volumes_cm3/mol"] = (MV[:, 0].tolist())
	const["component_names"] = spec_namelist
	const["factor_for_multiplying_ppb_to_get_molec/cm3"] = Cfactor
	with open(os.path.join(output_by_sim,'model_and_component_constants'),'w') as f:
		for key in const.keys():
			f.write("%s,%s\n"%(key, const[key]))
	
	# convert gas-phase concentrations from molecules/cc (air) into ppb
	# leaving any particle-phase concentrations as molecules/cc (air)
	y_mat[:, 0:num_speci] = y_mat[:, 0:num_speci]/Cfactor
	
	# prepare header for concentrations file
	y_header = str('')
	x2_header = str('') # prepare header for files relating to size bins
	for i in range(numsb+1):

		if i == 0:
			end = '_g'
		if i > 0 and i<numsb:
			end = str('_p'+str(i))
			if i!=numsb-1:
				x2_header = str(x2_header+str(i)+', ')
			else:
				x2_header = str(x2_header+str(i))
		if i == numsb:
			end = '_w'
		for ii in range(num_speci):
			if i  == 0 and ii == 0:
				start = ''
			else:
				start = ', '
			y_header = str(y_header+str(start+spec_namelist[ii])+end)
			
	# saving both gas- and particle-phase concentrations of species
	np.savetxt(os.path.join(output_by_sim, 'concentrations_all_components_all_times_gas_particle_wall'), y_mat, delimiter=',', header=str('time changes with rows which correspond to the time output file, species in columns, with _g representing gas phase (ppb), _pi representing particle phase where i is the size bin number (starting at 1) (molecules/cc (air)) and _w is the wall phase (molecules/cc (air))\n'+y_header)) 	
	# saving time of outputs
	np.savetxt(os.path.join(output_by_sim,'time'),t_out,delimiter=',',header='time (s), these correspond to the rows in the concentrations_all_components_all_times_gas_particle_wall, particle_number_concentration and size_bin_radius output files') 	
		

	if numsb>1: # if particles present
		
		np.savetxt(os.path.join(output_by_sim, 'particle_number_concentration'), Nresult,delimiter=',',
					header=('particle number concentration (#/cc (air)), with time changing with rows (corresponding times given in the time output file) and size bin changing with columns with size bin numbers given in the second row of the header\n'+x2_header))	
		
		np.savetxt(os.path.join(output_by_sim, 'size_bin_radius'), x2, delimiter=',',
					header= str('particle radii (um) per size_bin (including water contribution to size), with size bins represented by columns and their number (starting from 1) given in second line of header, per time step which is represented by rows and corresponding times given in the time output file \n'+x2_header))
		
		np.savetxt(os.path.join(output_by_sim, 'size_bin_bounds'), rbou.reshape(1, -1), delimiter=',',
					header=str('particle size bin bounds (um), with size bin number (starting at 1 and in line with the lower bound) given in second line of header\n'+x2_header))		
		
	return(output_by_sim)