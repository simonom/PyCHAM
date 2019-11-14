'''module to save PyCHAM results to files'''

import os
import sys
import numpy as np

def saving(filename, y_mat, t_out, Nresult, x2, numsb, y_mw, num_speci, 
			savefolder, rbou, Cfactor, MV, testf):
			
	#  -------------------------------
	# inputs:
	
	# filename - name of inputs file name
	# y_mat - species (columns) concentrations with time (rows) (molecules/cc (air))
	# Cfactor - conversion factor to change gas-phase concentrations from molecules/cc 
	# (air) into ppb
	# testf - flag to show whether in normal mode (0) or test mode (1)
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
	
	
	# constants
	
	const = np.zeros((3, num_speci))
	const[0,0] = numsb
	const[0,1] = num_speci
	const[1, :] = y_mw[:, 0]
	const[2, :] = MV[:, 0]
	np.savetxt(os.path.join(output_by_sim,'constants'),const,delimiter=',',
			header='constants')
	
	# convert gas-phase concentrations from molecules/cc (air) into ppb
	# leaving any particle-phase concentrations as molecules/cc (air)
	y_mat[:, 0:num_speci] = y_mat[:, 0:num_speci]/Cfactor
	# saving both gas- and particle-phase concentrations of species
	np.savetxt(os.path.join(output_by_sim,'y'),y_mat,delimiter=',',header='conc_vs_time') 	
	# saving time of outputs
	np.savetxt(os.path.join(output_by_sim,'t'),t_out,delimiter=',',header='time (s)') 	
		

	if numsb>1: # if particles present
		np.savetxt(os.path.join(output_by_sim,'N'),Nresult,delimiter=',',
					header='particle_number_conc')	
		np.savetxt(os.path.join(output_by_sim,'x'),x2,delimiter=',',
					header='particle_size')
		np.savetxt(os.path.join(output_by_sim,'sbb'),rbou,delimiter=',',
					header='particle_bin_bounds')		
		
	return(output_by_sim)