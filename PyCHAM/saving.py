'''module to save PyCHAM results to files'''

import os
import sys
import numpy as np

def saving(filename, y_mat, t_out, Nresult_dry, Nresult_wet, x2, numsb, y_mw, num_speci, 
			savefolder, rbou, Cfactor, MV, testf, dydt_vst, dydt_trak, spec_namelist, 
			rbou00):
			
	#  -------------------------------
	# inputs:
	
	# filename - name of inputs file name
	# y_mat - species (columns) concentrations with time (rows) (molecules/cc (air))
	# Nresult_dry  - number concentration of dry particles per size bin (#/cc (air))
	# Nresult_wet  - number concentration of dry particles per size bin (#/cc (air))
	# Cfactor - conversion factor to change gas-phase concentrations from molecules/cc 
	# (air) into ppb
	# testf - flag to show whether in normal mode (0) or test mode (1)
	# dydt_vst - tendency to change of user-specified components
	# dydt_trak - user-input names of components to track
	# spec_namelist - names of species given by the equation file
	# rbou00 - initial lowest radius bound (um)
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
	
	
	# constants dictionary
	const = {}
	const["num_sb"] = numsb
	const["num_speci"] = num_speci
	const["mw"] = (np.squeeze(y_mw[:, 0]).tolist())
	const["mv"] = (MV[:, 0].tolist())
	const["spec_namelist"] = spec_namelist
	const["Cfactor"] = Cfactor
	with open(os.path.join(output_by_sim,'constants'),'w') as f:
		for key in const.keys():
			f.write("%s,%s\n"%(key,const[key]))
	
	# convert gas-phase concentrations from molecules/cc (air) into ppb
	# leaving any particle-phase concentrations as molecules/cc (air)
	y_mat[:, 0:num_speci] = y_mat[:, 0:num_speci]
	# saving both gas- and particle-phase concentrations of species
	np.savetxt(os.path.join(output_by_sim,'y'),y_mat,delimiter=',',header='conc_vs_time') 	
	# saving time of outputs
	np.savetxt(os.path.join(output_by_sim,'t'),t_out,delimiter=',',header='time (s)') 	
	
	# if tracking of tendencies to change requested by user, loop through the components
	# and save the tendency record for each of these (%/hour)
	if len(dydt_vst)>0:
		compind = 0
		# loop through components to record the tendency of change \n')
		for compi in dydt_vst.get('comp_index'):
			# open relevant dictionary value, to get the 2D numpy array for saving
			dydt_rec = np.array(dydt_vst.get(compi))
			
			# get user-input name of this component
			comp_name = str(dydt_trak[compind] +'_dydt')
			# save
			np.savetxt(os.path.join(output_by_sim, comp_name), dydt_rec, delimiter=',', header='tendency to change, top row gives equation number (molecules/cc.s (air))')
			compind += 1
	
	
	if numsb>1: # if particles present (numsb includes wall)
		# ensure N has expected shape for plotting
		Nresult_dry = Nresult_dry.reshape(len(t_out), numsb-1)
		np.savetxt(os.path.join(output_by_sim,'N_dry'), Nresult_dry, delimiter=',',
					header='dry_particle_number_conc (#/cc (air))')
		Nresult_wet = Nresult_wet.reshape(len(t_out), numsb-1)
		np.savetxt(os.path.join(output_by_sim,'N_wet'), Nresult_wet, delimiter=',',
					header='wet_particle_number_conc (#/cc (air))')
		np.savetxt(os.path.join(output_by_sim,'x'),x2,delimiter=',',
					header='particle_size')
		# reverse change done in Size_distributions on uppermost boundary and change done
		# in pp_intro on lowermost boundary (um)
		rbou[-1] = rbou[-1]*1.0e-6
		rbou[0] = rbou00
		
		np.savetxt(os.path.join(output_by_sim,'sbb'),rbou,delimiter=',',
					header='particle_bin_bounds')		
		
	return(output_by_sim)