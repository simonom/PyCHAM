'''plotting function for showing the mole fraction of atoms/functional groups contained in each component'''
# mole fraction as a function of time

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import retr_out
import numpy as np

def plotter(caller, dir_path, atom_name, atom_num, self): # define function

	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0) or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# atom_name - SMILE string names of atom/functional group to target
	# atom_num - number of top contributors to plot
	# self - reference to GUI
	# --------------------------------------------------------------------------
	
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, SMILES, 
		y_mw, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, _, _, _, _, _) = retr_out.retr_out(dir_path)
	
	# empty lists to contain results
	cnt_list = []
	ind_list = []
	
	cn = 0 # count on components
	for ci in SMILES: # loop through component names
		at_cnt = ci.count(atom_name)
		if (at_cnt > 0): # if a contributor
			cnt_list.append(at_cnt)
			ind_list.append(int(cn))
			
		cn += 1
	
	# empty results array for contributing component index and number of occurrences
	res = np.zeros((len(cnt_list), 2))
	res[:, 0] = np.array((ind_list)) # index in first column
	res[:, 1] = np.array((cnt_list)) # count in second column
	res = res.astype(int)
	
	# reshape so that time in rows and components per size bin in columns
	yrec = yrec.reshape(len(timehr), num_comp*(num_sb+1))
	# convert gas-phase concentrations to molecules/cm3 from ppb
	yrec[:, 0:num_comp] = yrec[:, 0:num_comp]*(np.array((Cfac)).reshape(-1, 1))
	
	# extract just the relevant concentrations
	yrel = np.zeros((len(timehr), res.shape[0])) # molecules/cm3
	nam_rel = []
	cin = 0 # count on the relevant components
	for ci in res[:, 0]: # loop through relevant components
		
		yrel[:, cin] = np.sum(yrec[:, ci::num_comp], axis = 1)*res[cin, 1]
		nam_rel.append(comp_names[ci])
		cin += 1
	
	# identify the components to be plotted
	ytot = np.sum(yrel, axis=0)
	ytot_asc = np.sort(ytot)
	cutoff = ytot_asc[-atom_num] # the cutoff for the requested number of contributors
	yreln = yrel[:, ytot>=cutoff]
	nam_rel_indx = (np.where(ytot>=cutoff))[0]
	
	nam_up = [] # empty list
	for nri in nam_rel_indx:
		nam_up.append(nam_rel[nri])
	
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7)) # prepare figure
	
	for ci in range(len(yreln[0, :])): # loop through relevant components
		if (any(yreln[:, ci] > 0.)):
			ax0.semilogy(timehr, yreln[:, ci]/np.sum(yrel, axis=1), label = nam_up[ci])

	# details of plot
	ax0.set_ylabel(r'Fraction contribution', fontsize = 14)
	ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
	ax0.set_title(str('Fraction contribution to the atom/functional group ' + str(atom_name)), fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.legend(fontsize = 14, loc = 'lower right')
	plt.show()

	return()