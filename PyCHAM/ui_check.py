'''checking that user inputs are valid'''

import os
import sys
import numpy as np

def ui_check(sav_nam, sch_name, wall_on, caller, siz_str, num_sb, pmode, pconc, pconct, lowsize, std, mean_rad, new_partr, chamSA, chem_sch_mark, af_path, int_tol, update_stp, tot_time, RH, uman_up, light_stat, light_time, injectt, Ct, dens_comp, dens, seed_name, seedVr, seed_diss, partit_cutoff):

	# inputs: ------------------------------------------------------------
	# sav_nam - name of folder to save results to
	# sch_name - name of chemical scheme file
	# wall_on - marker for whether wall on or off
	# caller - marker for the calling module
	# siz_str - the size structure
	# num_sb - number of particle size bins
	# pmode - whether particle number concentrations expressed by mode or explicitly
	# pconc - number concentration of particles
	# pconct - times of particle injection (s)
	# lowsize - lower bound of particle sizes (um)
	# std - standard deviation for particle sizes
	# mean_rad - mean radius of particle size distributions (um)
	# new_partr - radius of newly nucleated particles (cm)
	# chamSA - chamber surface area (m2)
	# chem_sch_mark - markers in chemical scheme
	# af_path - actinic flux path
	# int_tol - integration tolerances
	# update_stp - time step for updating ODE initial conditions (s)
	# tot_time - total experiment time (s)
	# RH - relative humidity (0-1)
	# uman_up - marker for whether to update UManSysProp folder
	# light_stat - marker for whether lights on or off
	# light_time - time that lights attain status
	# injectt - times of instantaneous injections of gas-phase components (s)
	# Ct - concentrations of components instantaneously injected (ppb)
	# dens_comp - chemical scheme names of components with density 
	# manually assigned
	# dens - manually assigned densities (g/cc)
	# seed_name - name of component(s) comprising seed particles
	# seedVr - volume ratio of component(s) comprising seed particles
	# seed_diss - dissociation constant for seed particle component(s)
	# partit_cutoff - product of vapour pressure and activity coefficient
	# 	above which gas-particle partitioning assumed zero 
	# --------------------------------------------------------------------

	print('Checking user inputs')

	dir_path = os.getcwd() # current working directory
	output_root = 'PyCHAM/output'
	filename = os.path.basename(sch_name)
	filename = os.path.splitext(filename)[0]
	# one folder for one simulation
	output_by_sim = os.path.join(dir_path, output_root, filename, sav_nam)
	
	# constrain wall_on marker
	if (wall_on>=1):
		wall_on = 1
	else:
		wall_on = 0 

	if os.path.isdir(output_by_sim) == True and caller == 0:
		print('Error: results file name (' +output_by_sim+ ') already exists, please use an alternative')
		sys.exit()

	# ensure size structure marker is sensible
	if siz_str<0 or siz_str>1:
		siz_str = 0
	
	# consistency between number of particle size bins and particle number concentration
	if num_sb == 0 and (sum(pconc != 0) > 0):
		pconc[:] = 0.0
		print('Note that zero particle size bins registered (number_size_bins in model variables input file), however total particle number concentration (pconc in model variables input file) contains a non-zero value, therefore assuming no particle size bins')
	
	# if lower bound of particle sizes set to 0, this will induce error when taking log10
	# in pp_intro, so change to very small value (um)
	if lowsize == 0.:
		lowsize = 1.e-3

	if (num_sb > 0): # if size bins present

		if (pmode == 0): # number concentrations expressed as mode
			# if mean_rad inconsistent with pconc
			if (mean_rad.shape != pconc.shape):
				mean_rad = np.ones((pconc.shape))*-1.e6 # default flag
			if (std.shape != pconc.shape):
				std = np.ones((pconc.shape))*1.2 # default value  

	# if radius of newly nucleated particles empty, set to default value of 1 nm (cm)
	if new_partr == 0:
		new_partr = 2.e-7

	# convert chamber surface area (m2) to spherical equivalent radius (m)
	# (below eq. 2 in Charan (2018))
	chamR = (chamSA/(4.0*np.pi))**0.5

	if (len(chem_sch_mark) == 0): # if empty default to MCM kpp format
		chem_sch_mark = ['{', 'RO2', '+', 'C(ind_', ')', '', '&', '', '', ':', '}', ';']

	if not af_path: # actinic flux path
		af_path = str('no')

	if not int_tol: # integration tolerances
		int_tol = [1.e-3, 1.e-4]

	if not update_stp: # update step (s)
		update_stp = 6.e1

	if not tot_time: # total experiment time (s)
		tot_time = 3.6e3*4.

	if RH>1.0:
		print('Note: RH set above 1.0; simulation will be attempted, but please note that RH is interpreted as fraction, not a percentage, where 1 represents saturation of water vapour')

	# for UManSysProp if no update requested, check that there 
	# is an existing UManSysProp folder
	if uman_up == 0: # check for existing umansysprop folder
		if not os.path.isdir(os.getcwd() + '/umansysprop'):
			uman_up = 1
	
	if uman_up == 1: # test whether UManSysProp can be updated
		import urllib.request # module for checking internet connection
		# function for testing internet connection
		def connect(host='https://github.com/loftytopping/UManSysProp_public.git'):
			try:
				urllib.request.urlopen(host) #Python 3.x
				return True
			except:
				return False
		# test internet connection
		if connect():
			print('Internet connection confirmed and either user has requested cloning of UManSysProp via the model variables input file or no existing UManSysProp folder was found') 
		else:
			print('Error: either user has requested cloning of UManSysProp via the model variables input file or no existing UManSysProp folder has been found, but connection to the page failed, possibly due to no internet connection (UManSysProp repository site: https://github.com/loftytopping/UManSysProp_public.git)')
			sys.exit()

	# ensure that if multiple instantaneous injections of particles, that corresponding variables
	# have the correct shape, specifically that they cover the same number of times
	if (pconc.shape[1] != pconct.shape[1] or pconc.shape[1] != std.shape[1] or pconc.shape[1] != mean_rad.shape[1] or pconct.shape[1] != std.shape[1] or pconct.shape[1] != std.shape[1] or std.shape[1] != mean_rad.shape[1]):
		print(str('Error: inconsistent number of times for instantaneous injection of particles represented by model variable inputs (number of times represented in brackets) for: pconc ('+str(pconc.shape[1])+'), pconct ('+str(pconct.shape[1])+'), mean_rad ('+str(mean_rad.shape[1])+') and/or std ('+str(std.shape[1])+').  Please see README for guidance.'))
		sys.exit()
	
	if (len(light_stat) != len(light_time)):
		print('Error: length of input model variables light_status and light_time have different lengths and they should be the same, please see README for guidance.')
		sys.exit()

	
	if (light_time[0] != 0): # if no status provided at simulation start default to lights off
		light_stat = np.concatenate((np.zeros((1)).astype(int), light_stat))
		light_time = np.concatenate((np.zeros((1)), light_time))

	# check consistency between number of times components instantaneously injected and
	# the number of concentrations given
	if (len(injectt) != Ct.shape[1]):
		print('Error: number of times given for the instantaneous injection of gas-phase components (injectt in the model variables input file) is inconsistent with the number of concentrations of instantaneous injection provided (Ct in the model variables input file), please see the README for guidance')
		sys.exit()

	# check on consistency of manually assigned densities
	if (len(dens_comp) != len(dens)):
		print('Error: the number of chemical scheme names provided in the dens_Comp model variables input does not match the number of densities provided in the dens model variables input, please see README for details')
		sys.exit()

	# check on consistency of names of seed component(s) and their volume ratio
	if (len(seed_name) != len(seedVr)):
		if (len(seed_name) > 1) and (len(seedVr) == 1):
			seedVr = np.ones((len(seed_name)))
		else:
			print('Error: the number of seed particle component names (seed_name in model variables input file) and the number of seed particle component volume ratios (seedVr in model variables input file) are inconsistent, please see README for guidance')
			sys.exit()

	# check on consistency of names of seed component(s) and their dissociation constant
	if (len(seed_name) != len(seed_diss)):
		if (len(seed_name) > 1) and (len(seed_diss) == 1):
			seed_diss = np.ones((len(seed_name)))
		else:
			print('Error: the number of seed particle component names (seed_name in model variables input file) and the number of seed particle component dissociation constants (seed_diss in model variables input file) are inconsistent, please see README for guidance')
			sys.exit()

	if (len(partit_cutoff) > 1):
		print('Length of the model variables input partit_cutoff should have a maximum length of one, but is greater, please see README for guidance')
		sys.exit()

	return(wall_on, pconc, lowsize, std, mean_rad, new_partr, chamR, chem_sch_mark, af_path, int_tol, update_stp, tot_time, siz_str, light_stat, light_time, seedVr, seed_diss, uman_up)
