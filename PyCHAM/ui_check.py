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
'''checking that user inputs are valid and providing appropriate modification to the PyCHAM GUI'''
# once model variables are validated, this module also allows options in the GUI for continuing with simulations

import os
import sys
import numpy as np
import pickle
import chem_sch_SMILES
from PyQt5.QtWidgets import *
from PyQt5.QtGui import  *
from PyQt5.QtCore import *
import write_hyst_eq
import hyst_eq
import importlib

def ui_check(self):

	# inputs: ----------------------------------------------------------------------
	# self - reference to GUI
	# --------------------------------------------------------------------------------

	# path to store for variables
	input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
	with open(input_by_sim, 'rb') as pk:
		[sav_nam, update_stp, 
		tot_time, comp0, y0, RH, RHt, Press,
		Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
		save_step, Compt, injectt, Ct, seed_name,
		seed_mw, seed_diss, seed_dens, seedx,
		con_infl_t, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
		accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
		nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
		inflectk, chamSA, Rader, p_char, e_field, partit_cutoff, ser_H2O, 
		wat_hist, drh_str, erh_str, pcont, Vwat_inc, seed_eq_wat, z_prt_coeff, 
		chamV] = pickle.load(pk)
		pk.close()	

	# loaded variables: ------------------------------------------------------------
	# sav_nam - name of folder to save results to
	# self.sch_name - name of chemical scheme file
	# self.wall_on - marker for whether wall on or off
	# siz_stru - the size structure
	# num_sb - number of particle size bins
	# pmode - whether particle number concentrations expressed by mode or explicitly
	# pconc - number concentration of particles
	# pconct - times of particle injection (s)
	# lowsize - lower bound of particle sizes (um)
	# uppsize - upper bound of particle sizes (um)
	# std - standard deviation for particle sizes
	# mean_rad - mean radius of particle size distributions (um)
	# new_partr - radius of newly nucleated particles (cm)
	# chamSA - chamber surface area (m2)
	# self.chem_sch_mark - markers in chemical scheme
	# self.af_path - actinic flux path
	# int_tol - integration tolerances
	# update_stp - time step for updating ODE initial conditions (s)
	# tot_time - total experiment time (s)
	# RH - relative humidity (0-1)
	# uman_up - marker for whether to update UManSysProp folder
	# self.light_stat - marker for whether lights on or off
	# self.light_time - time that lights attain status
	# injectt - times of instantaneous injections of gas-phase components (s)
	# Ct - concentrations of components instantaneously injected (ppb)
	# dens_comp - chemical scheme names of components with density 
	# manually assigned
	# dens - manually assigned densities (g/cc)
	# seed_name - name of component(s) comprising seed particles
	# seedx - mole fraction of dry seed components
	# seed_diss - dissociation constant for seed particle component(s)
	# partit_cutoff - product of vapour pressure and activity coefficient
	# 	above which gas-particle partitioning assumed zero
	# wat_hist - flag for particle-phase history with respect to water 
	#	(0 on the deliquescence curve, 1 on the efflorescence curve)
	# con_infl_t - times of continuous influx of components
	# self.con_infl_nam - chemical scheme names of components with continuous influx
	# self.con_infl_C - influx rate of components with continuous influx (ppb/s)
	# Vwat_inc - flag for whether (1) or not (0) the number size 
	# distribution of seed particles includes the volume of water
	# seed_eq_wat - whether (1) or not (0) to allow water to equilibrate with
	# seeds before experiment starts
	# z_prt_coeff - fraction of total gas-particle partitioning coefficient
	#	below which partitioning treated as negligible, e.g. because a 
	#	size bin has relatively very small surface area
	# self.testf - flag for testing and plotting checks of inputs
	# --------------------------------------------------------------------
	
	# to begin assume no errors, so message is that simulation ready
	err_mess = 'Model variables fine - simulation ready'
	em_flag = 0 # flag that no problematic errors detected
	
	# saving path (copied from saving module) ------------
	dir_path = os.getcwd() # current working directory
	output_root = 'PyCHAM/output'
	filename = os.path.basename(self.sch_name)
	filename = os.path.splitext(filename)[0]
	# one folder for one simulation
	output_by_sim = os.path.join(dir_path, output_root, filename, sav_nam)
	# -------------------------------------------------------------------
	
	# incorrect seed molecular weight
	if seed_mw[0] == 'fail':
		if (em_flag < 2):
			err_mess = str('Error - the molecular weight of seed particle component (seed_mw in model variables input file) could not be interpreted, please check it adheres to the guidelines in README')
		em_flag = 2
	
	# incorrect wall_on marker
	if (self.wall_on != 1 and self.wall_on != 0):
		if (em_flag < 2):
			err_mess = str('Error - wall_on model variable must be either 0 for no wall or 1 for wall, please see the notes on the wall_on model variable in README')
			em_flag = 2
	
	if (os.path.isdir(output_by_sim) == True and em_flag < 2): # in case proposed results folder already proposed
		err_mess = str('Error - results folder (' +output_by_sim+ ') already exists, please use an alternative.  This can be changed by the res_file_name variable in the model variables file, as explained in README.')
		em_flag = 2
	
	if (output_by_sim in self.output_list and em_flag < 2 and self.chck_num == 1): # in case proposed results folder already on computer hard drive
		err_mess = str('Error - the proposed path for the results folder (' +output_by_sim+ ') has been taken by another simulation in the batch, please use an alternative.  This can be changed by the res_file_name variable in the model variables file, as explained in README.')
		em_flag = 2
	
	# let user know that results will be automatically deleted when 
	# default name used for folder to save to
	if (sav_nam == 'default_res_name' and em_flag < 2):
		err_mess_n = str('Note - default name for save folder used, therefore results folder will be automatically deleted at the end of the simulation to avoid future duplication\n')
		if (em_flag == 0):
			err_mess = err_mess_n
		else:
			err_mess = str(err_mess+err_mess_n)
		em_flag = 1	
	
	# check on particle number concentration inputs ----------------------------------
	
	# ensure size structure marker is sensible
	if (siz_stru < 0 or siz_stru > 1 and em_flag < 2):
		err_mess = str('Error - the size structure must be either 0 for moving-centre or 1 for full-moving, see the notes on the size_structure model variable in README')
		em_flag = 2
	
	# consistency between number of particle size bins and particle number concentration
	if (num_sb == 0 and (sum(sum(pconc != 0)) > 0) and em_flag < 2):
		err_mess = str('Error - zero particle size bins registered (number_size_bins in model variables input file), however total particle number concentration (pconc in model variables input file) contains a non-zero value, therefore please reconcile')
		em_flag = 2 # error message flag for error
	
	# consistency between length of total particle concentration, 
	# mean particle radius and standard deviation when in modal mode
	if (num_sb > 0 and em_flag < 2):
		if (pmode == 0): # particle number concentrations expressed as modes
			if (mean_rad.shape != pconc.shape and em_flag < 2):
				err_mess = str('Error - particle number concentration of seed particles has been detected in modal form (as colons separate values), however the length of the mean radius per mode does not match the length of the particle number concentration per mode, and it must, please see the pconc and mean_rad model variables in README.')
				em_flag = 2 # error message flag for error
			if (std.shape != pconc.shape and em_flag < 2):
				err_mess = str('Error - particle number concentration of seed particles has been detected in modal form (as colons separate values), however the length of the standard deviation per mode does not match the length of the particle number concentration per mode, and it must, please see the pconc and std model variables in README.')
				em_flag = 2 # error message flag for error
	
	# ensure that the mean radius of modes (mean_rad) is consistent with the particle size range
	if (num_sb > 0 and em_flag < 2):
		
		if (pmode == 0): # particle number concentrations expressed as modes
			if (sum(sum(mean_rad < lowsize)) > 0 and sum(sum(mean_rad == -1.e6)) == 0):
				err_mess = str('Error - A value for the mean radius (um) (determined by the mean_rad model variable) of particles is smaller than the particle size range (determined by the lower_part_size model variable).  The mean radius must be within bounds set by the lower_part_size and upper_part_size model variables.  Please see README for more guidance.')
				em_flag = 2 # error message flag for error
			if (sum(sum(mean_rad > uppsize)) > 0):
				err_mess = str('Error - A value for the mean radius (um) (determined by the mean_rad model variable) of particles is greater than the particle size range (determined by the upper_part_size model variable).  The mean radius must be within bounds set by the lower_part_size and upper_part_size model variables.  Please see README for more guidance.')
				em_flag = 2 # error message flag for error

		# if particle concentration per size bin supplied explicitly
		# check that number of time consisent across relevant variables
		if (pmode == 1):
			if (pconc.shape[1] != pconct.shape[1] or pconc.shape[1] != pcont.shape[1] or pconct.shape[1] != pcont.shape[1]):
				err_mess = str('Error: inconsistent number of times for injection of particles as represented by the following model variable inputs (number of times represented in brackets) for: pconc ('+str(pconc.shape[1])+'), pconct ('+str(pconct.shape[1])+') and/or pcont ('+str(pcont.shape[1])+').  Please see README for guidance.')
				em_flag = 2 # error message flag for error
			# in this mode ensure fillers for mean_rad and std match the same
			# length of first dimension as times
			mean_rad = np.ones((1, pconct.shape[1]))*-1.e6
			std = np.ones((1, pconct.shape[1]))*1.e20

		# if particle concentration per size bin supplied by modes
		# check that number of time consistent across relevant variables
		if (pmode == 0):
			if (pconct.shape[1] != std.shape[1] or pconct.shape[1] != mean_rad.shape[1] or pconct.shape[1] != pcont.shape[1] or pconct.shape[1] != pconc.shape[1]):
				err_mess = str('Error: inconsistent number of times for injection of particles as represented by the following model variable inputs (number of times represented in brackets) for: pconc ('+str(pconc.shape[1])+'), pconct ('+str(pconct.shape[1])+'), pcont ('+str(pcont.shape[1])+'), mean_rad ('+str(mean_rad.shape[1])+'), std ('+str(std.shape[1])+').  Please see README for guidance.')
				em_flag = 2 # error message flag for error

	# ensure only one particle concentration given at start of experiment
	if (sum(sum(pconct == 0.)) > 1):
		if (em_flag < 2):
			err_mess = str('Error: only one initial (pconct = 0.0 s) number size distribution is allowed, but for the input pconct variable, time = 0.0 s has been detected more than once.  If you wish to have injection of particles after experiment start but close to time = 0 s please use a pconct value that is greater than 0.0 s but small compared to the recording time step.')
			em_flag = 2 # error message flag for error
	
	# ensure std greater than 1.0
	if (em_flag < 2 and sum(sum(std <= 1.0)) > 0):
		err_mess = str('Error: a standard deviation for particle number size distribution modes (std model variable) less than or equal to 1.0 has been detected, but values must exceed 1.0.  See sdt in the model variables section of README for more details.')
		em_flag = 2 # error message flag for error

	# check on temperature inputs ----------------------------------------------
	
	if (em_flag<2):
		if (len(self.TEMP) != len(self.tempt)):
			err_mess = str('Error - length of temperature inputs does not match length of corresponding times for temperatures, and the two must match, please see the README notes on the temperature and tempt model variables for further guidance')
			em_flag = 2 # error message flag for a note
		if (len(self.TEMP) > 0 and self.tempt[0] != 0):
			err_mess = str('Error - no temperature given for the start of the experiment (0 seconds), and one must be provided, please see the README notes on the temperature and tempt model variables for further guidance')
			em_flag = 2 # error message flag for a note
		
	if (em_flag<2):
		# check on rate of temperature change 
		if (len(self.TEMP) == len(self.tempt) and len(self.TEMP)>1):
			# maximum change
			max_change = np.max(np.abs(np.diff(self.TEMP)))
			if (max_change > 1.):
				err_mess_n = str('Note - the maximum change of chamber temperature exceeds 1 K.  This may destabilise the solver for gas-particle partitioning of water.  It is recommended that changes to temperature given in the temperature model variable are limited to 1 K, with their corresponding times given in the tempt model variable (PyCHAM will automatically limit the update step according to the time differences given by the tempt model variable).  If the solver does become unstable, PyCHAM will assume this results from a change in chamber condition, such as temperature change, and will attempt to rectify the instability by attempting incrementally smaller temperature changes over incrementally smaller time steps.\n')
				if (em_flag == 0):
					err_mess = err_mess_n
				else:
					err_mess = str(err_mess+err_mess_n)
				em_flag = 1 # error message flag for a note

	# -------------------------------------------------------------------------------------

	# check on relative humidity inputs -----------------------------------------
	if (any(RH > 1.) and em_flag < 2):
		err_mess_n = str('Note - relative humidity set above 1.; simulation will be attempted, but please note that RH is interpreted as a fraction, not a percentage, where 1 represents saturation of water vapour\n')
		if (em_flag == 0):
			err_mess = err_mess_n
		else:
			err_mess = str(err_mess+err_mess_n)
		em_flag = 1 # error message flag for a note
	
	# check on big jumps in relative humidity
	if (em_flag<2):
		if (len(RH) > 1):
			if (np.max(np.abs(np.diff(RH))) > 0.01):
				err_mess_n = str('Note - the maximum change of chamber relative humidity exceeds 0.01.  This may destabilise the solver for gas-particle partitioning of water.  It is recommended that changes to relative humidity given in the rh model variable are limited to 0.01, with their corresponding times given in the rht model variable (PyCHAM will automatically limit the update step according to the time differences given by the rht model variable).  If the solver does become unstable, PyCHAM will assume this results from a change in chamber condition, such as relative humidity change, and will attempt to rectify the instability by attempting incrementally smaller relative humidity changes over incrementally smaller time steps.\n')
				if (em_flag == 0):
					err_mess = err_mess_n
				else:
					err_mess = str(err_mess+err_mess_n)
				em_flag = 1			
	
	if (len(RH) != len(RHt) and em_flag < 2):
		err_mess = str('Error - the number of relative humidities does not match the number of times through simulation at which relative humidities are reached, please refer to the notes on the rh and rht model variables in README.')
		em_flag = 2 # error message flag for error
	
	if (RHt[0] != 0 and em_flag < 2):
		err_mess = str('Error - the first time (seconds) through simulation at which a relative humidity time is given is not zero, which represents the start of the experiment, but the model requires the relative humidity at the start of the experiment, please refer to the notes on the rh and rht model variables in README.')
		em_flag = 2 # error message flag for error
	
	# check on UManSysProp ---------------------------------------------------

	# note that whilst conditionally setting uman_up below is required for this section to work, it does not
	# change the input to the model, rather this is done, if needed (and it's needed when umansysprop 
	# needs downloading even when the user hasn't specified this), by mod_var_read

	# for UManSysProp if no update requested, check that there 
	# is an existing UManSysProp folder
	if (uman_up == 0): # check for existing umansysprop folder
		if not os.path.isdir(os.getcwd() + '/umansysprop'): # if no existing folder then force update
			uman_up = 1
	
	if (uman_up == 1): # test whether UManSysProp can be updated
		import urllib.request # module for checking internet connection
		# function for testing internet connection
		def connect(host='https://github.com/loftytopping/UManSysProp_public.git'):
			try:
				urllib.request.urlopen(host) #Python 3.x
				return True
			except:
				return False
		# test internet connection
		if (connect() == False and em_flag < 2):
			err_mess = str('Error - either user has requested cloning of UManSysProp via the model variables input file or no existing UManSysProp folder has been found, but connection to the page failed, possibly due to no internet connection (UManSysProp repository site: https://github.com/loftytopping/UManSysProp_public.git)')
			em_flag = 2
	# ----------------------------------------------------------------------------

	
	
	if (len(self.light_stat) != len(self.light_time) and em_flag < 2):
		err_mess = str('Error - length of input model variables light_status and light_time have different lengths and they should be the same, please see README for guidance.')
		em_flag = 2
	
	# if no status provided at simulation start default to lights off, request one
	if (self.light_time[0] != 0 and em_flag < 2):		
		err_mess = str('Error - light status times have been supplied but not at experiment start (0 in the light_time).  Please see notes in README on the light_time and light_status model variables.')
		em_flag = 2
		
	# check consistency between number of times components instantaneously injected and
	# the number of concentrations given
	if (len(injectt) != Ct.shape[1] and em_flag < 2):
		err_mess = str('Error - number of times given for the instantaneous injection of gas-phase components (injectt in the model variables input file) is inconsistent with the number of concentrations of instantaneous injection provided (Ct in the model variables input file), please see the README for guidance')
		em_flag = 2
		
	# check on consistency of manually assigned densities
	if (len(dens_comp) != len(dens) and em_flag < 2):
		err_mess = str('Error: the number of chemical scheme names provided in the dens_Comp model variables input does not match the number of densities provided in the dens model variables input, please see README for details')
		em_flag = 2
	
	# check on consistency of names of seed components, number of size bins and the mole fraction of dry seed components, note components in rows and size bins in columns
	if (len(seed_name) != seedx.shape[0] and em_flag < 2):
		
		# if seedx is equal to the default value, then assume equal mole fractions of dry seed components
		if ((len(seed_name) > 1) and (seedx.shape[0] == 1) and (seedx[0, 0] == 1)):
			seedx = np.ones((len(seed_name), 1))
			seedx[:, 0] = 1./len(seed_name)
		else:
			err_mess = str('Error - the number of seed particle component names (seed_name in model variables input file) is inconsistent with the number of mole fractions of dry (excluding water) seed particle components (seedx in model variables input file), please see README for guidance')
			em_flag = 2
	
	# check on consistency of dry seed component mole fractions and number of size bins
	# note that mole fraction size bins vary with columns
	if (seedx.shape[1] > 1): # if size bins stated explicitly
		if (seedx.shape[1] != num_sb): # if inconsistent
			err_mess = str('Error - the number of size bins for which the dry (excluding water) seed particle component mole fractions (given by seedx in the model variables input file) is inconsistent with the number of size bins (number_size_bins in model variables input file), please see README for guidance.')
			em_flag = 2

	# check consistency between names of initial components and their concentrations ------------------
	if (len(comp0) != len(y0) and em_flag < 2):
		err_mess = str('Error - the number of gas-phase components present at simulation start (the Comp0 model variable) is different to the number of initial concentrations of these components (the C0 model variable), and they must be the same length.  Please see README for guidance')
		em_flag = 2	

	# check on consistency of names of seed component(s) and their dissociation constant --------------
	if (len(seed_name) != len(seed_diss) and em_flag < 2):
		if (len(seed_name) > 1) and (len(seed_diss) == 1):
			seed_diss = np.ones((len(seed_name)))
		else:
			err_mess = str('Error - the number of seed particle component names (seed_name in model variables input file) and the number of seed particle component dissociation constants (seed_diss in model variables input file) are inconsistent, please see README for guidance')
			em_flag = 2

	if (len(partit_cutoff) > 1 and em_flag < 2):
		err_mess = str('Error - length of the model variables input partit_cutoff should have a maximum length of one, but is greater, please see README for guidance')
		em_flag = 2
	
	# check on continuous injection of components inputs -----------------------------------------------
	if (em_flag < 2): # if no other errors
		try: # note, won't work on the empty con_infl_C default variable
			if (self.con_infl_C.shape[0] != len(self.con_infl_nam)):
				err_mess = str('Error - input for influx rate (ppb/s) of components with continuous influx does not correspond to the number of component names for continuous influx, note that for the continuous influx rate a semicolon should separate values for different components and commas should separate different times.  Please see README for more guidance.')
				em_flag = 2
		except:
			em_flag = em_flag
				
		try: # note, won't work on the empty con_infl_C default variable
			if (self.con_infl_C.shape[1] != len(con_infl_t)):
				err_mess = str('Error - input for influx rate (ppb/s) of components with continuous influx does not correspond to the number of times for continuous influx, note that for the continuous influx rate a comma should separate values for different times and semicolons should separate different components.  Please see README for more guidance.')
				em_flag = 2
		except:
			em_flag = em_flag

	# ----------------------------------------------------------
	# check on presence of actinic flux file for photolysis - note this has to 
	# be before chemical scheme check to stop that check crashing when actinic flux file has a problem
	if hasattr(self, 'af_path') and self.af_path != []  and self.af_path != 'no': # if file provided
		try: # try opening as in lamp_photo module
			f = open(self.af_path, 'r') # open file
			f.close() # close excel file
		except:
			err_mess = str('Error- actinic flux file ' + self.af_path + ' could not be found, please check file and/or the act_flux_file model variable in the model variable file.')
			em_flag = 2
	
	

	
	# chemical scheme check-------------------------------------
	if (em_flag < 2 and self.sch_name != 'Not found'): # if a chemical scheme has been identified
		# get chemical scheme names and SMILE strings of components present in chemical scheme file
		[comp_namelist, _, err_mess_new, 
		H2Oi] = chem_sch_SMILES.chem_scheme_SMILES_extr(self)
		
		if (err_mess_new == '' and em_flag < 2):
			
			# check for water presence in components instantaneously injected, note need to call function above to get H2Oi
			if (comp_namelist[H2Oi] in Compt):
				err_mess = str('Error - water is included in the components being instantaneously injected via the model variables Ct, Compt and injectt, however changes to water vapour content should be made using the rh and rht model variables.  Please see README for further guidance.')
				em_flag = 2
		
		# in case of note message
		if (err_mess_new[0:4] == 'Note' and em_flag < 2):
			err_mess = err_mess_new
			em_flag = 1

		# in case error message produced by function checking the chemical scheme
		if (err_mess_new[0:5] == 'Error' and em_flag < 2):
			err_mess = err_mess_new
			em_flag = 2
	
	if (em_flag < 2 and self.sch_name == 'Not found'): # if a chemical scheme has not been identified
	
		err_mess_n = str('Note - a new chemical scheme has not been identified automatically.  For PyCHAM to identify automatically the chemical scheme file must be inside the selected folder and its file name must contain the string \'chem\'.  Alternatively use the relevant button to select the chemical scheme individually.\n')
		if (em_flag == 0):
			err_mess = err_mess_n
		else:
			err_mess = str(err_mess+err_mess_n)
		em_flag = 1
	
	# -------------------------------------
	# check on water history
	if (em_flag < 2):
		if (wat_hist != 0 and wat_hist != 1):
			err_mess = str('Error - flag for history of particle-phase with respect to water (H2O_hist model variable) should be 0 or 1, please see notes on H2O_hist in README')
			em_flag = 2
	# ------------------------------------
	
	# -------------------------------------
	# check on water-partitioning hysteresis curves
	if (em_flag < 2):
		if (drh_str == -1):
			err_mess = str('Error - the expression for deliquescence relative humidity dependence on temperature could not be converted to a string, please see the README notes on the drh_ft model variable')
			em_flag = 2
		if (erh_str == -1):
			err_mess = str('Error - the expression for efflorescence relative humidity dependence on temperature could not be converted to a string, please see the README notes on the erh_ft model variable')
			em_flag = 2
	if (em_flag < 2):
		# write the module for hysteresis
		importlib.reload(write_hyst_eq)
		write_hyst_eq.write_hyst_eq(drh_str, erh_str)
		# test the module for hysteresis works
		importlib.reload(hyst_eq)
		try:
			drh = hyst_eq.drh(298.15)
		except:
			err_mess = str('Error - the expression for deliquescence relative humidity dependence on temperature was unsuccessfully transferred to a module, please see the README notes on the drh_ft model variable')
			em_flag = 2
		if (em_flag < 2):
			try:
				erh = hyst_eq.erh(298.15)
			except:
				err_mess = str('Error - the expression for efflorescence relative humidity dependence on temperature was unsuccessfully transferred to a module, please see the README notes on the erh_ft model variable')
				em_flag = 2
	# ------------------------------------
	# check on manually assigned vapour pressure of seed component

	if (seed_name.count('core') < len(seed_name)):
		for seedi in seed_name: # loop through seed components
			# check whether this seed component is not core 
			# and it has not been allocated a vapour pressure manually
			if (seedi != 'core' and vol_comp.count(seedi) == 0 and em_flag < 2):
				# add component to manually assigned vapour pressures
				vol_comp.append(seedi)
				# assume low volatility for this seed component (Pa)
				volP.append(1.e-20)
				err_mess = str('Note - a seed component was identified without a manually specified vapour pressure.  This risks evaporation of seed particles and instability for the ODE solver, therefore a vapour pressure of 1.e-20 Pa has been assumed.  To change this setting, specify the vapour pressure of the seed components using the vol_Comp and volP model variables.  Please see README for more guidance.')
				em_flag = 1
		
	# ------------------------------------
	# check on presence of observation file for tracking specific components
	if hasattr(self, 'obs_file') and self.obs_file != []: # if observation file provided for constraint
		try: # try opening as in eqn_pars module
			import openpyxl # require module
			self.obs_file = str(os.getcwd() + self.obs_file) # path to file
			wb = openpyxl.load_workbook(filename = self.obs_file) # open file
			wb.close() # close excel file
		except:
			err_mess = str('Error- observation file ' + self.obs_file + ' could not be found, please check observation file and/or the obs_file model variable in the model variable file.')
			em_flag = 2
	# ------------------------------------
	

	# store in pickle file
	list_vars = [sav_nam, update_stp, 
			tot_time, comp0, y0, RH, RHt, Press, 
			Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, 
			uppsize, space_mode, std, mean_rad, save_step, 
			Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, 
			seedx, con_infl_t, dens_comp, dens, vol_comp, 
			volP, act_comp, act_user, accom_comp, accom_val, uman_up, 	
			int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, 
			coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, 
			Rader, p_char, e_field, partit_cutoff, ser_H2O, 
			wat_hist, drh_str, erh_str, pcont, Vwat_inc, seed_eq_wat, 
			z_prt_coeff, chamV]

	input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
	with open(input_by_sim, 'wb') as pk: # the file to be used for pickling
		pickle.dump(list_vars, pk) # pickle
		pk.close() # close
	
	# update model variables message in GUI
	if (em_flag < 2): # if no error message
	
		self.l80.setText(str('Setup Status: \n' + err_mess))
		self.l80.setStyleSheet(0., '0px', 0., 0.) # remove any borders
		self.bd_st = 3 # change border status to ready for change
		
		# remember path to output, in it needs adding to list for batch
		self.output_by_sim = output_by_sim
		
		# if not (yet) in batch mode
		if (self.btch_no == 1):
		
			self.fab = 1 # remember that single simulation widgets are showing
			self.atb = 1 # remember that add to batch button showing
			
			# allow single simulation to start
			self.b81 = QPushButton('Start Single Simulation', self)
			self.b81.setToolTip('Button to run the simulation')
			self.b81.clicked.connect(self.on_click81sing)
			self.NSlayout.addWidget(self.b81, 5, self.mvpn, 1, 1)
		
			# add 'or' label
			self.l81 = QLabel(self)
			self.l81.setText('or')
			self.NSlayout.addWidget(self.l81, 5, self.mvpn+1)
		
			# allow adding to batch list
			self.b82 = QPushButton('Add To Batch', self)
			self.b82.setToolTip('Add this simulation setup to batch')
			self.b82.clicked.connect(self.on_click82)
			self.NSlayout.addWidget(self.b82, 5, self.mvpn+2)
			
		
		# if in batch mode
		else:
			self.atb = 1 # remember that add to batch button is showing
			# allow adding to batch list
			self.b82 = QPushButton('Add To Batch', self)
			self.b82.setToolTip('Add this simulation setup to batch')
			self.b82.clicked.connect(self.on_click82)
			self.NSlayout.addWidget(self.b82, 5, self.mvpn+2)
			
			# add 'or' label
			self.l81 = QLabel(self)
			self.l81.setText('or')
			self.NSlayout.addWidget(self.l81, 5, self.mvpn+1)
		
	else: # if there is a message
		
		# update error message
		self.l80.setText(str('Setup Status: \n' + err_mess))
		# change border accordingly
		if (self.bd_st == 1):
			self.l80.setStyleSheet(0., '2px dashed red', 0., 0.)
		if (self.bd_st >= 2):
			self.l80.setStyleSheet(0., '2px solid red', 0., 0.)

		self.bd_st += 2 # prepare for change to border status
		# change border status
		if (self.bd_st == 3):
			self.bd_st = 2
		if (self.bd_st >= 4):
			self.bd_st = 1
		
		if (self.fab == 1): # if showing, remove single simulation widgets
			self.b81.deleteLater()
			self.l81.deleteLater()
			self.fab = 0 # remember their removal
		if (self.atb == 1): # if showing remove add to batch button
			self.b82.deleteLater()
			self.atb = 0 # remember that add to batch button not showing
		
		
	import mod_var_up # update displayed model variables in case checking has modified any
	mod_var_up.mod_var_up(self)

	return()