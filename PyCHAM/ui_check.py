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
		[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, 
		tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on,
		Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
		save_step, const_comp, Compt, injectt, Ct, seed_name,
		seed_mw, seed_diss, seed_dens, seedVr,
		light_stat, light_time, daytime, lat, lon, af_path, 
		dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, 
		dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
		accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
		nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
		inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, 
		wat_hist, drh_str, erh_str, pcont] = pickle.load(pk)
		pk.close()	

	# loaded variables: ------------------------------------------------------------
	# sav_nam - name of folder to save results to
	# sch_name - name of chemical scheme file
	# wall_on - marker for whether wall on or off
	# siz_stru - the size structure
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
	# wat_hist - flag for particle-phase history with respect to water 
	#	(0 on the deliquescence curve, 1 on the efflorescence curve)
	# con_infl_t - times of continuous influx of components
	# con_infl_nam - chemical scheme names of components with continuous influx
	# con_infl_C - influx rate of components with continuous influx (ppb/s)
	# --------------------------------------------------------------------
	
	# to begin assume no errors, so message is that simulation ready
	err_mess = 'Model variables fine - simulation ready'
	em_flag = 0 # flag that no problematic errors detected

	# saving path (copied from saving module) ------------
	dir_path = os.getcwd() # current working directory
	output_root = 'PyCHAM/output'
	filename = os.path.basename(sch_name)
	filename = os.path.splitext(filename)[0]
	# one folder for one simulation
	output_by_sim = os.path.join(dir_path, output_root, filename, sav_nam)
	# -------------------------------------------------------------------
	
	# incorrect wall_on marker
	if (wall_on != 1 and wall_on != 0):
		if (em_flag < 2):
			err_mess = str('Error - wall_on model variable must be either 0 for no wall or 1 for wall, please see the notes on the wall_on model variable in README')
			em_flag = 2
	
	if (os.path.isdir(output_by_sim) == True and em_flag < 2): # in case proposed results folder already proposed
		err_mess = str('Error - results folder (' +output_by_sim+ ') already exists, please use an alternative.  This can be changed by the res_file_name variable in the model variables file, as explained in README.')
		em_flag = 2
	
	if (output_by_sim in self.output_list and em_flag < 2): # in case proposed results folder already on computer hard drive
		err_mess = str('Error - the proposed results folder path (' +output_by_sim+ ') has been taken by another simulation in the batch, please use an alternative.  This can be changed by the res_file_name variable in the model variables file, as explained in README.')
		em_flag = 2
	
	# let user know that results will be automatically deleted when 
	#default name used for folder to save to
	if (sav_nam == 'default_res_name' and em_flag < 2):
		err_mess_n = str('Note - default name for save folder used, therefore results folder will be automatically deleted at the end of the simulation to avoid future duplication\n')
		if (em_flag == 0):
			err_mess = err_mess_n
		else:
			err_mess = str(err_mess+err_mess_n)
		em_flag = 1	
	
	# check on particle number concentration inputs ----------------------------------
	
	# ensure size structure marker is sensible
	if (siz_stru<0 or siz_stru>1 and em_flag < 2):
		err_mess = str('Error - the size structure must be either 0 for moving-centre or 1 for full-moving, see the notes on the size_structure model variable in README')
		em_flag = 2
	
	# consistency between number of particle size bins and particle number concentration
	if (num_sb == 0 and (sum(pconc != 0) > 0) and em_flag < 2):
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
	
	# ensure that if multiple instantaneous injections of particles, that corresponding variables
	# have the correct shape, specifically that they cover the same number of times
	if (pconc.shape[1] != pconct.shape[1] or pconc.shape[1] != std.shape[1] or pconc.shape[1] != mean_rad.shape[1] or pconct.shape[1] != std.shape[1] or pconct.shape[1] != std.shape[1] or std.shape[1] != mean_rad.shape[1] or pcont.shape[1] != pconct.shape[1]):
		if (em_flag < 2):
			err_mess = str('Error: inconsistent number of times for instantaneous injection of particles represented by model variable inputs (number of times represented in brackets) for: pconc ('+str(pconc.shape[1])+'), pconct ('+str(pconct.shape[1])+'), mean_rad ('+str(mean_rad.shape[1])+'), ('+str(std.shape[1])+') and/or ('+str(pcont.shape[1])+').  Please see README for guidance.')

	# check on temperature inputs ----------------------------------------------
	
	if (em_flag<2):
		if (len(temp) != len(tempt)):
			err_mess = str('Error - length of temperature inputs does not match length of corresponding times for temperatures, and the two must match, please see the README notes on the temperature and tempt model variables for further guidance')
			em_flag = 2 # error message flag for a note
		if (len(temp) > 0 and tempt[0] != 0):
			err_mess = str('Error - no temperature given for the start of the experiment (0 seconds), and one must be provided, please see the README notes on the temperature and tempt model variables for further guidance')
			em_flag = 2 # error message flag for a note
		
	if (em_flag<2):
		# check on rate of temperature change 
		if (len(temp) == len(tempt) and len(temp)>1):
			# maximum change
			max_change = np.max(np.abs(np.diff(temp)))
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

	
	
	if (len(light_stat) != len(light_time) and em_flag < 2):
		err_mess = str('Error - length of input model variables light_status and light_time have different lengths and they should be the same, please see README for guidance.')
		em_flag = 2
	
	# if no status provided at simulation start default to lights off, request one
	if (light_time[0] != 0 and em_flag < 2):		
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
	
	# check on consistency of names of seed component(s) and their volume ratio
	if (len(seed_name) != len(seedVr) and em_flag < 2):
		if (len(seed_name) > 1) and (len(seedVr) == 1):
			seedVr = np.ones((len(seed_name)))
		else:
			err_mess = str('Error - the number of seed particle component names (seed_name in model variables input file) is inconsistent with the number of seed particle component volume ratios (seedVr in model variables input file), please see README for guidance')
			em_flag = 2

	# check on consistency of names of seed component(s) and their dissociation constant
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
			if (con_infl_C.shape[0] != len(con_infl_nam)):
				err_mess = str('Error - input for influx rate (ppb/s) of components with continuous influx does not correspond to the number of component names for continuous influx, note that for the continuous influx rate a semicolon should separate values for different components and commas should separate different times.  Please see README for more guidance.')
				em_flag = 2
		except:
			em_flag = em_flag
				
		try: # note, won't work on the empty con_infl_C default variable
			if (con_infl_C.shape[1] != len(con_infl_t)):
				err_mess = str('Error - input for influx rate (ppb/s) of components with continuous influx does not correspond to the number of times for continuous influx, note that for the continuous influx rate a comma should separate values for different times and semicolons should separate different components.  Please see README for more guidance.')
				em_flag = 2
		except:
			em_flag = em_flag
	# -------------------------------------
	# check on whether water included in instantaneously injected components, as it should not be, instead it should be
	# included in the rh model variable input
	
	if (em_flag < 2):
		# get chemical scheme names and SMILE strings of components present in chemical scheme file
		[comp_namelist, _, err_mess_new, H2Oi] = chem_sch_SMILES.chem_scheme_SMILES_extr(sch_name, xml_name, chem_sch_mark)
		if (err_mess_new == '' and em_flag < 2):
			
			# check for water presence in components instantaneously injected, note need to call function above to get H2Oi
			if (comp_namelist[H2Oi] in Compt):
				err_mess = str('Error - water is included in the components being instantaneously injected via the model variables Ct, Compt and injectt, however changes to water vapour content should be made using the rh and rht model variables.  Please see README for further guidance.')
				em_flag = 2
		else: # in case error message produced by function checking the chemical scheme
			err_mess = err_mess_new
			em_flag = 2
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
	

	# store in pickle file
	list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str, pcont]

	input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
	with open(input_by_sim, 'wb') as pk: # the file to be used for pickling
		pickle.dump(list_vars, pk) # pickle
		pk.close() # close
	
	# update model variables message in GUI
	if (em_flag<2): # if no error message
	
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
		self.l80.setText(str('Setup Status: \n' +err_mess))
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