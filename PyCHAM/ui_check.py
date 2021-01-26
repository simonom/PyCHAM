'''checking that user inputs are valid and providing appropriate modification to the PyCHAM GUI'''
# once model variables are validated, this module also allows options in the GUI for continuing with simulations

import os
import sys
import numpy as np
import pickle
from PyQt5.QtWidgets import *
from PyQt5.QtGui import  *
from PyQt5.QtCore import *

def ui_check(self):

	# inputs: ----------------------------------------------------------------------
	# self - reference to GUI
	# --------------------------------------------------------------------------------

	# path to store for variables
	input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
	with open(input_by_sim, 'rb') as pk:
		[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, 
		tot_time, comp0, y0, temp, tempt, RH, Press, wall_on,
		Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
		save_step, const_comp, Compt, injectt, Ct, seed_name,
		seed_mw, seed_diss, seed_dens, seedVr,
		light_stat, light_time, daytime, lat, lon, af_path, 
		dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, 
		dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
		accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
		nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
		inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O] = pickle.load(pk)
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
	# --------------------------------------------------------------------
	
	# to begin assume no errors, so message is that simulation ready
	err_mess = 'Model variables fine - simulation ready'

	# saving path (copied from saving module) ------------
	dir_path = os.getcwd() # current working directory
	output_root = 'PyCHAM/output'
	filename = os.path.basename(sch_name)
	filename = os.path.splitext(filename)[0]
	# one folder for one simulation
	output_by_sim = os.path.join(dir_path, output_root, filename, sav_nam)
	# -------------------------------------------------------------------
	
	# constrain wall_on marker
	if (wall_on>=1):
		wall_on = 1
	else:
		wall_on = 0
	
	if (os.path.isdir(output_by_sim) == True): # in case proposed results folder already proposed
		err_mess = str('Error - results folder (' +output_by_sim+ ') already exists, please use an alternative.  This can be changed by the res_file_name variable in the model variables file, as explained in README.')
	
	if (output_by_sim in self.output_list): # in case proposed results folder already on computer hard drive
		err_mess = str('Error - the proposed results folder path (' +output_by_sim+ ') has been taken by another simulation in the batch, please use an alternative.  This can be changed by the res_file_name variable in the model variables file, as explained in README.')
	
	# let user know that results will be automatically deleted when 
	#default name used for folder to save to
	if (sav_nam == 'default_res_name'):
		err_mess = str('Note - default name for save folder used, therefore results folder will be automatically deleted at the end of the simulation to avoid future duplication')
			
	# ensure size structure marker is sensible
	if (siz_stru<0 or siz_stru>1):
		siz_stru = 0
	
	# consistency between number of particle size bins and particle number concentration
	if num_sb == 0 and (sum(pconc != 0) > 0):
		pconc[:] = 0.0
		err_mess = str('Note that zero particle size bins registered (number_size_bins in model variables input file), however total particle number concentration (pconc in model variables input file) contains a non-zero value, therefore please reconcile')
	
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

	if (RH > 1.):
		err_mess = str('Note - RH set above 1.; simulation will be attempted, but please note that RH is interpreted as fraction, not a percentage, where 1 represents saturation of water vapour')

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
		if connect() == False:
			err_mess = str('Error - either user has requested cloning of UManSysProp via the model variables input file or no existing UManSysProp folder has been found, but connection to the page failed, possibly due to no internet connection (UManSysProp repository site: https://github.com/loftytopping/UManSysProp_public.git)')

	# ensure that if multiple instantaneous injections of particles, that corresponding variables
	# have the correct shape, specifically that they cover the same number of times
	if (pconc.shape[1] != pconct.shape[1] or pconc.shape[1] != std.shape[1] or pconc.shape[1] != mean_rad.shape[1] or pconct.shape[1] != std.shape[1] or pconct.shape[1] != std.shape[1] or std.shape[1] != mean_rad.shape[1]):
		err_mess = str('Error: inconsistent number of times for instantaneous injection of particles represented by model variable inputs (number of times represented in brackets) for: pconc ('+str(pconc.shape[1])+'), pconct ('+str(pconct.shape[1])+'), mean_rad ('+str(mean_rad.shape[1])+') and/or std ('+str(std.shape[1])+').  Please see README for guidance.')
	
	if (len(light_stat) != len(light_time)):
		err_mess = str('Error - length of input model variables light_status and light_time have different lengths and they should be the same, please see README for guidance.')

	
	if (light_time[0] != 0): # if no status provided at simulation start default to lights off
		light_stat = np.concatenate((np.zeros((1)).astype(int), light_stat))
		light_time = np.concatenate((np.zeros((1)), light_time))

	# check consistency between number of times components instantaneously injected and
	# the number of concentrations given
	if (len(injectt) != Ct.shape[1]):
		err_mess = str('Error - number of times given for the instantaneous injection of gas-phase components (injectt in the model variables input file) is inconsistent with the number of concentrations of instantaneous injection provided (Ct in the model variables input file), please see the README for guidance')

	# check on consistency of manually assigned densities
	if (len(dens_comp) != len(dens)):
		err_mess = str('Error: the number of chemical scheme names provided in the dens_Comp model variables input does not match the number of densities provided in the dens model variables input, please see README for details')

	# check on consistency of names of seed component(s) and their volume ratio
	if (len(seed_name) != len(seedVr)):
		if (len(seed_name) > 1) and (len(seedVr) == 1):
			seedVr = np.ones((len(seed_name)))
		else:
			err_mess = str('Error - the number of seed particle component names (seed_name in model variables input file) is inconsistent with the number of seed particle component volume ratios (seedVr in model variables input file), please see README for guidance')

	# check on consistency of names of seed component(s) and their dissociation constant
	if (len(seed_name) != len(seed_diss)):
		if (len(seed_name) > 1) and (len(seed_diss) == 1):
			seed_diss = np.ones((len(seed_name)))
		else:
			err_mess = str('Error - the number of seed particle component names (seed_name in model variables input file) and the number of seed particle component dissociation constants (seed_diss in model variables input file) are inconsistent, please see README for guidance')

	if (len(partit_cutoff) > 1):
		err_mess = str('Error - length of the model variables input partit_cutoff should have a maximum length of one, but is greater, please see README for guidance')

	# store in pickle file
	list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O]

	input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
	with open(input_by_sim, 'wb') as pk: # the file to be used for pickling
		pickle.dump(list_vars, pk) # pickle
		pk.close() # close


	# update model variables message in GUI
	if (err_mess[0:5] != 'Error'): # if no error message
	
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