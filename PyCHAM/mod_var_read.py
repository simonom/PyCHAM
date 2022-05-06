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
'''module to read and store (via pickle file and self parameter) model variables from file'''

import pickle
import numpy as np
import os

def mod_var_read(self):
	
	# inputs: ------------------------------------------
	# self - reference to PyCHAM
	# --------------------------------------------------
	
	def read():
	
		# prepare by opening existing model variables, ready for modification
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
		
		if (self.inname != 'Default' and self.inname != 'Not found'): # if not using defaults
			inputs = open(self.inname, mode= 'r' ) # open model variables file
			in_list = inputs.readlines() # read file and store everything into a list
			inputs.close() # close file
		else: # if using defaults
			in_list = []
		
		err_mess = '' # initial (blank) error message
		self.bd_st = 3 # change border/error message status to ready for change

		for i in range(len(in_list)): # loop through supplied model variables to interpret
			
			# ----------------------------------------------------
			# if commented out continue to next line
			if (in_list[i][0] == '#'):
				continue
			
			key, value = in_list[i].split('=') # split values from keys
			# model variable name - a string with bounding white space removed
			key = key.strip()
			# ----------------------------------------------------
			if key == 'res_file_name' and (value.strip()): # name of folder to save results in
				sav_nam = str(value.strip())
				
			if key == 'chem_scheme_markers' and (value.strip()): # formatting for chemical scheme
				self.chem_sch_mrk = [str(i).strip() for i in (value.split(','))]

			if key == 'update_step' and (value.strip()): # time step (s) for updating ODE initial conditions
				update_stp = float(value.strip())

			if key == 'total_model_time' and (value.strip()):
				try:
					tot_time = float(value.strip())
				except:
					err_mess = 'Could not convert string to float for total_model_time model variable, please check model variables file and see README for guidance'

			if key == 'Comp0' and (value.strip()): # names of components present at experiment start
				comp0 = [str(i).strip() for i in (value.split(','))]			

			if key == 'C0' and (value.strip()): # initial concentrations of components present at experiment start (ppb)
				try:
					y0 = [float(i) for i in (value.split(','))]
				except:
					err_mess = 'Could not read in the C0 model variable, please check the model variables file and see README for guidance'
			
			if key == 'temperature' and (value.strip()): # chamber temperature (K)
				self.TEMP = [float(i) for i in ((value.strip()).split(','))]

			if key == 'tempt' and (value.strip()): # times (s) that temperature values correspond to
				self.tempt = [float(i) for i in ((value.strip()).split(','))]

			if key == 'rh' and (value.strip()): # relative humidity in chamber (0-1)
				RH = np.array(([float(i) for i in ((value.strip()).split(','))]))
				
			if key == 'rht' and (value.strip()): # times through simulation (s) at which relative humidity reached
				RHt = np.array(([float(i) for i in ((value.strip()).split(','))]))

			if key == 'p_init' and (value.strip()): # pressure inside chamber
				Press = float(value.strip())

			if key == 'daytime_start' and (value.strip()): # time of day at experiment start (s)
				self.daytime = float(value.strip())

			if key == 'wall_on' and (value.strip()): # marker for whether or not to consider wall
				self.wall_on = int(value.strip())

			if key == 'eff_abs_wall_massC' and (value.strip()): # effective absorbing mass concentration of wall
				Cw = float(value.strip())

			if key == 'mass_trans_coeff' and (value.strip()): # mass transfer coefficient of vapours with wall
				kw_in = float(value.strip())
				try:
					kw = np.zeros((len(kw_in)))
					kw[:] = kw_in
				except:
					kw = np.zeros((1))
					kw[0] = kw_in

			if key == 'chamSA' and (value.strip()): # chamber surface area (m2)
				chamSA = float(value.strip())

			if key == 'chamV' and (value.strip()): # chamber volume (m3)
				chamV = float(value.strip())

			if key == 'size_structure' and (value.strip()): # the size structure
				siz_stru = int(value.strip())

			if key == 'number_size_bins' and (value.strip()): # number of particle size bins
				num_sb = int(value.strip())

			if key == 'pconc' and (value.strip()): # seed particle number concentrations (#/cc)
				time_cnt = 1 # track number of times
				sb_cnt = 1 # track number of size bins
				mode_cnt = 1 # track number of modes
			
				for i in value:
					if i == ';': # semi-colon represents a time difference
						time_cnt += 1 # increase time count
					if (time_cnt == 1 and i == ','):
						sb_cnt += 1 # increase size bin count
						pmode = 1 # explicitly stated particle concentrations
					if (time_cnt == 1 and i == ':'):
						mode_cnt += 1 # mode count
						pmode = 0 # particle concentrations expressed as modes
				# a possible situation where there is just one size bin and 
				# therefore only one concentration given for that size bin
				if (sb_cnt ==1 and mode_cnt ==1):
					 pmode = 1 # explicitly stated particle concentrations

				# if number concentration per size bin given explicitly
				if (sb_cnt > 1):
					pconc = np.zeros((sb_cnt, time_cnt))
					for i in range(time_cnt):
						pconc[:, i] = [float(ii.strip()) for ii in ((value.split(';')[i]).split(','))]
				else: # mode quantities provided
					pconc = np.zeros((mode_cnt, time_cnt))
					for i in range(time_cnt):
						pconc[:, i] = [float(ii.strip()) for ii in ((value.split(';')[i]).split(':'))]
				
			
			if key == 'pconct' and (value.strip()): # seed particle input times (s)
				time_cnt = 1 # track number of times
				for i in value:
					if (i == ';'):
						time_cnt += 1 # increase time count
					
				# times in columns
				pconct = np.zeros((1, time_cnt))
				pconct[0, :] = [float(i) for i in ((value.strip()).split(';'))]
			
			# whether seed particle injection continuous or instantaneous
			if key == 'pcont' and (value.strip()):
				time_cnt = 1 # track number of times
				for i in value:
					if (i == ';'):
						time_cnt += 1 # increase time count
					
				# times in columns
				pcont = np.zeros((1, time_cnt), dtype=int)
				pcont[0, :] = [int(i) for i in ((value.strip()).split(';'))]
			

			if key == 'lower_part_size' and (value.strip()): # lowest size bin bound
				lowsize = float(value.strip())
			
			if key == 'upper_part_size' and (value.strip()): # radius of uppermost size bin boundary
				uppsize = float(value.strip())

			if key == 'space_mode' and (value.strip()): # method of spacing size bin sizes
				space_mode = str(value.strip())

			if key == 'std' and (value.strip()): # seed particle standard deviation of number size distribution
				time_cnt = 1 # track of number of times
				mode_cnt =1 # track number of modes 
				for i in value:
					if (i==';'):
						time_cnt += 1 # increase time count
						# reset mode count
						mode_cnt = 1 # track number of modes

					if (i==':'):
						mode_cnt += 1 # increase mode count

				std = np.zeros((mode_cnt, time_cnt))
				for i in range(time_cnt):
					std[:, i] = [float(ii.strip()) for ii in ((value.split(';')[i]).split(':'))]
			
			if key == 'mean_rad' and (value.strip()): # seed particle mean radius (um)
				time_cnt = 1 # track of number of times
				mode_cnt = 1 # track number of modes
				for i in value:
					if (i==';'):
						time_cnt += 1 # increase time count
						# reset mode count
						mode_cnt = 1 # track number of modes
					if (i==':'):
						mode_cnt += 1 # increase mode count
				
				mean_rad = np.zeros((mode_cnt, time_cnt))
				for i in range(time_cnt):
					mean_rad[:, i] = [float(ii.strip()) for ii in ((value.split(';')[i]).split(':'))]
				
			if key == 'recording_time_step' and (value.strip()): # frequency (s) of storing results
				save_step = float(value.strip())

			if key == 'const_comp' and (value.strip()): # names of components with later continuous injections
				self.const_comp = [str(i).strip() for i in (value.split(','))]

			if key == 'obs_file' and (value.strip()): # name of file containing observations to constrain by
				self.obs_file =str(value.strip())
			
			# names of components with instantaneous gas-phase injections
			if key == 'Compt' and (value.strip()):
				Compt = [str(i).strip() for i in (value.split(','))]

			# times of later gas-phase instantaneous injections (s)
			if key == 'injectt' and (value.strip()):
				injectt = [float(i.strip()) for i in (value.split(','))]
				injectt = np.array((injectt))

			# concentration of components with instantanetous injection (ppb)
			if key == 'Ct' and (value.strip()):
				comp_count = 1 # count number of components
				time_count = 1 # track number of times
				for i in value:
					if i==';':
						comp_count += 1 # record number of components
						time_count = 1 # reset time count
					if i==',':
						time_count += 1
				Ct = np.zeros((comp_count, time_count))
				for i in range(comp_count):
					Ct[i, :] = [float(ii.strip()) for ii in ((value.split(';')[i]).split(','))]


			if key == 'seed_name' and value.strip(): # name(s) of component(s) comprising seed particles
				seed_name = [str(i).strip() for i in (value.split(','))]
				
			if key == 'seed_mw' and value.strip():
				try:
					seed_mw =  [float(i) for i in (value.split(','))]
					seed_mw =  np.array((seed_mw))
				except:
					seed_mw = ['fail']
		
			if key == 'seed_diss' and value.strip(): # dissociation constant of seed material
				seed_diss = [float(i) for i in (value.split(','))]

			if key == 'seed_dens' and value.strip():
				seed_dens = [float(i) for i in (value.split(','))]
				seed_dens = np.array((seed_dens))

			if key == 'seedx' and value.strip(): # mole fraction of components in dry seed particles
				comp_count = 1 # count number of components
				sb_count = 1 # track number of size bins
				for i in value:
					if i==';':
						comp_count += 1 # record number of components
						sb_count = 1 # reset size bin count
					if i==',':
						sb_count += 1 # size bin count
				seedx = np.zeros((comp_count, sb_count))
				for i in range(comp_count):
					seedx[i, :] = [float(ii.strip()) for ii in ((value.split(';')[i]).split(','))]
			
			if key == 'Vwat_inc' and value.strip(): # whether number size distribution includes volume of water
				Vwat_inc = int(value.strip())
			
			# whether to allow water to equilibrate with seed prior to experiment
			if key == 'seed_eq_wat' and value.strip():
				seed_eq_wat = int(value.strip())

			# fraction below which gas-particle partitioning coefficient treated as zero,
			# e.g. because size bin has relatively very small surface area
			if key == 'z_prt_coeff'  and value.strip():
				z_prt_coeff = float(value.strip())			

			if key == 'light_status' and value.strip(): # status of lights (on or off)
				light_stat = [int(i) for i in (value.split(','))]
				self.light_stat = np.array((light_stat))
				
			if key == 'light_time' and value.strip(): # times (s) corresponding to light status
				light_time = [float(i) for i in (value.split(','))]
				self.light_time = np.array((light_time))
				
			if key == 'lat': # latitude (degrees)
				if (value.strip()): self.lat = float(value.strip())

			if key == 'lon': # longitude (degrees)
				if (value.strip()): self.lon = float(value.strip())

			if key == 'act_flux_file' and (value.strip()): # for indoor actinic flux
				self.af_path = str(os.getcwd() + '/PyCHAM/photofiles/' + value.strip())

			if key == 'DayOfYear' and (value.strip()):
				self.dayOfYear = int(value.strip())		
			
			# name of file with wavelength-dependent absorption cross-section and quantum yield 
			# calculations
			if key == 'photo_par_file' and (value.strip()):
				self.photo_path = str(os.getcwd() + '/PyCHAM/photofiles/' + value.strip())
			
			if key == 'trans_fac' and (value.strip()): # transmission factor for natural light
				self.tf = float(value.strip())
				
			if key == 'tf_UVC' and (value.strip()): # transmission factors for 254 nm light
				self.tf_UVC = [float(i.strip()) for i in (value.split(','))]
				
			if key == 'tf_UVCt' and (value.strip()): # transmission factor times for 254 nm light
				self.tf_UVCt = np.array(([float(i.strip()) for i in (value.split(','))]))

			if key == 'light_adapt' and (value.strip()): # whether to adapt time step to naturally varying light intensity
				self.light_ad = int(value.strip())

			if key == 'const_infl' and (value.strip()): # names of components with continuous influx
				self.con_infl_nam = [str(i).strip() for i in (value.split(','))]

			if key == 'const_infl_t' and (value.strip()): # times of continuous influxes (s)
				con_infl_t = [float(i.strip()) for i in (value.split(','))]
				con_infl_t = np.array((con_infl_t))

			if key == 'Cinfl' and (value.strip()): # influx rate of components with continuous influx (ppb/s)
				comp_count = 1 # count number of components
				time_count = 1 # track number of times
				for i in value:
					if (i==';'):
						comp_count += 1 # record number of components
					if (i==',' and comp_count == 1):
						time_count += 1 # record number of times
				self.con_infl_C = np.zeros((comp_count, time_count))
				
				try:	
					for i in range(comp_count): # loop through components
						for ii in range(time_count): # loop through times
							self.con_infl_C[i, ii] = float((((value.split(';')[i]).split(',')))[ii].strip())
							
				# in case semicolons and commas messed up on input, note this will invoke an 
				# error message from the user input check module
				except:
					self.con_infl_C = np.empty(0)
				
			if key == 'tracked_comp' and (value.strip()): # names of components whose tendency to change will be tracked
				self.dydt_trak = [str(i).strip() for i in (value.split(','))]

			if key == 'dens_Comp' and (value.strip()):
				dens_comp = [str(i).strip() for i in (value.split(','))]

			if key == 'dens' and (value.strip()):
				dens = [float(i) for i in (value.split(','))]

			if key == 'vol_Comp' and (value.strip()):
				vol_comp = [str(i).strip() for i in (value.split(','))]

			if key == 'volP' and (value.strip()):
				volP = [float(i) for i in (value.split(','))]

			if key == 'act_comp' and (value.strip()): # names of components with specified activity coefficients
				act_comp = [i for i in (((value.strip()).split(',')))]

			if key == 'act_user' and (value.strip()): # activity coefficients (dimensionless set by user)
				act_user = [i for i in (((value.strip()).split(',')))]

			if key == 'accom_coeff_comp' and (value.strip()): # names of componenets with specified accommodation coefficients
				accom_comp = [i for i in (((value.strip()).split(',')))]

			if key == 'accom_coeff_user' and (value.strip()): # value(s) of accommodation coefficient set by user
				accom_val = [i for i in (((value.strip()).split(',')))]

			if key == 'partit_cutoff' and (value.strip()): # value(s) of the gas-particle partitioning cutoff
				partit_cutoff = [float(i) for i in (((value.strip()).split(',')))]

			if key == 'umansysprop_update' and (value.strip()): # marker for whether to clone latest version of UManSysProp from web
				uman_up = int(value.strip())

			if key == 'int_tol' and (value.strip()): # tolerances for integration
				int_tol = [float(i) for i in (value.split(','))]

			if key == 'new_partr' and (value.strip()): # radius of newly nucleated particles (cm)
				new_partr = float(value.strip())

			if key == 'nucv1' and (value.strip()): # first parameter in the nucleation equation
				nucv1 = float(value.strip())

			if key == 'nucv2' and (value.strip()): # second nucleation parameter (onset)
				nucv2 = float(value.strip())

			if key == 'nucv3' and (value.strip()): # third nucleation parameter (duration)
				nucv3 = float(value.strip())

			if key == 'nuc_comp' and (value.strip()): # chemical scheme name of nucleating component
				nuc_comp = [str(i).strip() for i in (value.split(','))]

			if key == 'nuc_adapt' and (value.strip()): # marker for whether to adapt time interval to nucleation
				nuc_ad = int(value.strip())

			if key == 'coag_on' and (value.strip()): # marker for whether to model coagulation
				coag_on = int(value.strip())

			if key == 'inflectDp' and (value.strip()): # diameter at which wall deposition of particles inflection occurs
				inflectDp = float(value.strip())

			if key == 'Grad_pre_inflect' and (value.strip()): # gradient of wall deposition of particles before inflection
				pwl_xpre = float(value.strip())

			if key == 'Grad_post_inflect' and (value.strip()): # gradient of particle deposition to wall after inflection
				pwl_xpro = float(value.strip())

			if key == 'Rate_at_inflect' and (value.strip()): # particle deposition to wall rate at inflection
				inflectk = float(value.strip())

			if key == 'ChamSA' and (value.strip()): # chamber surface area (m2) used for particle loss to walls
				ChamSA = float(value.strip())

			if key == 'McMurry_flag' and (value.strip()): # marker for whether to use the McMurry model for particle deposition to wall
				Rader = int(value.strip())

			if key == 'part_charge_num' and (value.strip()): # average number of charges per particle
				p_char = float(value.strip())
			
			if key == 'elec_field' and (value.strip()): # electric field in chamber
				e_field = float(value.strip())
			
			if key == 'dil_fac' and (value.strip()): # dilution factor rate
				self.dil_fac = float(value)
				
			if key == 'ser_H2O' and (value.strip()): # whether to serialise water gas-particle partitioning
				ser_H2O = int(value)
				
			if (key == 'H2O_hist' and (value.strip())): # history of particle-phase with respect to particle-phase water
				try:
					wat_hist = int(value)
				except:
					wat_hist = -1 # will cause error message
			if (key == 'drh_ft' and (value.strip())): # deliquescence RH dependence on temperature
				try:
					drh_str = str(value)
				except:
					drh_str = -1 # will cause error message
			if (key == 'erh_ft' and (value.strip())): # efflorescence RH dependence on temperature
				try:
					erh_str = str(value)
				except:
					erh_str = -1 # will cause error message
		
		
		# UManSysProp check ----------------------------------
		# for UManSysProp if no update requested, check that there 
		# is an existing UManSysProp folder
		if (uman_up == 0): # check for existing umansysprop folder
			if not os.path.isdir(os.getcwd() + '/umansysprop'): # if no existing folder then force update
				uman_up = 1
		# -------------------------------------------
		
		# update model variables message in GUI
		if (err_mess != ''): # if error message occurs
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

		# prepare for pickling
		list_vars = [sav_nam, update_stp, 
				tot_time, comp0, y0, RH, RHt, Press, 
				Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, 
				uppsize, space_mode, std, mean_rad, save_step, 
				Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, 
				seedx, con_infl_t, dens_comp, dens, vol_comp, volP, 
				act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, 
				new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, 
				inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, 
				e_field, partit_cutoff, ser_H2O, wat_hist, drh_str, 
				erh_str, pcont, Vwat_inc, seed_eq_wat, z_prt_coeff, 
				chamV]

		input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
		with open(input_by_sim, 'wb') as pk: # the file to be used for pickling
			pickle.dump(list_vars, pk) # pickle
			pk.close() # close
		
			
	read() # call on function to read the model variables