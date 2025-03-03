########################################################################
#								       #
# Copyright (C) 2018-2025					       #
# Simon O'Meara : simon.omeara@manchester.ac.uk			       #
#								       #
# All Rights Reserved.                                                 #
# This file is part of PyCHAM                                          #
#                                                                      #
# PyCHAM is free software: you can redistribute it and/or modify it    #
# under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or    #
# (at  your option) any later version.                                 #
#                                                                      #
# PyCHAM is distributed in the hope that it will be useful, but        #
# WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# General Public License for more details.                             #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with PyCHAM.  If not, see <http://www.gnu.org/licenses/>.      #
#                                                                      #
########################################################################
'''module to read and store (via pickle file and self parameter) model 
variables from file'''

import pickle
import numpy as np
import os

def mod_var_read(self):
	
	# inputs: ------------------------------------------
	# self - reference to PyCHAM
	# --------------------------------------------------
	
	def read(self):
	
		err_mess = '' # initial (blank) error message

		# prepare by opening existing model variables, ready for modification
		input_by_sim = str(self.PyCHAM_path + '/PyCHAM/pickle.pkl')
		
		with open(input_by_sim, 'rb') as pk:
			[y0_gas, siz_stru, num_sb, 
			lowsize, uppsize, 
			std, Compt, injectt, Ct,
			seed_mw, seed_dens,
			dens_comp, dens, vol_comp, volP, act_comp, 
			act_user, accom_comp, accom_val, uman_up, 
			int_tol, new_partr, coag_on, 
			inflectDp, pwl_xpre, pwl_xpro, 
			inflectk, chamSA, Rader, p_char, e_field, 
			ser_H2O, 
			wat_hist, drh_str, erh_str, z_prt_coeff, 
			chamV] = pickle.load(pk)
		pk.close()
		
		# if not using defaults	
		if (self.inname != 'Default' and 
			self.inname != 'Not found'): 
			# open model variables file
			inputs = open(self.inname, mode= 'r' ) 
			try: # in case reading in model variables fine
				# read file and store everything into 
				# a list
				in_list = inputs.readlines() 
				inputs.close() # close file
			
			# in case problem with reading in model 
			# variables
			except:
				while True:
					line = file.readline()
					if not line:
						err_mess = str('''Error:
						 could not interpret the 
						following line in the 
						model variables file: 
						''' + str(line))
						# remove old progress 
						# message	
						self.l81b.setText('')
						self.l81b.setText(
						err_mess)
						# change border 
						# accordingly
						if (self.bd_st == 1):
							self.l80.\
							setStyleSheet(
							0., 
							'''2px dashed 
							red''', 0., 
							0.)
						if (self.bd_st >= 2):
							self.l80.\
							setStyleSheet(
							0., '''2px 
							solid red''', 
							0., 0.)
						
						# prepare for change to 
						# border status
						self.bd_st += 2
						# change border status
						if (self.bd_st == 3):
							self.bd_st = 2
						if (self.bd_st >= 4):
							self.bd_st = 1
						return()
		else: # if using defaults
			in_list = []

		# if parameters automatically supplied via the 
		# automated_setup_and_call module
		if (type(self.param_const) == dict and 
			self.inname == 'Not found'):
			in_list = [] # prepare list
			# loop through supplied parameters
			for key, value in self.param_const.items():
				# if no editing needed, then just append
				in_list.append(str(str(key) + ' = ' + 
					str(value)))
		
		# change border/error message status to ready for change
		self.bd_st = 3
		# default value for number of modes represented by 
		# particle number concentration
		self.pmode_cnt = 0
		
		# loop through supplied model variables to interpret
		for i in range(len(in_list)):
			
			# ----------------------------
			# if commented out continue to next line
			if (in_list[i][0] == '#'):
				continue
			try:
				# split values from keys
				key, value = in_list[i].split('=')
				
			except:
				err_mess = str('Did not see an \'=\' symbol in ' +
					'a line in the model variables file, ' +
					'please check model variables file and ' +
					'ensure that any lines supposed to be ' +
					'comments begin with #, ' +
					'and see README for guidance')
				continue
			# model variable name - a string with bounding white 
			# space removed
			key = key.strip()
			
			# ----------------------------------
			# name of folder to save results in
			if (key == 'res_file_name' and (value.strip())):
				self.sav_nam = str(value.strip())
			
			# formatting for chemical scheme
			if (key == 'chem_scheme_markers' and (value.strip())):
				self.chem_sch_mrk = [str(i).strip() for i 
				in (value.split(','))]

			# whether to define components by atomic composition
			# given in the chemical scheme file, rather than by
			# SMILES
			if (key == 'ac_by_cs' and (value.strip())):
				self.ac_by_cs = int(value.strip())

			# path to chemical scheme
			if (key == 'chem_sch_name' and (value.strip())): 
				self.sch_name = str(value.strip())

			# path to chemical scheme rate constants
			if (key == 'rate_constant_file' and (value.strip())): 
				self.rate_cons_name = str(value.strip())

			# path to xml file
			if (key == 'xml_name' and (value.strip())):
				self.xml_name = str(value.strip())
			
			# time step (s) for updating ODE initial 
			# conditions
			if (key == 'update_step' and (value.strip())): 
				self.update_stp = float(value.strip())

			if (key == 'total_model_time' and 
				(value.strip())):
				try:
					self.tot_time = float(
						value.strip())
				except:
					err_mess = str('Could not \
					convert \
					string to float for \
					total_model_time model \
					variable, \
					please check the model \
					variables file and see README \
					for guidance')

			# whether or not to skip parsing of the 
			# chemical scheme
			if (key == 'pars_skip' and (value.strip())):
				try:
					# in case a numerical flag
					self.pars_skip = int(value.strip())
					# if 1, then check that local 
					# variables are stored, ready for use
					if (self.pars_skip == 1):
						try:
							self.Psat = self.Psat_rec0
						except:
							self.pars_skip = 0
				except:
					# in case a path to saved variables
					self.pars_skip_path = str(value.strip())
					self.pars_skip = 2

			if (key == 'spin_up' and (value.strip())):
			
				try:
					self.spin_up = int(value.strip())
				except:
					err_mess = 'Could not convert string to integer for spin_up variable, please check model variables file and see README for guidance'

			# names of components present at 
			# experiment start
			if (key == 'Comp0' and (value.strip())):
				self.comp0 = [str(i).strip() for i in 
				(value.split(','))]			

			# initial concentrations of components present 
			# at experiment start (ppb)
			if (key == 'C0' and (value.strip())):
				
				# treat as path to file
				if '/' in value or '\\' in value:
					self.path_to_C0 = \
						str(value.strip())
					[y0_gas, self] = C0_open(self)
					err_mess = self.err_mess
					if err_mess[0:5] == 'Error':
						break
				else:

					try:
						y0_gas = [float(i) for i 
						in (value.split(','))]
					except:
					
						err_mess = '''Error - 
						could not read in the 
						C0 model variable, please check the model 
						variables file and see README for guidance'''
			
			# chamber temperature (K)
			if (key == 'temperature' and (value.strip())):
				self.TEMP = [float(i) for i in ((value.strip()).split(','))]
			# times (s) that temperature values correspond to
			if (key == 'tempt' and (value.strip())): 
				self.tempt = [float(i) for i in ((value.strip()).split(','))]

			# relative humidity in chamber (0-1)
			if (key == 'rh' and (value.strip())):
				self.RH = np.array(([float(i) for i in
					((value.strip()).split(','))]))
				
			# times through simulation (s) at which 
			# relative humidity reached
			if (key == 'rht' and (value.strip())):
				self.RHt = np.array(([float(i) for i 
					in ((value.strip()).split(','))]))

			# pressure inside chamber
			if (key == 'p_init' and (value.strip())):
				self.Press = np.array(([float(i) for i in
					((value.strip()).split(','))]))

			# time of day at experiment start (s)
			if key == 'daytime_start' and (value.strip()):
				self.daytime = float(value.strip())

			# marker for whether or not to consider wall
			if key == 'wall_on' and (value.strip()):
				if value.strip()[-1] == '.':
					wall_on_value = value.strip()[0:-1]
				else:
					wall_on_value = value.strip()
				if (int(wall_on_value) == 0): # in case wall to be turned off
					self.wall_on = int(wall_on_value)
				if (int(wall_on_value) == 1): # in case wall to be turned on
					# turn on wall if it has been turned off prior
					if (self.wall_on < 1):
						self.wall_on = int(wall_on_value)
					# in case the number of wall bins is already 
					# registered, then we don't need to 
					# adjust the self.wall_on value
					else:
						continue
								
			if key == 'number_wall_bins' and (value.strip()): # number of wall bins
				
				# note, that if wall_on already set to zero, then this 
				# overrides any number_wall_bins setting
				if 'self.wall_on' in locals():
					if (self.wall_on == 0):
						continue
					else:
						self.wall_on = int(
						value.strip())
				else:
					self.wall_on = int(
					value.strip())
				
			# effective absorbing mass concentration of wall
			if (key == 'eff_abs_wall_massC' and 
				(value.strip())): 
				self.Cw = np.array([float(i) for i in 
				((value.strip()).split(','))])
			
			# mass transfer coefficient of vapours with wall
			# set to -1 for Huang et al. 2018 
			# (Eq. 2 and accompanying text), 
			# https://doi.org/10.1021/acs.est.7b05575
			if (key == 'mass_trans_coeff' and 
				(value.strip())):

				# check if this is a path to a file 
				# containing values, if not treat as a list of
				# component names
				# treat as path to file containing mass 
				# transfer coefficients
				if '/' in value or '\\' in value:
					if 'D_ig' not in value:
						self.mtc_path = str(
						value.strip())
						value = mass_trans_coeff_open(
						self)
						err_mess = self.err_mess
						if err_mess[0:5] == 'Error':
							break
					
				max_comp = 1
				max_of_all = []
				# prepare to read in string between 
				# semi-colons
				reader = ''
				# list of prescribed mass transfer 
				# coefficients
				prescribes = []
				# number of walls mentioned in the 
				# mass transfer coefficient
				num_wall_mtfs = 1
				# count maximum number of components for 
				# a single wall
				for i in value.strip():
					if i == ';': # new component
						max_comp += 1
						# remember old component 
						# and wall and 
						# mass transfer 		
						# coefficient
						prescribes.append(
						reader)
						# prepare to read in 
						# string between 
						# semi-colons
						reader = ''
						
					if (i == ','): # new wall
						max_of_all.append(max_comp)
						max_comp = 1
						num_wall_mtfs += 1
						prescribes.append(reader)
						# prepare to read in string between semi-colons
						reader = ''

					if (i != ';' and i != ',' and i != ' '):
						reader = str(reader + i)
				
				# ensure final wall is registered
				prescribes.append(reader)
				
				# ensure final wall accounted for
				max_of_all.append(max_comp)

				# remember the components that are 
				# prescribed, their wall number 
				# and mass transfer coefficient
				self.wmtc_names = []
				self.wmtc_wn = []
				self.wmtc = []
				# count on number of prescribed mass 
				# transfer coefficients
				pre_cnt = 0

				# prepare holding array for rate 
				# coefficient information
				self.kw = np.ones((num_wall_mtfs, 
				max(max_of_all)))*-1.e-6
				
				wall_num = 0 # count on walls
				# number of entries for each wall
				ind_wall_cnt = 1
				for prei in prescribes:

					# if this gives the coefficient 
					# for a specific component
					if prei.count('_') == 2 and 'D_ig' not in prei:
						
						# get component name
						self.wmtc_names.append(prei[0:prei.index('_')])
						# get wall number (note spreadsheet values 
						# start at 1 for first wall (not zero))
						self.wmtc_wn.append(int(prei[prei.index(
						'_')+5:prei.index('_')+5+prei[prei.index(
						'_')+5::].index('_')])-1)
						# get mass transfer coefficient
						self.wmtc.append(float(prei[prei.index('_')+1 
						+ prei[prei.index('_')+1::].index('_')+1::]))
						self.kw[int(self.wmtc_wn[-1]), 
							ind_wall_cnt-1] = -1.e-7
					# if just a generic mass transfer coefficient for 
					# this wall
					else:
						if 'D_ig' not in prei:
							self.wmtc_names.append(
								'all_other_components')
							self.wmtc.append(float(prei))
							self.kw[wall_num, 
								ind_wall_cnt-1] = -1.e-7

					# if this is an equation, e.g. eq. 1 of Methods 
					# section in Roldin et al. 2019
					# (https://doi.org/10.1038/s41467-019-12338-8)
					if 'D_ig' in prei:
						import write_mass_trans_coeff
						write_mass_trans_coeff.write(str(
						prei), self)
						# tell partit_var_prep and cham_up that
						# mass transfer coefficients are provided
						# in equation form
						self.mtc_calc_flag = 1
					
					# count on number of prescribed mass transfer 
					# coefficients
					pre_cnt += 1
					ind_wall_cnt += 1 # number of entries for each wall
					# if moving up a wall
					if pre_cnt == sum(max_of_all[0:wall_num+1]):
						wall_num += 1
						# number of entries for each wall
						ind_wall_cnt = 1
			
			if key == 'chamSA' and (value.strip()): # chamber surface area (m2)
				chamSA = float(value.strip())

			if key == 'chamV' and (value.strip()): # chamber volume (m3)
				chamV = float(value.strip())

			if key == 'size_structure' and (value.strip()): # the size structure
				siz_stru = int(value.strip())

			# number of particle size bins
			if (key == 'number_size_bins' and (value.strip())):
				num_sb = int(value.strip())
				self.num_asb = int(value.strip())

			# seed particle number concentrations (# particles/cm3)
			if (key == 'pconc' and (value.strip())):
			
				time_cnt = 1 # track number of times
				sb_cnt = 1 # track number of size bins
				self.pmode_cnt = 1 # track number of modes
			
				for i in value:
					if i == ';': # semi-colon represents a time difference
						time_cnt += 1 # increase time count
					if (time_cnt == 1 and i == ','):
						sb_cnt += 1 # increase size bin count
						# explicitly stated particle 
						# concentrations
						self.pmode = 1
						self.pmode_cnt = 0 # no modes
					if (time_cnt == 1 and i == ':'):
						self.pmode_cnt += 1 # mode count
						# particle concentrations 
						# expressed as modes
						self.pmode = 0

				# a possible situation where there is just one size bin and 
				# therefore only one concentration given for that size bin
				if (sb_cnt == 1 and self.pmode_cnt == 1):
					# explicitly stated particle concentrations
					self.pmode = 1 
						
				# see below model input loop for final determination of pmode (
				# specifically what happens if more 
				# than one size bin, but just one mode given)
					
				# if number concentration per size bin given explicitly
				if (self.pmode == 1):
					self.pconc = np.zeros((sb_cnt, time_cnt))
					for i in range(time_cnt):
						self.pconc[:, i] = [float(ii.strip()) for ii 
						in ((value.split(';')[i]).split(','))]
				else: # mode quantities provided
					self.pconc = np.zeros((self.pmode_cnt, time_cnt))
					for i in range(time_cnt):
						self.pconc[:, i] = [float(ii.strip()) for ii 
						in ((value.split(';')[i]).split(':'))]

			# seed particle input times (s)
			if (key == 'pconct' and (value.strip())):
				time_cnt = 1 # track number of times
				for i in value:
					if (i == ';'):
						time_cnt += 1 # increase time count
					
				# times in columns
				self.pconct = np.zeros((1, time_cnt))
				self.pconct[0, :] = [float(i) 
					for i in ((value.strip()).split(';'))]

			# whether to repeat particle influxes every 24 hours
			if (key == 'pconctf' and (value.strip())):
				self.pconctf = float(value.strip())

			# whether seed particle injection continuous or 
			# instantaneous
			if (key == 'pcont' and (value.strip())):
				time_cnt = 1 # track number of times
				for i in value:
					if (i == ';'):
						time_cnt += 1 # increase time count
					
				# times in columns
				self.pcont = np.zeros((1, time_cnt), dtype=int)
				self.pcont[0, :] = [int(i) for i in 
				((value.strip()).split(';'))]
			

			# lowest size bin bound
			if (key == 'lower_part_size' and (value.strip())):
				lowsize = str(value.strip())
				# if a list, representing manually 
				# set radius (um) bounds
				if ',' in lowsize:
					self.manual_rbounds = [float(i) for i in 
						((value.strip()).split(','))]
					lowsize = self.manual_rbounds[0]
				else:
					lowsize = float(value.strip())
			
			# radius of uppermost size bin boundary
			if key == 'upper_part_size' and (value.strip()):
				uppsize = float(value.strip())
			
			# method of spacing size bin sizes
			if (key == 'space_mode' and (value.strip())):
				self.space_mode = str(value.strip())

			# seed particle standard deviation of 
			# number size distribution
			if key == 'std' and (value.strip()):
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
					std[:, i] = [float(ii.strip()) for ii in 
						((value.split(';')[i]).split(':'))]
			
			# seed particle mean radius (um)
			if (key == 'mean_rad' and (value.strip())):
				time_cnt = 1 # track of number of times
				mode_cnt = 1 # track number of modes
				for i in value:
					if (i==';'):
						time_cnt += 1 # increase time count
						# reset mode count
						mode_cnt = 1 # track number of modes
					if (i==':'):
						mode_cnt += 1 # increase mode count
				
				self.mean_rad = np.zeros((mode_cnt, time_cnt))
				for i in range(time_cnt):
					self.mean_rad[:, i] = [float(ii.strip()) for ii in 
						((value.split(';')[i]).split(':'))]
				
			# whether (1) or not (0) to treat gas-particle paritioning as 
			# equilibrium (i.e. not dynamic)
			if (key == 'equi_gtop_partit'and (value.strip())):
				self.equi_gtop_partit = int(value.strip())

			# frequency (s) of storing results	
			if (key == 'recording_time_step' and 
				(value.strip())): 
				self.save_step = float(value.strip())

			# names of components with constant abundances
			if (key == 'const_comp' and (value.strip())): 
				value = str(value).replace(' ', '')
				comp_count = 1 # count number of components
				time_count = 1 # track number of times

				# prepare to hold component names
				const_comp_list = []
				ic = 0 # keep count on character in value
				for i in value:
					if i==';':
						time_count += 1 # record number of times
						self.const_comp = np.concatenate(
							(self.const_comp, 
							self.const_comp[:, 0].reshape(-1, 1))
							, axis=1)
						self.const_comp[:, -1] = ''

						try:
							const_comp = value[
							ic+1:ic+1+value[ic+1::].index(',')]
							# in case change of time sign included
							if (';' in const_comp):
								const_comp = const_comp[
								0:const_comp.index(';')]
							const_comp = const_comp.replace(
								'\n', '')
						except: # in case at final component
							const_comp = value[ic+1::]
						
						indx = np.where(self.const_comp[
							:, 0] == const_comp)[0][0]
						self.const_comp[indx, -1] = const_comp
						
					if i==',':
						try:
							const_comp = value[
							ic+1:ic+1+value[ic+1::].index(',')]
							# in case change of time sign included
							if (';' in const_comp):
								const_comp = const_comp[
								0:const_comp.index(';')]
						except: # in case at final component
							const_comp = value[ic+1::]
							const_comp = const_comp.replace(
								'\n', '')
						if (time_count == 1): # if on first time
							self.const_comp = np.append(
							self.const_comp, 
							np.zeros((1,1)), axis=0)
							self.const_comp[-1, 0] = const_comp
						else: # if on later times
							indx = np.where(self.const_comp[
							:, 0] == const_comp)[0][0]
							self.const_comp[indx, -1] = const_comp

					# if on first character
					if (ic == 0):
						try:
							const_comp = (value[
								ic:value.index(',')])
						except: # if just one component given
							const_comp = (value[:])
							const_comp = const_comp.replace(
								'\n', '')
						# prepare to hold component names (rows) 
						# at the required
						# times (columns)
						self.const_comp = np.zeros((1, 1)).astype('str')
						self.const_comp[0, 0] = const_comp
						
					ic += 1 # keep count on character in value
				
			# times that constant components are for
			if (key == 'const_compt' and (value.strip())):
				self.const_compt = np.array(([str(i).strip() for i in 
					(value.split(';'))])).astype('float')

			# name of file containing observations to 
			# constrain by
			if (key == 'obs_file' and (value.strip())):
				self.obs_file = str(value.strip())
				# get observed values
				from obs_file_open import obs_file_open
				self = obs_file_open(self)
			
			# name of file to save calculated continuous
			# influx rates to
			if (key == 'sim_cont_infl_file' and 
				(value.strip())):
				self.sim_ci_file = str(value.strip())

			# names of components with instantaneous 
			# gas-phase injections
			if (key == 'Compt' and (value.strip())):
				Compt = [str(i).strip() for i in 
					(value.split(','))]

			# times of later gas-phase instantaneous 
			# injections (s)
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
					Ct[i, :] = [float(ii.strip()) for ii 
					in ((value.split(';')[i]).split(','))]

			# name(s) of component(s) comprising seed particles
			if (key == 'seed_name' and value.strip()):
				self.seed_name = [str(i).strip() 
					for i in (value.split(','))]
				
			if (key == 'seed_mw' and value.strip()):
				try:
					seed_mw =  [float(i) for i in (value.split(','))]
					seed_mw =  np.array((seed_mw))
				except:
					seed_mw = ['fail']
		
			# dissociation constant(s) of seed component(s)
			if (key == 'seed_diss' and value.strip()):
				self.core_diss = [float(i) for i in 
				(value.split(','))]

			# dissociation constant(s) of seed component(s) with
			# respect to water
			if (key == 'seed_diss_wrtw' and value.strip()):
				self.core_diss_wrtw = [float(i) for i in 
				(value.split(','))]

			# dissociation constant(s) of nonseed component(s) with
			# respect to water
			if (key == 'nonseed_diss_wrtw' and value.strip()):
				self.noncore_diss_wrtw = [float(i) for i in 
				(value.split(','))]

			# dissociation constant of water with respect to organics
			if (key == 'H2O_diss_wrtorg' and value.strip()):
				self.H2O_diss_const_org = [float(i) for i in 
				(value.split(','))]

			if (key == 'seed_dens' and value.strip()):
				seed_dens = [float(i) for i in 
				(value.split(','))]
				seed_dens = np.array((seed_dens))

			# mole fraction of components in dry seed 
			# particles
			if (key == 'seedx' and value.strip()):
				# count number of components
				comp_count = 1
				sb_count = 1 # track number of size bins
				ti_count = 1 # track number of times
				for i in value:
					if (i == ':'):
						# count number of times
						ti_count += 1
						# reset number of 
						# components
						comp_count = 1
						# reset size bin count
						sb_count = 1
					if (i == ';'):
						# record number of 
						# components
						comp_count += 1
						# reset size bin count
						sb_count = 1 
					if (i == ','):
						# size bin count
						sb_count += 1
				self.seedx = np.zeros((comp_count, sb_count, 
					ti_count))
				# count on seed mole fractions through 
				# time
				self.seedx_tcnt = 0
				# split by times
				seedx_split = value.split(':')
				# loop through times
				for i0 in range(ti_count):
					seedxn = seedx_split[i0]
					# loop through components
					for i in range(comp_count):
						self.seedx[i, :, i0] = [
						float(ii.strip()) for ii
						 in ((seedxn.split(';')
						[i]).split(','))]
				
			# whether number size distribution includes 
			# volume of water
			if (key == 'Vwat_inc' and value.strip()):
				self.Vwat_inc = int(value.strip())
			
			# whether to allow water to equilibrate with 
			# seed prior to experiment
			if (key == 'seed_eq_wat' and value.strip()):
				self.seed_eq_wat = int(value.strip())

			# fraction of total partitioning coefficient
			# below which gas-particle partitioning rate
			# is set to zero, e.g. because of tiny
			# surface area in a particle size bin
			if (key == 'z_prt_coeff' and value.strip()):
				z_prt_coeff = float(value.strip())

			# fraction of total (excluding water and seed)
			# gas- and particle-phase molecular concentration
			# below which partitioning to particles and wall
			# is set to 0
			if (key == 'z_prt_coeff_loC' and value.strip()):
				self.z_prt_coeff_loC = float(value.strip())			

			# status of lights (on or off)
			if (key == 'light_status' and value.strip()):

				# check if a path to file containing photolysis rates
				try:
					from J_value_file_open import J_value_file_open
					self = J_value_file_open(str(value.strip()), self)
					self.light_stat = (np.ones((
						len(self.light_time)))*3).astype('int')
				except:
					light_stat = [int(i) for i in (value.split(','))]
					self.light_stat = np.array((light_stat))
					
					# checl for flag to use pre-calculated transmission
					# factor for solar radiation (Hayman method)
					# through clear glass (see photolysisRates
					# for more information)
					if (sum(np.array((self.light_stat)) == 3) > 0):
						# ensure lights on
						self.light_stat = np.ones((1)).astype('int')
						# set a tf_range value that tells
						# photolysisRates to apply
						# a pre-calcalated transmission factor 
						self.tf_range = 2
					
			
			# times (s) corresponding to light status
			if (key == 'light_time' and value.strip()):
				light_time = [float(i) for i in (value.split(','))]
				self.light_time = np.array((light_time))
				
			if (key == 'lat'): # latitude (degrees)
				if (value.strip()): self.lat = float(value.strip())

			if (key == 'lon'): # longitude (degrees)
				if (value.strip()): self.lon = float(value.strip())

			# for specified actinic flux, note 'act_flux_file' 
			# is old version
			if (key == 'act_flux_path' or key == 'act_flux_file'): 
				if (value.strip()): 					
	   				self.af_path = str(value.strip())

			if key == 'DayOfYear' and (value.strip()):
				self.dayOfYear = int(value.strip())		
			
			# name of file with wavelength-dependent absorption 
			# cross-section and quantum yield 
			# calculations
			if key == 'photo_par_file' and (value.strip()):
				self.photo_path = str(os.getcwd() + 
				'/PyCHAM/photofiles/' + value.strip())
			
			# transmission factor for light
			if (key == 'trans_fac' and (value.strip())):

			
				# if wavelength-dependent transmission factors supplied
				if '_' in value.strip():
					self.tf = [str(i.strip()) for i in (value.split(','))]
					# flag for wavelength-dependency
					self.tf_range = 1
				else: # if a single transmission factor supplied
					self.tf = (([float(i.strip()) for i in
					 (value.split(';'))]))
	
					# so long as tf_range not set to
					# 2 in light_status intepretation above
					if (self.tf_range != 2):
						# flag for no wavelength-dependency
						self.tf_range = 0

			# time through simulation (s) transmission factor 
			# for light applies to
			if (key == 'trans_fact' and (value.strip())):
				
				self.tft = np.array(([float(i.strip()) for i 
				in (value.split(';'))]))

			# transmission factors for 254 nm light
			if key == 'tf_UVC' and (value.strip()):
				self.tf_UVC = [float(i.strip()) for i in (value.split(','))]
				
			if key == 'tf_UVCt' and (value.strip()): # transmission factor times for 254 nm light
				self.tf_UVCt = np.array(([float(i.strip()) for i in (value.split(','))]))

			# whether to adapt time step to naturally 
			# varying light intensity
			if key == 'light_adapt' and (value.strip()):
				self.light_ad = int(value.strip())

			# to keep solar intensity constant
			if key == 'secx' and (value.strip()):
				self.secx = float(value.strip())

			# to keep solar intensity constant
			if key == 'cosx' and (value.strip()):
				self.cosx = float(value.strip())

			# names of components with continuous influx
			if (key == 'const_infl' or key == 'cont_infl'):
				if (value.strip()):				
					# start by assuming 
					# that constant influx unit is 
					# ppb
					self.abun_unit = 'ppb'

					# check if this is a path to a 
					# file containing continuous 
					# influxes, if not treat as a 
					# list of component names
					# attempt as path to file 
					# containing continuous influxes
					try: 
						from cont_infl_file_open import cont_infl_open
						self.const_infl_path = str(
						value.strip())
						self = cont_infl_open(
							self)
						
					# treat as list of components
					except:
						self.con_infl_nam = np.array((
						[str(i).strip() for i in (value.split(','))]))

					if ('not in a file' in self.con_infl_nam):
						self.con_infl_nam = np.array((
						[str(i).strip() for i in (value.split(','))]))	


			if (key == 'const_infl_t' or key == 'cont_infl_t'):
				if (value.strip()): # times of continuous influxes (s)
					self.con_infl_t = [float(i.strip()) for i in (value.split(','))]
					self.con_infl_t = np.array((self.con_infl_t))

			# flag for how to treat times of continuous influxes
			if (key == 'cont_infl_tf' and value.strip()):
				self.con_infl_tf = float(value.strip())
				
				# in case continuous influxes and times
				# allocated to variables prior to this
				if (self.con_infl_tf == 1):
						
					if hasattr(self, 'con_infl_t'):
						clim = sum(self.con_infl_t < 24.*3.6e3)
						self.con_infl_t = self.con_infl_t[0:clim]
						self.con_infl_C = self.con_infl_C[:, 0:clim]

			# influx rate of components with continuous influx (ppb/s)
			if (key == 'Cinfl' and (value.strip())):
				
				comp_count = 1 # count number of components
				time_count = 1 # track number of times
				for i in value:
					if (i==';'):
						comp_count += 1 # record number of components
					if (i==',' and comp_count == 1):
						time_count += 1 # record number of times
				self.con_infl_C = np.zeros((comp_count, time_count))
				
				try:
					# loop through components
					for i in range(comp_count):
						# loop through times
						for ii in range(time_count):
							self.con_infl_C[i, ii] = float(
							(((value.split(';')[i]).split(
							',')))[ii].strip())
							
				# in case semicolons and commas messed 
				# up on input, note this will invoke an 
				# error message from the user input 
				# check module
				except:
					self.con_infl_C = np.empty(0)

			# names of components whose tendency to change 
			# will be tracked
			if (key == 'tracked_comp' and (value.strip())):
				self.dydt_trak = [str(i).strip() for i 
					in (value.split(','))]
			
			if (key == 'dens_Comp' and (value.strip())):
				dens_comp = [str(i).strip() for i in 
				(value.split(','))]

			if (key == 'dens' and (value.strip())):
				dens = [float(i) for i in 
					(value.split(','))]

			if (key == 'vol_Comp' and (value.strip())):
				vol_comp = [str(i).strip() for i in 
				(value.split(','))]

			if key == 'volP' and (value.strip()):
				volP = [float(i) for i in (value.split(','))]

			if key == 'volP' and (value.strip()):
				volP = [float(i) for i in (value.split(','))]

			# vapour pressure estimation method for non-HOMs
			if (key == 'nonHOMs_vp_method' and (value.strip())):
				self.nonHOMs_vp =  str(value.strip())

			# vapour pressure estimation method for HOMs
			if (key == 'HOMs_vp_method' and (value.strip())):
				self.HOMs_vp =  str(value.strip())

			# names of components with specified activity coefficients
			if (key == 'inorg_part_flag' and (value.strip())):
				self.inorg_part_flag = float(value.strip())

			# names of components with specified activity coefficients
			if (key == 'act_comp' and (value.strip())):
				act_comp = [i for i in 
				(((value.strip()).split(',')))]

			# activity coefficients (dimensionless set by user)
			if (key == 'act_user' and (value.strip())):
				act_user = [i for i in 
				(((value.strip()).split(',')))]

			# names of componenets with specified 
			# accommodation coefficients
			if (key == 'accom_coeff_comp' and (value.strip())):
				accom_comp = [i for i in (((value.strip()).split(',')))]

			# value(s) of accommodation coefficient set by user
			if (key == 'accom_coeff_user' and (value.strip())):
				accom_val = [i for i in (((value.strip()).split(',')))]

			# value(s) of the gas-particle
			# partitioning cutoff (Pa)
			if (key == 'ppartit_cutoff' and (value.strip())):
				self.ppartit_cutoff = [float(i) for 
				i in (((value.strip()).split(',')))]

			# value(s) of the gas-wall
			# partitioning cutoff (Pa)
			if (key == 'wpartit_cutoff' and (value.strip())):
				self.wpartit_cutoff = [float(i) for 
				i in (((value.strip()).split(',')))]
			
			# marker for whether to clone latest version 
			# of UManSysProp from web
			if key == 'umansysprop_update' and (value.strip()): 
				uman_up = int(value.strip())

			if key == 'int_tol' and (value.strip()): # tolerances for integration
				int_tol = [float(i) for i in (value.split(','))]

			# radius of newly nucleated particles (cm)
			if key == 'new_partr' and (value.strip()):
				new_partr = float(value.strip())

			# first parameter in the nucleation equation
			if (key == 'nucv1' and (value.strip())): 
				self.nucv1 = float(value.strip())

			# second nucleation parameter (onset)
			if (key == 'nucv2' and (value.strip())): 
				self.nucv2 = float(value.strip())

			# third nucleation parameter (duration)
			if (key == 'nucv3' and (value.strip())): 
				self.nucv3 = float(value.strip())

			# chemical scheme name of nucleating component
			if (key == 'nuc_comp' and (value.strip())): 
				self.nuc_comp = [str(i).strip() for i in
 					(value.split(','))]

			# marker for whether to adapt time interval to 
			# nucleation
			if (key == 'nuc_adapt' and (value.strip())): 
				self.nuc_ad = int(value.strip())

			# marker for whether to model coagulation
			if key == 'coag_on' and (value.strip()):
				coag_on = int(value.strip())

			# diameter at which wall deposition of particles inflection occurs
			if key == 'inflectDp' and (value.strip()):
				inflectDp = float(value.strip())

			# gradient of wall deposition of particles before inflection
			if key == 'Grad_pre_inflect' and (value.strip()):
				pwl_xpre = float(value.strip())

			# gradient of particle deposition to wall after inflection
			if key == 'Grad_post_inflect' and (value.strip()):
				pwl_xpro = float(value.strip())

			# particle deposition to wall rate at inflection
			if key == 'Rate_at_inflect' and (value.strip()):
				inflectk = float(value.strip())

			# chamber surface area (m2) used for particle loss to walls
			if key == 'ChamSA' and (value.strip()):
				ChamSA = float(value.strip())

			# marker for whether to use the McMurry model 
			# for particle deposition to wall
			if key == 'McMurry_flag' and (value.strip()):
				Rader = int(value.strip())

			# average number of charges per particle
			if key == 'part_charge_num' and (value.strip()):
				p_char = float(value.strip())

			# electric field in chamber
			if key == 'elec_field' and (value.strip()):
				e_field = float(value.strip())
			
			# dilution factor rate
			if key == 'dil_fac' and (value.strip()):
				self.dil_fac = np.array(([float(i) for i in 
				(((value.strip()).split(',')))]))
			
			# dilution factor rate times through experiment (s)
			if key == 'dil_fact' and (value.strip()):
				self.dil_fact = np.array(([float(i) for i in 
				(((value.strip()).split(',')))]))
			
			# whether to serialise water gas-particle partitioning
			if key == 'ser_H2O' and (value.strip()):
				ser_H2O = int(value)
			
			# history of particle-phase with respect to particle-phase water
			if (key == 'H2O_hist' and (value.strip())):
				try:
					wat_hist = int(value)
				except:
					wat_hist = -1 # will cause error message
			# deliquescence RH dependence on temperature
			if (key == 'drh_ft' and (value.strip())):
				try:
					drh_str = str(value)
				except:
					drh_str = -1 # will cause error message
			# efflorescence RH dependence on temperature
			if (key == 'erh_ft' and (value.strip())):
				try:
					erh_str = str(value)
				except:
					erh_str = -1 # will cause error message

			# whether to remove influxes of components that aren't seen in
			# chemical scheme
			if (key == 'remove_influx_not_in_scheme' and (value.strip())):
				self.remove_influx_not_in_scheme = int(value)

			# name of component who's abundance needs to be nudged to give
			# the prescribed RO2 pool abundance, followed by the prescribed
			# RO2 abundance (molecules/cm3)
			if (key == 'comp_nudge_RO2' and (value.strip())):
				self.comp_nudge_RO2 = [str(i).strip() for i in
 					(value.split(','))][0]
				self.RO2_nudge_target = float([str(i).strip() for i in
 					(value.split(','))][1])
			
			# user-defined outputs (molecules/cm3)
			if (key == 'user_output' and (value.strip())):
				# begin with the unconditional variables
				self.user_output = ['unconditional_variables']
				# append user-defined variables
				user_output_list = [str(i).strip() for i in
 					(value.split(','))]
				self.user_output += user_output_list
		
		# UManSysProp check ----------------------------------
		# for UManSysProp if no update requested, check that there 
		# is an existing UManSysProp folder
		if (uman_up == 0): # check for existing umansysprop folder
			# if no existing folder then force update
			if not os.path.isdir(self.PyCHAM_path + '/umansysprop'):
				uman_up = 1
		# -------------------------------------------
		
		# update model variables message in GUI
		# if error message occurs
		if (err_mess != '' or self.err_mess != ''):
			# update error message
			if (err_mess != ''):
				self.l80.setText(str('Setup Status: \n' + err_mess))
			if (self.err_mess != ''):
				self.l80.setText(str('Setup Status: \n' + self.err_mess))
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
			
			()

		# a possible situation where there are multiple particle size bins,
		# but just one mode given for the particle number size distribution
		if (num_sb > 1 and self.pmode_cnt == 1):
			self.pmode = 0 # modal particle concentrations
		
		# prepare for pickling
		list_vars = [y0_gas, 
				siz_stru, num_sb, lowsize, 
				uppsize, std, 
				Compt, injectt, Ct, seed_mw, 
				seed_dens, 
				dens_comp, dens, vol_comp, volP, 
				act_comp, act_user, accom_comp, 
				accom_val, uman_up, int_tol, 
				new_partr, coag_on, 
				inflectDp, pwl_xpre, pwl_xpro, inflectk,
				 chamSA, Rader, p_char, 
				e_field, ser_H2O, 
				wat_hist, drh_str, 
				erh_str, z_prt_coeff, chamV]

		# the file to be used for pickling
		with open(input_by_sim, 'wb') as pk: 
			pickle.dump(list_vars, pk) # pickle
			pk.close() # close
		
		return() # end function	
		
	read(self) # call on function to read the model variables

	return()

# function for converting excel spreasheet of surface depositions 
# into a string that commands above can interpret
def mass_trans_coeff_open(self):

	import openpyxl
	import os

	try: # try to open the file at the user-supplied path
		wb = openpyxl.load_workbook(filename = self.mtc_path)
	except: # if file not found tell user
		self.err_mess = str('Error: file path provided by user in model variables file for gas-wall mass transfer coefficient of components was not found, file path attempted was: ' + self.mtc_path)
		return(self)

	sheet = wb['mtc']
	# component names are in second column, mass transfer coefficients are in first column		
	ic = -1 # count on row iteration
	
	# prepare to store component names and mass transfer coefficient
	value = ''
	
	
	for i in sheet.iter_rows(values_only=True): # loop through rows
		ic += 1 # count on row iteration
		if (ic == 0): # do not need header
			continue				

		# get names of components (matching chemical scheme names) 
		# and their emission rates (abundance unit given above/s)
		else:
			if (i[1] is None): # reached end of contiguous components
				break
			# if all components then just need mass transfer coefficient value
			if i[1] == 'all_other_components':
				if (value == ''): # first installment
					walln_0 = i[3]	
				value = str(value + str(i[2]))

			else: # if component specified then get name, coefficient and wall number (in spreadsheet first wall = 1)
				if (value == ''): # if first installment
					value = str(value + i[1] + '_wall' + str(i[3]) + '_' + str(i[2]))
					walln_0 = i[3]
				else:
					# get wall number
					walln_now = i[3]
					if (walln_now == walln_0):
						value = str(value + ';'  + i[1] + '_wall' + str(i[3]) + '_' + str(i[2]))
					else:	
						value = str(value + ',' + i[1] + '_wall' + str(i[3]) + '_' + str(i[2]))
					walln_0 = i[3]	
	wb.close() # close excel file
	
	return(value)

# function for converting csv of initial concentrations
def C0_open(self):

	import os
	import ast
	from textwrap import dedent

	try: # try to open the file at the user-supplied path
		fname = '''\
		/concentrations_all_components_all
		_times_gas_particle_wall
		'''
		fname = dedent(fname)
		fname = fname.replace('\n', '')
		fname = str(self.path_to_C0 + fname)
		
		y0 = np.loadtxt(fname, delimiter=',', skiprows=1)
		# keep just the final time step concentrations (ppb)
		y0 = y0[-1, :]

	except: # if file not found tell user
		self.err_mess = str(''''Error: file path provided by
		 user in model variables file for initial concentrations
		 of components not found, file path attempted was: ''' 
		+ fname)
		# bypass error, e.g. when you know that the file will
		# become available following completion of preceding
		# simulations
		self.err_mess = ''
		return([], self)

	# now get names of components that concentrations correspond to
	# get chemical scheme names

	# name of file where experiment constants saved
	fname = str(self.path_to_C0 + '/model_and_component_constants')
	
	const_in = open(fname)
	self.comp0 = [] # prepare to store chemical scheme names
	# get chemical scheme names of components
	for line in const_in.readlines():
		if (str(line.split(',')[0]) == 'chem_scheme_names'):
			# find index of first [ and index of last ]
			icnt = 0 # count on characters
			for i in line:
				if i == '[':
					st_indx = icnt
					break
				icnt += 1 # count on characters

			for cnt in range(10):
				if line[-cnt] == ']':
					fi_indx = -cnt+1
					break

			self.comp0 = ast.literal_eval(line[st_indx:fi_indx])

	if (self.comp0 == []): # if new style of saving results used
		fname = str(self.path_to_C0 + '/comp_namelist.npy')
		self.comp0 = (np.load(fname, 
			allow_pickle=True)).tolist()
	
	# if particle results are saved from previous simulation
	if os.path.exists(str(self.path_to_C0 + '/particle_number_concentration_wet')):

		# set starting particle properties
		# withdraw number-size distributions (# particles/cm3 (air))
		fname = str(self.path_to_C0 + '/particle_number_concentration_wet')
		self.N_perbin0_prev_sim = np.loadtxt(fname, delimiter=',', skiprows=1)
		# keep just final time
		self.N_perbin0_prev_sim = self.N_perbin0_prev_sim[-1, :].reshape(-1, 1)	

		# particle radii (um)
		fname = str(self.path_to_C0 + '/size_bin_radius')
		# skiprows=1 omits header
		self.x0_prev_sim = np.loadtxt(fname, delimiter=',', skiprows=1)
		# keep just final time
		self.x0_prev_sim = self.x0_prev_sim[-1, :]

		# particle volumes (um3)
		self.Varr0_prev_sim = (4./3.)*np.pi*self.x0_prev_sim**3.
	

	# starting concentrations (ppb) of just the gas-phase
	y0_gas = y0[0:len(self.comp0)]

	# get concentrations (molecules/cm3) of any other
	# phases at this final time step
	if (y0.shape[0]>len(self.comp0)):
		self.y0_other_phase = y0[len(self.comp0)::]

	# now just keep the concentrations for components with
	# concentrations above zero
	self.comp0 = np.array((self.comp0))[y0_gas > 0.]
	y0_gas = y0_gas[y0_gas > 0.]
	# remove water
	y0_gas = (y0_gas[self.comp0 != 'H2O']).tolist()
	self.comp0 = (self.comp0[self.comp0 != 'H2O']).tolist()
	
	return(y0_gas, self)
