########################################################################
#								       #
# Copyright (C) 2018-2024					       #
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
'''The module that generates takes the inputs for 
PyCHAM_lite, and connects with the core PyCHAM model'''

# first module called when PyCHAM started from the terminal/command 
# window, takes inputs
# and sends to model modules, also calls the saving module

import pickle # for storing inputs
import sys
import os
import def_mod_var
import numpy as np
import re
import importlib
import mod_var_read
import time
		
class PyCHAM_lite():
	
	# set the reference (self) to the current instance of the PyCHAM class,
	# note that by setting a default value for param_const, a user doesn't need
	# to provide an input for param_const.  But, having param_const as an
	# optional argument allows the automated setup and call 
	# (automated_setup_and_call.py) to call this module
	def __init__(self, param_const=0):
		
		super().__init__()

		# store information about simulations in self
		self.param_const = param_const

		# get path to PyCHAM folder
		self.PyCHAM_path = os.path.dirname(__file__)	
		
		self.err_mess = '' # begin with no error message
		
		# default variables for all required input model 
		# variables -------------------------
		[y0, Press, siz_stru, num_sb, 
		lowsize, uppsize, std, 
		Compt, injectt, Ct, seed_mw, seed_diss, seed_dens, 
		dens_comp, dens, vol_comp, volP, act_comp, 
		act_user, accom_comp, accom_val, uman_up, int_tol, 
		new_partr, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
		inflectk, chamSA, Rader, 
		p_char, e_field, ser_H2O, wat_hist, 
		drh_str, erh_str, 
		z_prt_coeff, chamV, 
		self] = def_mod_var.def_mod_var(0, self)
		
		# if parameters have been provided automatically (without using GUI)
		# then automatically setup and run the simulation, e.g. using the
		# automated_setup_and_call.py
		if (type(param_const) == dict):
			self.autorun() # call on automatic setup and run
			
		return()
	
	def on_click81sing(self): # when single simulation button pressed
			
		# action to simulate a single simulation
		err_mess = self.on_click81b()
			
		# once all simulations done in single simulation mode, tidy up
		# clear the list of file combinations for simulations in batch
		self.output_list = [] # reset list of output paths
		
	def on_click81b(self): # when button to run simulation pressed
			
		# path where results are being saved to
		output_by_sim = os.path.join(self.sav_nam)
			
		# call on function to simulate
		err_mess = self.act_81(output_by_sim, 0)
		
		# if no error message then return with no error message
		return(err_mess)
		
	# start the simulation
	def act_81(self, output_by_sim, sim_num):
		
		from middle import middle # prepare to communicate with main program		

		note_messf = 0 # cancel note message flag
		
		for prog in middle(self): # call on modules to solve problem
			
			if (isinstance(prog, str)): # check if it's a message
				mess = prog
				# if it's an error message
				if (mess[0:5] == 'Error'):
					
					return(mess)
					
					# remove the directory that was due to hold results
					import shutil	
					shutil.rmtree(self.output_by_sim)
				else:
					# flag that a note message has been generated
					note_messf = 1 
				
			else: # if no message from model
			
				# if there was previously a note message
				if (note_messf == 1):
					note_messf = 0 # cancel note message flag

		# if this point reached then no error message generated
		mess = ''
		return(mess)
	
		
	
	def autorun(self): # function to automatically run PyCHAM
		# note that default model variables are already set 
		# before autorun function called	

		# get the string for checking on change in chemistry
		wo_str = self.param_const['wo_str']

		ri = 0 # count on simulations
	
		# loop through simulations
		for inname in self.param_const['mod_var_name']:
			
			# default variables for all required input model 
			# variables -------------------------
			[y0, Press, siz_stru, num_sb, 
			lowsize, uppsize, std, 
			Compt, injectt, Ct, seed_mw, seed_diss, seed_dens, 
			dens_comp, dens, vol_comp, volP, act_comp, 
			act_user, accom_comp, accom_val, uman_up, int_tol, 
			new_partr, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
			inflectk, chamSA, Rader, 
			p_char, e_field, ser_H2O, wat_hist, 
			drh_str, erh_str, 
			z_prt_coeff, chamV, 
			self] = def_mod_var.def_mod_var(0, self)

			# store model variable name
			self.inname = inname

			# get index of final /
			file_name_indx = -1*self.inname[::-1].index('/')
			# get just file name
			file_name = self.inname[file_name_indx::]

			# establish parameters provided by user by calling mod_var_read
			mod_var_read.mod_var_read(self)
			
			# if past the first simulation
			if (inname != self.param_const['mod_var_name'][0]):
				# check whether chemistry has changed and therefore
				# can't skip parsing the reactions
				if (wo_str in sim_name0 and 
					wo_str not in file_name):
					print('new chemistry')
					self.pars_skip = 0
				if (wo_str not in sim_name0 
					and wo_str in file_name):
					print('new chemistry')
					self.pars_skip = 0
			

			# run simulation
			self.on_click81sing() # run simulation
			print('completed ', file_name, ri)
			# remember this name
			sim_name0 = file_name
			# pause
			time.sleep(0.5)
			ri += 1 # count on simulations

		return()