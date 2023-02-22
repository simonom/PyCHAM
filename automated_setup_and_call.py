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
'''code to automatically setup a series of simulations and run PyCHAM without using the PyCHAM graphical user interface'''
# useful for users that want to setup and run simulations in an automated manner, e.g. because manual selection 
# of simulations via the graphical user interface is not practical 


# import requirements
import numpy as np
import sys
import os
from PyQt5.QtWidgets import QApplication, QWidget


# dictionary to hold parameter ranges
param_range = {}
# dictionary to hold parameters that will be held constant
param_const = {}

param_const['sim_num'] = 2 # number of simulations

# state path to chemical scheme
param_const['sch_name'] = 'C:\\Users\\Psymo\\Desktop\\PyCHAM\\PyCHAM\\PyCHAM\\input\\auto_call_test\\apinene_ch4_mcm_PRAM_schem.kpp'

# state path to xml file
param_const['xml_name'] = 'C:\\Users\\Psymo\\Desktop\\PyCHAM\\PyCHAM\\PyCHAM\\input\\auto_call_test\\apinene_ch4_mcm_xml.xml'

# state parameter ranges
param_const['res_file_name'] = 'ambient_run_num'
param_const['total_model_time'] = 1.2e3
param_const['update_step'] = 6.e2
param_const['recording_time_step'] = 6.e2
param_const['light_status'] = 1
param_const['light_time'] = 0.

# need to check that MAC actinic flux*7 is consistent with natural solar flux in Andalusia (Spain) on June 21st
param_const['act_flux_path'] = 'MAC_Actinic_Flux_Spectrum_wUVCx7p0.csv'
# need to check which wavelengths this input affects
param_range['trans_fac'] = [0., 1.] 

param_range['temperature'] = [273.15, 323.15]
param_const['tempt'] = 0.
param_const['p_init'] = 101325.
param_range['rh'] = [0.05, 0.95]
param_const['wall_on'] = 0

# check this input
#param_const['const_comp'] = 'APINENE, BENZENE, CH4, CO, NO2, SO2'
param_const['const_comp'] = 'APINENE, CO, NO2, SO2'

#param_const['Comp0'] = 'APINENE, BENZENE, CH4, CO, NO2, SO2'
param_const['Comp0'] = 'APINENE, CO, NO2, SO2'
#param_range['C0'] = [[0., 1.e1], [0., 2.e1], [3.e2, 2.e3], [0., 2.e3], [0., 2.e2], [0., 1.e1]]
param_range['C0'] = [[0., 1.e1], [0., 2.e3], [0., 2.e2], [0., 1.e1]]
param_const['number_size_bins'] = 7
param_const['space_mode'] = 'log'
param_const['coag_on'] = 1
param_const['pconct'] = 0.
param_const['pconc'] = 1.e4
param_const['pcont'] = 0

# setup dictionary items to hold chosen values
param_const['trans_fac'] = 0.
param_const['temperature'] = 0.
param_const['rh'] = 0.

# save ranges inside constants, for use inside gui
param_const['param_ranges'] = param_range

def auto_setup_and_call(param_const): # define the function

	# note, to test that steady state being reached, 
	# it is recommended to run just two simulations to 
	# begin and checking on these whether steady state reached
		
	# pass chosen parameters to PyCHAM for running the simulation ------
	# ensure gui module can be seen
	sys.path.append(str(os.getcwd() + '/PyCHAM'))

	import gui
	app = QApplication(sys.argv)

	ex = gui.PyCHAM(param_const)
	
	sys.exit() # fully disconnect from QApplication
	
	return() # end call to this module
	# ------------------------------------------------------------------
		

# call the function
auto_setup_and_call(param_const)