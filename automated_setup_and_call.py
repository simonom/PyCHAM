##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2023 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
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
'''code to automatically setup a series of simulations 
and run PyCHAM without using the PyCHAM graphical user interface'''
# useful for users that want to setup and run simulations in 
# an automated manner, e.g. because manual selection 
# of simulations via the graphical user interface is not practical 


# import required modules
import numpy as np
import sys
import os
from PyQt5.QtWidgets import QApplication, QWidget
import ast  # converting path text to list


# dictionary to hold parameter possibilities
param_range = {}
# dictionary to hold parameters that will be held constant
param_const = {}

# number of simulations (in which case provide a number), or whether to work through given values (in which case state 'set')
param_const['sim_num'] = 'set'
# either 'starter' or 'finisher'
param_const['sim_type'] = 'finisher'

# markers for interpreting chemical scheme - note that all chemical scheme markers should be 
# included inside quotation marks, e.g.: 
# param_const['chem_scheme_markers'] = '{, RO2, +, C(ind_, ), , &, , , :, }, ;,'
param_const['chem_scheme_markers'] = '{, RO2, +, C(ind_, ), , &, , , :, }, ;,'
param_const['pars_skip'] = 0 # need to parse equations and estimate properties on first go

# state path to chemical scheme and xml files
if sys.platform == 'win32':
	param_const['sch_name'] = 'C:\\Users\\Psymo\\Desktop\\PyCHAM\\PyCHAM\\PyCHAM\\input\\auto_call_test\\AP_BZ_MCM_PRAMAP_autoAPRAMBZ_scheme.dat'
	param_const['xml_name'] = 'C:\\Users\\Psymo\\Desktop\\PyCHAM\\PyCHAM\\PyCHAM\\input\\auto_call_test\\MCM_PRAM_xml.xml'

if sys.platform == 'darwin':
	param_const['sch_name'] = '/Users/user/Documents/GitHub/PyCHAM/PyCHAM/input/auto_call_test/AP_BZ_MCM_PRAMAP_autoAPRAMBZ_scheme.dat'
	param_const['xml_name'] = '/Users/user/Documents/GitHub/PyCHAM/PyCHAM/input/auto_call_test/MCM_PRAM_xml.xml'

# state parameter ranges

if (param_const['sim_type'] == 'starter'): # 24 hours (8.64e4s) spin up
	# note that results file name for finisher runs set below
	param_const['res_file_name'] = 'ambient_run_num'
	param_const['total_model_time'] = 8.64e4
	param_const['update_step'] = 6.e2
	param_const['recording_time_step'] = 3.6e3
	param_const['light_status'] = 1 # set to 2 for sinusoidal light, 1 for constant light

# at 38o North (Athens, Greece), http://www.csgnetwork.com/degreelenllavcalc.html says that
# 1o longitude is 9e4 m and 1o latitude is 1e5 m, 
# and file://data.cas.manchester.ac.uk/database6/Chamber/2021_Got-Jul-Man/Model/EMEP/EMEP_factsheet_Feb2020.pdf 
# says that EMEP has a lon-lat grid resolution of 0.125ox0.0625o, which for Athens,
# works out to 1.1250e4mx6.250e3m lon-lat. Taking the greater distance, then taking
# the European average wind speed of 3 m/s, the average time in an EMEP grid cell is:
# (1.125e4m/3m/s) = 3.75e3 s (1.04 hours)
if (param_const['sim_type'] == 'finisher'):
	param_const['total_model_time'] = 4.32e4
	param_const['update_step'] = 4.5e2
	param_const['recording_time_step'] = 4.5e2
	param_const['light_status'] = 1 # constant light intensity

param_const['light_time'] = 0.

# 12km altitude cloud-free solar actinic flux spectrum from Greece in June: https://doi.org/10.1029/2001JD900142
param_const['act_flux_path'] = 'Greece_obs_doi_10dot10292001JD900142.csv'
# affects all wavelengths, note that we want 0 so that darkness is represented
param_range['trans_fac'] = [0.0, 0.5, 1.]

# minimum temperature and relative humidity ranges given by Porter et al. 2021 (doi.org/10.1021/acsearthspacechem.1c00090)
param_range['temperature'] = [273.15, 293.15, 313.15]
param_const['tempt'] = 0.
param_const['p_init'] = 101325.
param_const['rh'] = 0.50
# allow wall to be on so that particles can be lost as they would be in the real atmosphere
param_const['wall_on'] = 0
#param_const['McMurry_flag'] = 1
#param_const['ChamSA'] = 1
# allow gas-wall partitioning
#param_const['eff_abs_wall_massC'] = 0.0
#param_const['mass_trans_coeff'] = 1.e-4
param_const['tracked_comp'] = 'NO2, APINENE, BENZENE, O3, NO, OH, HO2, CH3O2, HNO3'

# for the SOAPRA project, benzene is used as a proxy for the OH-reactivity-equivalent of: benzene, toluene, ethylbenzene, xylene (BTEX)
# the minimum concentration for benzene is then 0., and the maximum is based on the maximum observed concentrations of these components:
# benzene: 300 ppb (doi.org/10.1016/S0048-9697(01)00986-X); toluene 1000 ppb (doi.org/10.1016/S0048-9697(01)00986-X); 
# ethylbenzene 200 ppb (since ethyl-benzene is convincingly 1/5th the value of xylene in doi.org/10.1016/j.scitotenv.2022.158873, and the 
# xylene maximum is cited next in this comment); xylene 1000 ppb (doi.org/10.1016/S0048-9697(01)00986-X)
# OH reactivity of BTEX components according to MCM v3.3.1:
# benzene: 2.3D-12*EXP(-190/TEMP)*0.352
# toluene: 1.8D-12*EXP(340/TEMP)*0.07
# ethylbenzene: 7.00D-12*0.07
# xylene: 1.36D-11*0.55
# so, assuming the maximum temperature allowed here, we can get the total OH reactivity of the maximum concentration of 
# all BTEX components:
# first find the factor for converting ppb into # molecules/cm3 assuming a pressure of 101325 Pa (R has units cm3.Pa/K.mol) (2.3e10)
import scipy.constants as si
Cfac = 101325.*(si.N_A/((si.R*1.e6)*max(param_range["temperature"])))*1e-9
# reactivity with OH: 
OHreac = (2.3e-12*np.exp(-190./max(param_range["temperature"]))*0.352*(300.*Cfac)) + (1.8e-12*np.exp(340./max(param_range["temperature"]))*0.07*(1000.*Cfac)) + (7.00e-12*0.07*(200.*Cfac)) + (1.36e-11*0.55*(1000.*Cfac))

# benzene concentration (ppb) equivalent:
benzC = OHreac/(2.3e-12*np.exp(-190/max(param_range["temperature"]))*0.352)/Cfac

# the components to keep at constant concentration throughout the simulation
#param_const['const_comp'] = 'CH4, CO'

# the range of starting concentrations for starter simulations 
if (param_const['sim_type'] == 'starter'):
	# the components present at the start of the simulation, note do not leave white space
	param_const['Comp0'] = 'APINENE,BENZENE,CH4,CO,NO2,NO,SO2'
	param_range['C0'] = [[0.], [0.], [3.e2,1.*10**2.9,2.e3], [1.5,1.*10.**1.2,1.5e2], [1.e-2,1.e0,1.e2], [1.e-2,1.e0,1.e2], [0.]]

# if finishing a simulation, we want to use the end 
# result of the starter simulation as initial concentration
if (param_const['sim_type'] == 'finisher'):
	
	# look up the available starter simulations
	# get current working directory (in standard PyCHAM, 
	# this is the PyCHAM home directory)
	cwd = os.getcwd()
	# path to automated run results
	if sys.platform == 'win32': # windows
		init_conc_path = str('C:\\Users\\Psymo\\OneDrive - The University of Manchester\\PyCHAM\\outputs\\interact\\AP_BZ_MCM_PRAMAP_autoAPRAMBZ_scheme\\')
	if sys.platform == 'darwin': # mac
		init_conc_path = str('/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/PyCHAM/outputs/interact/AP_BZ_MCM_PRAMAP_autoAPRAMBZ_scheme/')
	# prepare to hold required variables
	y_arrays = []
	starter_path_name = []
	temp_arrays = []
	j_arrays = []
	rh_arrays = []
	press_arrays = []
	#N_perbin_arrays = []
	
	# get directories in this top-level directory
	dirlist = [item for item in os.listdir(init_conc_path)]
	
	for path in dirlist: # loop through objects in this directory
		
		if (path[-9::] == '.DS_Store' or path[-7::] == 'lab_sim' or path[-9::] == 'filter.py'): # ignore any irrelevant folders 	
			continue
		# extract initial concentrations in gas-phase (ppb)
		# withdraw concentrations (ppb in gas, # molecules/cm3 in particle and wall)
		fname = str(init_conc_path + path + '/concentrations_all_components_all_times_gas_particle_wall')
		# note, just keep results from final time
		
		y = (np.loadtxt(fname, delimiter=',', skiprows=1))[-1, :]

		# convert gas-phase concentration from ppb to # molecules/cm3
		fname = str(init_conc_path + '/' + path + '/model_and_component_constants')
		const_in = open(fname)
		for line in const_in.readlines():
	
			if str(line.split(',')[0]) == 'factor_for_multiplying_ppb_to_get_molec/cm3_with_time':

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

				# conversion factor to change gas-phase concentrations from # molecules/cm3 
				# (air) into ppb
				Cfactor = (ast.literal_eval(line[st_indx:fi_indx]))[-1]

			# get number of components
			for i in line.split(',')[1::]:
				if str(line.split(',')[0]) == 'number_of_components':
					num_comp = (int(i))
				else:
					continue
		
		# withdraw final time step number-size distributions (# particles/cm3 (air))
		#fname = str(init_conc_path + '/' + path + '/particle_number_concentration_dry')
		#N = (np.loadtxt(fname, delimiter=',', skiprows=1))[-1, :]

		# convert gas-phase concentration from ppb to # molecules/cm3
		y[0:num_comp] = y[0:num_comp]*Cfactor

		y_arrays.append(y) # store this concentration			
		starter_path_name.append(path) # store path to starter file
		#N_perbin_arrays.append(N)

		# withdraw environmental conditions (s)
		fname = str(init_conc_path + '/' + path + '/chamber_environmental_conditions')
		cham_env = np.loadtxt(fname, delimiter=',', skiprows=1)
		temp_arrays.append(cham_env[-1, 0]) # store temperature (K)
		press_arrays.append(cham_env[-1, 1]) # store pressure (Pa)
		rh_arrays.append(cham_env[-1, 2]) # store relative humidity (0-1)
		j_arrays.append(cham_env[-1, 3])

	# store required outputs
	param_range['ys'] = y_arrays
	param_range['starter_paths'] = starter_path_name
	#param_range['Ns'] = N_perbin_arrays
	param_range['temps'] = temp_arrays
	param_range['js'] = j_arrays
	param_range['rhs'] = rh_arrays
	param_range['pressures'] = press_arrays
# the components with constant influx, note, do not leave whitespace
param_const['const_infl'] = 'APINENE,BENZENE,CH4,CO,NO2,NO,SO2,O3'

# range of concentration of components present throughout simulation  - see above in this module for provenance of benzene,
# a maximum of alpha-pinene of 10 ppb is given by doi.org/10.1016/j.scitotenv.2020.144129, the range in CH4 is from a minimum of 400 ppb,
# which is from ice-core data (https://data.ess-dive.lbl.gov/view/doi:10.3334/CDIAC/ATG.030) and a maximum of 2000 ppb (2 ppm), which is
# from NOAA (https://gml.noaa.gov/ccgg/trends_ch4/). for carbon monoxide minimum is from doi.org/10.3402/tellusb.v50i3.16101, maximum is from https://scied.ucar.edu/learning-zone/air-quality/carbon-monoxide and https://earthobservatory.nasa.gov/global-maps/MOP_CO_M, for NO2, and the maximum is from doi.org/10.1007/s41810-023-00175-8, which sees a maximum NOx of 150 ug/m3 in urban India, which equates to 150*1e-12/32g/mol*si.N_A/Cfac  = 124 ppb, whilst the minimum for NOx is likely below the detection limit of instruments, as indicated by this paper: 10.5194/acp-22-12025-2022. For SO2, the minimum is from this paper: doi.org/10.1007/s10874-011-9185-2, maximum from Fig. 4 of doi.org/10.1016/j.partic.2012.09.005. For O3, we apply a constant influx to match the constant influx observed by measurements during Boreal winter in the Arctic (which gives around 30 ppb in Alert during winter in Figure 2A of doi.org/10.1016/j.atmosenv.2006.09.053), even though there is no photochemical production in winter, therefore background influx is the only source. Note that the loJ simulations for starter runs show that the influx stated below are sufficient to reproduce the observed O3 in high latitude Boreal winter. 
# note that we divide the absolute concentration (ppb) ranges by the time influx occurs over to get the emission rate (ppb/s)
#param_range['Cinfl'] = [[1.e-4/3.6e3, 1.e1/3.6e3], [1.e-4/3.6e3, benzC/3.6e3], [4.e2/3.6e3, 2.e3/3.6e3], [4.e1/3.6e3, 2.e4/3.6e3], [5.e-5/3.6e3, 1.e2/3.6e3], [5.e-5/3.6e3, 1.e2/3.6e3], [1.e-1/3.6e3, 1.e2/3.6e3]]


param_range['Cinfl'] = [[1.e-6, 5.7e-5, 3.3e-3], [benzC*5.e-8, benzC*4.e-6, benzC*3.e-4], [1.e-7, 5.e-7, 1.e-6], [1.e-12, 1.e-11, 1.e-10], [1.e-6, 1.*10.**-4.3, 5.e-4], [1.e-6, 1.*10.**-4.3, 5.e-4], [0.], [5.e-5, 9.e-5, 3.e-4]] 

# time over which influx of components occurs
if (param_const['sim_type'] == 'starter'):
	param_const['const_infl_t'] = '0.'
	param_const['pconct'] = '0.1; 3.6e3'
	param_const['number_size_bins'] = 0
if (param_const['sim_type'] == 'finisher'):
	param_const['const_infl_t'] = '0.'
	#param_const['pconct'] = '0.1; 3.75e3'
	#param_const['number_size_bins'] = 3

#param_const['space_mode'] = 'log'
#param_const['coag_on'] = 1

# minimum particulate mass concentration from doi.org/10.1021/acsearthspacechem.1c00090, maximum from https://indianexpress.com/article/cities/delhi/delhi-pm2-5-pm10-levels-shoot-through-the-roof-morning-after-diwali-7608039/
#param_range['pconc'] = [1.e-1, 1.e2]
#param_const['pcont'] = '1; 0' # ensure continuous influx

# setup dictionary items to hold chosen values
param_const['trans_fac'] = 0.
param_const['temperature'] = 0.
param_const['Cinfl'] = 0.

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
	print('the automated_setup_and_call module has completed its execution')
	return() # end call to this module
	# ------------------------------------------------------------------
		

# call the function
auto_setup_and_call(param_const)
