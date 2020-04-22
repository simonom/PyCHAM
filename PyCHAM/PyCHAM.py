'''The module that generates the Graphical User Interface for PyCHAM, and connects that GUI with the core PyCHAM model'''

import sys
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QFileDialog, QLabel
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot, Qt
import pickle # for storing inputs
import threading # for calling on further PyCHAM functions
import os
import numpy as np

class PyCHAM(QWidget):

	def __init__(self):
		super().__init__()
		self.title = 'PyCHAM'
		self.left = 10
		self.top = 10
		self.width = 450
		self.height = 300
		self.initUI()
    
	def initUI(self):
		self.setWindowTitle(self.title)
		self.setGeometry(self.left, self.top, self.width, self.height)
		
		label = QLabel(self)
		label.setText("Welcome to PyCHAM.  Please see the README file for guidance.")
		label.move(20,35)
		label.show()
		
		
		button = QPushButton('Chemical Scheme .txt File', self)
		button.setToolTip('Select the .txt file containing the desired chemical scheme')
		button.move(100,70)
		button.clicked.connect(self.on_click1)
		
		button = QPushButton('Chemical Scheme .xml File', self)
		button.setToolTip('Select the .xml file containing the desired conversion file')
		button.move(100,100)
		button.clicked.connect(self.on_click2)
		
		button = QPushButton('Model Variables .txt File', self)
		button.setToolTip('Select the desired file containing the model variables')
		button.move(100,130)
		button.clicked.connect(self.on_click3)
		
		button = QPushButton('Run Model', self)
		button.setToolTip('Start the simulation')
		button.move(100,160)
		button.clicked.connect(self.on_click4)
		
		button = QPushButton('Plot Results', self)
		button.setToolTip('Plot output from simulation')
		button.move(100,190)
		button.clicked.connect(self.on_click5)
		
		self.show()
		

	@pyqtSlot()
	def on_click1(self):
		
		dirpath = os.getcwd() # get current path
		fname = dirpath+'/PyCHAM/inputs/limonene_MCM_PRAM.txt' # hard-code chemical scheme input
# 		fname = self.openFileNameDialog() # ask for location of input chemical scheme file
		with open(dirpath+'/fname.txt','w') as f:
			f.write(fname)
		f.close()
	
	@pyqtSlot()
	def on_click2(self):
		
		dirpath = os.getcwd() # get current path
		xmlname = dirpath+'/PyCHAM/inputs/Example_Run_xml.xml' # hard-code xml input
# 		xmlname = self.openFileNameDialog()
		with open(dirpath+'/xmlname.txt','w') as f:
			f.write(xmlname)
		f.close()
		
	@pyqtSlot()
	def on_click3(self):
		dirpath = os.getcwd() # get current path
		inname = dirpath+'/PyCHAM/inputs/limonene_inputs.txt' # hard-code model variables input
# 		inname = self.openFileNameDialog() # name of model variables inputs file
		
		# open the file
		inputs = open(inname, mode='r')
		
		# read the file and store everything into a list
		in_list = inputs.readlines()
		inputs.close()
		
		# check on whether correct number of inputs supplied
		input_len = 60
		if len(in_list) != input_len:
			print(('Error: The number of variables in the model variables file is incorrect, should be ' + str(input_len) + ', but is ' + str(len(in_list))))
			sys.exit()
			
		for i in range(len(in_list)):
			key, value = in_list[i].split('=')
			key = key.strip() # a string with bounding white space removed
			if key == 'Res_file_name':
				if value.split(',')==['\n']:
					print('Error: no requested results file name detected in model inputs file, please supply')
					sys.exit()
				resfname = str(value.strip())
			if key == 'Total_model_time':
				if value.split(',')==['\n']:
					print('Error: no total model run time detected in model inputs file, please supply')
					sys.exit()
				else:
					end_sim_time = float(value.strip())
			if key == 'Time_step':
				if value.split(',')==['\n']:
					print('Notice: No model time step detected in model inputs file, defaulting to 60s')
					tstep_len = float(60.0)
				else:
					tstep_len = float(value.strip())
			if key == 'Recording_time_step':
				if value.split(',')==['\n']:
					print('Notice: no recording time step detected in model inputs file, defaulting to 60s')
					save_step = float(60.0)
				else:
					save_step = float(value.strip())
			if key == 'Number_size_bins':
				if value.split(',')==['\n']:
					print('Notice: no number of size bins detected in model inputs, defaulting to zero')
					num_sb = int(0)
				else:
					num_sb = int(value.strip())
			if key == 'lower_part_size':
				if value.split(',')==['\n']:
					lowersize = float(0.0)
				else:
					lowersize = float(value.strip())
			if key == 'upper_part_size':
				if value.split(',')==['\n']:
					uppersize = float(0.0)
				else:
					uppersize = float(value.strip())
			if key == 'space_mode':
				if value.split(',')==['\n']:
					space_mode = str('none')
				else:
					space_mode = str(value.strip())
			if key == 'mass_trans_coeff':
				if value.split(',')==['\n']:
					kgwt = float(0.0)
				else:
					kgwt = float(value.strip())
			if key == 'eff_abs_wall_massC':
				if value.split(',')==['\n']:
					Cw = float(0.0)
				else:
					Cw = float(value.strip())
			if key == 'Temperature':
				if value.split(',')==['\n']:
					print('Error: no air temperature detected in model inputs file')
					sys.exit()
				else:
					TEMP = float(value.strip())
			if key == 'PInit':
				if value.split(',')==['\n']:
					print('Error: no air pressure detected in model inputs file')
					sys.exit()
				else:
					PInit = float(value.strip())
			if key == 'RH':
				if value.split(',')==['\n']:
					print('Error: no relative humidity detected in model variables input file')
					sys.exit()
				else:
					RH = float(value.strip())
					if RH>1.0:
						print('Note: RH set above 1.0, code will attempt to run but please note that RH is interpreted as fraction, not a percentage, where 1 represents saturation of water vapour')
			if key == 'lat':
				if value.split(',')==['\n']:
					lat = float(0.0)
				else:
					lat = float(value.strip())
			if key == 'lon':
				if value.split(',')==['\n']:
					lon = float(0.0)
				else:
					lon = float(value.strip())
			if key == 'DayOfYear':
				if value.split(',')==['\n']:
					DayOfYear = int(1)
				else:
					DayOfYear = int(value.strip())		
			if key == 'daytime_start': # for outdoor actinic flux equation
				if value.split(',')==['\n']:
					dt_start = float(0.0)
				else:
					dt_start = float(value.strip())
			if key == 'act_flux_file': # for indoor actinic flux
				cwd = os.getcwd() # address of current working directory
				if (value.strip()).split(',')==['']:
					act_flux_path = str('no')
				else:
					act_flux_path = str(cwd + '/PyCHAM/photofiles/' + value.strip())
			# name of file with wavelength-dependent absorption cross-section and quantum yield 
			# calculations
			if key == 'photo_par_file':
				cwd = os.getcwd() # address of current working directory
				if (value.strip()).split(',')==['']:
					photo_par_file = str(cwd + '/PyCHAM/photofiles/' + 'MCMv3.2')
				else:
					photo_par_file = str(cwd + '/PyCHAM/photofiles/' + value.strip())
			if key == 'ChamSA': # chamber surface area used for particle loss to walls
				if value.split(',')==['\n']:
					ChamSA = float(0.0)
				else:
					ChamSA = float(value.strip())
			if key == 'nucv1': # first parameter in the nucleation equation
				if value.split(',')==['\n']:
					nucv1 = float(0.0)
				else:
					nucv1 = float(value.strip())
			if key == 'nucv2':
				if value.split(',')==['\n']:
					nucv2 = float(0.0)
				else:
					nucv2 = float(value.strip())
			if key == 'nucv3':
				if value.split(',')==['\n']:
					nucv3 = float(0.0)
				else:
					nucv3 = float(value.strip())
			if key == 'nuc_comp': # name of nucleating component
				if (value.strip()).split(',')==['']:
					nuc_comp = [] # empty list
				else:
					nuc_comp = [str(i).strip() for i in (value.split(','))]
			if key == 'new_partr': # radius of newly nucleated particles (cm)
				if (value.strip()).split(',')==['']:
					new_partr = float(2.0e-7)
				else:
					new_partr = float(value.strip())
			if key == 'inflectDp':
				if (value.strip()).split(',')==['']:
					inflectDp = float(0.0)
				else:
					inflectDp = float(value.strip())
			if key == 'Grad_pre_inflect':
				if (value.strip()).split(',')==['']:
					pwl_xpre = float(0.0)
				else:
					pwl_xpre = float(value.strip())
			if key == 'Grad_post_inflect':
				if (value.strip()).split(',')==['']:
					pwl_xpro = float(0.0)
				else:
					pwl_xpro = float(value.strip())
			if key == 'Rate_at_inflect':
				if (value.strip()).split(',')==['']:
					inflectk = float(0.0)
				else:
					inflectk = float(value.strip())
			if key == 'part_charge_num':
				if (value.strip()).split(',')==['']:
					p_char = -1.0e-6
				else:
					p_char = float(value.strip())
			if key == 'elec_field':
				if (value.strip()).split(',')==['']:
					e_field = -1.0e6
				else:
					e_field = float(value.strip())
			
			if key == 'McMurry_flag':
				if value.split(',')==['\n']:
					Rader = int(-1)
				else:
					Rader = int(value.strip())
					# exit if no chamber surface area supplied and McMurry and Rader 
					# (1985) asked for
					if Rader == 1 and ChamSA == 0.0:
						print('Error: inside model variables input file, the McMurry flag is set to 1 indicating that the McMurry and Rader (1985) model for particle wall deposition be used, however the chamber surface area (ChamSA) is set to zero.  Please either increase surface area to use the McMurry and Rader model, or set the McMurry_flag to -1 to turn off particle deposition to walls')
						sys.exit()
					if Rader == 1 and p_char == -1.0e6:
						print('Error: inside model variables input file, the McMurry flag is set to 1 indicating that the McMurry and Rader (1985) model for particle wall deposition be used, however the average number of charges per particle (part_charge_num) is empty.  Please either state average charge per particle to use the McMurry and Rader model, or set the McMurry_flag to -1 to turn off particle deposition to walls.  A setting of zero is acceptable and implies no electrostatic effects.')
						sys.exit()
					if Rader == 1 and e_field == -1.0e6:
						print('Error: inside model variables input file, the McMurry flag is set to 1 indicating that the McMurry and Rader (1985) model for particle wall deposition be used, however the average electric field inside the chamber (elec_field) is empty.  Please either state average electric field inside the chamber to use the McMurry and Rader model, or set the McMurry_flag to -1 to turn off particle deposition to walls.  A setting of zero is acceptable and implies no electrostatic effects.')
						sys.exit()
			if key == 'C0': # initial concentrations of components given in Comp0
				if (value.strip()).split(',')==['']:
					C0 = np.empty(0)
				else:
					C0 = [float(i) for i in (value.split(','))]
			if key == 'Comp0': # strip removes white space, names of components for C0
				if (value.strip()).split(',')==['']:
					Comp0 = [] # empty list
				else:
					Comp0 = [str(i).strip() for i in (value.split(','))]
			
			if key == 'injectt': # times of later injections (s)
				if (value.strip()).split(',')==['']:
					injectt = np.empty(0) # empty list
				else:
					injectt = [float(i.strip()) for i in (value.split(','))]
					injectt = np.array((injectt))
			
			if key == 'Compt': # names of components with later injections
				if (value.strip()).split(',')==['']:
					Compt = [] # empty list
				else:
					Compt = [str(i).strip() for i in (value.split(','))]
			
			if key == 'Ct': # concentration of component injections (ppb)
				if (value.strip()).split(',')==['']:
					Ct = [] # empty list
				else:
					Ct = np.zeros((len(Compt), len(injectt)))
					# check if more than one component being injected
					comp_count = 1
					# keep track of number of injection concentrations given for an 
					# individual component
					time_count = 1
					for i in value:
						if i==';':
							comp_count += 1 # record number of components
							# stop if number of submitted concentrations does not match number of injection times
							if time_count != len(injectt):
								print('Error noticed in Ct input inside Model Variables input file, number of injection concentrations given does not match number of injection times given by the injectt input, please see README for help')
								sys.exit()
							# reset time count
							time_count = 1
						if i==',':
							time_count += 1
					if time_count != len(injectt):
						print('Error noticed in Ct input inside Model Variables input file, number of injection concentrations given does not match number of injection times given by the injectt input, please see README for help')
						sys.exit()
					if comp_count!= len(Compt):
						print('Error noticed in Ct input inside Model Variables input file, number of component concentrations given does not match number of component given by the Compt input, please see README for help')	
						sys.exit()
					for i in range(comp_count):
						Ct[i, :] = [float(ii.strip()) for ii in ((value.split(';')[i]).split(','))]
			
			if key == 'const_comp': # names of components with later injections
				if (value.strip()).split(',')==['']:
					const_comp = [] # empty list
				else:
					const_comp = [str(i).strip() for i in (value.split(','))]
			if key == 'const_infl': # names of components with later injections
				if (value.strip()).split(',')==['']:
					const_infl = np.array(([])) # empty numpy array
				else:
					const_infl = [str(i).strip() for i in (value.split(','))]
					const_infl = np.squeeze(np.array(const_infl))
					
			if key == 'const_infl_t': # times of constant influxes (s)
				if (value.strip()).split(',')==['']:
					const_infl_t = np.empty(0) # empty list
				else:
					const_infl_t = [float(i.strip()) for i in (value.split(','))]
					const_infl_t = np.array((const_infl_t))
					
			if key == 'Cinfl': # constant influx concentrations of components (ppb/s)
				if (value.strip()).split(',')==['']:
					Cinfl = [] # empty list
				else: # fill constant influx rates
					# count on number of components with constant influx
					comp_count = 1
					# keep track of number of injection concentrations given for an 
					# individual component
					time_count = 1
					for i in value:
						if i==';':
							comp_count += 1 # record number of components
							time_count = 1 # reset time count
						if i==',':
							time_count += 1
					# components represented on rows, times in columns
					Cinfl = np.zeros((comp_count, time_count))
					for i in range(comp_count):
						Cinfl[i, :] = [float(ii.strip()) for ii in ((value.split(';')[i]).split(','))]
			
			
			if key == 'vol_Comp':
				if (value.strip()).split(',')==['']:
					vol_Comp = [] # empty list
				else:
					vol_Comp = [str(i).strip() for i in (value.split(','))]	
			if key == 'volP':
				if (value.strip()).split(',')==['']:
					volP = np.empty(0)
				else:
					volP = [float(i) for i in (value.split(','))]
			if key == 'act_wi':
				if (value.strip()).split(',')==['']:
					act_wi = np.empty(0)
				else:
					act_wi = [int(i) for i in (value.split(','))]	
			if key == 'act_w':
				if (value.strip()).split(',')==['']:
					act_w = np.empty(0)
				else:
					act_w = [float(i) for i in (value.split(','))]
					
			if key == 'pconct': # seed particles input times (s)
				if (value.strip()) == ['']:
					pconct = []
				else:
					# keep track of number of times given
					time_count = 1
					for i in value:
						if i==';':
							time_count += 1 # increase time count
					# times in columns
					pconct = np.zeros((1, time_count))
					pconct[0, :] = [float(i) for i in ((value.strip()).split(';'))]
					
			if key == 'pconc': # seed particle number concentrations (#/cc)
				if (value.strip()).split(',')==['']:
					pconc = []
				else:
					# keep track of number of times given
					time_count = 1
					# keep track of number of size bins given
					sb_count = 1
					for i in value:
						if i==';':
							time_count += 1 # increase time count
						if i==',':
							sb_count += 1 # increase size bin count
							
					# size bins in rows, times in columns
					pconc = np.zeros((sb_count, time_count))
					
					for i in range(time_count):
						pconc[:, i] = [float(ii.strip()) for ii in ((value.split(';')[i]).split(','))]
				
			if key == 'seed_name':
				if (value.strip()).split(',')==['']:
					seed_name = 'core'
				else:
					seed_name = str(value.strip())
			if key == 'seed_mw':
				if (value.strip()).split(',')==['']:
					seed_mw = 132.14
				else:
					seed_mw = float(value.strip())
			if key == 'seed_dens':
				if (value.strip()).split(',')==['']:
					seed_dens = 1.00
				else:
					seed_dens = float(value.strip())
					
			if key == 'mean_rad': # seed particle mean radius (um)
				if (value.strip()).split(',')==['']:
					mean_rad = np.zeros((1,1))
					mean_rad[0, 0] = -1.0e6
				else:
					# keep track of number of times given
					time_count = 1
					for i in value:
						if i==';':
							time_count += 1 # increase time count
					# times in columns
					mean_rad = np.zeros((1, time_count))
					mean_rad[0, :] = [float(i) for i in ((value.strip()).split(';'))]
			
			if key == 'std':
				if value.split(',')==['\n']:
					std = float(1.1)
				else:
					# keep track of number of times given
					time_count = 1
					for i in value:
						if i==';':
							time_count += 1 # increase time count
					# times in columns
					std = np.zeros((1, time_count))
					std[0, :] = [float(i) for i in ((value.strip()).split(';'))]
			
			if key == 'core_diss':
				if value.split(',')==['\n']:
					core_diss = float(1.0)
				else:
					core_diss = float(value.strip())
			if key == 'light_time':
				if (value.strip()).split(',')==['']:
					light_time = np.empty(0)
				else:
					light_time = [float(i) for i in (value.split(','))]
			if key == 'light_stat':
				if (value.strip()).split(',')==['']:
					light_stat = np.empty(0)
				else:
					light_stat = [int(i) for i in (value.split(','))]
				if sum(light_stat)== 0:
					print('Note, light_stat input variable (set in Model Variables .txt file) set to zero for entire simulation therefore photochemical reaction rates will be zero')
			# names of components whose tendency to change will be tracked
			if key == 'tracked_comp':
				if (value.strip()).split(',')==['']:
					dydt_trak = []
				else:
					dydt_trak = [str(i).strip() for i in (value.split(','))]
			# flag to say whether or not to clone latest version of UManSysProp from web
			# defaults to zero (no update)
			if key == 'umansysprop_update':
				if (value.strip()).split(',')==['']:
					umansysprop_update = int(0)
				else:
					umansysprop_update = int(i)
					
				# if no update requested, check that there is an existing UManSysProp
				# folder
				if umansysprop_update == 0:
					cwd = os.getcwd() # address of current working directory
					# check if there is an existing umansysprop folder
					if os.path.isdir(cwd + '/umansysprop'):
						continue
					else:
						print('Note, no download of UManSysProp requested by user via model variables input file, but no existing UManSysProp module found, so will try to update via internet')
						umansysprop_update = 1
				
				if umansysprop_update == 1: # test whether UManSysProp can be updated
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
						print('Internet connection confirmed and user has requested cloning of UManSysProp via the model variables input file') 
					else:
						print('Error: user has requested cloning of UManSysProp via the model variables input file but connection to the page failed, possibly due to no internet connection (UManSysProp repository site: https://github.com/loftytopping/UManSysProp_public.git)')
						sys.exit()
			if key == 'chem_scheme_markers': # formatting for chemical scheme
				if (value.strip()).split(',')==['']:
					# default to MCM FACSIMILE inputs
					chem_scheme_markers = ['* Reaction definitions. ;', '%', '(.*) End (.*)', '* Generic Rate Coefficients ;', ';', '\*\*\*\*', 'RO2', '+', '*;']
				else:
					chem_scheme_markers = [str(i).strip() for i in (value.split(','))]
			if key == 'int_tol': # tolerances for integration
				if (value.strip()).split(',')==['']:
					# default to minimum tolerances
					int_tol = np.array(([1.0e-3, 1.0e-4]))
				else:
					int_tol = np.array(([float(i) for i in (value.split(','))]))
			if key == 'dil_fac': # dilution factor rate
				if (value.strip()).split(',')==['']:
					dil_fac = float(0.0)
				else:
					dil_fac = float(value)
		# --------------------------------------------------------------------------------
		# checks on inputs
		
		# can't have a recording time step (s) less than the ode solver time step (s),
		# so force to be equal and tell user
		if save_step<tstep_len:
			print('Recording time step cannot be less than the ode solver time step, so increasing recording time step to be equal to input ode solver time step')
			save_step = tstep_len
		
		# if initial particle concentration is an array, its elements must align with the
		# number of size bins (pp_intro.py is where initial concentrations are sorted)
		if len(pconc)>1 and len(pconc)!=num_sb:
			print('If pconc (set in Model Variables .txt file) is an array, its length must equal the number of particle size bins, but currently length of pconc is '+str(len(pconc))+' and number of size bins is '+str(num_sb))
			sys.exit()
		
		# overide any particle number concentration inputs if number of size bins left empty or set to zero
		if num_sb == 0 and pconc != 0.0:
			print('Notice, since no size bins detected in Model Variables .txt file, pconc variable will be set to empty (even if values supplied)')
			pconc = []
		
		if max(pconct.shape) != pconc.shape[1]: # pconct should have dimensions 1, number of times pconc given for
			print('Error: number of times given for pconct input variable in model variables input file does not match number of times particle number concentration given for in the pconc input variable in the same file.  Please see README for guidance.')
			sys.exit()
		# components with constant influx for set periods of time
		if len(Cinfl)>0:
			if len(const_infl)!=(Cinfl.shape[0]):
				print('Error: the number of components given for constant influx by the const_infl variable inside the model variables input file does not match the number of components with constant influx concentrations provided by the Cinfl variable of that file, please see the README for guidance.')
				sys.exit()
			if len(const_infl_t)!=(Cinfl.shape[1]):
				print('Error: the number of times given for constant influx by the const_infl_t variable inside the model variables input file does not match the number of times with constant influx concentrations provided by the Cinfl variable of that file, please see the README for guidance.')
				sys.exit()
				
		if len(chem_scheme_markers)!=9:
			print('Error: length of chem_scheme_markers (specified in model variables input file) is not 8 and should be, please see README for guidance')
		# components with assigned vapour pressures
		if len(vol_Comp)!=len(volP):
			print('Error: the number of components with assigned vapour pressures does not equal the number of assigned vapour pressures (vol_Comp and volP variables, respectively, in the model variables input folder), please see the README for guidance')
			sys.exit()
		
		# nucleation inputs
		if nucv1>0.0 and nuc_comp == []:
			print('Error: the nucleation parameter nucv1 set in the model variables input folder is greater than zero, but no nucleating component is recognised from the nuc_comp variable, so nucleation cannot proceed.  Please see README for guidance.')
			sys.exit()
		
		# check actinic flux and photochemical parameter files in expected place
		if (photo_par_file != str(cwd + '/PyCHAM/photofiles/MCMv3.2')): # photochemical parameters
			
			if os.path.isfile(photo_par_file):
				pass
			else:
				print('Error: file given for the photo_par_file input in the model variables input file has not been found, please check README for guidance')
				sys.exit()
			
		if act_flux_path != 'no': # actinic flux file
			
			if os.path.isfile(act_flux_path):
				pass
			else:
				print('Error: path name given for the act_flux_path input in the model variables input file has not been found, please check README for guidance')
				sys.exit()
		
		# --------------------------------------------------------------------------------
		# get names of chemical scheme and xml files
		# set path prefix based on whether this is a test or not (i.e. whether test flag 
		# file exists)
		dirpath = os.getcwd() # get current path
			
		f = open(dirpath+'/fname.txt','r')
		content = f.readlines()
		f.close()
		fname = str(content[0])
		f = open(dirpath+'/xmlname.txt','r')
		content = f.readlines()
		f.close()
		xmlname = str(content[0])
		# remove temporary files
		os.remove(dirpath+'/fname.txt')
		os.remove(dirpath+'/xmlname.txt')
		
		# write variable values to pickle file
		list_vars = [fname, num_sb, lowersize, uppersize, end_sim_time, 
		resfname, tstep_len, TEMP, PInit, RH, lat, lon, DayOfYear, dt_start, 
		act_flux_path, Cw,  
		save_step, ChamSA, nucv1, nucv2, nucv3, nuc_comp, new_partr,   
		inflectDp, pwl_xpre, pwl_xpro, inflectk, Rader, xmlname, C0, Comp0, 
		vol_Comp, volP, pconc, std, mean_rad, core_diss, light_stat, light_time,
		kgwt, dydt_trak, space_mode, Ct, Compt, injectt, seed_name, const_comp,
		const_infl, Cinfl, act_wi, act_w, seed_mw, umansysprop_update, seed_dens, p_char, 
		e_field, const_infl_t, chem_scheme_markers, int_tol, photo_par_file, dil_fac, 
		pconct]
		
		if os.path.isfile(dirpath+'/testf.txt'):
			print('Model input buttons work successfully')
			with open('test_var_store.pkl','wb') as f:
				pickle.dump(list_vars,f)
				print('Pickle file dumped successfully')
				os.remove('test_var_store.pkl')
			f.close()
		else:
			with open('PyCHAM/var_store.pkl','wb') as f:
				pickle.dump(list_vars,f)
			f.close()
	
	@pyqtSlot()	
	def on_click4(self):
		import front as model
		dirpath = os.getcwd() # get current path
		if os.path.isfile(dirpath+'/testf.txt'):
			testf=2	
		else:
			testf=0
		# call on model to run
		t = threading.Thread(target=model.run(testf))
		t.daemon = False
		t.start()
	
	@pyqtSlot()
	def on_click5(self):
		import res_plot_super as plotter
		
		dirpath = os.getcwd() # get current path
		if os.path.isfile(dirpath+'/testf.txt'):
			testf=1		
		else:
			testf=0
		
		# pass the name of the folder where results are saved
		t = threading.Thread(target=plotter.run(testf))
		t.daemon = False
		t.start()
		
		if testf==1:
			# remove the test flag file
			# remove test file
			os.remove(dirpath+'/testf.txt')
			print('test_PyCHAM complete, please close the user interface')
		
		
	def openFileNameDialog(self):
		options = QFileDialog.Options()
		fname, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
		return fname