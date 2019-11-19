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
		# set path prefix based on whether this is a test or not (i.e. whether test flag 
		# file exists)
		dirpath = os.getcwd() # get current path
		
		fname = self.openFileNameDialog() # ask for location of input file

		with open(dirpath+'/fname.txt','w') as f:
			f.write(fname)
		f.close()
	
	@pyqtSlot()
	def on_click2(self):
		# set path prefix based on whether this is a test or not (i.e. whether test flag 
		# file exists)
		dirpath = os.getcwd() # get current path
		
		xmlname = self.openFileNameDialog()
		
		with open(dirpath+'/xmlname.txt','w') as f:
			f.write(xmlname)
		f.close()
		
	@pyqtSlot()
	def on_click3(self):
		inname = self.openFileNameDialog()
		
		# open the file
		inputs = open(inname, mode='r')
		
		# read the file and store everything into a list
		in_list = inputs.readlines()
		inputs.close()
		
		if len(in_list) != 36:
			print('Error: The number of model inputs is incorrect, should be 36, but is ' + str(len(in_list)) )
			sys.exit()
		for i in range(len(in_list)):
			key, value = in_list[i].split('=')
			key = key.strip() # a string
			if key == 'Res_file_name':
				resfname = str(value.strip())
			if key == 'Total_model_time':
				end_sim_time = float(value.strip())
			if key == 'Time_step':
				tstep_len = float(value.strip())
			if key == 'Recording_time_step':
				save_step = float(value.strip())
			if key == 'Number_size_bins':
				num_sb = int(value.strip())
			if key == 'lower_part_size':
				lowersize = float(value.strip())
			if key == 'upper_part_size':
				uppersize = float(value.strip())
			if key == 'mass_trans_coeff':
				kgwt = float(value.strip())
			if key == 'eff_abs_wall_massC':
				Cw = float(value.strip())
			if key == 'Temperature':
				TEMP = float(value.strip())
			if key == 'PInit':
				PInit = float(value.strip())
			if key == 'RH':
				RH = float(value.strip())
			if key == 'lat':
				lat = float(value.strip())
			if key == 'lon':
				lon = float(value.strip())
			if key == 'daytime_start':
				dt_start = float(value.strip())
			if key == 'ChamSA':
				ChamSA = float(value.strip())
			if key == 'nucv1':
				nucv1 = float(value.strip())
			if key == 'nucv2':
				nucv2 = float(value.strip())
			if key == 'nucv3':
				nucv3 = float(value.strip())
			if key == 'nuc_comp':
				nuc_comp = int(value.strip())
			if key == 'inflectDp':
				inflectDp = float(value.strip())
			if key == 'Grad_pre_inflect':
				pwl_xpre = float(value.strip())
			if key == 'Grad_post_inflect':
				pwl_xpro = float(value.strip())
			if key == 'Kern_at_inflect':
				inflectk = float(value.strip())
			if key == 'Rader_flag':
				Rader = int(value.strip())
			if key == 'C0':
				C0 = [float(i) for i in (value.split(','))]
			if key == 'Comp0': # strip removes white space
				Comp0 = [str(i).strip() for i in (value.split(','))]
			if key == 'voli':
				
				if value.split(',')==['\n']:
					voli = np.empty(0)
				else:
					voli = [int(i) for i in (value.split(','))]	
			if key == 'volP':
				if value.split(',')==['\n']:
					volP = np.empty(0)
				else:
					volP = [float(i) for i in (value.split(','))]
			if key == 'pconc':
				pconc = float(value.strip())
			if key == 'std':
				std = float(value.strip())
			if key == 'loc':
				loc = float(value.strip())
			if key == 'scale':
				scale = float(value.strip())
			if key == 'core_diss':
				core_diss = float(value.strip())
			if key == 'light_time':
				light_time = [float(i) for i in (value.split(','))]
			if key == 'light_stat':
				light_stat = [int(i) for i in (value.split(','))]
		
		# can't have a recording time step (s) less than the ode solver time step (s),
		# so force to be equal and tell user
		if save_step<tstep_len:
			print('Recording time step cannot be less than the ode solver time step, so increasing recording time step to be equal to input ode solver time step')
			save_step = tstep_len
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
		resfname, tstep_len, TEMP, PInit, RH, lat, lon, dt_start, Cw,  
		save_step, ChamSA, nucv1, nucv2, nucv3, nuc_comp,   
		inflectDp, pwl_xpre, pwl_xpro, inflectk, Rader, xmlname, C0, Comp0, 
		voli, volP, pconc, std, loc, scale, core_diss, light_stat, light_time,
		kgwt]
			
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