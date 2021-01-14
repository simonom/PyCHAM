'''The module that generates the Graphical User Interface for PyCHAM, and connects that GUI with the core PyCHAM model'''
# first module called when PyCHAM started from the terminal/command window, takes inputs
# and sends to model modules, also calls the saving module

from PyQt5.QtWidgets import *
from PyQt5.QtGui import  *
from PyQt5.QtCore import *
import pickle # for storing inputs
import sys
import os
import def_mod_var
import numpy as np

# class for scrollable label 
class ScrollLabel(QScrollArea): 
  
	def __init__(self, *args, **kwargs): 
		QScrollArea.__init__(self, *args, **kwargs) 
  
		# making widget resizable 
		self.setWidgetResizable(True) 
  
		# making qwidget object 
		content = QWidget(self) 
		self.setWidget(content) 
  
		# horizontal box layout 
		lay = QHBoxLayout(content) 
  
		# creating label 
		self.label = QLabel(content) 
  
		# setting alignment to the text 
		self.label.setAlignment(Qt.AlignLeft | Qt.AlignTop) 
  
		# making label multi-line 
		#self.label.setWordWrap(True) 
  
		# adding label to the layout 
		lay.addWidget(self.label)
		
	def setText(self, text): # the setText method
		self.label.setText(text) # setting text to the label

	def clear(self): # the clear method
		self.label.clear() # clearing text

class PyCHAM(QWidget):

	def __init__(self):
		super().__init__()
		self.title = 'PyCHAM'
		self.left = 10
		self.top = 10
		self.width = 800
		self.height = 530
		self.initUI()
		
		
		
		return
    
	def initUI(self):
	
		self.setWindowTitle(self.title)
		self.setGeometry(self.left, self.top, self.width, self.height)
		
		label = QLabel(self)
		label.setOpenExternalLinks(True)
		label.setText("Welcome to PyCHAM - please see <a href=\"http://www.github.com/simonom/PyCHAM\">README</a> for guidance")
		label.setFont(QFont('Arial', 20))
		label.move(150, 35)
		label.show()
		
		# horizontal line
		self.hline = QFrame()
		self.hline.setFrameShape(QFrame.HLine)
		self.hline.setFrameShadow(QFrame.Sunken)
		self.hline.setStyleSheet("background-color: rgb(200, 0, 0)")
		self.hline.setGeometry(100, 100, 100, 100)
		vbox = QVBoxLayout()
		vbox.addWidget(self.hline)
		self.setLayout(vbox)
		
	
		
		# default variables for all required input model variables 
		# stored to pickle file and output here
		[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O] = def_mod_var.def_mod_var(0)
		
		
		button = QPushButton('Folder Containing Input Files', self)
		button.setToolTip('Select the folder containing the required input files')
		button.move(20, 70)
		button.clicked.connect(self.on_click1)
		
		button = QPushButton(r'Check Model Variables', self)
		button.setToolTip('Ensure variables are suitable')
		button.move(350, 70)
		button.clicked.connect(self.on_click5)
		
		button = QPushButton(r'Update Model Variables', self)
		button.setToolTip('Reset model variables to those presented')
		button.move(350, 480)
		button.clicked.connect(self.on_click6)
		
		button = QPushButton(r'Run Model', self)
		button.setToolTip('Start the simulation')
		button.move(700, 420)
		button.clicked.connect(self.on_click7)
		
		button = QPushButton('Plot Results', self)
		button.setToolTip('Plot output from simulation')
		button.move(700, 450)
		button.clicked.connect(self.on_click8)

		button = QPushButton('Quit', self)
		button.setToolTip('Finish with PyCHAM and close this window')
		button.move(700, 480)		
		button.clicked.connect(self.on_click9)
		
		# file selection -------------------------------------------------------------
		
		l2 = QLabel(self)
		l2.setText("The following files have been found: ")
		l2.move(20, 100)
		l2.show()
		
		# creating scroll label
		l3a = QLabel(self)
		l3a.setText("Chemical scheme file: ")
		l3a.move(20, 130)
		l3a.show()
			
		self.l3 = ScrollLabel(self)
		self.l3.setText(sch_name)
		self.l3.setGeometry(20, 150, 300, 50)
		self.l3.show()
			
		button1 = QPushButton('Select new file', self)
		button1.setToolTip('Select the file containing the required chemical scheme file')
		button1.move(20, 195)
		button1.clicked.connect(self.on_click2)
		button1.show()
			
		l4a = QLabel(self)
		l4a.setText("xml file: ")
		l4a.move(20, 230)
		l4a.show()
			
		self.l4 = ScrollLabel(self)
		self.l4.setText(xml_name)
		self.l4.setGeometry(20, 250, 300, 50)
		self.l4.show()
			
		button2 = QPushButton('Select new file', self)
		button2.setToolTip('Select the file containing the required xml file')
		button2.move(20, 295)
		button2.clicked.connect(self.on_click3)
		button2.show()
			
		l5a = QLabel(self)
		l5a.setText("Model variables file: ")
		l5a.move(20, 330)
		l5a.show()
			
		self.l5 = ScrollLabel(self)
		self.l5.setText(inname)
		self.l5.setGeometry(20, 350, 300, 50)
		self.l5.show()
			
		button3 = QPushButton('Select new file', self)
		button3.setToolTip('Select the file containing the required model variables file')
		button3.move(20, 395)
		button3.clicked.connect(self.on_click4)
		button3.show()
		
		# model variables ----------------------------------------------------------
		
		l6 = QLabel(self)
		l6.setText("The following model variables have been found: ")
		l6.move(350, 100)
		l6.show()
			
		l7 = QLabel(self)
		l7.setText('Folder to save to: ')
		l7.move(350, 130)
		l7.show()
		self.e7 = QLineEdit(self)
		self.e7.setText(sav_nam)
		self.e7.move(460, 130)
		self.e7.show()
		
		l8 = QLabel(self)
		l8.setText('Chemical scheme markers: ')
		l8.move(350, 160)
		l8.show()
		self.e8 = QLineEdit(self)
		self.e8.setText((str(chem_sch_mark)).replace('\'', '').replace(' ', '')[1:-1])
		self.e8.move(520, 160)
		self.e8.show()
		
		l9 = QLabel(self)
		l9.setText('Total experiment time (s): ')
		l9.move(350, 190)
		l9.show()
		self.e9 = QLineEdit(self)
		self.e9.setText((str(tot_time)).replace('\'', '').replace(' ', ''))
		self.e9.move(520, 190)
		self.e9.show()
		
		l10 = QLabel(self)
		l10.setText('Update time interval (s): ')
		l10.move(350, 220)
		l10.show()
		self.e10 = QLineEdit(self)
		self.e10.setText((str(update_stp)).replace('\'', '').replace(' ', ''))
		self.e10.move(520, 220)
		self.e10.show()
		
		l11 = QLabel(self)
		l11.setText('Recording time interval (s): ')
		l11.move(350, 250)
		l11.show()
		self.e11 = QLineEdit(self)
		self.e11.setText((str(save_step)).replace('\'', '').replace(' ', ''))
		self.e11.move(520, 250)
		self.e11.show()
		
		
		self.show()
		return
		

	@pyqtSlot()
	def on_click1(self):
	
	
		@pyqtSlot()
		def click1_up(self): # update labels following selection of directory
		
			# creating scroll label
			self.l3 = ScrollLabel(self)
			self.l3.setText(sch_name)
			self.l3.setGeometry(20, 150, 300, 50)
			self.l3.show()
			
			self.l4 = ScrollLabel(self)
			self.l4.setText(xml_name)
			self.l4.setGeometry(20, 250, 300, 50)
			self.l4.show()
			
			self.l5 = ScrollLabel(self)
			self.l5.setText(inname)
			self.l5.setGeometry(20, 350, 300, 50)
			self.l5.show()
			
			self.show()
			return()
		
		# prepare by opening default inputs, ready for modification
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
		
		# button to get relevant files
		options = QFileDialog.Options()
		fol_nme = QFileDialog.getExistingDirectory(self, "Select Folder Containing Required Input Files", "./PyCHAM/input/")
		
		# list of files (and any folders) here
		dir_con = os.listdir(fol_nme)
		
		# unknown names
		sch_name = xml_name = inname = 'Not found'
		
		# look for corresponding files here
		for i in dir_con:
			if ('chem' in i): # chemical scheme file
				sch_name = str(fol_nme+'/'+i)
			if ('xml' in i): # xml file
				xml_name = str(fol_nme+'/'+i)
			if ('var' in i): # model variables file
				inname = str(fol_nme+'/'+i)
		
		# pickle with new file names	
		list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O]
		with open(input_by_sim, 'wb') as pk:
			pickle.dump(list_vars, pk) # pickle
			pk.close() # close
			
		# update label and buttons to show found files and allow choice of different file
		click1_up(self)
		
		return()

	@pyqtSlot()
	def on_click2(self): # when different chemical scheme file requires selection
	
		# prepare by opening default inputs, ready for modification
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
	
		sch_name, _ = QFileDialog.getOpenFileName(self, "Select Chemical Scheme File", "./PyCHAM/input/") # get path of file
		self.l3.clear() # clear old label
		self.l3.setText(str("Chemical scheme file: " + sch_name))
		self.l3.show()

		# pickle with new chemical scheme name	
		list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O]
		with open(input_by_sim, 'wb') as pk:
			pickle.dump(list_vars, pk) # pickle
			pk.close() # close
			
	@pyqtSlot()
	def on_click3(self): # when different xml file requires selection
	
		# prepare by opening default inputs, ready for modification
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
	
		xml_name, _ = QFileDialog.getOpenFileName(self, "Select xml File", "./PyCHAM/input/") # get path of file
		self.l4.clear() # clear old label
		self.l4.setText(str("xml file: " + xml_name))
		self.l4.show()
		
		# pickle with new chemical scheme name	
		list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O]
		with open(input_by_sim, 'wb') as pk:
			pickle.dump(list_vars, pk) # pickle
			pk.close() # close
			
	@pyqtSlot()
	def on_click4(self): # when different model variables file requires selection

		inname, _ = QFileDialog.getOpenFileName(self, "Select Model Variables File", "./PyCHAM/input/") # get path of file
		
		self.l5.clear() # clear old label
		self.l5.setText(str("Model Variables file: " + inname))
		self.l5.show()
		
		return()
		

	@pyqtSlot()
	def on_click5(self): # interrogating the model variables file to check inputs
	
		@pyqtSlot()
		def click5_up(self): # update model variables following click of check button
			
			self.e7 = QLineEdit(self)
			self.e7.setText(sav_nam)
			self.e7.move(460, 130)
			self.e7.show()
			
			self.e8 = QLineEdit(self)
			self.e8.setText((str(chem_sch_mark)).replace('\'', '').replace(' ', '')[1:-1])
			self.e8.move(520, 160)
			self.e8.show()
			
			self.e9 = QLineEdit(self)
			self.e9.setText((str(tot_time)).replace('\'', '').replace(' ', ''))
			self.e9.move(520, 190)
			self.e9.show()
			
			self.e10 = QLineEdit(self)
			self.e10.setText((str(update_stp)).replace('\'', '').replace(' ', ''))
			self.e10.move(520, 220)
			self.e10.show()
			
			self.e11 = QLineEdit(self)
			self.e11.setText((str(save_step)).replace('\'', '').replace(' ', ''))
			self.e11.move(520, 250)
			self.e11.show()
			
			self.show()
			return()

		# prepare by opening default inputs, ready for modification
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
			pk.close() # close pickle file
		
		if (inname != ''): # if not using defaults
			inputs = open(inname, mode= 'r' ) # open model variables file
			in_list = inputs.readlines() # read file and store everything into a list
			inputs.close() # close file
		else: # if using defaults
			in_list = []
			
		for i in range(len(in_list)): # loop through supplied model variables to interpret

			key, value = in_list[i].split('=') # split values from keys
			# model variable name - a string with bounding white space removed
			key = key.strip()

			if key == 'res_file_name' and (value.strip()): # name of folder to save results in
				sav_nam = str(value.strip())

			if key == 'chem_scheme_markers' and (value.strip()): # formatting for chemical scheme
				chem_sch_mark = [str(i).strip() for i in (value.split(','))]

			if key == 'update_step' and (value.strip()): # time step (s) for updating ODE initial conditions
				update_stp = float(value.strip())

			if key == 'total_model_time' and (value.strip()):
				tot_time = float(value.strip())

			if key == 'Comp0' and (value.strip()): # names of components present at experiment start
				comp0 = [str(i).strip() for i in (value.split(','))]			

			if key == 'C0' and (value.strip()): # initial concentrations of components present at experiment start (ppb)
				y0 = [float(i) for i in (value.split(','))]

			if key == 'temperature' and (value.strip()): # chamber temperature (K)
				temp = [float(i) for i in ((value.strip()).split(','))]

			if key == 'tempt' and (value.strip()): # times (s) that temperature values correspond to
				tempt = [float(i) for i in ((value.strip()).split(','))]

			if key == 'rh' and (value.strip()): # relative humidity in chamber (0-1)
				RH = float(value.strip())

			if key == 'p_init' and (value.strip()): # pressure inside chamber
				Press = float(value.strip())

			if key == 'daytime_start' and (value.strip()): # time of day at experiment start (s)
				daytime = float(value.strip())

			if key == 'wall_on' and (value.strip()): # marker for whether or not to consider wall
				wall_on = int(value.strip())

			if key == 'eff_abs_wall_massC' and (value.strip()): # effective absorbing mass concentration of wall
				Cw = float(value.strip())

			if key == 'mass_trans_coeff' and (value.strip()): # mass transfer coefficient of vapours with wall
				kw = float(value.strip())

			if key == 'size_structure' and (value.strip()): # the size structure
				siz_stru = int(value.strip())

			if key == 'number_size_bins' and (value.strip()): # number of particle size bins
				num_sb = int(value.strip())

			if key == 'pconc' and (value.strip()): # seed particle number concentrations (#/cc)
				time_cnt = 1 # track number of times
				sb_cnt = 1 # track number of size bins
				mode_cnt = 1 # track number of modes
				for i in value:
					if i==';':
						time_cnt += 1 # increase time count
					if (time_cnt == 1 and i==','):
						sb_cnt += 1 # increase size bin count
						pmode = 1 # explicitly stated particle concentrations
					if (time_cnt == 1 and i == ':'):
						mode_cnt += 1 # mode count
						pmode = 0 # particle concentrations expressed as modes

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
					if i==';':
						time_cnt += 1 # increase time count
					if i==':':
						mode_cnt += 1 # increase mode count
				std = np.zeros((mode_cnt, time_cnt))
				for i in range(time_cnt):
					std[:, i] = [float(ii.strip()) for ii in ((value.split(';')[i]).split(':'))]
			
			if key == 'mean_rad' and (value.strip()): # seed particle mean radius (um)
				time_cnt = 1 # track of number of times
				mode_cnt = 1 # track number of modes
				for i in value:
					if i==';':
						time_cnt += 1 # increase time count
					if i==':':
						mode_cnt += 1 # increase mode count
				
				mean_rad = np.zeros((mode_cnt, time_cnt))
				for i in range(time_cnt):
					mean_rad[:, i] = [float(ii.strip()) for ii in ((value.split(';')[i]).split(':'))]
				
			if key == 'recording_time_step' and (value.strip()): # frequency (s) of storing results
				save_step = float(value.strip())

			if key == 'const_comp' and (value.strip()): # names of components with later continuous injections
				const_comp = [str(i).strip() for i in (value.split(','))]
			
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
				seed_mw = float(value.strip())
		
			if key == 'seed_diss' and value.strip(): # dissociation constant of seed material
				seed_diss = [float(i) for i in (value.split(','))]

			if key == 'seed_dens' and value.strip():
				seed_dens = float(value.strip())

			if key == 'seedVr' and value.strip(): # volume ratio of components in seed particles
				seedVr = [float(i) for i in (value.split(','))]

			if key == 'light_status' and value.strip(): # status of lights (on or off)
				light_stat = [int(i) for i in (value.split(','))]
				light_stat = np.array((light_stat))
				
			if key == 'light_time' and value.strip(): # times (s) corresponding to light status
				light_time = [float(i) for i in (value.split(','))]
				light_time = np.array((light_time))
				
			if key == 'lat': # latitude (degrees)
				if (value.strip()): lat = float(value.strip())

			if key == 'lon': # longitude (degrees)
				if (value.strip()): lon = float(value.strip())

			if key == 'act_flux_file' and (value.strip()): # for indoor actinic flux
				af_path = str(os.getcwd() + '/PyCHAM/photofiles/' + value.strip())

			if key == 'DayOfYear' and (value.strip()):
				dayOfYear = int(value.strip())		
			
			# name of file with wavelength-dependent absorption cross-section and quantum yield 
			# calculations
			if key == 'photo_par_file' and (value.strip()):
				photo_path = str(os.getcwd() + '/PyCHAM/photofiles/' + value.strip())
			
			if key == 'trans_fac' and (value.strip()): # transmission factor for natural light
				tf = float(value.strip())

			if key == 'light_adapt' and (value.strip()): # whether to adapt time step to naturally varying light intensity
				light_ad = int(value.strip())

			if key == 'const_infl' and (value.strip()): # names of components with continuous influx
				con_infl_nam = [str(i).strip() for i in (value.split(','))]

			if key == 'const_infl_t' and (value.strip()): # times of continuous influxes (s)
				con_infl_t = [float(i.strip()) for i in (value.split(','))]
				con_infl_t = np.array((con_infl_t))

			if key == 'Cinfl' and (value.strip()): # influx rate of components with continuous influx (ppb/s)
				comp_count = 1 # count number of components
				time_count = 1 # track number of times
				for i in value:
					if i==';':
						comp_count += 1 # record number of components
						time_count = 1 # reset time count
					if i==',':
						time_count += 1
				con_infl_C = np.zeros((comp_count, time_count))
				for i in range(comp_count):
					con_infl_C[i, :] = [float(ii.strip()) for ii in ((value.split(';')[i]).split(','))]

			if key == 'tracked_comp' and (value.strip()): # names of components whose tendency to change will be tracked
				dydt_trak = [str(i).strip() for i in (value.split(','))]

			if key == 'dens_Comp' and (value.strip()):
				dens_comp = [str(i).strip() for i in (value.split(','))]

			if key == 'dens' and (value.strip()):
				dens = [float(i) for i in (value.split(','))]

			if key == 'vol_Comp' and (value.strip()):
				vol_comp = [str(i).strip() for i in (value.split(','))]

			if key == 'volP' and (value.strip()):
				volP = [float(i) for i in (value.split(','))]

			if key == 'act_comp' and (value.strip()): # names of componentes with specified activity coefficients
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
				dil_fac = float(value)
				
			if key == 'ser_H2O' and (value.strip()): # whether to serialise water gas-particle partitioning
				ser_H2O = int(value)
		
		# prepare for pickling
		list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O]

		input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
		with open(input_by_sim, 'wb') as pk: # the file to be used for pickling
			pickle.dump(list_vars, pk) # pickle
			pk.close() # close
	
		import ui_check # check on inputs module
		ui_check.ui_check() # check on inputs
		
		# update checked inputs
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
			pk.close() # close pickle file
			
		# display inputs
		click5_up(self)
	
	@pyqtSlot()	
	def on_click6(self): # button to reset model variables to those presented
		
		sav_nam = self.e7.text()
		chem_sch_mark = self.e8.text().split(',')
		tot_time = float(self.e9.text())
		update_stp = float(self.e10.text())
		save_step = float(self.e11.text())
		
	
	@pyqtSlot()	
	def on_click7(self): # button to run simulation
		import middle # prepare to communicate with main programme
		# call on modules to solve problem
		middle.middle()
	
	@pyqtSlot() # button to plot results graphically
	def on_click8(self):	
		import plotter
		print('Plotting and saving standard results graph')
		plotter.plotter(0) # plot results

	@pyqtSlot() # button to quit software
	def on_click9(self):
		QWidget.close(self)
		
	def openFileNameDialog(self): # allows opening of system's directory navigator
		options = QFileDialog.Options()
		fname, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
		return(fname)
