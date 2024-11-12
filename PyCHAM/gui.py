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
'''The module that generates the Graphical User Interface (GUI) for 
PyCHAM, and connects that GUI with the core PyCHAM model'''

# first module called when PyCHAM started from the terminal/command 
# window, takes inputs
# and sends to model modules, also calls the saving module

from PyQt5.QtWidgets import *
from PyQt5.QtGui import  *
from PyQt5.QtCore import *
import pickle # for storing inputs
import sys
import os
import def_mod_var
import numpy as np
import re
import vol_contr_analys
import importlib
		
class PyCHAM(QWidget):
	
	# set the reference (self) to the current instance of the PyCHAM class,
	# note that by setting a default value for param_const, a user doesn't need
	# to provide an input for param_const.  But, having param_const as an
	# optional argument allows the automated setup and call 
	# (automated_setup_and_call.py) to call gui
	def __init__(self, param_const=0):
		
		super().__init__()

		self.param_const = param_const

		# get path to PyCHAM folder
		self.PyCHAM_path = os.path.dirname(__file__)[0:-7]	
		
		self.title = 'PyCHAM'
		self.left = 10
		self.top = 10
		self.width = 800
		self.height = 530
		self.initUI() # call on initialisation function to fill window	

		# changing the background color
		#self.setStyleSheet("background-color: darkCyan;")
  
		# set the title
		#self.setWindowTitle('PyCHAM')		

		self.err_mess = '' # begin with no error message
		
		# if parameters have been provided automatically (without using GUI)
		# then automatically setup and run the simulation, e.g. using the
		# automated_setup_and_call.py
		if (type(param_const) == dict):
			self.autorun() # call on automatic setup and run
			
		return
		
	def initUI(self):
	
		# general window setup --------------------------------------------------------
		self.setWindowTitle(self.title)
		self.setGeometry(self.left, self.top, self.width, self.height)
		
		# define grid layout
		grid = QGridLayout() 
		self.setLayout(grid)
		
		# title ------------------------------------------------------------------------
		# PyCHAM logo
		l00 = QLabel(self)
		pixmap = QPixmap('PyCHAM/PyCHAM_logo_transparent.png')
		l00.setPixmap(pixmap.scaled(400, 100, transformMode = Qt.SmoothTransformation))
		PyCHlogo_hindx = 0 # horizontal position on grid
		grid.addWidget(l00, 0, PyCHlogo_hindx, 1, 1)
		# link to PyCHAM website
		bn00 = QPushButton('', self)
		bn00.setStyleSheet('background : transparent; border : 0px; color : white')
		bn00.setToolTip('Visit the PyCHAM README (online)')
		bn00.clicked.connect(self.on_clickn00)
		grid.addWidget(bn00, 0, PyCHlogo_hindx, 1, 1)

		#l0 = QLabel(self)
		#l0.setOpenExternalLinks(True)
		#l0.setText(str('Welcome to PyCHAM - please see ' + 
		#'<a href=\"\">README</a> for guidance'))
		#l0.setFont(QFont('Arial', 20))
		#grid.addWidget(l0, 0, 0, 1, 1)
		
		# EUROCHAMP logo
		l0b= QLabel(self)
		pixmap = QPixmap('PyCHAM/logos_eurochamp2020-orange-noir.png')
		l0b.setPixmap(pixmap.scaled(100, 100, transformMode = Qt.SmoothTransformation))
		EUROlogo_hindx = 3 # horizontal position on grid
		grid.addWidget(l0b, 0, EUROlogo_hindx)
		# link to EUROCHAMP website
		bn1 = QPushButton('', self)
		bn1.setStyleSheet('background : transparent; border : 0px; color : white')
		bn1.setToolTip('Visit the EUROCHAMP website')
		bn1.clicked.connect(self.on_clickn1)
		grid.addWidget(bn1, 0, EUROlogo_hindx)
		
		# NCAS logo
		l0c= QLabel(self)
		pixmap = QPixmap('PyCHAM/NCAS_national_centre_logo_transparent.png')
		l0c.setPixmap(pixmap.scaled(160, 40, transformMode = Qt.SmoothTransformation))
		NCASlogo_hindx = EUROlogo_hindx+1 # horizontal position on grid
		grid.addWidget(l0c, 0, NCASlogo_hindx)
		# link to NCAS website
		bn1a = QPushButton('', self)
		bn1a.setStyleSheet('background : transparent; border : 0px; color : white')
		bn1a.setToolTip('Visit the NCAS website')
		bn1a.clicked.connect(self.on_clickn1a)
		grid.addWidget(bn1a, 0, NCASlogo_hindx)
		
		# University of Manchester logo
		l0d= QLabel(self)
		pixmap = QPixmap('PyCHAM/TAB_col_background.png')
		l0d.setPixmap(pixmap.scaled(96, 32, transformMode = Qt.SmoothTransformation))
		UoMlogo_hindx = EUROlogo_hindx+2 # horizontal position on grid
		grid.addWidget(l0d, 0, UoMlogo_hindx)
		# link to EUROCHAMP website
		bn1b = QPushButton('', self)
		bn1b.setStyleSheet('background : transparent; border : 0px; color : white')
		bn1b.setToolTip('Visit the University of Manchester website')
		bn1b.clicked.connect(self.on_clickn1b)
		grid.addWidget(bn1b, 0, UoMlogo_hindx)
		
		# add tabs --------------------------------------------------------------------
		
		tabs = QTabWidget()
		tabs.addTab(self.NStab(), "Simulate")
		tabs.addTab(self.PLtab(), "Plot")
		grid.addWidget(tabs, 1, 0, 1, UoMlogo_hindx)

		# README button ------------------------------------------
		b89p = QPushButton('README', self)
		b89p.setToolTip('Go to the online PyCHAM README page (alternatively see README in your local PyCHAM repository)')
		b89p.clicked.connect(self.on_clickn00)
		grid.addWidget(b89p, 0, PyCHlogo_hindx+1)

		# Quit pane ----------------------------------------------
		b89 = QPushButton('Quit', self)
		b89.setToolTip('Finish with PyCHAM and close this window')
		b89.clicked.connect(self.on_click89)
		grid.addWidget(b89, 1, UoMlogo_hindx)
		
		
		self.show()
		return
	
	# simulate tab - note that this called only once to
	# set uo the simulation tab when PyCHAM started
	def NStab(self):
	
		NSTab = QWidget()
		self.NSlayout = QGridLayout() 
		NSTab.setLayout(self.NSlayout)
		
		# folder and file starting column number
		ffscn = 0
		
		# folder selection for new simulation -----------------------------------------
		b0 = QPushButton('Select Folder Containing Input Files', self)
		b0.setToolTip('Select the folder containing the required input files')
		b0.clicked.connect(self.on_click1)
		self.NSlayout.addWidget(b0, 0, ffscn, 1, 2)
		
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
		
		# listing input files ----------------------------------
		l1 = QLabel(self)
		l1.setText("The following files have been found: ")
		self.NSlayout.addWidget(l1, 1, ffscn, 1, 2)
		
		l2 = QLabel(self)
		l2.setText("Chemical \nScheme: ")
		self.NSlayout.addWidget(l2, 2, ffscn, 1, 1)
			
		self.l3 = ScrollLabel(self)
		self.l3.setText(self.sch_name)
		self.NSlayout.addWidget(self.l3, 2, ffscn+1,  1, 1)
		
		b1 = QPushButton('Select Different Chemical Scheme', self)
		b1.setToolTip('Select the file containing the required chemical scheme file')
		b1.clicked.connect(self.on_click2)
		self.NSlayout.addWidget(b1, 3, ffscn, 1, 2)
		
		l4 = QLabel(self)
		l4.setText("XML: ")
		self.NSlayout.addWidget(l4, 4, ffscn, 1, 1)
			
		self.l5 = ScrollLabel(self)
		self.l5.setText(self.xml_name)
		self.NSlayout.addWidget(self.l5, 4, ffscn+1, 1, 1)
			
		b2 = QPushButton('Select Different XML', self)
		b2.setToolTip('Select the file containing the required xml file')
		b2.clicked.connect(self.on_click3)
		self.NSlayout.addWidget(b2, 5, ffscn, 1, 2)
		
		l6 = QLabel(self)
		l6.setText("Model \nVariables: ")
		self.NSlayout.addWidget(l6, 6, ffscn, 1, 1)
			
		self.l7 = ScrollLabel(self)
		self.l7.setText(self.inname)
		self.NSlayout.addWidget(self.l7, 6, ffscn+1, 1, 1)
			
		b3 = QPushButton('Select Different Model Variables', self)
		b3.setToolTip('Select the file containing the required model variables file')
		b3.clicked.connect(self.on_click4)
		self.NSlayout.addWidget(b3, 7, ffscn, 1, 2)
		
		
		# vertical separator line -------------------------------
		self.separatorLine = QFrame()
		self.separatorLine.setFrameShape(QFrame.VLine)
		self.separatorLine.setFrameShadow(QFrame.Raised)
		self.NSlayout.addWidget(self.separatorLine, 0, ffscn+2, 10, 1)
		self.separatorLine.show()
		
		# model variables list ---------------------------------------------
		self.mvpn = ffscn+3 # column number for model variables list
		
		l8 = QLabel(self)
		l8.setText("The following model variables have been found: ")
		self.NSlayout.addWidget(l8, 0, self.mvpn, 1, 3)
		
		# begin a scroll area to contain model variables, note that
		# setting the location of the scroll area is done below its contents
		self.scroll = QScrollArea()
		self.scroll.setWidgetResizable(True)
		self.scrollwidget = QWidget() 
		self.varbox = QGridLayout()
		
		
		# contents of model variables scroll area ---------------
		
		# General ------------------------------
		
		l9a = QLabel(self)
		l9a.setText('General')
		self.varbox.addWidget(l9a, 0, 0)
		l9a.setFont(QFont("Arial", 13, QFont.Bold))
		gen_row = 0
		
		l9 = QLabel(self)
		l9.setText('Folder to save to: ')
		self.varbox.addWidget(l9, gen_row+1, 0)
		self.l9a = QLabel(self)
		self.l9a.setText(self.sav_nam)
		self.varbox.addWidget(self.l9a, gen_row+1, 1)
		
		l10 = QLabel(self)
		l10.setText('Chemical scheme markers: ')
		self.varbox.addWidget(l10, gen_row+2, 0)
		self.l10a = QLabel(self)
		self.l10a.setText((str(self.chem_sch_mrk)).replace(
		'\'', '').replace(' ', '')[1:-1])
		self.varbox.addWidget(self.l10a, gen_row+2, 1)
		
		l11 = QLabel(self)
		l11.setText('Total experiment time (s): ')
		self.varbox.addWidget(l11, gen_row+3, 0)
		self.l11a = QLabel(self)
		self.l11a.setText((str(self.tot_time)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.l11a, gen_row+3, 1)
		
		l12 = QLabel(self)
		l12.setText('Update time interval (s): ')
		self.varbox.addWidget(l12, gen_row+4, 0)
		self.l12a = QLabel(self)
		self.l12a.setText((str(self.update_stp)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.l12a, gen_row+4, 1)
		
		l13 = QLabel(self)
		l13.setText('Recording time interval (s): ')
		self.varbox.addWidget(l13, gen_row+5, 0)
		self.l13a = QLabel(self)
		self.l13a.setText((str(self.save_step)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.l13a, gen_row+5, 1)
		
		l13_1 = QLabel(self)
		l13_1.setText('Whether UManSysProp should (1) \nor should not (0) be updated: ')
		self.varbox.addWidget(l13_1, gen_row+6, 0)
		self.l13_1a = QLabel(self)
		self.l13_1a.setText((str(uman_up)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.l13_1a, gen_row+6, 1)
		
		l13_2 = QLabel(self)
		l13_2.setText('Absolute and relative integration \ntolerances: ')
		self.varbox.addWidget(l13_2, gen_row+7, 0)
		self.l13_2a = QLabel(self)
		self.l13_2a.setText((str(int_tol)).replace(
			'\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l13_2a, gen_row+7, 1)
		
		l13_2 = QLabel(self)
		l13_2.setText('Absolute and relative integration \ntolerances: ')
		self.varbox.addWidget(l13_2, gen_row+7, 0)
		self.l13_2a = QLabel(self)
		self.l13_2a.setText((str(int_tol)).replace(
			'\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l13_2a, gen_row+7, 1)
		
		# Chamber Environment ----------
		
		l14a = QLabel(self)
		l14a.setText('Chamber Environment')
		l14a.setFont(QFont("Arial", 13, QFont.Bold))
		env_row = gen_row+8
		self.varbox.addWidget(l14a, env_row+0, 0)
		
		l14 = QLabel(self)
		l14.setText('Chamber temperature(s) (K): ')
		self.varbox.addWidget(l14, env_row+1, 0)
		self.l14a = QLabel(self)
		self.l14a.setText((str(self.TEMP)).replace('\'', '').replace(' ', '')[1:-1])
		self.varbox.addWidget(self.l14a, env_row+1, 1)
		
		l15 = QLabel(self)
		l15.setText('Time(s) for chamber \ntemperature(s) (s): ')
		self.varbox.addWidget(l15, env_row+2, 0)
		self.l15a = QLabel(self)
		self.l15a.setText((str(self.tempt)).replace('\'', '').replace(' ', '')[1:-1])
		self.varbox.addWidget(self.l15a, env_row+2, 1)
		
		l16 = QLabel(self)
		l16.setText('Relative humidity (0-1): ')
		self.varbox.addWidget(l16, env_row+3, 0)
		self.l16a = QLabel(self)
		self.l16a.setText((str(self.RH)).replace('\'', '').replace(' ', ','))
		self.varbox.addWidget(self.l16a, env_row+3, 1)
		
		l16b = QLabel(self)
		l16b.setText('Relative humidity times (s): ')
		self.varbox.addWidget(l16b, env_row+4, 0)
		self.l16c = QLabel(self)
		self.l16c.setText((str(self.RHt)).replace('\'', '').replace(' ', ','))
		self.varbox.addWidget(self.l16c, env_row+4, 1)
		
		l17 = QLabel(self)
		l17.setText('Chamber air pressure (Pa): ')
		self.varbox.addWidget(l17, env_row+5, 0)
		self.l17a = QLabel(self)
		self.l17a.setText((str(Press)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.l17a, env_row+5, 1)
		
		l17_1 = QLabel(self)
		l17_1.setText('Dilution factor (volume fraction \nof chamber per second): ')
		self.varbox.addWidget(l17_1, env_row+6, 0)
		self.l17_1a = QLabel(self)
		self.l17_1a.setText((str(self.dil_fac)).replace(
			'\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l17_1a, env_row+6, 1)
		
		# particle properties ----------------
		
		l18a = QLabel(self)
		l18a.setText('Particle Properties')
		l18a.setFont(QFont("Arial", 13, QFont.Bold))
		par_row = env_row+7
		
		self.varbox.addWidget(l18a, par_row+0, 0)
		
		l18 = QLabel(self)
		l18.setText('Size structure (MC=0, FM=1): ')
		self.varbox.addWidget(l18, par_row+1, 0)
		self.l18a = QLabel(self)
		self.l18a.setText((str(siz_stru)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.l18a, par_row+1, 1)
		
		l19 = QLabel(self)
		l19.setText('Number of particle size bins: ')
		self.varbox.addWidget(l19, par_row+2, 0)
		self.l19a = QLabel(self)
		self.l19a.setText((str(num_sb)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.l19a, par_row+2, 1)
		
		l20 = QLabel(self)
		l20.setText('Particle concentrations (#/cc): ')
		self.varbox.addWidget(l20, par_row+3, 0)
		self.l20a = QLabel(self)
		self.l20a.setText((str(self.pconc)).replace(
			'\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l20a, par_row+3, 1)
		
		l21 = QLabel(self)
		l21.setText('Particle concentration times (s): ')
		self.varbox.addWidget(l21, par_row+4, 0)
		self.l21a = QLabel(self)
		self.l21a.setText((str(self.pconct)).replace(
			'\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l21a, par_row+4, 1)
		
		l22 = QLabel(self)
		l22.setText('Molecular weight of seed particle \ncomponent (g/mol): ')
		self.varbox.addWidget(l22, par_row+5, 0)
		self.l22a = QLabel(self)
		self.l22a.setText((str(seed_mw)).replace('\'', '').replace(
			' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l22a, par_row+5, 1)
		
		l23 = QLabel(self)
		l23.setText('Dissociation constant(s) of seed \ncomponent(s): ')
		self.varbox.addWidget(l23, par_row+6, 0)
		self.l23a = QLabel(self)
		self.l23a.setText((str(seed_diss)).replace('\'', '').replace(
			' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l23a, par_row+6, 1)
		
		l24 = QLabel(self)
		l24.setText('Density of seed particles (g/cc): ')
		self.varbox.addWidget(l24, par_row+7, 0)
		self.l24a = QLabel(self)
		self.l24a.setText((str(seed_dens)).replace('\'', '').replace(
			' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l24a, par_row+7, 1)
		
		l25 = QLabel(self)
		l25.setText('Name of seed component: ')
		self.varbox.addWidget(l25, par_row+8, 0)
		self.l25a = QLabel(self)
		self.l25a.setText((str(self.seed_name)).replace('\'', '').replace(
			' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l25a, par_row+8, 1)
		
		l26 = QLabel(self)
		l26.setText('Mole fraction of non-water \ncomponents in \ndry seed particles: ')
		self.varbox.addWidget(l26, par_row+9, 0)
		self.l26a = QLabel(self)
		self.l26a.setText((str(self.seedx)).replace('\'', '').replace(
			' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l26a, par_row+9, 1)
		
		l26b = QLabel(self)
		l26b.setText(str('Whether (1) or not (0) \nvolume of water ' +
			'included \nin seed particle \nnumber size distribution: '))
		self.varbox.addWidget(l26b, par_row+10, 0)
		self.l26b = QLabel(self)
		self.l26b.setText((str(self.Vwat_inc)).replace('\'', '').replace(
			' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l26b, par_row+10, 1)
		
		l26c = QLabel(self)
		l26c.setText(str('Whether (1) or not (0) \nwater to be ' +
			'equilibrated \nwith seed particle \nprior to experiment start: '))
		self.varbox.addWidget(l26c, par_row+11, 0)
		self.l26c = QLabel(self)
		self.l26c.setText((str(self.seed_eq_wat)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l26c, par_row+11, 1)
		
		l27 = QLabel(self)
		l27.setText('Smallest radius size bin \nboundary (um): ')
		self.varbox.addWidget(l27, par_row+12, 0)
		self.l27a = QLabel(self)
		self.l27a.setText((str(lowsize)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l27a, par_row+12, 1)
		
		l28 = QLabel(self)
		l28.setText('Largest radius size bin \nboundary (um): ')
		self.varbox.addWidget(l28, par_row+13, 0)
		self.l28a = QLabel(self)
		self.l28a.setText((str(uppsize)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l28a, par_row+13, 1)
		
		l29 = QLabel(self)
		l29.setText('Method for spacing size bins (lin) \nor (log): ')
		self.varbox.addWidget(l29, par_row+14, 0)
		self.l29a = QLabel(self)
		self.l29a.setText((str(self.space_mode)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l29a, par_row+14, 1)
		
		l30 = QLabel(self)
		l30.setText('Standard deviation for particle \nnumber size distribution: ')
		self.varbox.addWidget(l30, par_row+15, 0)
		self.l30a = QLabel(self)
		self.l30a.setText((str(std)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l30a, par_row+15, 1)
		
		l31 = QLabel(self)
		l31.setText('Mean radius (um) for particle \nnumber size distribution: ')
		self.varbox.addWidget(l31, par_row+16, 0)
		self.l31a = QLabel(self)
		self.l31a.setText((str(self.mean_rad)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l31a, par_row+16, 1)
		
		l32 = QLabel(self)
		l32.setText('Radius (cm) of newly nucleated \nparticles: ')
		self.varbox.addWidget(l32, par_row+17, 0)
		self.l32a = QLabel(self)
		self.l32a.setText((str(new_partr)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l32a, par_row+17, 1)
		
		l33 = QLabel(self)
		l33.setText('First nucleation \nparameterisation ' 
			'parameter: ')
		self.varbox.addWidget(l33, par_row+18, 0)
		self.l33a = QLabel(self)
		self.l33a.setText((str(self.nucv1)).replace('\'', 
			'').replace(' ', '').replace('[', 
			'').replace(']', ''))
		self.varbox.addWidget(self.l33a, par_row+18, 1)
		
		l34 = QLabel(self)
		l34.setText('Second nucleation \nparameterisation ' 
			'parameter: ')
		self.varbox.addWidget(l34, par_row+19, 0)
		self.l34a = QLabel(self)
		self.l34a.setText((str(self.nucv2)).replace('\'', 
		'').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l34a, par_row+19, 1)
		
		l35 = QLabel(self)
		l35.setText('Third nucleation \nparameterisation '
			'parameter: ')
		self.varbox.addWidget(l35, par_row+20, 0)
		self.l35a = QLabel(self)
		self.l35a.setText((str(self.nucv3)).replace('\'', 
		'').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l35a, par_row+20, 1)
		
		l36 = QLabel(self)
		l36.setText('Chemical scheme name of \nnucleating '
			'component: ')
		self.varbox.addWidget(l36, par_row+21, 0)
		self.l36a = QLabel(self)
		self.l36a.setText((str(self.nuc_comp)).replace('\'', 
		'').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l36a, par_row+21, 1)
		
		l37 = QLabel(self)
		l37.setText('Whether (1) or not (0) to adapt \n'
			'integration time interval \nto nucleation: ')
		self.varbox.addWidget(l37, par_row+22, 0)
		self.l37a = QLabel(self)
		self.l37a.setText((str(self.nuc_ad)).replace('\'', 
		'').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l37a, par_row+22, 1)
		
		l38 = QLabel(self)
		l38.setText('Whether (1) or not (0) to serialise \ngas-particle partitioning of water: ')
		self.varbox.addWidget(l38, par_row+23, 0)
		self.l38a = QLabel(self)
		self.l38a.setText((str(ser_H2O)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l38a, par_row+23, 1)
		
		l38_1 = QLabel(self)
		l38_1.setText('Whether (1) or not (0) to model \nparticle coagulation: ')
		self.varbox.addWidget(l38_1, par_row+24, 0)
		self.l38_1a = QLabel(self)
		self.l38_1a.setText((str(coag_on)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l38_1a, par_row+24, 1)
		
		l38_2 = QLabel(self)
		l38_2.setText('Particle-phase history with \nrespect to water partitioning (0 or 1): ')
		self.varbox.addWidget(l38_2, par_row+25, 0)
		self.l38_2a = QLabel(self)
		self.l38_2a.setText((str(wat_hist)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l38_2a, par_row+25, 1)
		
		l38_3 = QLabel(self)
		l38_3.setText('Deliquescence relative humidity \nas a function of temperature: ')
		self.varbox.addWidget(l38_3, par_row+26, 0)
		self.l38_3a = QLabel(self)
		self.l38_3a.setText((str(drh_str)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l38_3a, par_row+26, 1)
		
		l38_4 = QLabel(self)
		l38_4.setText('Efflorescence relative humidity \nas a function of temperature: ')
		self.varbox.addWidget(l38_4, par_row+27, 0)
		self.l38_4a = QLabel(self)
		self.l38_4a.setText((str(erh_str)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l38_4a, par_row+27, 1)
		
		l38_5 = QLabel(self)
		l38_5.setText('Whether injection of seed particle \nis instantaneous or continuous: ')
		self.varbox.addWidget(l38_5, par_row+28, 0)
		self.l38_5a = QLabel(self)
		self.l38_5a.setText((str(self.pcont)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l38_5a, par_row+28, 1)

		l38_6 = QLabel(self)
		l38_6.setText('Fraction of total gas-particle \npartitioning coefficient below \nwhich partitioning is zero: ')
		self.varbox.addWidget(l38_6, par_row+29, 0)
		self.l38_6a = QLabel(self)
		self.l38_6a.setText((str(z_prt_coeff)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l38_6a, par_row+29, 1)
		
		# gas inputs ----------------
		
		l39a = QLabel(self)
		l39a.setText('Gas')
		l39a.setFont(QFont("Arial", 13, QFont.Bold))
		gas_row = par_row+30
		self.varbox.addWidget(l39a, gas_row+0, 0)
		
		l40 = QLabel(self)
		l40.setText('Chemical scheme name of \ncomponents present at \nexperiment start: ')
		self.varbox.addWidget(l40, gas_row+1, 0)
		self.l40a = QLabel(self)
		self.l40a.setText((str(self.comp0)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l40a, gas_row+1, 1)
		
		l41 = QLabel(self)
		l41.setText('Concentrations of \ncomponents present at \nexperiment start (ppb): ')
		self.varbox.addWidget(l41, gas_row+2, 0)
		self.l41a = QLabel(self)
		self.l41a.setText((str(y0)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l41a, gas_row+2, 1)
		
		l42 = QLabel(self)
		l42.setText('Chemical scheme names of \n components with continuous \n influx: ')
		self.varbox.addWidget(l42, gas_row+3, 0)
		self.l42a = QLabel(self)
		self.l42a.setText((str(self.con_infl_nam)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l42a, gas_row+3, 1)
		
		l43 = QLabel(self)
		l43.setText('Influx rate of components with \ncontinuous influx (ppb/s): ')
		self.varbox.addWidget(l43, gas_row+4, 0)
		self.l43a = QLabel(self)
		self.l43a.setText((str(self.con_infl_C)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l43a, gas_row+4, 1)
		
		l44 = QLabel(self)
		l44.setText('Times of component influx (s): ')
		self.varbox.addWidget(l44, gas_row+5, 0)
		self.l44a = QLabel(self)
		self.l44a.setText((str(self.con_infl_t)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l44a, gas_row+5, 1)
		
		l45 = QLabel(self)
		l45.setText('Chemical scheme names of components \nwith constant concentration: ')
		self.varbox.addWidget(l45, gas_row+6, 0)
		self.l45a = QLabel(self)
		self.l45a.setText((str(self.const_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l45a, gas_row+6, 1)
		
		l46 = QLabel(self)
		l46.setText('Chemical scheme names of components \ninjected instantaneously after start \nof experiment: ')
		self.varbox.addWidget(l46, gas_row+7, 0)
		self.l46a = QLabel(self)
		self.l46a.setText((str(Compt)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l46a, gas_row+7, 1)

		l47 = QLabel(self)
		l47.setText('Times through experiment at which \ninstantaneous injection of components \noccur (s): ')
		self.varbox.addWidget(l47, gas_row+8, 0)
		self.l47a = QLabel(self)
		self.l47a.setText((str(injectt)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l47a, gas_row+8, 1)
		
		l48 = QLabel(self)
		l48.setText('Concentrations of components injected \ninstantaneously (ppb): ')
		self.varbox.addWidget(l48, gas_row+9, 0)
		self.l48a = QLabel(self)
		self.l48a.setText((str(Ct)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l48a, gas_row+9, 1)
		
		# light inputs ---------------------------------------

		l49b = QLabel(self)
		l49b.setText('Lights')
		l49b.setFont(QFont("Arial", 13, QFont.Bold))
		light_row = gas_row+10
		self.varbox.addWidget(l49b, light_row+0, 0)
		
		l49 = QLabel(self)
		l49.setText('Whether lights off (0) or on (1): ')
		self.varbox.addWidget(l49, light_row+1, 0)
		self.l49a = QLabel(self)
		self.l49a.setText((str(self.light_stat)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l49a, light_row+1, 1)
		
		l50 = QLabel(self)
		l50.setText('Time through simulation that \nlight status attained (s): ')
		self.varbox.addWidget(l50, light_row+2, 0)
		self.l50a = QLabel(self)
		self.l50a.setText((str(self.light_time)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l50a, light_row+2, 1)
		
		l51 = QLabel(self)
		l51.setText('Time of day experiment starts \n(s since 00:00): ')
		self.varbox.addWidget(l51, light_row+3, 0)
		self.l51a = QLabel(self)
		self.l51a.setText((str(self.daytime)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l51a, light_row+3, 1)
		
		l52 = QLabel(self)
		l52.setText('Latitude of experiment (degrees): ')
		self.varbox.addWidget(l52, light_row+4, 0)
		self.l52a = QLabel(self)
		self.l52a.setText((str(self.lat)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l52a, light_row+4, 1)
		
		l53 = QLabel(self)
		l53.setText('Longitude of experiment (degrees): ')
		self.varbox.addWidget(l53, light_row+5, 0)
		self.l53a = QLabel(self)
		self.l53a.setText((str(self.lon)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l53a, light_row+5, 1)
		
		l54 = QLabel(self)
		l54.setText('Path to customised (non-MCM) \nactinic flux file: ')
		self.varbox.addWidget(l54, light_row+6, 0)
		self.l54a = QLabel(self)
		self.l54a.setText((str(self.af_path)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l54a, light_row+6, 1)
		
		l55 = QLabel(self)
		l55.setText('Path to file containing absorption \ncross-sections and quantum yields: ')
		self.varbox.addWidget(l55, light_row+7, 0)
		self.l55a = QLabel(self)
		self.l55a.setText((str(self.photo_path)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l55a, light_row+7, 1)
		
		l56 = QLabel(self)
		l56.setText('Day number of experiment \n(number of days since the \npreceding 31st December): ')
		self.varbox.addWidget(l56, light_row+8, 0)
		self.l56a = QLabel(self)
		self.l56a.setText((str(self.dayOfYear)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l56a, light_row+8, 1)
		
		l57 = QLabel(self)
		l57.setText('Transmission factor for natural \nsunlight (0-1 fraction): ')
		self.varbox.addWidget(l57, light_row+9, 0)
		self.l57a = QLabel(self)
		self.l57a.setText((str(self.tf)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l57a, light_row+9, 1)
		
		l58 = QLabel(self)
		l58.setText('Whether to (1) or not to (0) adapt \nintegration time interval and initial \ncondition update to changing \nnatural light intensity: ')
		self.varbox.addWidget(l58, light_row+10, 0)
		self.l58a = QLabel(self)
		self.l58a.setText((str(self.light_ad)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l58a, light_row+10, 1)
		
		l58b = QLabel(self)
		l58b.setText('Transmission factor for \n254 nm wavelength light: ')
		self.varbox.addWidget(l58b, light_row+11, 0)
		self.l58bb = QLabel(self)
		self.l58bb.setText((str(self.tf_UVC)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l58bb, light_row+11, 1)

		l58c = QLabel(self)
		l58c.setText('Transmission factor for \n254 nm wavelength light times (s): ')
		self.varbox.addWidget(l58c, light_row+12, 0)
		self.l58cc = QLabel(self)
		self.l58cc.setText((str(self.tf_UVCt)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l58cc, light_row+12, 1)
		
		
		
		# wall inputs ---------------------------
		
		l59b = QLabel(self)
		l59b.setText('Walls')
		l59b.setFont(QFont("Arial", 13, QFont.Bold))
		wall_row = light_row+13
		self.varbox.addWidget(l59b, wall_row+0, 0)
		
		l59 = QLabel(self)
		l59.setText('Whether wall off (0) or on (1): ')
		self.varbox.addWidget(l59, wall_row+1, 0)
		self.l59a = QLabel(self)
		self.l59a.setText((str(self.wall_on)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l59a, wall_row+1, 1)
		
		l60 = QLabel(self)
		l60.setText('Effective absorbing mass of \nwall (g/m3 (air)): ')
		self.varbox.addWidget(l60, wall_row+2, 0)
		self.l60a = QLabel(self)
		self.l60a.setText((str(self.Cw)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l60a, wall_row+2, 1)
		
		l61 = QLabel(self)
		l61.setText('Gas-wall mass transfer \ncoefficient (/s): ')
		self.varbox.addWidget(l61, wall_row+3, 0)
		self.l61a = QLabel(self)
		self.l61a.setText((str(self.kw)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l61a, wall_row+3, 1)
		
		l62 = QLabel(self)
		l62.setText('Diameter of deposition function \ninflection (m): ')
		self.varbox.addWidget(l62, wall_row+4, 0)
		self.l62a = QLabel(self)
		self.l62a.setText((str(inflectDp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l62a, wall_row+4, 1)
		
		l63 = QLabel(self)
		l63.setText('Gradient of the logarithm of \nparticle deposition rate vs. logarithm \nof particle diameter at diameters \nbelow the inflection: ')
		self.varbox.addWidget(l63, wall_row+5, 0)
		self.l63a = QLabel(self)
		self.l63a.setText((str(pwl_xpre)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l63a, wall_row+5, 1)
		
		l64 = QLabel(self)
		l64.setText('Gradient of the logarithm of \nparticle deposition rate \nvs. logarithm \nof particle diameter at diameters \nabove the inflection: ')
		self.varbox.addWidget(l64, wall_row+6, 0)
		self.l64a = QLabel(self)
		self.l64a.setText((str(pwl_xpro)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l64a, wall_row+6, 1)
		
		l65 = QLabel(self)
		l65.setText('Particle deposition rate at the \ninflection (/s): ')
		self.varbox.addWidget(l65, wall_row+7, 0)
		self.l65a = QLabel(self)
		self.l65a.setText((str(inflectk)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l65a, wall_row+7, 1)
		
		l66 = QLabel(self)
		l66.setText('Chamber surface area (m2) : ')
		self.varbox.addWidget(l66, wall_row+8, 0)
		self.l66a = QLabel(self)
		self.l66a.setText((str(chamSA)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l66a, wall_row+8, 1)

		l66aa = QLabel(self)
		l66aa.setText('Chamber volume (m3) : ')
		self.varbox.addWidget(l66aa, wall_row+9, 0)
		self.l66aaa = QLabel(self)
		self.l66aaa.setText((str(chamV)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l66aaa, wall_row+9, 1)
		
		l67 = QLabel(self)
		l67.setText(str('Whether particle deposition to \nwall ' + 
		'treated by Rader and \nMcMurry (1) or customised (0): '))
		self.varbox.addWidget(l67, wall_row+10, 0)
		self.l67a = QLabel(self)
		self.l67a.setText((str(Rader)).replace('\'', '').replace(
			' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l67a, wall_row+10, 1)
		
		l68 = QLabel(self)
		l68.setText('Average number of charges per \nparticle (/particle): ')
		self.varbox.addWidget(l68, wall_row+11, 0)
		self.l68a = QLabel(self)
		self.l68a.setText((str(p_char)).replace('\'', '').replace(
			' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l68a, wall_row+11, 1)
		
		l69 = QLabel(self)
		l69.setText('Average electric field inside chamber \n(g.m/A.s3): ')
		self.varbox.addWidget(l69, wall_row+12, 0)
		self.l69a = QLabel(self)
		self.l69a.setText((str(e_field)).replace('\'', '').replace(
			' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l69a, wall_row+12, 1)
		
		# specific component properties ----------------------
		
		l70b = QLabel(self)
		l70b.setText('Specific Component Properties')
		l70b.setFont(QFont("Arial", 13, QFont.Bold))
		scp_row = wall_row+12
		self.varbox.addWidget(l70b, scp_row+0, 0)
		
		l70 = QLabel(self)
		l70.setText('''Chemical scheme name of components \nto 
			track process tendencies: ''')
		self.varbox.addWidget(l70, scp_row+1, 0)
		self.l70a = QLabel(self)
		self.l70a.setText((str(self.dydt_trak)).replace('\'', 
			'').replace(' ', '').replace('[', 
			'').replace(']', ''))
		self.varbox.addWidget(self.l70a, scp_row+1, 1)
		
		l71 = QLabel(self)
		l71.setText('''Chemical scheme name of components 
			\nwith specified densities: ''')
		self.varbox.addWidget(l71, scp_row+2, 0)
		self.l71a = QLabel(self)
		self.l71a.setText((str(dens_comp)).replace('\'', 
			'').replace(' ', '').replace('[', 
			'').replace(']', ''))
		self.varbox.addWidget(self.l71a, scp_row+2, 1)
		
		l72 = QLabel(self)
		l72.setText('Specified densities of \ncomponents (g/cc): ')
		self.varbox.addWidget(l72, scp_row+3, 0)
		self.l72a = QLabel(self)
		self.l72a.setText((str(dens)).replace('\'', '').replace(
			' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l72a, scp_row+3, 1)
		
		l73 = QLabel(self)
		l73.setText(str('Chemical scheme name of components ' +
		'\nwith specified vapour pressures: '))
		self.varbox.addWidget(l73, scp_row+4, 0)
		self.l73a = QLabel(self)
		self.l73a.setText((str(vol_comp)).replace(
			'\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l73a, scp_row+4, 1)
		
		l74 = QLabel(self)
		l74.setText('Specified vapour pressures of \ncomponents (Pa): ')
		self.varbox.addWidget(l74, scp_row+5, 0)
		self.l74a = QLabel(self)
		self.l74a.setText((str(volP)).replace('\'', '').replace(
			' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l74a, scp_row+5, 1)
		
		l75 = QLabel(self)
		l75.setText(str('Chemical scheme name of components ' +
		'\nwith specified activity coefficients: '))
		self.varbox.addWidget(l75, scp_row+6, 0)
		self.l75a = QLabel(self)
		self.l75a.setText((str(act_comp)).replace(
		'\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l75a, scp_row+6, 1)
		
		l76 = QLabel(self)
		l76.setText('Specified activity coefficients of \ncomponents:  ')
		self.varbox.addWidget(l76, scp_row+7, 0)
		self.l76a = QLabel(self)
		self.l76a.setText((str(act_user)).replace(
		'\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l76a, scp_row+7, 1)

		l77 = QLabel(self)
		l77.setText(str('Chemical scheme name of components \nwith ' + 
			'specified accommodation \ncoefficients: '))
		self.varbox.addWidget(l77, scp_row+8, 0)
		self.l77a = QLabel(self)
		self.l77a.setText((str(accom_comp)).replace(
			'\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l77a, scp_row+8, 1)
		
		l78 = QLabel(self)
		l78.setText('Specified accommodation coefficients of \ncomponents:  ')
		self.varbox.addWidget(l78, scp_row+9, 0)
		self.l78a = QLabel(self)
		self.l78a.setText((str(accom_val)).replace(
			'\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l78a, scp_row+9, 1)
		
		
		# properties of model variables scroll area ----------------
		self.scrollwidget.setLayout(self.varbox)
		self.scroll.setWidget(self.scrollwidget)
		self.NSlayout.addWidget(self.scroll, 1, self.mvpn, 3, 3)
		
		# --------------------------------------------------------------------
		# label to let user know preparedness of simulation - displayed text 
		# updated in ui_check module below
		self.l80 = ScrollLabel(self)
		self.NSlayout.addWidget(self.l80, 4, self.mvpn, 1, 3)
		self.bd_st = 2 # border status

		# --------------------------------------------------------------------
		# drop down button to let users view a variety of variables 
		# determined by model variables
		self.b80s = QComboBox(self)
		self.b80s.addItem('Photolysis Rates')
		self.b80s.addItem('Particle Number Size Distributions')
		self.b80s.addItem('Gas-phase Diffusion Coefficients')
		self.b80s.addItem('Gas-phase Mean Thermal Speeds')
		self.b80s.addItem('Molar Masses')
		self.b80s.addItem('Vapour Pressures')
		self.b80s.addItem('Nucleation Function')
		self.NSlayout.addWidget(self.b80s, 7, self.mvpn+0, 1, 2)
		
		# button to run checks on variables selected in drop down button
		self.b80 = QPushButton('Check Values', self)
		self.b80.setToolTip(str('See the values of the variables ' +
		'selected in the drop down button to the left'))
		self.b80.clicked.connect(self.on_clickb80)
		self.NSlayout.addWidget(self.b80, 7, self.mvpn+2, 1, 1)

		# -------------------------------------------------------------------

		# label to let users know file combinations included in batch
		self.btch_str = str('File combinations included in batch: ' +
		'chemical scheme, xml, model variables\n')
		# begin count on number of simulations in batch
		self.btch_no = 1
		
		# label to let users know progress through simulation
		self.l81b = ScrollLabel(self)
		self.NSlayout.addWidget(self.l81b, 6, self.mvpn, 1, 3)
		self.l81b.setText('')
		
		self.fab = 0 # to begin single simulation widgets turned off
		self.atb = 0 # to begin add to batch button turned off
		self.ssb = 0 # to begin series simulation button turned off
		
		self.output_list = [] # begin list of output paths
		
		# running check on default model variables ---------------------------------
		# let checking module know this is the first call
		# from gui
		self.chck_num = 1

		import ui_check # module for checking on model variables
		# check on inputs - note this loads the last saved pickle file and saves 
		# any change to this pickle file
		ui_check.ui_check(self)
		
		# finished check on model variables -----------------------------------------
		
		# relative stretching (width-wise) of each column in Simuate tab
		self.NSlayout.setColumnStretch(0, 1)
		self.NSlayout.setColumnStretch(1, 1)
		self.NSlayout.setColumnStretch(2, 1)
		self.NSlayout.setColumnStretch(3, 1)
		self.NSlayout.setColumnStretch(4, 1)
		self.NSlayout.setColumnStretch(5, 1)

		return(NSTab)
	
	def PLtab(self): # Plotting tab
	
		PLTab = QWidget()
		self.PLlayout = QGridLayout() 
		PLTab.setLayout(self.PLlayout)
		
		# results folder and dialogue row --------------
		
		self.l201 = ScrollLabel(self)
		cwd = os.getcwd() # current working directory
		path = str('Select results folder using Select New Folder')
		self.l201.setText(path)
		self.PLlayout.addWidget(self.l201, 0, 0, 1, 1)
		
		b202 = QPushButton('Select New Folder', self)
		b202.setToolTip('Select the folder containing the result files to plot')
		b202.clicked.connect(self.on_click202)
		self.PLlayout.addWidget(b202, 0, 1)
		
		# status for border around plotting messages
		self.bd_pl = 3
		
		# plotting dialogue box:
		self.l203a = ScrollLabel(self)
		self.l203a.setText('No message currently')
		self.PLlayout.addWidget(self.l203a, 0, 2, 1, 2)
		
		# vertical separator line -------------------------------
		#cil = 1# column index for line
		
		#self.separatorLine = QFrame()
		#self.separatorLine.setFrameShape(QFrame.VLine)
		#self.separatorLine.setFrameShadow(QFrame.Raised)
		#self.PLlayout.addWidget(self.separatorLine, 0, cil, 13, 1)
		#self.separatorLine.show()
		
		# ---------------------------------------------------------
		
		# include tabs
		
		PLtabs = QTabWidget()
		PLtabs.addTab(self.PRIMtab(), "Quick")
		PLtabs.addTab(self.SECtab(), "Flux")
		PLtabs.addTab(self.VOLtab(), "Volatility")
		PLtabs.addTab(self.PARtab(), "Particle") # for particle-phase things
		PLtabs.addTab(self.PYIELDtab(), "Yield") # for particle-phase things
		PLtabs.addTab(self.RADtab(), "Radicals") # for radical chemicals
		PLtabs.addTab(self.INSTRtab(), "Convolution")
		PLtabs.addTab(self.OBStab(), "Observed")
		PLtabs.addTab(self.EXPtab(), "Exp. Prep.")
		self.PLlayout.addWidget(PLtabs, 1, 0, 1, 5)
		
		# relative stretching (width-wise) of each column in Plot tab
		self.PLlayout.setColumnStretch(0, 1)
		self.PLlayout.setColumnStretch(1, 2)
		self.PLlayout.setColumnStretch(2, 1)
		self.PLlayout.setColumnStretch(3, 1)
		self.PLlayout.setColumnStretch(4, 2)
		
		return(PLTab)
		
		
	def PRIMtab(self): # basic plotting tab definition
	
		PRIMTab = QWidget()
		self.PRIMlayout = QGridLayout() 
		PRIMTab.setLayout(self.PRIMlayout)
	
		self.b203 = QPushButton('Standard Results Plot', self)
		self.b203.setToolTip('Create standard results plot (gas-phase concentration temporal profiles of components with specified initial concentrations and particle properties)')
		self.b203.clicked.connect(self.on_click203)
		self.PRIMlayout.addWidget(self.b203, 0, 0)

		# drop down button for units
		self.b203a = QComboBox(self)
		self.b203a.addItem('ppb')
		self.b203a.addItem(str(u'\u03BC' + 'g/m' +u'\u00B3'))
		self.b203a.addItem(str(u'\u0023' + ' molecules/cm' +u'\u00B3'))
		self.PRIMlayout.addWidget(self.b203a, 0, 1, 1, 1)

		# label for input line that receives component names to plot temporal profiles of
		self.l203 = QLabel(self)
		self.l203.setText('Write the chemical scheme names of components you want plotted in the box below (use H2O for water and RO2 for organic peroxy radicals)')
		self.l203.setWordWrap(True)
		self.PRIMlayout.addWidget(self.l203, 1, 0, 1, 1)

		# input bar for names of components to plot temporal profiles of
		self.e205 = QLineEdit(self)
		# show left most point first
		self.e205.setStyleSheet('qproperty-cursorPosition : 0')
		self.PRIMlayout.addWidget(self.e205, 2, 0)

		# drop down button for whether to display abundances of single 
		# components, or the sum of their abundances
		self.b205 = QComboBox(self)
		self.b205.setToolTip('Chose whether to display abundances of single components, or the sum of their abundances')
		self.b205.addItem('Individual Values')
		self.b205.addItem('Sum Value')
		self.PRIMlayout.addWidget(self.b205, 2, 1, 1, 1)

		# gas-phase concentrations temporal profiles -------------
		
		# button to plot temporal profile of gas-phase concentrations
		self.b206 = QPushButton(str('Gas-phase concentrations'), self)
		self.b206.setToolTip('Plot gas-phase concentration temporal profile for the specified components')
		self.b206.clicked.connect(self.on_click206)
		self.PRIMlayout.addWidget(self.b206, 3, 0)
		
		# drop down button for units
		self.b206b = QComboBox(self)
		self.b206b.addItem('ppb linear')
		self.b206b.addItem(str(u'\u03BC' + 'g/m' + u'\u00B3' + ' linear'))
		self.b206b.addItem(str(u'\u0023' + ' molecules/cm' + u'\u00B3' + ' linear'))
		self.b206b.addItem('ppb log.')
		self.b206b.addItem(str(u'\u03BC' + 'g/m' + u'\u00B3' + ' log.'))
		self.b206b.addItem(str(u'\u0023' + ' molecules/cm' + u'\u00B3' + ' log.'))
		self.PRIMlayout.addWidget(self.b206b, 3, 1, 1, 1)

		# button to plot ozone isopleth as relevant to the used 
		# chemical scheme and the oberved
		# range of VOC and NOx
		self.b206c = QPushButton(str('Ozone isopleth'), self)
		self.b206c.setToolTip('Plot equilibrium ozone concentrations over the range of simulated concentrations of NOx and VOC')
		self.b206c.clicked.connect(self.on_click206c)
		self.PRIMlayout.addWidget(self.b206c, 3, 2, 1, 1)

		# particle-phase concentrations temporal profiles -------------

		# button to plot temporal profile of total particle-phase concentration of supplied components
		self.b209 = QPushButton(str('Particle-phase concentration \n(summed over size bins) ('+u'\u03BC'+'g/m'+u'\u00B3'+')'), self)
		self.b209.setToolTip('Plot particle-phase concentration temporal profile of these components')
		self.b209.clicked.connect(self.on_click209)
		self.PRIMlayout.addWidget(self.b209, 4, 0)

		# button to plot temporal profile of total particle-phase 
		# concentration excluding seed and water 
		self.b209a = QPushButton(str('Particle-phase concentrations (summed over components and size bins) \n('+u'\u03BC'+'g/m'+u'\u00B3'+') excluding seed and water'), self)
		self.b209a.setToolTip('Plot total particle-phase concentration of all components except for seed and water')
		self.b209a.clicked.connect(self.on_click209a)
		#self.b209a.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.PRIMlayout.addWidget(self.b209a, 4, 1, 1, 2)
		
		# button to plot temporal profile of air quality index, using the DEFRA definition (see reference inside plotter.py)
		self.b209b = QPushButton(str('Air Quality Index'), self)
		self.b209b.setToolTip('Plot air quality index as a function of time')
		self.b209b.clicked.connect(self.on_click209b)
		self.PRIMlayout.addWidget(self.b209b, 5, 1, 1, 1)

		# button to plot temporal profile of total VOCs (including methane) in gas-phase
		self.b209c = QPushButton(str('TVOCS in gas phase'), self)
		self.b209c.setToolTip('Plot total volatile organic compounds (including methane) as a function of time')
		self.b209c.clicked.connect(self.on_click209c)
		self.PRIMlayout.addWidget(self.b209c, 5, 2, 1, 1)
		
		# wall (from gas-wall partitioning) concentrations temporal profiles -------------
		
		# button to plot temporal profile of total particle-phase concentrations
		self.b212 = QPushButton('Wall concentrations (excluding from particle deposition to wall) ('+u'\u03BC'+'g/m'+u'\u00B3'+')', self)
		self.b212.setToolTip('Plot the temporal profile of wall concentration (from gas-wall partitioning) for the specified components')
		self.b212.clicked.connect(self.on_click212)
		self.PRIMlayout.addWidget(self.b212, 5, 0)
		
		# wall (from particle deposition to wall) concentrations temporal profiles -------------
		
		# button to plot temporal profile of wall concentrations (from particle deposition to wall)
		self.b215 = QPushButton('Wall concentrations (from particle deposition to wall) ('+u'\u03BC'+'g/m'+u'\u00B3'+')', self)
		self.b215.setToolTip('Plot the temporal profile of the wall concentration (from particle deposition to wall) for the specified components')
		self.b215.clicked.connect(self.on_click215)
		self.PRIMlayout.addWidget(self.b215, 6, 0)
		
		# chamber conditions temporal profiles -------------
		
		# button to plot temporal profile of chamber environmental variables
		self.b215_a = QPushButton('Physical Conditions (T, P, RH)', self)
		self.b215_a.setToolTip('''Plot the temporal profile of chamber variables (temperature, pressure, relative humidity)''')
		self.b215_a.clicked.connect(self.on_click215_a)
		self.PRIMlayout.addWidget(self.b215_a, 0, 2)
		
		return(PRIMTab)
	
	def SECtab(self): # more detailed plotting tab definition

		SECTab = QWidget()
		self.SEClayout = QGridLayout() 
		SECTab.setLayout(self.SEClayout)
		
		# input bar for names of components to plot change tendencies
		self.e217 = QLineEdit(self)
		self.e217.setText('Provide the chemical scheme names of components for plotting change tendencies')
		self.e217.setStyleSheet('qproperty-cursorPosition : 0')
		self.SEClayout.addWidget(self.e217, 0, 0, 1, 1)
		
		# input bar for top number of reactions to consider
		self.e217a = QLineEdit(self)
		self.e217a.setText('Provide the number of chemical reactions to plot (arranged in descending order)')
		self.e217a.setStyleSheet('qproperty-cursorPosition : 0')
		self.SEClayout.addWidget(self.e217a, 0, 2, 1, 1)
		
		# button to plot temporal profile of change tendencies
		self.b218 = QPushButton('Plot change tendencies', self)
		self.b218.setToolTip('Plot the rate of change due to relevant processes')
		self.b218.clicked.connect(self.on_click218)
		self.SEClayout.addWidget(self.b218, 1, 0)
		
		# button to plot temporal profiles of individual chemical reaction change tendencies
		self.b218aa = QPushButton('Plot change tendency due to each chemical reaction', self)
		self.b218aa.setToolTip('Plot the rate of change of this component due to individual chemical reactions')
		self.b218aa.clicked.connect(self.on_click218aa)
		self.SEClayout.addWidget(self.b218aa, 1, 2)

		# drop down button to select units for change tendencies
		self.b218aaa = QComboBox(self)
		self.b218aaa.addItem('ppb/s')
		self.b218aaa.addItem(str(u'\u03BC' + 'g/m' +u'\u00B3' + '/s'))
		self.b218aaa.addItem(str(u'\u0023' + ' molecules/cm' +u'\u00B3' + '/s'))
		self.SEClayout.addWidget(self.b218aaa, 1, 3, 1, 1)

		# input bar for time period to consider production over
		self.e217aa = QLineEdit(self)
		self.e217aa.setText('Time to start production calculation at (hours), time to finish production calculation at (hours)')
		self.e217aa.setStyleSheet('qproperty-cursorPosition : 0')
		self.SEClayout.addWidget(self.e217aa, 2, 0, 1, 1)

		# button to estimate production integrated over stated period 
		self.b218ab = QPushButton('Production over this period', self)
		self.b218ab.setToolTip('Show production integrated over stated time period in the message box above')
		self.b218ab.clicked.connect(self.on_click218ab)
		self.SEClayout.addWidget(self.b218ab, 2, 2)

		# horizontal separator line -------------------------------
		self.separatorLine3 = QFrame()
		self.separatorLine3.setFrameShape(QFrame.HLine)
		self.separatorLine3.setFrameShadow(QFrame.Raised)
		self.SEClayout.addWidget(self.separatorLine3, 3, 0, 1, 4)
		self.separatorLine3.show()
		
		# vertical separator line -------------------------------
		self.separatorLine4 = QFrame()
		self.separatorLine4.setFrameShape(QFrame.VLine)
		self.separatorLine4.setFrameShadow(QFrame.Raised)
		self.SEClayout.addWidget(self.separatorLine4, 3, 1, 5, 1)
		self.separatorLine4.show()

		# ---------------------------------------------------------
		# input bar for reaction numbers (starting from 1), to get ratios for
		self.e220a = QLineEdit(self)
		self.e220a.setText('Reaction numbers (starting from 1) for numerator/reaction numbers for denominator (separate reactions with a comma and the numerator from denominator with a /)')
		self.e220a.setStyleSheet('qproperty-cursorPosition : 0')
		self.SEClayout.addWidget(self.e220a, 4, 0, 1, 1)

		# button to calculate and plot reaction rate ratios
		self.b220b = QPushButton('Reaction ratios', self)
		self.b220b.setToolTip('')
		self.b220b.clicked.connect(self.on_click220b)
		self.SEClayout.addWidget(self.b220b, 5, 0, 1, 1)

		# ---------------------------------------------------------
		
		# input bar for atom or functional group to plot contributions from
		self.e218a = QLineEdit(self)
		self.e218a.setText('Provide the SMILES names of atoms or functional groups for plotting component contributions (use RO2 for organic peroxy radicals)')
		self.e218a.setStyleSheet('qproperty-cursorPosition : 0')
		self.SEClayout.addWidget(self.e218a, 4, 2, 1, 1)
		
		# input bar for top number of components containing the relevant atom or functional groups
		self.e218b = QLineEdit(self)
		self.e218b.setText('Provide the number of components (in descending order) containing the atom/functional group to plot')
		self.e218b.setStyleSheet('qproperty-cursorPosition : 0')
		self.SEClayout.addWidget(self.e218b, 5, 2, 1, 1)
		
		# button to plot temporal profile of component contributions
		self.b218b = QPushButton('Plot component contributions (mole fraction)', self)
		self.b218b.setToolTip('Plot the contributions to this atom/functional group by component (mole fraction)')
		self.b218b.clicked.connect(self.on_click218b)
		self.SEClayout.addWidget(self.b218b, 6, 2, 1, 1)

		# button to plot carbon reservoirs with time
		self.b220c = QPushButton('Carbon reservoirs', self)
		self.b220c.setToolTip('View carbon flux accumulation with time')
		self.b220c.clicked.connect(self.on_click220c)
		self.SEClayout.addWidget(self.b220c, 6, 0, 1, 1)

		# horizontal separator line -------------------------------
		self.separatorLine3a = QFrame()
		self.separatorLine3a.setFrameShape(QFrame.HLine)
		self.separatorLine3a.setFrameShadow(QFrame.Raised)
		self.SEClayout.addWidget(self.separatorLine3a, 7, 0, 1, 4)
		self.separatorLine3a.show()

		# input bar for name of component to view molar mass (g/mol) of
		self.e218c = QLineEdit(self)
		self.e218c.setText('Provide the chemical scheme name of component for displaying the property selected below')
		self.e218c.setStyleSheet('qproperty-cursorPosition : 0')
		self.SEClayout.addWidget(self.e218c, 8, 0, 1, 1)

		# selection button for property of component to display
		# drop down button to select units for change tendencies
		self.b220e = QComboBox(self)
		self.b220e.addItem('Molar Mass (g/mol)')
		self.b220e.addItem('Pure component saturation vapour pressure at starting temperature of simulation (Pa)')
		self.b220e.addItem('Pure component saturation vapour pressure at 298.15 K (Pa)')
		self.SEClayout.addWidget(self.b220e, 8, 2, 1, 1)

		# button to plot carbon reservoirs with time
		self.b220d = QPushButton('Display Property', self)
		self.b220d.setToolTip('View selected property in dialogue box above')
		self.b220d.clicked.connect(self.on_click220d)
		self.SEClayout.addWidget(self.b220d, 8, 3, 1, 1)

		# column and row relative lengths---------------------------------
		
		# relative stretching (height-wise) of each row in Plot tab
		#self.PLlayout.setRowStretch(0, 1)
		#self.PLlayout.setRowStretch(1, 1)
		#self.PLlayout.setRowStretch(2, 1)
		#self.PLlayout.setRowStretch(3, 1)
		#self.PLlayout.setRowStretch(4, 1)
		#self.PLlayout.setRowStretch(5, 1)
		#self.PLlayout.setRowStretch(6, 1)
		#self.PLlayout.setRowStretch(7, 1)
		#self.PLlayout.setRowStretch(8, 1)
		#self.PLlayout.setRowStretch(9, 1)
		#self.PLlayout.setRowStretch(10, 1)
		#self.PLlayout.setRowStretch(11, 1)
		#self.PLlayout.setRowStretch(12, 1)
		#self.PLlayout.setRowStretch(13, 1)
		
		return(SECTab)

	def VOLtab(self): # more detailed plotting tab definition

		VOLTab = QWidget()
		self.VOLlayout = QGridLayout() 
		VOLTab.setLayout(self.VOLlayout)

		# drop down button for which phase to consider for volatility plotting
		self.b221a = QComboBox(self)
		self.b221a.addItem('Particle Phase')
		self.b221a.addItem('Gas Phase')
		self.b221a.addItem('Particle and Gas Phase Combined')
		self.b221a.addItem('Particle Phase Excluding Seed and Water')
		self.b221a.addItem('Gas Phase Only C>1, O>0')
		self.VOLlayout.addWidget(self.b221a, 0, 1, 1, 1)

		# volatility basis set ------------------
		
		
		# button to plot temporal profile of volatility basis set mass fractions without water
		self.b221 = QPushButton('Volatility Basis Set', self)
		self.b221.setToolTip('Plot the temporal profile of volatility basis set mass fractions')
		self.b221.clicked.connect(self.on_click221)
		#self.b221.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.VOLlayout.addWidget(self.b221, 0, 0)
		
		# two-dimensional volatility basis set ------------------
		
		# input bar for time through experiment to plot the 2D VBS for
		self.e222 = QLineEdit(self)
		self.e222.setText('Provide the time (seconds) through experiment at which to plot the options below - the closest recorded time to this will be used')
		self.e222.setStyleSheet('qproperty-cursorPosition : 0')
		self.VOLlayout.addWidget(self.e222, 2, 0)
		
		# button to plot 2D VBS
		self.b223 = QPushButton('2D Volatility Basis Set (excluding seed and water)', self)
		self.b223.setToolTip('Plot the two-dimensional volatility basis set (O:C ratio and vapour pressures) at the time through experiment specified above')
		self.b223.clicked.connect(self.on_click223)
		#self.b223.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.VOLlayout.addWidget(self.b223, 3, 0)
		
		# button to plot gas-phase concentration ordered by vapour pressure at the time through experiment defined by user
		self.b223a = QPushButton('Gas-phase Concentrations Ordered by Volatility', self)
		self.b223a.setToolTip('Plot gas-phase concentrations (ppb) of all components (ordered by volatility) at the time through experiment specified above')
		self.b223a.clicked.connect(self.on_click223a)
		self.VOLlayout.addWidget(self.b223a, 4, 0)
		
		# button to plot pie chart of mass fraction by top n components at n time after experiment start, with vapour pressure given
		self.b2231 = QPushButton('Pie Chart of Top Contributors', self)
		self.b2231.setToolTip('Plot a pie chart of the top (number specified to the right) mass contributors to the phase specified above at the time through experiment specified above')
		self.b2231.clicked.connect(self.on_click2231)
		self.VOLlayout.addWidget(self.b2231, 5, 0)
		
		# input bar for number of components to include in pie chart
		self.e2231 = QLineEdit(self)
		self.e2231.setText('Provide the number of components to consider in the pie chart')
		self.e2231.setStyleSheet('qproperty-cursorPosition : 0')
		self.VOLlayout.addWidget(self.e2231, 5, 1)

		return(VOLTab)

	# more detailed particle-phase plotting tab definition
	def PARtab(self): 		
		PARTab = QWidget()
		self.PARlayout = QGridLayout() 
		PARTab.setLayout(self.PARlayout)
		
		# input bar for diameter (um) limits for particle-phase
		# concentration properties
		self.e303p = QLineEdit(self)
		self.e303p.setText('Diameter (um) limits for options below, separate multiple values with a comma, e.g. 1.0, 2.5 for PM1 and PM2.5')
		self.e303p.setStyleSheet('qproperty-cursorPosition : 0')
		self.PARlayout.addWidget(self.e303p, 0, 0, 1, 3)
		
		# button to plot cumulative particle-phase mass 
		# concentration by different sizes of particle
		self.b303 = QPushButton(str('''Cumulative particle ''') 		+ str('''mass \nconcentration without water'''))
		self.b303.setToolTip('See the time series of particle mass grouped by upper diameter limits (mass and particle size excludes water)')
		self.b303.clicked.connect(self.on_click303)
		self.b303.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.PARlayout.addWidget(self.b303, 1, 0, 1, 1)
		
		# button to plot cumulative particle-phase number concentration by different 
		# sizes of particle
		self.b303q = QPushButton(str('Cumulative particle number \nconcentration without water'))
		self.b303q.setToolTip('See the time series of particle number grouped by upper diameter limits (particle size excludes water)')
		self.b303q.clicked.connect(self.on_click303q)
		self.b303q.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.PARlayout.addWidget(self.b303q, 1, 1, 1, 1)
		
		# button to plot cumulative particle-phase surface area concentration by different 
		# sizes of particle
		self.b303r = QPushButton(str('''Cumulative particle ''')
		+ str('''surface area \nconcentration without water'''))
		self.b303r.setToolTip(str('''See the time series of ''')
		+ str('''particle surface area grouped by upper ''') +
		str('''diameter limits (particle size excludes water)'''		))
		self.b303r.clicked.connect(self.on_click303r)
		self.b303r.setStyleSheet('''background-color : white; 
		border-width : 1px; border-radius : 7px; 
		border-color: silver; padding: 2px; 
		border-style : solid''')
		self.PARlayout.addWidget(self.b303r, 1, 2, 1, 1)
		
		# horizontal separator line ---------------------------
		self.separatorLine3 = QFrame()
		self.separatorLine3.setFrameShape(QFrame.HLine)
		self.separatorLine3.setFrameShadow(QFrame.Raised)
		self.PARlayout.addWidget(self.separatorLine3, 2, 0, 1, 
		3)
		self.separatorLine3.show()
		
		# input bar for number of components contributing
		# to particle-phase concentration
		self.e300 = QLineEdit(self)
		self.e300.setText(str('''Provide the top number of ''') 
		+ str('''components contributing to particle-phase'''))
		self.e300.setStyleSheet('qproperty-cursorPosition : 0')
		self.PARlayout.addWidget(self.e300, 3, 0, 1, 1)
	
		# button to plot particle-phase contributions
		self.b300 = QPushButton(str('Particle-phase contributions (%)'))
		self.b300.setToolTip('Show the contribution to particle-phase by the top contributors (number given in box above)')
		self.b300.clicked.connect(self.on_click300)
		#self.b300.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.PARlayout.addWidget(self.b300, 3, 2)

		# drop down button for type of contributors
		self.b300a = QComboBox(self)
		self.b300a.addItem('All components')
		self.b300a.addItem('Excluding Seed and Water')
		self.PARlayout.addWidget(self.b300a, 3, 1, 1, 1)

		# button to plot particle-phase surface concentration
		self.b301 = QPushButton(str('Particle surface area'))
		self.b301.setToolTip('Graph the temporal profile of particle surface area')
		self.b301.clicked.connect(self.on_click301)
		#self.b301.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.PARlayout.addWidget(self.b301, 4, 2)

		# drop down button for type of surface area
		self.b301a = QComboBox(self)
		self.b301a.addItem('All components')
		self.b301a.addItem('Seed Only')
		self.PARlayout.addWidget(self.b301a, 4, 1, 1, 1)

		# button to plot particle-phase mass contribution by 
		# different generations of oxidised organic molecules
		self.b302 = QPushButton(str('''Organic molecule ''') + 
		str('''contribution by generation'''))
		self.b302.setToolTip(str('''See the particle-phase ''') 
		+ str('''mass contribution by different generations ''')
		+ str('''of oxidised organic molecules'''))
		self.b302.clicked.connect(self.on_click302)
		self.b302.setStyleSheet('''background-color : 
		white; border-width : 1px; border-radius : 7px; 
		border-color: silver; padding: 2px; 
		border-style : solid''')
		self.PARlayout.addWidget(self.b302, 5, 0, 1, 3)

		# input bar for names of individual or groups of 
		# components contributing to particle mass
		self.e304 = QLineEdit(self)
		self.e304.setText(str('''Provide the names of ''') 
		+ str('''individual components or component groups ''')
		+ str('''for plotting mass contribution through time''')
		)
		self.e304.setStyleSheet('qproperty-cursorPosition : 0')
		self.PARlayout.addWidget(self.e304, 6, 0, 1, 1)

		# button to plot particle-phase mass contribution by 
		# different components
		self.b304 = QPushButton(str('''Contribution by ''') + 
		str('''mass per component with time'''))
		self.b304.setToolTip(str('''See the particle-phase ''') 
		+ str('''mass contribution by the provided ''')
		+ str('''components'''))
		self.b304.clicked.connect(self.on_click304)
		self.b304.setStyleSheet('''background-color : 
		white; border-width : 1px; border-radius : 7px; 
		border-color: silver; padding: 2px; 
		border-style : solid''')
		self.PARlayout.addWidget(self.b304, 6, 1, 1, 1)	

		# button to plot particle-phase risk contribution by 
		# different components
		self.b305 = QPushButton(str('''Contribution by ''') + 
		str('''risk per component with time'''))
		self.b305.setToolTip(str('''See the particle-phase ''') 
		+ str('''risk contribution by the provided ''')
		+ str('''components'''))
		self.b305.clicked.connect(self.on_click305)
		self.b305.setStyleSheet('''background-color : 
		white; border-width : 1px; border-radius : 7px; 
		border-color: silver; padding: 2px; 
		border-style : solid''')
		self.PARlayout.addWidget(self.b305, 6, 2, 1, 1)

		# button to plot particle-phase carbon oxidation state 
		# contribution by 
		# different components
		self.b306 = QPushButton(str('''Contribution by ''') + 
		str('''carbon oxidation state per component with ''') +
		str(''' time'''))
		self.b306.setToolTip(str('''See the particle-phase ''') 
		+ str('''carbon oxidation state contribution by the ''')
		+ str('''provided components'''))
		self.b306.clicked.connect(self.on_click306)
		self.b306.setStyleSheet('''background-color : 
		white; border-width : 1px; border-radius : 7px; 
		border-color: silver; padding: 2px; 
		border-style : solid''')
		self.PARlayout.addWidget(self.b306, 7, 1, 1, 1)
	
		return(PARTab)

	def PYIELDtab(self): # information on consumption and SOA yield
		
		PYIELDTab = QWidget()
		self.PYIELDlayout = QGridLayout() 
		PYIELDTab.setLayout(self.PYIELDlayout)
		
		# section for consumption and yield calculations -------
		# input bar for component to estimate consumption for
		self.e224 = QLineEdit(self)
		self.e224.setText('Provide the chemical scheme names of components to view consumption/yield for (result displayed in message box above, separate chemical names by a comma, e.g. APINENE, BENZENE)')
		self.e224.setStyleSheet('qproperty-cursorPosition : 0')
		self.PYIELDlayout.addWidget(self.e224, 0, 0, 1, 2)

		# input bar for starting time to estimate consumption for
		self.e224a = QLineEdit(self)
		self.e224a.setText('Provide the starting time to calculate consumption for (hours)')
		self.e224a.setStyleSheet('qproperty-cursorPosition : 0')
		self.PYIELDlayout.addWidget(self.e224a, 1, 0, 1, 2)

		# input bar for finshing time to estimate consumption for
		self.e224b = QLineEdit(self)
		self.e224b.setText('Provide the finishing time to calculate consumption for (hours)')
		self.e224b.setStyleSheet('qproperty-cursorPosition : 0')
		self.PYIELDlayout.addWidget(self.e224b, 2, 0, 1, 2)

		# button to estimate consumption
		self.b224 = QPushButton(str('Consumption by chemical reaction (' + u'\u03BC' + 'g/m' + u'\u00B3' +')'))
		self.b224.setToolTip('For the component specified above show the mass concentration consumed throughout the whole simulation in the message box above')
		self.b224.clicked.connect(self.on_click224)
		#self.b224.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.PYIELDlayout.addWidget(self.b224, 3, 0, 1, 1)

		# button to estimate SOA yield
		self.b225 = QPushButton(str('SOA yield by chemical reaction (fraction 0-1)'))
		self.b225.setToolTip('In the message box above show the SOA yield from the component given in the box above between the times stated in the boxes above.')
		self.b225.clicked.connect(self.on_click225)
		#self.b225.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.PYIELDlayout.addWidget(self.b225, 4, 1, 1, 1)
		
		return(PYIELDTab)

	def RADtab(self): # more detailed plotting tab for radical chemicals

		RADTab = QWidget()
		self.RADlayout = QGridLayout() 
		RADTab.setLayout(self.RADlayout)
		
		# input bar for number of components contributing
		# to radical population
		self.e360 = QLineEdit(self)
		self.e360.setText('Provide the top number of components to plot')
		self.e360.setStyleSheet('qproperty-cursorPosition : 0')
		self.RADlayout.addWidget(self.e360, 0, 0, 1, 2)

		# drop down button for type of radical
		self.b361a = QComboBox(self)
		self.b361a.addItem('RO2 (alkyl peroxy radical) Concentration (ppb)')
		self.b361a.addItem('RO (alkoxy radical) Concentration (ppb)')
		self.b361a.addItem('RO2 (alkyl peroxy radical) Flux (molecules/cm3/s)')
		self.b361a.addItem('RO (alkoxy radical) Flux (molecules/cm3/s)')
		self.RADlayout.addWidget(self.b361a, 1, 0, 1, 2)
		
		# input bar for carbon number threshold
		self.e361 = QLineEdit(self)
		self.e361.setText('Provide the minimum carbon number of a radical to be considered (defaults to 1)')
		self.e361.setStyleSheet('qproperty-cursorPosition : 0')
		self.RADlayout.addWidget(self.e361, 2, 0, 1, 1)
	
		# button to plot radical output
		self.b362 = QPushButton(str('Plot'))
		self.b362.setToolTip('Show the selected output by the top number of components (number given in box above)')
		self.b362.clicked.connect(self.on_click362)
		#self.b362.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.RADlayout.addWidget(self.b362, 2, 1)

		return(RADTab)

	def INSTRtab(self): # instrument comparison plotting tab definition
	
		INSTRTab = QWidget()
		self.INSTRlayout = QGridLayout() 
		INSTRTab.setLayout(self.INSTRlayout)
		
		# label to explain what happens on this instrument comparison tab
		l240 = QLabel(self)
		l240.setText('Use the tabs on either side to convolve model output into a form comparable with the corresponding instruments')
		l240.setWordWrap(True)
		self.INSTRlayout.addWidget(l240, 0, 0, 1, 1)
		
		INSTRtabs = QTabWidget()
		INSTRtabs.addTab(self.CPCtab(), "CPC")
		INSTRtabs.addTab(self.SMPStab(), "SMPS")
		INSTRtabs.addTab(self.CIMStab(), "CIMS")
		self.INSTRlayout.addWidget(INSTRtabs, 1, 0, 1, 1)
		INSTRtabs.setTabPosition(2)
		
		# -------------------------------------------------------------------------
	
		return(INSTRTab)
	
	def CPCtab(self): # instrument comparison plotting tab definition
	
		self.CPCTab = QWidget() # define tab widget
		self.CPClayout = QGridLayout() 
		self.CPCTab.setLayout(self.CPClayout)
		
		# create scrollable widget inside tab
		self.scroll = QScrollArea()
		self.scroll.setWidgetResizable(True)
		self.scrollwidget = QWidget()
		self.CPCscrolllayout  = QGridLayout()
		
		# text explaining purpose of CPC tab
		# label to explain what happens on this instrument comparison tab
		l220 = QLabel(self)
		l220.setText('The Condensation Particle Counter instrument and its associated software may have corrected for coincidence (the default setting here assumes this is the case).  Other settings here are essential to compare simulation results with instrument results.  Although all settings have a default, this should be checked against the value for the relevant instrument and its operation.')
		l220.setWordWrap(True)
		self.CPCscrolllayout.addWidget(l220, 0, 0, 1, 10)
	
		# input for equilibrium humidity on reaching the condensing 
		# section of the condensation particle counter
		self.e220 = QTextEdit(self)
		self.e220.setText('Relative humidity (fraction (0-1)) on reaching CPC condensing unit (particles assumed to equilibrate to this) (defaults to 0.65)')
		self.CPCscrolllayout.addWidget(self.e220, 1, 0, 1, 3)
		
		# input for minimum particle concentration (false background counts) detectable
		self.e220_a = QTextEdit(self)
		self.e220_a.setText('False background counts (used as the minimum detectable particle concentration) (# particles cm<sup>-3</sup>, defaults to 1.e-2)')
		self.CPCscrolllayout.addWidget(self.e220_a, 2, 0, 1, 3)
		
		# input for maximum particle concentration detectable
		self.e220_b = QTextEdit(self)
		self.e220_b.setText('Maximum detectable particle concentration (# particles cm<sup>-3</sup> (e.g. 3.e5), defaults to -1, which implies no maximum)')
		self.CPCscrolllayout.addWidget(self.e220_b, 3, 0, 1, 3)
		
		# input for coincidence correction
		self.e220_j = QTextEdit(self)
		self.e220_j.setText('Coincidence inputs: volumetric flow rate through counting unit (cm<sup>3</sup> s<sup>-1</sup>), instrument dead time (s) and upper limit of actual particle concentration (# particles cm<sup>-3</sup>) this can be applied to; should be three numbers separated by a comma (e.g. 5., 2.e-6, 3.e5).  Defaults to -1, -1, -1, which indicates no convolution needed for coincidence (e.g. because coincidence already corrected for by instrument)')
		self.CPCscrolllayout.addWidget(self.e220_j, 4, 0, 1, 3)
		
		# input for detection efficiency curve of counter
		self.e220_c = QTextEdit(self)
		self.e220_c.setText('Particle diameter (nm) at 50 % detection efficiency (>0 nm), factor for detection efficiency dependency on particle size (>0) (determines the range of particle sizes affected by reduced detection efficiency and assumes a sigmoid function.  Defaults to 5, 0.5)')
		self.CPCscrolllayout.addWidget(self.e220_c, 1, 3, 1, 3)
		
		# input for maximum detectable size of particle
		self.e220_d = QTextEdit(self)
		self.e220_d.setText('Maximum detectable particle diameter (nm), e.g. 3.e5.  Defaults to -1 which implies no maximum')
		self.CPCscrolllayout.addWidget(self.e220_d, 2, 3, 1, 3)
		
		# input for uncertainty in total particle number concentration (%)
		self.e220_e = QTextEdit(self)
		self.e220_e.setText('Uncertainty (%) around total particle number concentration, defaults to 10')
		self.CPCscrolllayout.addWidget(self.e220_e, 3, 3, 1, 3)
		
		# input for response time function
		self.e220_f = QTextEdit(self)
		self.e220_f.setText('Inputs for accounting for instrument response time and any mixing of particles from different times entering inlet (five inputs in total, all separated by a comma: i) shortest delay (s) in particles reaching counting unit, ii) delay at which weighting of particles at a maximum (s), iii) function of weighting against delay time (s) (use t for time and np for numpy) for particles between the shortest delay (i) and the delay at which weighting at maximum (ii), iv) longest delay (s) in particles reaching counting unit, v) function of weighting against delay time (s)  (use t for time and np for numpy) for particles between the delay at which weighting at maximum (ii) and the longest delay (iv).  For example: 0.1, 1.3, np.exp(t), 2.5, np.exp(np.flip(t-1.3)).  Note that weighting is normalised by the integral so that the final integral is one.  Defaults to: 1., 1., 1.*t, 1., 1.*t, which represents a response time of 1 s with no mixing of particles of different ages.')
		self.CPCscrolllayout.addWidget(self.e220_f, 4, 3, 1, 3)
		
		# input for frequency of instrument output
		self.e220_g = QTextEdit(self)
		self.e220_g.setText('Frequency of instrument output (Hz), defaults to 1.')
		self.CPCscrolllayout.addWidget(self.e220_g, 1, 7, 1, 3)
		
		# input for particle loss during inlet passage
		self.e220_h = QTextEdit(self)
		self.e220_h.setText('Loss rate (fraction s<sup>-1</sup>) as a function of particle size (um) (using Dp for diameter (um), np for numpy functions and python math symbols for math functions); time of passage through inlet (s).  E.g.: np.append(10.**(-5.5-0.5*np.log10(Dp[Dp &lt;= 1.e-1])), 10.**(-4.2+0.8*np.log10(Dp[Dp &gt; 1.e-1]))); 5..  Defaults to 0.; 0., which implies no particle losses in inlet')
		self.CPCscrolllayout.addWidget(self.e220_h, 2, 7, 1, 3)
		
		# input for averaging interval (s)
		self.e220_i = QTextEdit(self)
		self.e220_i.setText('Averaging interval (s).  Defaults to 1.')
		self.CPCscrolllayout.addWidget(self.e220_i, 3, 7, 1, 3)
		
		# button to plot counting efficiency dependence on particle size 
		self.b220_0 = QPushButton('Counting \nefficiency \ndependence \non size', self)
		self.b220_0.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.b220_0.setToolTip('Counting efficiency dependence on particle size')
		self.b220_0.clicked.connect(self.on_click234_a)
		self.CPCscrolllayout.addWidget(self.b220_0, 1, 6)
		
		# button to plot weighting as a function of response time
		self.b220_f = QPushButton('Weighting \ndependency \non response \ntime', self)
		self.b220_f.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.b220_f.setToolTip('Plot the weighting of particles by age due to the response time function')
		self.b220_f.clicked.connect(self.on_click222_b)
		self.CPCscrolllayout.addWidget(self.b220_f, 4, 6)
		
		# button to plot inlet loss rate as a function of particle diameter
		self.b220_h = QPushButton('Inlet loss \nrate with \nparticle \ndiameter', self)
		self.b220_h.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.b220_h.setToolTip('Plot the inlet loss rate as a function of particle diameter')
		self.b220_h.clicked.connect(self.on_click222_h)
		self.CPCscrolllayout.addWidget(self.b220_h, 2, 10)
		
		# button to plot temporal profile of total particle number concentration
		self.b220_a = QPushButton('CPC \nobservations', self)
		self.b220_a.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.b220_a.setToolTip('Plot the total particle number concentration as observed by a condensation particle counter')
		self.b220_a.clicked.connect(self.on_click222)
		self.CPCscrolllayout.addWidget(self.b220_a, 4, 10)
		
		# properties of CPC scroll area ----------------
		self.scrollwidget.setLayout(self.CPCscrolllayout)
		self.scroll.setWidget(self.scrollwidget)
		self.CPClayout.addWidget(self.scroll, 0, 0, 1, 1)

		return(self.CPCTab)
	
	def SMPStab(self): # instrument comparison plotting tab definition
	
		SMPSTab = QWidget()
		self.SMPSlayout = QGridLayout() 
		SMPSTab.setLayout(self.SMPSlayout)
	
		# create scrollable widget inside tab
		self.scroll = QScrollArea()
		self.scroll.setWidgetResizable(True)
		self.scrollwidget = QWidget()
		self.SMPSscrolllayout  = QGridLayout()
		
		# text explaining purpose of SMPS tab
		# label to explain what happens on this instrument comparison tab
		l230 = QLabel(self)
		l230.setText('The Scanning Mobility Particle Spectrometer instrument and its associated software may have corrected for the instrument characteristics provided here.  Although all settings have a default, this should be checked against the value for the relevant instrument and its operation.')
		l230.setWordWrap(True)
		self.SMPSscrolllayout.addWidget(l230, 0, 0, 1, 10)
		
		# input for equilibrium humidity on reaching the condensing 
		# section of the condensation particle counter
		self.e230_b = QTextEdit(self)
		self.e230_b.setText('Relative humidity (fraction (0-1)) on reaching sizing unit (particles assumed to equilibrate to this) (defaults to 0.65)')
		self.SMPSscrolllayout.addWidget(self.e230_b, 1, 0, 1, 3)
		
		# input for minimum particle concentration (false background counts) detectable
		self.e230_c = QTextEdit(self)
		self.e230_c.setText('False background counts (used as the minimum detectable particle concentration) (# particles cm<sup>-3</sup>, defaults to 1.e-2)')
		self.SMPSscrolllayout.addWidget(self.e230_c, 2, 0, 1, 3)
		
		# input for maximum particle concentration detectable
		self.e230_d = QTextEdit(self)
		self.e230_d.setText('Maximum detectable particle concentration (# particles cm<sup>-3</sup> (e.g. 3.e5), defaults to -1, which implies no maximum)')
		self.SMPSscrolllayout.addWidget(self.e230_d, 3, 0, 1, 3)
		
		# input for coincidence correction
		self.e230_e = QTextEdit(self)
		self.e230_e.setText('Coincidence inputs: volumetric flow rate through counting unit (cm<sup>3</sup> s<sup>-1</sup>), instrument dead time (s) and upper limit of actual particle concentration (# particles cm<sup>-3</sup>) this can be applied to; should be three numbers separated by a comma (e.g. 5., 2.e-6, 3.e5).  Defaults to -1, -1, -1, which indicates no convolution needed for coincidence (e.g. because coincidence already corrected for by instrument)')
		self.SMPSscrolllayout.addWidget(self.e230_e, 4, 0, 1, 3)
		
		# input for detection efficiency curve of counter
		self.e230_f = QTextEdit(self)
		self.e230_f.setText('Particle diameter (nm) at 50 % detection efficiency (>0 nm), factor for detection efficiency dependency on particle size (>0) (determines the range of particle sizes affected by reduced detection efficiency and assumes a sigmoid function.  Defaults to 5, 0.5)')
		self.SMPSscrolllayout.addWidget(self.e230_f, 1, 3, 1, 3)
		
		# input for maximum detectable size of particle
		self.e230_g = QTextEdit(self)
		self.e230_g.setText('Maximum detectable particle diameter (nm), e.g. 3.e5.  Defaults to -1 which implies no maximum')
		self.SMPSscrolllayout.addWidget(self.e230_g, 2, 3, 1, 3)
		
		# input for uncertainty in total particle number concentration (%)
		self.e230_h = QTextEdit(self)
		self.e230_h.setText('Uncertainty (%) around total particle number concentration, defaults to 10')
		self.SMPSscrolllayout.addWidget(self.e230_h, 3, 3, 1, 3)
		
		# input for response time function
		self.e230_i = QTextEdit(self)
		self.e230_i.setText('Inputs for accounting for instrument response time and any mixing of particles from different times entering inlet (five inputs in total, all separated by a comma: i) shortest delay (s) in particles reaching counting unit, ii) delay at which weighting of particles at a maximum (s), iii) function of weighting against delay time (s) (use t for time and np for numpy) for particles between the shortest delay (i) and the delay at which weighting at maximum (ii), iv) longest delay (s) in particles reaching counting unit, v) function of weighting against delay time (s)  (use t for time and np for numpy) for particles between the delay at which weighting at maximum (ii) and the longest delay (iv).  For example: 0.1, 1.3, np.exp(t), 2.5, np.exp(np.flip(t-1.3)).  Note that weighting is normalised by the integral so that the final integral is one.  Defaults to: 1., 1., 1.*t, 1., 1.*t, which represents a response time of 1 s with no mixing of particles of different ages.')
		self.SMPSscrolllayout.addWidget(self.e230_i, 4, 3, 1, 3)
		
		# input for frequency of instrument output
		self.e230_j = QTextEdit(self)
		self.e230_j.setText('Frequency of instrument output (Hz), defaults to 1.')
		self.SMPSscrolllayout.addWidget(self.e230_j, 1, 7, 1, 3)
		
		# input for particle loss during passage through instrument
		self.e230_k = QTextEdit(self)
		self.e230_k.setText('Loss rate (fraction s<sup>-1</sup>) as a function of particle size (um) (using Dp for diameter (um), np for numpy functions and python math symbols for math functions); time of passage through inlet (s).  E.g.: np.append(10.**(-5.5-0.5*np.log10(Dp[Dp &lt;= 1.e-1])), 10.**(-4.2+0.8*np.log10(Dp[Dp &gt; 1.e-1]))); 5..  Defaults to 0., 0., which implies no particle losses in inlet')
		self.SMPSscrolllayout.addWidget(self.e230_k, 2, 7, 1, 3)
		
		# input for averaging interval (s)
		self.e230_l = QTextEdit(self)
		self.e230_l.setText('Averaging interval (s).  Defaults to 1.')
		self.SMPSscrolllayout.addWidget(self.e230_l, 3, 7, 1, 3)
		
		# input for number of channels per decade (channels means size bins)
		self.e230_m = QTextEdit(self)
		self.e230_m.setText('Number of channels per decade of particle size.  Defaults to 128.')
		self.SMPSscrolllayout.addWidget(self.e230_m, 4, 7, 1, 3)
		
		# button to plot counting efficiency dependence on particle size 
		self.b230_i = QPushButton('Counting \nefficiency \ndependence \non size', self)
		self.b230_i.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.b230_i.setToolTip('Counting efficiency dependence on particle size')
		self.b230_i.clicked.connect(self.on_click230_i)
		self.SMPSscrolllayout.addWidget(self.b230_i, 1, 6)
		
		# button to plot weighting as a function of response time
		self.b230_j = QPushButton('Weighting \ndependency \non response \ntime', self)
		self.b230_j.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.b230_j.setToolTip('Plot the weighting of particles by age due to the response time function')
		self.b230_j.clicked.connect(self.on_click230_j)
		self.SMPSscrolllayout.addWidget(self.b230_j, 4, 6)
		
		# button to plot inlet loss rate as a function of particle diameter
		self.b230_k = QPushButton('Inlet loss \nrate with \nparticle \ndiameter', self)
		self.b230_k.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.b230_k.setToolTip('Plot the inlet loss rate as a function of particle diameter')
		self.b230_k.clicked.connect(self.on_click230_k)
		self.SMPSscrolllayout.addWidget(self.b230_k, 2, 10)
		
		# button to plot temporal profile of number size distribution
		self.b230_m = QPushButton('SMPS observations', self)
		self.b230_m.setToolTip('Plot the number size distribution as observed by a particle counter')
		self.b230_m.clicked.connect(self.on_click230_m)
		self.SMPSscrolllayout.addWidget(self.b230_m, 4, 10)
		
		# properties of SMPS scroll area ----------------
		self.scrollwidget.setLayout(self.SMPSscrolllayout)
		self.scroll.setWidget(self.scrollwidget)
		self.SMPSlayout.addWidget(self.scroll, 0, 0, 1, 1)
		
		return(SMPSTab)
	
	def CIMStab(self): # instrument comparison plotting tab definition
	
		CIMSTab = QWidget()
		self.CIMSlayout = QGridLayout() 
		CIMSTab.setLayout(self.CIMSlayout)
	
		# create scrollable widget inside tab
		self.scroll = QScrollArea()
		self.scroll.setWidgetResizable(True)
		self.scrollwidget = QWidget()
		self.CIMSscrolllayout  = QGridLayout()
	
		# input for mass:charge resolution
		self.e280 = QTextEdit(self)
		self.e280.setText(str('Mass:charge resolution determinants: ' +
		'interval of probability peaks (m:z), width of distribution ' +
		'per peak.  Defaults to 1., 0.3'))
		self.CIMSscrolllayout.addWidget(self.e280, 0, 0)
		
		# input for time to show mass spectrum for
		self.e281 = QTextEdit(self)
		self.e281.setText(str('Time through experiment to show mass ' + 
		'spectrum for (s), set to \'all times\' to integrate over all times'))
		self.CIMSscrolllayout.addWidget(self.e281, 1, 0)
		
		# type of ionisation source
		self.e282 = QTextEdit(self)
		self.e282.setText(str('Ionisation source (I for iodide, ' +
		'N for nitrate, defaults to I), whether to add molar mass ' +
		'to mass of ionised components (1 for yes, 0 for no), ' +
		'defaults to 0 for no'))
		self.CIMSscrolllayout.addWidget(self.e282, 2, 0)
		
		# sensitivity dependence on molar mass
		self.e283 = QTextEdit(self)
		self.e283.setText(str('Sensitivity (instrument count/real count) ' +
		'dependence on molar mass (g/mol), use y_MM to denote molar ' +
		'mass (g/mol) of components.  Defaults to 1.0 which implies ' +
		'no dependency on molar mass. To zero ranges of m/z use ' +
		'inequalities, e.g. <200. sets all m/z less than 200 to 0.' +
		'To zero inorganics, follow any (or blank) molar mass functions ' +
		'with a comma and the words organics only, e.g.: <73., organics only'))
		self.CIMSscrolllayout.addWidget(self.e283, 3, 0)
		
		# button to plot probability distribution function 
		# demonstrating mass:charge resolution
		self.b290_aa = QPushButton('PDF graph', self)
		self.b290_aa.setToolTip(str('Plot the probability ' +
		'distribution function causing mass:charge resolution'))
		self.b290_aa.clicked.connect(self.on_click290_aa)
		self.CIMSscrolllayout.addWidget(self.b290_aa, 0, 1)	

		# button to plot sensitivity dependence on molar mass
		self.b290_a = QPushButton('Sensitivity', self)
		self.b290_a.setToolTip('Plot the sensitivity to molar mass')
		self.b290_a.clicked.connect(self.on_click290_a)
		self.CIMSscrolllayout.addWidget(self.b290_a, 3, 1)
	
		# drop-down button to select linear or logarithmic y-axis	
		self.b290_ab = QComboBox(self)
		self.b290_ab.addItem('Linear y')
		self.b290_ab.addItem('Logarithmic y')	
		self.b290_ab.setToolTip(str('Select whether to have a linear ' +
		'or logarithmic spacing on the CIMS y-axis'))	
		self.CIMSscrolllayout.addWidget(self.b290_ab, 4, 0)	

		# drop-down button to select histogram or markers	
		self.b290_abb = QComboBox(self)
		self.b290_abb.addItem('Stem (molecules/cm3)')
		self.b290_abb.addItem('Markers (molecules/cm3)')
		self.b290_abb.addItem('Stem (normalised)')
		self.b290_abb.addItem('Markers (normalised)')	
		self.b290_abb.setToolTip(str('Select whether to have ' +
		'bars or markers layout'))	
		self.CIMSscrolllayout.addWidget(self.b290_abb, 4, 1)

		# drop-down button to select which phase is plotted	
		self.b290_abc = QComboBox(self)
		self.b290_abc.addItem('Gas')
		self.b290_abc.addItem('Particle')
		self.b290_abc.addItem('Gas and Particle')	
		self.b290_abc.setToolTip(str('Select the phase(s) to plot'))	
		self.CIMSscrolllayout.addWidget(self.b290_abc, 5, 0)

		# button to plot mass spectrum
		self.b290 = QPushButton('CIMS mass spectrum', self)
		self.b290.setToolTip(str('Plot the mass spectrum as observed ' +
		'by a chemical ionisation mass spectrometer'))
		self.b290.clicked.connect(self.on_click290)
		self.CIMSscrolllayout.addWidget(self.b290, 6, 0)
		
		# button to save output in CIMS form
		self.b290_aaa = QPushButton('Save in CIMS form', self)
		self.b290_aaa.setToolTip('Save the results in CIMS format at all time steps')
		self.b290_aaa.clicked.connect(self.on_click290_aaa)
		self.CIMSscrolllayout.addWidget(self.b290_aaa, 6, 1)

		# properties of CIMS scroll area ----------------
		self.scrollwidget.setLayout(self.CIMSscrolllayout)
		self.scroll.setWidget(self.scrollwidget)
		self.CIMSlayout.addWidget(self.scroll, 0, 0, 1, 1)
	
		return(CIMSTab)
	
	def OBStab(self): # comparing with observations

		OBSTab = QWidget()
		self.OBSlayout = QGridLayout() 
		OBSTab.setLayout(self.OBSlayout)
	
		# create scrollable widget inside tab
		self.scroll = QScrollArea()
		self.scroll.setWidgetResizable(True)
		self.scrollwidget = QWidget()
		self.OBSscrolllayout  = QGridLayout()

		# show path to observations file ----------------------------
		self.l401 = ScrollLabel(self)
		cwd = os.getcwd() # current working directory
		path = str(cwd + '/PyCHAM/output/26_10.xlsx')
		self.l401.setText(path)
		self.OBSscrolllayout.addWidget(self.l401, 0, 0, 1, 1)

		# select observations file -----------------------------------
		b400 = QPushButton('Select file containing observations', self)
		b400.setToolTip('Select the file containing the required observations')
		b400.clicked.connect(self.on_click400)
		self.OBSscrolllayout.addWidget(b400, 0, 1, 1, 1)

		# drop down button for what to plot
		self.b402 = QComboBox(self)
		self.b402.addItem('Gas-phase Concentrations (ppb)')
		self.b402.addItem(str('Particle Concentrations (' + u'\u03BC' + 'g/m' + u'\u00B3' + ')  with Standard Results Plot'))
		self.b402.addItem(str('Particle Concentrations (' + u'\u03BC' + 'g/m' + u'\u00B3' +') with Cumulative particle mass concentration without water plot'))
		self.b402.addItem('Van Krevelen (Time Profile, Averaged Over All Hydrocarbons)')
		self.b402.addItem('Van Krevelen (Time Profile, Averaged Over Non-methane Hydrocarbons)')
		self.b402.addItem('Van Krevelen (Time Profile, Averaged Over _ Extension Hydrocarbons)')
		self.b402.addItem('Van Krevelen (All Individual Hydrocarbons at The Time Through Experiment (s) Provided Below)')
		self.b402.addItem('Van Krevelen (Non-methane Individual Hydrocarbons at The Time Through Experiment (s) Provided Below)')
		self.b402.addItem('Van Krevelen (_ Extension Individual Hydrocarbons at The Time Through Experiment (s) Provided Below)')
		self.b402.addItem('Mass Defect of All Components in Chemical Scheme')
		self.b402.addItem('Mass Defect of All Hydrocarbons Scaled to Concentrations at Time Through Experiment (s) Provided Below')
		self.b402.addItem('CIMS mass spectrum (as defined in Plot/Convolution/CIMS)')
		self.OBSscrolllayout.addWidget(self.b402, 1, 0, 1, 1)

		# input for observed and modelled plotting options
		self.e403 = QTextEdit(self)
		self.e403.setText(str('Provide Plotting Options as ' +
		'Described in the Selected Drop-down Button Option Above'))
		self.OBSscrolllayout.addWidget(self.e403, 2, 0, 1, 1)

		# plot observations and model results
		# select observations file ---------------------------------
		b401 = QPushButton('Observed && Modelled', self)
		b401.setToolTip('Plot Observations and Model Results')
		b401.clicked.connect(self.on_click401)
		self.OBSscrolllayout.addWidget(b401, 3, 0, 1, 1)
		
		# properties of observations scroll area ----------------
		self.scrollwidget.setLayout(self.OBSscrolllayout)
		self.scroll.setWidget(self.scrollwidget)
		self.OBSlayout.addWidget(self.scroll, 0, 0, 1, 1)

		return(OBSTab)

	def EXPtab(self): # preparing for experiments

		EXPTab = QWidget()
		self.EXPlayout = QGridLayout() 
		EXPTab.setLayout(self.EXPlayout)

		# select observations file --------------------------------------------------------------
		b500 = QPushButton('Output', self)
		b500.setToolTip('Estimates the outputs typically useful for experiment preparation and save in excel file inside same folder as selected results stored in')
		b500.clicked.connect(self.on_click500)
		self.EXPlayout.addWidget(b500, 0, 0, 1, 1)
		
		return(EXPTab)

	@pyqtSlot() # PyCHAM website
	def on_clickn00(self): # EUROCHAMP website
		import webbrowser
		webbrowser.open('http://www.github.com/simonom/PyCHAM')
		return()

	@pyqtSlot() # eurochamp website under development 10/02/2021
	def on_clickn1(self): # EUROCHAMP website
		import webbrowser
		webbrowser.open('https://www.eurochamp.org')
		return()
		
	@pyqtSlot() 
	def on_clickn1a(self): # NCAS website
		import webbrowser
		webbrowser.open('https://ncas.ac.uk')
		return()
	
	@pyqtSlot() 
	def on_clickn1b(self): # University of Manchester website
		import webbrowser
		webbrowser.open('https://www.manchester.ac.uk')
		return()
	
	@pyqtSlot()
	def on_click1(self): # selecting folder containing input files
	
		
		if (self.fab == 1): # if showing, remove single simulation widgets
			self.b81.deleteLater()
			self.l81.deleteLater()
			# remember that single simulation widgets not showing
			self.fab = 0
		if (self.atb == 1): # if showing then remove add to batch button
			self.b82.deleteLater()
			self.l81.deleteLater() # 'or' label
			self.atb = 0 # remember that add to batch button not showing
		# remove any old message from previous run
		if (len(self.l81b.text()) > 0):
			if (self.l81b.text()[0:17] != 'File combinations'):
				self.l81b.setText('')
				self.l81b.setStyleSheet(0., '0px', 0., 0.)
		
		# prepare by enforcing default variables
		# default variables for all required input model 
		# variables -------------------------------------------- 
		[y0, Press, siz_stru, num_sb, 
		lowsize, uppsize, std, 
		Compt, 
		injectt, Ct, seed_mw, seed_diss, seed_dens, 
		dens_comp, dens, vol_comp, volP, act_comp, 
		act_user, accom_comp, accom_val, uman_up, int_tol, 
		new_partr, 
		coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, 
		chamSA, Rader, 
		p_char, e_field, ser_H2O, wat_hist, 
		drh_str, erh_str, 
		z_prt_coeff, chamV, self] = def_mod_var.def_mod_var(0, 
			self)
		
		# then open default variables, ready for modification
		input_by_sim = str(self.PyCHAM_path + 
			'/PyCHAM/pickle.pkl')
		
		with open(input_by_sim, 'rb') as pk:
			[y0, Press,
			siz_stru, num_sb, 
			lowsize, uppsize, std, 
			Compt, injectt, Ct,
			seed_mw, seed_diss, seed_dens,
			dens_comp, dens, vol_comp, volP, act_comp, 
			act_user, 
			accom_comp, accom_val, uman_up, int_tol, 
			new_partr, coag_on, inflectDp, pwl_xpre, 
			pwl_xpro, 
			inflectk, chamSA, Rader, p_char, e_field, 
			ser_H2O, 
			wat_hist, drh_str, erh_str, z_prt_coeff, 
			chamV] = pickle.load(pk)
			pk.close()
		
		# button to get path to folder/file containing 
		# relevant files
		options = QFileDialog.Options()
		#QFileDialog.getOpenFileName
		fol_nme = getExistingFilesAndDirs(self, "Select Input", "./PyCHAM/input/")
		
		if (fol_nme == []): # if no folder selected (e.g. because selection cancelled)
			return()
		fol_nme = fol_nme[0]
		
		# unknown names
		self.sch_name = self.xml_name = self.inname = 'Not found'
		
		try: # in case this is a directory
			# list of files (and any folders) here
			dir_con = os.listdir(fol_nme)
		
			# look for corresponding files here
			for i in dir_con:
				if ('chem' in i): # chemical scheme file
					self.sch_name = str(fol_nme+'/'
						+i)
				if ('xml' in i): # xml file
					self.xml_name = str(fol_nme+'/'
						+i)
				if ('var' in i): # model variables file
					self.inname = str(fol_nme+'/'+i)
		
		# if not a directory then assume it is a model 
		# variable file
		except:
			self.inname = fol_nme

		# read in model variables of this model variables file 
		# and store to pickle
		# Note that doing this call after searching for 
		# files with names containing 'chem' means that any 
		# chemical 
		# scheme files specified in the model variables overide
		# any found by automatic search
		import mod_var_read
		mod_var_read.mod_var_read(self)
		
		# end this function if an error thrown by reading of 
		# model variables
		if (self.bd_st == 1 or self.bd_st == 2):
			return()
			
		# updating scroll labels showing path to files
		self.l3.setText(self.sch_name)
		self.l5.setText(self.xml_name)
		self.l7.setText(self.inname)
			
		self.show()
		
		# update displayed model variables to those of the selected 
		# model variables file
		import mod_var_up
		mod_var_up.mod_var_up(self)
		
		# running check on model variables -----------------------------------------
		import ui_check # module for checking on model variables
		
		# let checking module know that this is the call once user-defined
		# model variables have been read in (in contrast to the first call
		# to checking module, which is just for default variables)
		self.chck_num = 2

		# check on inputs - note this loads the last saved pickle 
		# file and saves any change 
		# to this pickle file
		
		ui_check.ui_check(self)

		# move check number up by one
		self.chck_num = 3
		# finished check on model variables -------------------------------------------
		
		return()
		
	@pyqtSlot()
	def on_click2(self): # when different chemical scheme file requires selection
	
		if type(self.param_const) != dict: # if not called from automatic script (automated_setup_and_call.py)
			self.sch_name, _ = QFileDialog.getOpenFileName(self, "Select Chemical Scheme File", "./PyCHAM/input/") # get path of file
		else: # if called from automatic script (automated_setup_and_call.py)
			self.sch_name = self.param_const['sch_name']
		self.l3.clear() # clear old label
		self.l3.setText(str(self.sch_name))
		self.l3.show()
	
		# running check on model variables -------------------------------------------
		import ui_check # module for checking on model variables
		# check on inputs - note this loads the last saved pickle 
		# file and saves any change to this pickle file
		ui_check.ui_check(self)
		# finished check on model variables ------------------------------------------
	
	@pyqtSlot()
	def on_click3(self): # when different xml file requires selection
		
		if type(self.param_const) != dict: # if not called from automatic script (automated_setup_and_call.py)
			self.xml_name, _ = QFileDialog.getOpenFileName(self, "Select xml File", "./PyCHAM/input/") # get path of file
		else: # if called from automatic script (automated_setup_and_call.py)
			self.xml_name = self.param_const['xml_name']

		self.l5.clear() # clear old label
		self.l5.setText(str(self.xml_name))
		self.l5.show()
			
	@pyqtSlot()
	# when different model variables file requires selection
	def on_click4(self):
	
		if (self.fab == 1): # if showing, remove single simulation widgets
			self.b81.deleteLater()
			self.l81.deleteLater()
			self.fab = 0 # remember that single simulation widgets not showing
		if (self.atb == 1): # if showing then remove add to batch button
			self.b82.deleteLater()
			self.l81.deleteLater() # 'or' label
			self.atb = 0 # remember that add to batch button not showing
		# remove any old 'Simulation complete' message from previous run
		if ((self.l81b.text() == 'Simulation complete') or 
		(self.l81b.text() == 'Simulations complete')):
			self.l81b.setText('')

		# if not called from automatic script (automated_setup_and_call.py)
		if type(self.param_const) != dict:
			# user chooses path of file to model variables file
			self.inname, _ = QFileDialog.getOpenFileName(self, "Select Model Variables File", "./PyCHAM/input/")
		# if run from automatic script (automated_setup_and_call.py)
		if (type(self.param_const) == dict and self.inname == 'Not found'):
			self.inname = 'Automatically Set'

		# if no file selected, e.g. because selection was cancelled
		if (self.inname == ''):
			return()
		
		self.l7.clear() # clear old label
		self.l7.setText(self.inname)
		self.l7.show()
	
		# read in model variables of this model variables file and store to pickle
		import mod_var_read
		mod_var_read.mod_var_read(self)
		
		# end this function if an error thrown by reading of model variables
		if (self.bd_st == 1 or self.bd_st == 2):
			return()
		# update displayed model variables to those of the 
		# selected model variables file
		import mod_var_up
		mod_var_up.mod_var_up(self)
		
		# running check on model variables -------------------------------------------
		import ui_check # module for checking on model variables

		# let checking module know that this call comes after user-defined
		# model variables have been read in
		self.chck_num = 2

		# check on inputs - note this loads the last saved pickle file 
		# and saves any change
		# to this pickle file
		ui_check.ui_check(self)
		# finished check on model variables ---------------------------------

		return()
			
	
	@pyqtSlot()	
	def on_click81b(self): # when button to run simulation pressed
		
		for sim_num in range(self.btch_no):
			
			# --------------------------------------------
			# reach end of list in batch mode
			if ((self.btch_no > 1) and (sim_num == self.btch_no-1)):
				# once all simulations done tidy up
				# clear the list of file combinations for simulations in batch
				self.output_list = [] # reset list of output paths
				# return to single simulation mode
				self.btch_no = 1
				self.btch_str = str('File combinations included ' +
				'in batch: chemical scheme, xml, model variables\n')
				
				# remove old progress message
				self.l81b.setText('')
				# tell user that simulations finished
				self.l81b.setText(str('Simulations complete'))
				return(err_mess)
			
			# --------------------------------------------
		
			# reset to default variables to allow any new 
			# variables to arise
			# from the current model variables file only
			[y0, Press, 
			siz_stru, num_sb, 
			lowsize, uppsize, 
			std, Compt, injectt, 
			Ct, seed_mw, seed_diss, seed_dens,  
			dens_comp, dens, vol_comp, volP, act_comp, 
			act_user, accom_comp, 
			accom_val, uman_up, int_tol, new_partr, 
			coag_on, inflectDp, pwl_xpre, pwl_xpro, 
			inflectk, chamSA, Rader, p_char, e_field, 
			ser_H2O, wat_hist, drh_str, erh_str, 
			z_prt_coeff, chamV, 
			self] = def_mod_var.def_mod_var(0, self)
			
			# get text from batch list label
			btch_list = self.btch_str
			if ((sim_num+1) < (self.btch_no-1)): # if prior to final item in list
				# get location of text relevant to this simulation
				txt_st = (btch_list.find(str(str(sim_num+1) + '.' + 
				'\n')) + len(str(str(sim_num+1)+'.'+'\n')))
				txt_fi = btch_list.find(str(str(sim_num+2)+'.'+'\n'))
				txtn = btch_list[txt_st:txt_fi]
			else: # if final item in list
				# get location of text relevant to this simulation
				txt_st = (btch_list.find(str(str(sim_num+1) + '.' + 
				'\n'))+len(str(str(sim_num+1)+'.'+'\n')))
				txtn = btch_list[txt_st::]
					
			# split into separate file paths
			txtn = txtn.split('\n')
				
			# update file names
			self.sch_name = txtn[0]
			self.xml_name = txtn[1]
			self.inname = txtn[2]
			
			# read in model variables of this model variables file 
			# (as identified by self.inname) and store to pickle
			import mod_var_read
			mod_var_read.mod_var_read(self)
			
			# end this function if an error thrown by reading of model variables
			if (self.bd_st == 1 or self.bd_st == 2):
				return()

			# get the save path name variables
			input_by_sim = str(self.PyCHAM_path 
				+ '/PyCHAM/pickle.pkl')
			with open(input_by_sim, 'rb') as pk:
				[y0, Press,
				siz_stru, num_sb,
				lowsize, uppsize, std,
				Compt, injectt, Ct,
				seed_mw, seed_diss, seed_dens,
				dens_comp, dens, vol_comp, volP, 
				act_comp, act_user, 
				accom_comp, accom_val, uman_up, int_tol, 				new_partr, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, 
				p_char, e_field, 
				ser_H2O, wat_hist, drh_str, erh_str, 
				z_prt_coeff, 
				chamV] = pickle.load(pk)
				pk.close() # close pickle file
			
			# let check know this is a second call
			self.chck_num = 2
			
			# run another check on inputs - means any changes 
			# made by default are set
			import ui_check; ui_check.ui_check(self)

			# reset check number to a first call
			self.chck_num = 1

			# saving path - copied from user_input.py module
			dir_path = os.getcwd() # current working directory
			output_root = 'PyCHAM/output'
			filename = os.path.basename(self.sch_name)
			filename = os.path.splitext(filename)[0]
			# if path and folder name already stated for saving to
			if (('/' in self.sav_nam) or ('\\' in self.sav_nam)): 
				# one folder for one simulation
				output_by_sim = os.path.join(self.sav_nam)
			else: # if no path given for saving to
				# one folder for one simulation
				output_by_sim = os.path.join(dir_path, output_root, 
				filename, self.sav_nam)
			
			# display the progress label, including the name of the saving path
			if (self.btch_no == 1): # single run mode
				self.l81b.setText('')
				self.l81b.setText(str('Progress through simulation that will save to: \n' + str(output_by_sim) + '\n' + str(sim_num+1) + ' of ' + str(self.btch_no)))
			if (self.btch_no > 1): # batch run mode
				self.l81b.setText('')
				self.l81b.setText(str('Progress through simulation that will save to: \n' + str(output_by_sim) + '\n' + str(sim_num+1) + ' of ' + str(self.btch_no-1)))
			
			# if showing, disable then remove start simulation button
			if (self.fab == 1):
				self.b81.setEnabled(False)
				self.b81.deleteLater()
				self.fab = 0
			
			# if showing, disable then remove add to batch button
			if (self.atb == 1): 
				self.b82.setEnabled(False)
				self.b82.deleteLater()
				# remove or label
				self.l81.deleteLater()
				self.atb = 0
			
			# hide checking inputs buttons
			self.b80.hide()
			self.b80s.hide()
			
			# path to error log
			err_log = str(self.PyCHAM_path + '/PyCHAM/err_log.txt')
			if (sim_num == 0): # delete any existing error log and create new log
				# list upcoming simulation in the error log
				with open(err_log, 'w') as el:
					el.write(str(output_by_sim+'\n'))
					el.close
			else: # append to existing log
				# list upcoming simulation in the error log
				with open(err_log, 'a') as el:
					el.write(str(output_by_sim+'\n'))
					el.close
			
			# tell numpy what to do if error observed
			import err_log_code
			log = err_log_code
			saved_handler = np.seterrcall(log)
			# note that error messages saved to the log file are less verbose than 
			# those printed to the command line, to get the command line version 
			# (which outputs to the command line)
			# change 'log' below to 'warn'
			save_err = np.seterr(all='log')
			
			# call on function to simulate
			err_mess = self.act_81(output_by_sim, sim_num)

			if (err_mess != ''): # state error message if any generated
				self.l81b.setText('') # remove old progress message
				if (self.btch_no > 1): # in batch mode
					self.l81b.setText(str('See error message below generated during simulation saving to: \n' + str(output_by_sim) + '\n' + str(sim_num+1) + ' of ' + str(self.btch_no-1) + '\n' + err_mess))
				if (self.btch_no == 1): # in single simulation mode
					self.l81b.setText(str('See error message below generated during simulation saving to: \n' + str(output_by_sim) + '\n' + str(sim_num+1) + ' of ' + str(1) + '\n' + err_mess))
				self.l81b.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.output_list = [] # reset list of output paths
				self.btch_no = 1 # reset to single simulation run mode (rather than batch) 
				return(err_mess)
		
		# if no error message then return with no error message
		return(err_mess)
				
	
	@pyqtSlot()
	# start the simulation
	def act_81(self, output_by_sim, sim_num):
		
		from middle import middle # prepare to communicate with main program
		
		note_messf = 0 # cancel note message flag
		
		for prog in middle(self): # call on modules to solve problem
			
		
			if (isinstance(prog, str)): # check if it's a message
				mess = prog
				# if it's an error message
				if (mess[0:5] == 'Error'):
					# remove the progress bar
					self.progress.deleteLater()
					return(mess)
					
					# remove the directory that was due to hold results
					import shutil	
					shutil.rmtree(self.output_by_sim)
				else:
					self.l81b.setText(mess)
					# flag that a note message has been generated
					note_messf = 1 
			
			else: # if no message from model
			
				# if there was previously a note message
				if (note_messf == 1):
					note_messf = 0 # cancel note message flag
				
					if (self.btch_no == 1): # single run mode
						self.l81b.setText('')
						self.l81b.setText(str('Progress through simulation that will save to: \n' + str(output_by_sim) + '\n' + str(sim_num+1) + ' of ' + str(self.btch_no)))
					if (self.btch_no > 1): # batch run mode
						self.l81b.setText('')
						self.l81b.setText(str('Progress through simulation that will save to: \n' + str(output_by_sim) + '\n' + str(sim_num+1) + ' of ' + str(self.btch_no-1)))
			
				self.progress.setValue(int(prog)) # get progress
				
			QApplication.processEvents() # allow progress bar/message panel to update
		

		# remove the progress bar after each simulation
		self.progress.deleteLater()

		# set the path to folder to plot results to the latest simulation results
		self.l201.setText(output_by_sim)
		
		# tell user that output ready to plot
		self.l203a.setText('Output ready to plot')
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.bd_pl = 3
		
		# if this point reached then no error message generated
		mess = ''
		return(mess)
	
	@pyqtSlot()
	def on_click81sing(self): # when single simulation button pressed
		
		# note it can be useful to process events now, in case PyCHAM working in automatic mode
		QApplication.processEvents() # allow progress bar/message panel to update

		cs_file = self.l3.text()
		xml_file = self.l5.text()
		mv_file = self.l7.text()
		self.btch_no = 1 # ensure just one simulation when single simulation button pressed
		self.btch_str = (str(str(self.btch_no)+'.\n' + cs_file + '\n' + xml_file + '\n' + mv_file + '\n'))
		
		# if called from autorun function and if later simulations 
		# being conducted (i.e. not the first)  
		if (type(self.param_const) == dict and 'progress' in self.__dict__):
			self.progress.setValue(0) # reset to zero
		else:
			# display progress bar
			self.progress = QProgressBar(self)
		self.NSlayout.addWidget(self.progress, 7, self.mvpn, 1, 3)
		
		
		if (self.fab == 1): # if showing, remove single simulation widgets
			self.b81.setEnabled(False)
			self.b81.deleteLater()
			self.fab = 0 # remember that single simulation widgets not showing
		if (self.atb == 1): # if showing, disable and remove add to batch button
			self.b82.setEnabled(False)
			self.b82.deleteLater()
			self.atb = 0 # remember that add to batch button not showing
			# remove or label
			self.l81.deleteLater()
		if (self.ssb == 1):# if showing, remove series of simulations button
			self.b83.setEnabled(False)
			self.b83.deleteLater()
			self.ssb = 0 # remember this button removed
		
		self.l81b.setText('') # clear progress message
			
		# action to simulate a single simulation
		err_mess = self.on_click81b()
			
		# once all simulations done in single simulation mode, tidy up
		# clear the list of file combinations for simulations in batch
		self.output_list = [] # reset list of output paths
		# return to single simulation mode
		self.btch_no = 1
		self.btch_str = 'File combinations included in batch: chemical scheme, xml, model variables\n'
		
		# show checking inputs buttons
		self.b80.show()
		self.b80s.show()

		if (err_mess == ''): # if no error message generated
			# tell user that simulations finished
			self.l81b.setText('')
			self.l81b.setText(str('Simulation complete'))
	
	@pyqtSlot()
	def on_click81(self): # when 'start series of simulation' button pressed
		
		
		# display progress bar
		self.progress = QProgressBar(self)
		self.NSlayout.addWidget(self.progress, 7, self.mvpn, 1, 3)
		
		if (self.fab == 1): # if showing, remove single simulation widgets
			self.l81.deleteLater()
			self.fab = 0 # remember that single simulation widgets not showing
		if (self.atb == 1): # if showing, remove add to batch button
			self.b82.deleteLater()
			self.atb = 0 # remember that add to batch button not showing
		if (self.ssb == 1):# if showing, remove series of simulations button
			self.b83.deleteLater()
			self.ssb = 0 # remember this button removed
		self.l81b.setText('') # clear any old progress message from previous run
		
		# action to simulate series of simulations
		err_mess = self.on_click81b()
		
		# show checking inputs buttons
		self.b80.show()
		self.b80s.show()
	
	@pyqtSlot()
	def on_click82(self): # when button to add to batch pressed
		
		if (self.fab == 1): # if showing, remove single simulation widgets
			self.b81.deleteLater()
			self.l81.deleteLater() # or label
			self.fab = 0 # remember that single simulation widgets not showing
		if (self.atb == 1): # if showing remove add to batch button
			self.b82.deleteLater()
			self.l81.deleteLater() # or label
			self.atb = 0 # remember that add to batch button not showing
		
		self.l81b.setText('') # remove any old progress message from previous run
			
		if (self.btch_no == 1): # if first time adding to batch
			
			# add 'Start Series Of Simulations' button
			self.b81 = QPushButton('Start Series Of Simulations', self)
			self.b81.setToolTip('Start the series of simulations')
			self.b81.clicked.connect(self.on_click81)
			self.NSlayout.addWidget(self.b81, 5, self.mvpn, 1, 1)
		
		cs_file = self.l3.text()
		xml_file = self.l5.text()
		mv_file = self.l7.text()
		
		self.btch_str = (str(self.btch_str + str(self.btch_no)+'.\n' + cs_file + '\n' + xml_file + '\n' + mv_file + '\n'))
		self.l81b.setText(self.btch_str)
		
		self.btch_no += 1 # increase number in batch
		
		# append output path to list of output paths
		self.output_list.append(self.output_by_sim)
		
	@pyqtSlot() # button to quit software
	def on_click89(self):
		QWidget.close(self)
		sys.exit() # end program and release all memory
		
	# plot functions -----------------------------------------
	
	@pyqtSlot()
	def on_click202(self): # when model results folder requires selection

		# button to get path to folder containing relevant files
		options = QFileDialog.Options()
		fol_nme = QFileDialog.getExistingDirectory(self, "Select Folder Containing Required Input Files", "./PyCHAM/output/")
		
		# remember path of directory
		self.dir_path = fol_nme
		
		self.l201.clear() # clear old label
		self.l201.setText(fol_nme)

		# let user know data is being prepared
		self.l203a.setText(str('Progress through loading data' ))
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.bd_pl = 3
		
		# display progress bar
		self.progress = QProgressBar(self)
		self.PLlayout.addWidget(self.progress, 0, 4)

		QApplication.processEvents() # allow message panel to update
		
		# check whether required files present here
		try:
			
			# name of file where experiment constants saved
			fname = str(self.dir_path + '/model_and_component_constants')
			const_in = open(fname)
			const_in.close()
			
			# enable plotting buttons
			# disable plotting buttons
			self.b203.setEnabled(True)
			self.b206.setEnabled(True)
			self.b209.setEnabled(True)
			self.b212.setEnabled(True)
			self.b215.setEnabled(True)
			self.b218.setEnabled(True)
			self.b221.setEnabled(True)
			self.b223.setEnabled(True)
			
			# retrieve all outputs and include in self
			from retr_out import retr_out
			
			for prog in retr_out(self): # call on modules to solve problem
				
				if (isinstance(prog, str)): # check if it's a message
					mess = prog
					
					self.l203a.setText(mess)
					self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
					self.bd_pl = 3

				self.progress.setValue(int(prog)) # display progress	
				QApplication.processEvents() # allow message panel to update

			# remove any old progress message from previous run
			self.l203a.setText('Results loaded')
			# remove the progress bar after each simulation
			self.progress.deleteLater()
			QApplication.processEvents() # allow message panel to update
		except:
			# remove the progress bar after each simulation
			self.progress.deleteLater()
			self.l203a.setText('The required output files cannot be found at this path, please ensure the folder immediately above the output is selected')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			
			# disable plotting buttons
			self.b203.setEnabled(False)
			self.b206.setEnabled(False)
			self.b209.setEnabled(False)
			self.b212.setEnabled(False)
			self.b215.setEnabled(False)
			self.b218.setEnabled(False)
			self.b221.setEnabled(False)
			self.b223.setEnabled(False)
		QApplication.processEvents() # allow message panel to update
		return()
		
	@pyqtSlot() # button to plot standard results graphically
	def on_click203(self):	
		import plotter
		
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		dir_path = self.l201.text() # name of folder containing results
		gp_units = self.b203a.currentText() # gas-phase concentration units
		
		# convert units into number option
		if (gp_units[0] == 'p'):
			uc = 0
		if (gp_units[1] == 'g'):
			uc = 1
		if (gp_units[2] == 'm'):
			uc = 2

		# test whether upload has happened
		try:
			a_test = self.ro_obj.wf
		except:
			self.l203a.setText(str('Ensure that output is loaded using the Load Outputs button'))
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()

		plotter.plotter(0, dir_path, uc, self) # plot results
	
	@pyqtSlot() # button to plot gas-phase concentration temporal profile
	def on_click206(self):	
		
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		# test whether upload has happened
		try:
			a_test = self.ro_obj.wf
		except:
			self.l203a.setText(str('Ensure that output is loaded using the Load Outputs button'))
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()
		
		# get names of components to plot
		comp_names = [str(i) for i in self.e205.text(). split(',')]
		
		gp_units = self.b206b.currentText() # gas-phase concentration units
		
		# convert units into number option
		if (gp_units[0] == 'p' and gp_units[-1] == 'r'):
			caller = 1
		if (gp_units[1] == 'g' and gp_units[-1] == 'r'):
			caller = 0
		if (gp_units[2] == 'm' and gp_units[-1] == 'r'):
			caller = 3
		if (gp_units[0] == 'p' and gp_units[-1] == '.'):
			caller = 4
		if (gp_units[1] == 'g' and gp_units[-1] == '.'):
			caller = 5
		if (gp_units[2] == 'm' and gp_units[-1] == '.'):
			caller = 6

		import plotter_gp
		dir_path = self.l201.text() # name of folder with model results

		if (dir_path[-4::] != '.nc'):
			plotter_gp.plotter(caller, dir_path, comp_names, self) # plot results
		if (dir_path[-3::] == '.nc'):
			plotter_gp.plotter_noncsv(caller, dir_path, comp_names, self) # plot results
	
	@pyqtSlot() # button to plot ozone isopleth
	def on_click206c(self):

		# test whether upload has happened
		try:
			a_test = self.ro_obj.wf
		except:
			self.l203a.setText(str('Ensure that output is loaded using the Load Outputs button'))
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()

		import plotter_gp
		self.dir_path = self.l201.text() # name of folder with model results
	
		plotter_gp.O3_iso(self)

	# button to plot total particle-phase abundance for 
	# components' temporal profile
	@pyqtSlot()
	def on_click209(self):	
		
		# test whether upload has happened
		try:
			a_test = self.ro_obj.wf
		except:
			self.l203a.setText(str('Ensure that output is loaded using the Load Outputs button'))
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()
		
		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		# get names of components to plot
		comp_names = [str(i) for i in self.e205.text(). split(',')]
		
		# get whether to sum components or not
		sum_ornot = self.b205.currentText()
		if (sum_ornot[0] == 'I'):
			self.sum_ornot_flag = 0
		if (sum_ornot[0] == 'S'):
			self.sum_ornot_flag = 1

		import plotter_pp
		dir_path = self.l201.text() # name of folder with results
		plotter_pp.plotter(0, dir_path, comp_names, self) # plot results

	# button to plot temporal profile of total particle-phase concentration
	# excluding seed and water
	@pyqtSlot()
	def on_click209a(self):	
		
		# to ensure no confusion with plotting single components
		self.sum_ornot_flag = 0

		# test whether upload has happened
		try:
			a_test = self.ro_obj.wf
		except:
			self.l203a.setText(str('Ensure that output is loaded using the Select New Folder button'))
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()

		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		# get names of components to plot
		comp_names = []
		
		import plotter_pp
		# name of folder with results
		dir_path = self.l201.text()
		# plot results 
		plotter_pp.plotter(3, dir_path, comp_names, self) 

		return()
		
	# button to plot air quality index
	@pyqtSlot()
	def on_click209b(self):	
		
		# test whether upload has happened
		try:
			a_test = self.ro_obj.wf
		except:
			self.l203a.setText(str('Ensure that output is loaded using the Load Outputs button'))
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()
		
		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		import plotter
		plotter.aqi_calc(self) # plot results

	# button to plot total volatile organic compounds (including methane) against time
	@pyqtSlot()
	def on_click209c(self):

		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')

		import plotter
		plotter.tvoc_calc(self) # plot results

	# button to plot temporal profile of particle-phase contribution 
	# by top contributors to particle-phase concentration
	@pyqtSlot()
	def on_click300(self):	
		
		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		# test whether upload has happened
		try:
			a_test = self.ro_obj.wf
		except:
			self.l203a.setText(str('Ensure that output is loaded using the Load Outputs button'))
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()
		
		# get names of components (filler)
		comp_names = []
	
		# get number of components to plot
		try:
			self.e300r = int(self.e300.text())
		except: # default
			self.e300r = 10

		comp2cons = self.b300a.currentText() # components to consider
		# convert units into number option
		if (comp2cons[0] == 'A'):
			caller = 4
		if (comp2cons[0] == 'E'):
			caller = 7

		import plotter_pp
		dir_path = self.l201.text() # name of folder with results
		plotter_pp.plotter(caller, dir_path, comp_names, self) # plot results

	# graph the surface area concentration (m2/m3) of particles
	@pyqtSlot()
	def on_click301(self):	
		
		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')

		# test whether upload has happened
		try:
			a_test = self.ro_obj.wf
		except:
			self.l203a.setText(str('Ensure that output is loaded using the Load Outputs button'))
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()

		sa_type = self.b301a.currentText() # surface area to consider
		
		# convert units into number option
		if (sa_type[0] == 'A'):
			caller = 5
		if (sa_type[0] == 'S'):
			caller = 6
		
		# names of components (filler)
		comp_names = []

		import plotter_pp
		dir_path = self.l201.text() # name of folder with results
		plotter_pp.plotter(caller, dir_path, comp_names, self) # plot results

	# graph the particle-phase mass contribution from different generations of
	# oxidised organic molecules
	@pyqtSlot()
	def on_click302(self):	
		
		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')

		# test whether upload has happened
		try:
			a_test = self.ro_obj.wf
		except:
			self.l203a.setText(str('Ensure that output is loaded using the Load Outputs button'))
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()

		caller = 8
		
		# names of components (filler)
		comp_names = []

		import plotter_pp
		dir_path = self.l201.text() # name of folder with results
		plotter_pp.plotter(caller, dir_path, comp_names, self) # plot results
	
	# graph the contribution by user-supplied component by particle
	# mass against time 
	@pyqtSlot()
	def on_click304(self):	
		
		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')

		# test whether upload has happened
		try:
			a_test = self.ro_obj.wf
		except:
			self.l203a.setText(str(('''Ensure that ''') +
			str('''output is loaded using the Select ''') +
			str('''New Folder button''')))
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., 
				'2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., 
				'2px solid red', 0., 0.)
				self.bd_pl = 1
			return()

		import plotter_pp
		# name of folder with results
		self.dir_path = self.l201.text() 
		# names of components to plot
		self.comp_names_part_mass_vs_time = [str(i.strip()) for
			 i in self.e304.text().split(',')] 	
		# plot results
		plotter_pp.comp_part_mass_vs_time(self)	

	
	# graph the contribution by user-supplied component by particle
	# risk against time 
	@pyqtSlot()
	def on_click305(self):	
		
		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')

		# test whether upload has happened
		try:
			a_test = self.ro_obj.wf
		except:
			self.l203a.setText(str(('''Ensure that ''') +
			str('''output is loaded using the Select ''') +
			str('''New Folder button''')))
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., 
				'2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., 
				'2px solid red', 0., 0.)
				self.bd_pl = 1
			return()

		import plotter_pp
		# name of folder with results
		self.dir_path = self.l201.text() 
		# names of components to plot
		self.comp_names_part_mass_vs_time = [str(i.strip()) for
			 i in self.e304.text().split(',')] 	
		# plot results
		plotter_pp.comp_part_risk_vs_time(self)	

	# graph the contribution by user-supplied component by particle
	# carbon oxidation state against time 
	@pyqtSlot()
	def on_click306(self):	
		
		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')

		# test whether upload has happened
		try:
			a_test = self.ro_obj.wf
		except:
			self.l203a.setText(str(('''Ensure that ''') +
			str('''output is loaded using the Select ''') +
			str('''New Folder button''')))
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., 
				'2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., 
				'2px solid red', 0., 0.)
				self.bd_pl = 1
			return()

		import plotter_pp
		# name of folder with results
		self.dir_path = self.l201.text() 
		# names of components to plot
		self.comp_names_part_mass_vs_time = [str(i.strip()) for
			 i in self.e304.text().split(',')] 	
	
		# plot results
		plotter_pp.comp_part_cos_vs_time(self)
		

	# graph the particle-phase mass (excluding water) contribution 
	# from different upper limits of particle radius against time
	@pyqtSlot()
	def on_click303(self):	
		
		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')

		# test whether upload has happened
		try:
			a_test = self.ro_obj.wf
		except:
			self.l203a.setText(str('Ensure that output is loaded using the Load Outputs button'))
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()

		import plotter_pp
		self.dir_path = self.l201.text() # name of folder with results
		plotter_pp.part_mass_vs_time_sizeseg(self) # plot results
	
	# graph the particle-phase number (excluding water) contribution from 
	# different upper limits of particle radius against time
	@pyqtSlot()
	def on_click303q(self):	
		
		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')

		# test whether upload has happened
		try:
			a_test = self.ro_obj.wf
		except:
			self.l203a.setText(str('Ensure that output is loaded using the Load Outputs button'))
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()

		import plotter_pp
		self.dir_path = self.l201.text() # name of folder with results
		plotter_pp.part_num_vs_time_sizeseg(self) # plot results
	
	# graph the particle-phase surface area (excluding water) contribution from 
	# different upper limits of particle radius against time
	@pyqtSlot()
	def on_click303r(self):
		
		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')

		# test whether upload has happened
		try:
			a_test = self.ro_obj.wf
		except:
			self.l203a.setText(str('Ensure that output is loaded using the Load Outputs button'))
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()

		import plotter_pp
		self.dir_path = self.l201.text() # name of folder with results
		plotter_pp.part_area_vs_time_sizeseg(self) # plot results

	@pyqtSlot() # button to plot temporal profile of total concentration of 
	# components that have gas-wall partitioned to wall
	def on_click212(self):	
		
		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		# test whether upload has happened
		try:
			a_test = self.ro_obj.wf
		except:
			self.l203a.setText(str('Ensure that output is loaded using the Load Outputs button'))
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()
		
		# get names of components to plot
		comp_names = [str(i) for i in self.e205.text(). split(',')]
		
		import plotter_wp
		self.dir_path = self.l201.text() # name of folder with results
		plotter_wp.plotter(0, comp_names, self) # plot results
	
	@pyqtSlot() # button to plot temporal profile of total concentration of 
	# components on wall due to particle deposition to wall
	def on_click215(self):	

		# clear dialogue message
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')

		# get names of components to plot
		comp_names = [str(i) for i in self.e205.text(). split(',')]
		
		import plotter_wp_part
		dir_path = self.l201.text() # name of folder with results
		plotter_wp_part.plotter(0, dir_path, comp_names, self) # plot results
	
	@pyqtSlot() # button to plot temporal profile of chamber conditions
	# (temperature, pressure, relative humidity)
	def on_click215_a(self):	

		# clear dialogue message
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		import plotter_cham_env
		dir_path = self.l201.text() # name of folder with results
		plotter_cham_env.plotter(0, dir_path, self) # plot results
		
	@pyqtSlot() # button to plot change tendencies
	def on_click218(self):

		# clear dialogue message
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
			
		# no priority message to give to plotter module
		self.pre_mess = 0

		# get names of components to plot
		comp_names = [str(i) for i in self.e217.text().split(',')]

		import plotter_ct
		dir_path = self.l201.text() # name of folder with results
		plotter_ct.plotter(0, dir_path, comp_names, self) # plot results
	
	@pyqtSlot() # button to plot change tendencies due to individual chemical reactions
	def on_click218aa(self):

		# clear dialogue message
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
			
		# get names of components to plot
		comp_names = [str(i) for i in self.e217.text().split(',')]
		
		# get top number of chemical reactions to plot
		try:
			top_num = [int(i) for i in self.e217a.text().split(',')]
			self.pre_mess = 0
		except:
			# let plotter_ct know there is a 
			# priority message to user
			self.pre_mess = 1
			top_num = [10]

			mess = str('Please note that a user-defined number of ' +
			'reactions to plot was not seen, therefore the default of ' + 
			str(top_num[0]) + ' has been set.')
			self.l203a.setText(mess)

		ct_units = self.b218aaa.currentText() # change tendency units
		# convert units into number option
		if (ct_units[0] == 'p'):
			uc = 0
		if (ct_units[1] == 'g'):
			uc = 1
		if (ct_units[2] == 'm'):
			uc = 2
		
		import plotter_ct	
		# name of folder with results
		dir_path = self.l201.text()
		plotter_ct.plotter_ind(0, dir_path, comp_names, 
			top_num, uc, self) # plot results
	
	@pyqtSlot() # button to show production
	def on_click218ab(self):

		# clear dialogue message
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
			
		# get names of relevant components
		comp_names = [str(i) for i in self.e217.text().split(',')]
		
		# get times to integrate between (hours)
		tp = [float(i) for i in self.e217aa.text().split(',')]

		ct_units = self.b218aaa.currentText() # change tendency units
		# convert units into number option
		if (ct_units[0] == 'p'):
			uc = 0
		if (ct_units[1] == 'g'):
			uc = 1
		if (ct_units[2] == 'm'):
			uc = 2
		
		import plotter_ct
		dir_path = self.l201.text() # name of folder with results
		plotter_ct.plotter_prod(0, dir_path, comp_names, tp, uc, self) # plot results
	

	@pyqtSlot() # button to plot component contributions
	def on_click218b(self):

		# clear dialogue message
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
			
		# get SMILE string names of atom/functional group
		atom_name = str(self.e218a.text())
		
		# get number of components to plot
		atom_num = int(self.e218b.text())
		
		import plotter_atom_frac
		dir_path = self.l201.text() # name of folder with results
		plotter_atom_frac.plotter(0, dir_path, atom_name, atom_num, self) # plot results
	
	@pyqtSlot() # button to plot reaction ratios
	def on_click220b(self):
		
		try: # in case nunmerator numbers supplied
			# get numerator reaction numbers
			self.num_reac_num = np.array(([float(i) for i in self.e220a.text().split('/')[0].split(',')]))
		except:
			self.num_reac_num = []

		try: # in case denominator numbers supplied
			# get denominator reaction numbers
			self.den_reac_num = np.array(([float(i) for i in self.e220a.text().split('/')[1].split(',')]))
		except:
			self.den_reac_num = []
		self.dir_path = self.l201.text() # name of folder with results

		# call on plotter for reaction rate ratios
		import plotter_ct
		plotter_ct.plotter_reac_ratios(self)

	@pyqtSlot() # button to plot carbon reservoirs
	def on_click220c(self):

		self.dir_path = self.l201.text() # name of folder with results

		# call on plotter for reaction rate ratios
		import plotter_ct
		plotter_ct.plotter_carb_res(self)

	@pyqtSlot() # button to display molar mass of a single component
	def on_click220d(self):

		self.dir_path = self.l201.text() # name of folder with results

		# get name of single component
		self.mm_comp_name = self.e218c.text().strip()

		# get property to display
		self.single_comp_prop = self.b220e.currentText()

		# call on code for property display
		import plotter_ct
		plotter_ct.plotter_individ_prop(self)
		
	
	# button to plot volatility basis set mass fractions without water
	@pyqtSlot()
	def on_click221(self):

		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
			
		# get the phase to consider
		self.phase4vol = self.b221a.currentText()

		if 'Seed' in self.phase4vol and 'Water' in self.phase4vol:
			now = 1
		
		import vol_contr_analys
		dir_path = self.l201.text() # name of folder with results
		vol_contr_analys.plotter_wiw(0, dir_path, self, now) # plot results
	
	# button to plot two-dimensional volatility basis set mass fractions (VBS and O:C)
	@pyqtSlot()
	def on_click223(self):
		
		self.l203a.setStyleSheet(0., '0px dashed magenta', 0., 0.)
		self.l203a.setText('')
		
		try: # get time through experiment (s) at which to plot
			t_thro = float(self.e222.text())
		except: # give error message
			self.l203a.setText('Error - time through experiment (seconds) must be a single number')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1

			return()
		
		import vol_contr_analys	
		dir_path = self.l201.text() # name of folder with results
		vol_contr_analys.plotter_2DVBS(0, dir_path, self, t_thro) # plot results

		return()	
	# button to plot concentrations ordered by volatility
	@pyqtSlot()
	def on_click223a(self):
		
		self.l203a.setStyleSheet(0., '0px dashed magenta', 0., 0.)
		self.l203a.setText('')
		
		try: # get time through experiment (s) at which to plot
			t_thro = float(self.e222.text())
		except: # give error message
			self.l203a.setText('Error - time through experiment (seconds) must be a single number')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1

			return()
			
		import vol_contr_analys
		dir_path = self.l201.text() # name of folder with results
		vol_contr_analys.plotter_gpc(0, dir_path, self, t_thro) # plot results
		return()
		
	def on_click2231(self): # plot pie chart of top contributors to phase
		
		self.l203a.setStyleSheet(0., '0px dashed magenta', 0., 0.)
		self.l203a.setText('')
		
		# get the phase to consider
		self.phase4vol = self.b221a.currentText()
		
		try: # get time through experiment (s) at which to plot
			self.t_thro = float(self.e222.text())
		except: # give error message
			self.l203a.setText('Error - time through experiment (seconds) must be a single number')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1

		try: # get number of components to plot
			self.num_pie_comp = int(self.e2231.text())
		except: # give error message
			self.l203a.setText('Error - number of components for pie chart must be a single integer')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1

			return()
			
		import vol_contr_analys
		self.dir_path = self.l201.text() # name of folder with results
		vol_contr_analys.plotter_pie_top_n(self) # plot results
		return()
	
	@pyqtSlot() # button to plot total particle number concentration replication of particle counter
	def on_click222(self):
	
		# reset error message
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		try:
			dryf = float(self.e220.toPlainText()) # equilibrium humidity (fraction (0-1))
			
		except: # give error message
			self.l203a.setText('Note - relative humidity (fraction (0-1)) on reaching CPC condensing unit should be a single number, defaulting to 0.65')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
				
			dryf = 0.65 # default

		# false background counts ----------------------------------------------------------------
		try:
			cdt = float(self.e220_a.toPlainText()) # concentration detection limit (particles/cm3)
			
		except: # give error message
			self.l203a.setText('Note - false background particle number concentration counts (# particles cm<sup>-3</sup>) of counter should be a single number, defaulting to 1.e-2')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
				
			cdt = 1.e-2 # default
		# ----------------------------------------------------------------------------------------------
		
		# maximum detectable particle concentration (# particles/cm3) ------------
		try:
			max_dt = float(self.e220_b.toPlainText()) # concentration detection limit (# particles/cm3)
			
		except: # give error message
			self.l203a.setText('Note - maximum detectable particle number concentration (# particles cm<sup>-3</sup>) of counter should be a single number, defaulting to -1, which implies no maximum')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			max_dt = -1 # default

		# ----------------------------------------------------------------------------------------------
		
		# diameter at 50 % counting efficiency and width factor for sigmoidal 
		# curve of counting efficiency against size ---------------------------------
		try:
			# get particle diameter at 50 % counting efficiency and width 
			# factor for counting efficiency dependence on particle size
			sdt = ((self.e220_c.toPlainText()).split(','))
			sdt = [float(i) for i in sdt] 
			
			
		except: # give error message
			self.l203a.setText('Note - particle diameter (nm) at 50 % detection efficiency and factor for detection efficiency dependence on particle size should be two numbers separated by a comma, e.g.: 5, 1, defaulting to 5, 0.5')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			sdt = [5., 0.5]
		
		# -----------------------------------------------------------------------------------------
		
		# maximum particle size -------------------------------------------------------
		try:
			#  maximum particle size of counter
			max_size = float((self.e220_d.toPlainText()))
			
		except: # give error message
			self.l203a.setText('Note - maximum detectable particle diameter (nm) should be a single number, e.g. 3e5, defaulting to -1, which implies no maximum')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			max_size = -1
		# -------------------------------------------------------------------------------------
		
		# counter uncertainty (%) -------------------------------------------------------
		try:
			#  counter uncertainty (%)
			uncert = float((self.e220_e.toPlainText()))
			
		except: # give error message
			self.l203a.setText('Note - counter uncertainty (%) around particle number concentration should be a single number, e.g. 10, defaulting to 10')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			uncert = 10.
		# -------------------------------------------------------------------------------------
		
		# response time -----------------------------------------------------------------
		try:
			# get relevant inputs in list form
			ins = ((self.e220_f.toPlainText()).split(','))
			
			# obtain the delay times (s)
			delays = np.array((float(ins[0]), float(ins[1]), float(ins[3])))
			
			# obtain the weighting functions with delay time
			wfuncs = [str(ins[2]), str(ins[4])]
			
		except: # give error message
			self.l203a.setText('Note - there must be five inputs for weighting dependency on instrument response, all separated by commas.  The first, second and fourth must be numbers that represent the shortest response time at which particles counted, the response time at which greatest weighting is given to particle counts and the longest response time at which particles counted, respectively, with all these in seconds.  The third and fifth inputs should be the right hand side of the equations for weighting as a function of response time, with the third input providing this function between the shortest response time and the response time at maximum weighting and the fifth input providing this function between the response time at maximum weighting and the longest response time.  These functions should use t to represent response time, np for any numpy functions and python symbols for math functions such as multiply, e.g. 0.1, 1.3, np.exp(t), 2.5, np.exp(np.flip(t-1.3)).  Defaulting to 1., 1., 1.*t, 1., 1.*t which represents a response time of 1s with no mixing of particles of different ages.')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			delays = np.array((1., 1., 1.))
			wfuncs = ['1.*t', '1.*t']
		# -----------------------------------------------------------------------------------
		
		# inlet loss function and time ------------------------------------------------------------
		try:
			# get relevant inputs in list form
			ins = ((self.e220_h.toPlainText()).split(';'))
			
			# loss rate (fraction/s) as a function os particle size (um)
			loss_func_str = str(ins[0])
			losst = float(ins[1]) # time loss rate applies over (s)
			
		except: # give error message
			self.l203a.setText('Note - Loss rate (fraction s<sup>-1</sup>) as a function of particle size (um) and time of passage through inlet (s), should be a string (using Dp for diameter (um), np for numpy functions and python math symbols for math functions) followed by a  number, with the two separated by a semicolon, e.g.: np.append(10.**(-5.5-0.5*np.log10(Dp[Dp &lt;= 1.e-1])), 10.**(-4.2+0.8*np.log10(Dp[Dp &gt; 1.e-1]))); 5..  Defaulting to 0., 0. which implies no particle losses in inlet.')
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			loss_func_str = str('0.')
			losst = 0.
		# -----------------------------------------------------------------------------------
		
		# output frequency ------------------------------------------------------------
		try:
			# get relevant inputs in list form
			Hz = float((self.e220_g.toPlainText()))
			
		except: # give error message
			self.l203a.setText('Note - Output frequency of instrument (Hz) should be a single number, e.g. 0.1, defaulting to 1.')
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			Hz = 1.
		# -----------------------------------------------------------------------------------
		
		# averaging interval (s) ------------------------------------------------------
		try:
			# get relevant inputs in list form
			av_int = float((self.e220_i.toPlainText()))
			
		except: # give error message
			self.l203a.setText('Note - Averaging interval (s) should be a single number.  Defaulting to 1.')
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			av_int = 1.
		
		# averaging interval (s) ------------------------------------------------------
		try:
			# get relevant inputs in list form
			ins = ((self.e220_j.toPlainText()).split(','))
			Q = float(ins[0])
			tau = float(ins[1])
			coi_maxDp = float(ins[2])
			
		except: # give error message
			self.l203a.setText('Note - Coincidence inputs: volumetric flow rate through counting unit (cm<sup>3</sup> s<sup>-1</sup>), instrument dead time (s) and maximum actual particle concentration this can be applied to (# particles cm<sup>-3</sup>); should be three numbers separated by a comma (e.g. 5., 2.e-6, 3.e5).  Defaulting to -1, -1, -1, which indicates no convolution needed for coincidence (e.g. because coincidence already corrected for)')
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			Q = -1
			tau = -1
			coi_maxDp = -1
		# ---------------------------------------------------------------------------------
			
		# plot total particle number concentration (# particles/cm3)
		import plotter_counters
		importlib.reload(plotter_counters) # ensure latest version uploaded
		dir_path = self.l201.text() # name of folder with results
		plotter_counters.cpc_plotter(0, dir_path, self, dryf, cdt, max_dt, sdt, max_size, uncert, 
			delays, wfuncs, Hz, loss_func_str, losst, av_int, Q, tau, coi_maxDp) # plot results

		return()
	
	@pyqtSlot() # button to plot number size distribution replication of SMPS
	def on_click230_m(self):
	
		# reset error message
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		try:
			dryf = float(self.e230_b.toPlainText()) # equilibrium humidity (fraction (0-1))
			
		except: # give error message
			self.l203a.setText('Note - relative humidity (fraction (0-1)) on reaching CPC condensing unit should be a single number, defaulting to 0.65')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
				
			dryf = 0.65 # default

		# false background counts (minimum detection limit) ---------------------------
		try:
			cdt = float(self.e230_c.toPlainText()) # concentration detection limit (particles/cm3)
			
		except: # give error message
			self.l203a.setText('Note - false background particle number concentration counts (# particles cm<sup>-3</sup>) of counter should be a single number, defaulting to 1.e-2')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
				
			cdt = 1.e-2 # default
		
		# maximum detectable particle concentration (# particles/cm3) ------------
		try:
			max_dt = float(self.e230_d.toPlainText()) # concentration detection limit (# particles/cm3)
			
		except: # give error message
			self.l203a.setText('Note - maximum detectable particle number concentration (# particles cm<sup>-3</sup>) of counter should be a single number, defaulting to -1, which implies no maximum')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			max_dt = -1 # default

		# conincidence inputs (s) ------------------------------------------------------
		try:
			# get relevant inputs in list form
			ins = ((self.e220_e.toPlainText()).split(','))
			Q = float(ins[0])
			tau = float(ins[1])
			coi_maxD = float(ins[2])
			
		except: # give error message
			self.l203a.setText('Note - Coincidence inputs: volumetric flow rate through counting unit (cm<sup>3</sup> s<sup>-1</sup>), instrument dead time (s) and maximum actual particle concentration this can be applied to (# particles cm<sup>-3</sup>); should be three numbers separated by a comma (e.g. 5., 2.e-6, 3.e5).  Defaulting to -1, -1, -1, which indicates no convolution needed for coincidence (e.g. because coincidence already corrected for)')
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			Q = -1
			tau = -1
			coi_maxD = -1

		# diameter at 50 % counting efficiency and width factor for sigmoidal 
		# curve of counting efficiency against size ---------------------------------
		try:
			# get particle diameter at 50 % counting efficiency and width 
			# factor for counting efficiency dependence on particle size
			sdt = ((self.e230_f.toPlainText()).split(','))
			sdt = [float(i) for i in sdt] 
			
			
		except: # give error message
			self.l203a.setText('Note - particle diameter (nm) at 50 % detection efficiency and factor for detection efficiency dependence on particle size should be two numbers separated by a comma, e.g.: 5, 1, defaulting to 5, 0.5')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			sdt = [5., 0.5]
		
		# -----------------------------------------------------------------------------------------
		
		# maximum and minimum particle size -------------------------------------------------------
		try:
			#  maximum and minimum particle sizes
			max_size = ((self.e230_g.toPlainText()).split(','))
			max_size = [float(i) for i in max_size] 
			
		except: # give error message
			self.l203a.setText('Note - maximum and minimum detectable particle diameter (nm) should be two numbers separated by a comman, e.g. 1.0, 3e5, defaulting to -1, -1, which implies no minimum or maximum')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			max_size = [-1, -1]
		# -------------------------------------------------------------------------------------
		
		# counter uncertainty (%) -------------------------------------------------------
		try:
			#  counter uncertainty (%)
			uncert = float((self.e230_h.toPlainText()))
			
		except: # give error message
			self.l203a.setText('Note - counter uncertainty (%) around particle number concentration should be a single number, e.g. 10, defaulting to 10')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			uncert = 10.
		# -------------------------------------------------------------------------------------
		
		# response time -----------------------------------------------------------------
		try:
			# get relevant inputs in list form
			ins = ((self.e230_i.toPlainText()).split(','))
			
			# obtain the delay times (s)
			delays = np.array((float(ins[0]), float(ins[1]), float(ins[3])))
			
			# obtain the weighting functions with delay time
			wfuncs = [str(ins[2]), str(ins[4])]
			
		except: # give error message
			self.l203a.setText('Note - there must be five inputs for weighting dependency on instrument response, all separated by commas.  The first, second and fourth must be numbers that represent the shortest response time at which particles counted, the response time at which greatest weighting is given to particle counts and the longest response time at which particles counted, respectively, with all these in seconds.  The third and fifth inputs should be the right hand side of the equations for weighting as a function of response time, with the third input providing this function between the shortest response time and the response time at maximum weighting and the fifth input providing this function between the response time at maximum weighting and the longest response time.  These functions should use t to represent response time, np for any numpy functions and python symbols for math functions such as multiply, e.g. 0.1, 1.3, np.exp(t), 2.5, np.exp(np.flip(t-1.3)).  Defaulting to 1., 1., 1.*t, 1., 1.*t which represents a response time of 1s with no mixing of particles of different ages.')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			delays = np.array((1., 1., 1.))
			wfuncs = ['1.*t', '1.*t']
		# -----------------------------------------------------------------------------------
		
		# output frequency ------------------------------------------------------------
		try:
			# get relevant inputs in list form
			Hz = float((self.e230_j.toPlainText()))
			
		except: # give error message
			self.l203a.setText('Note - Output frequency of instrument (Hz) should be a single number, e.g. 0.1, defaulting to 1.')
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			Hz = 1.
		
		# inlet loss function and time ------------------------------------------------------------
		try:
			# get relevant inputs in list form
			ins = ((self.e230_k.toPlainText()).split(';'))
			
			# loss rate (fraction/s) as a function os particle size (um)
			loss_func_str = str(ins[0])
			losst = float(ins[1]) # time loss rate applies over (s)
			
		except: # give error message
			self.l203a.setText('Note - Loss rate (fraction s<sup>-1</sup>) as a function of particle size (um) and time of passage through inlet (s), should be a string (using Dp for diameter (um), np for numpy functions and python math symbols for math functions) followed by a  number, with the two separated by a semi-colon, e.g.: np.append(10.**(-5.5-0.5*np.log10(Dp[Dp &lt;= 1.e-1])), 10.**(-4.2+0.8*np.log10(Dp[Dp &gt; 1.e-1]))); 5..  Defaulting to 0., 0. which implies no particle losses in inlet.')
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			loss_func_str = str('0.')
			losst = 0.
		# -----------------------------------------------------------------------------------
		
		# averaging interval (s) ------------------------------------------------------
		try:
			# get relevant inputs in list form
			av_int = float((self.e230_l.toPlainText()))
			
		except: # give error message
			self.l203a.setText('Note - Averaging interval (s) should be a single number.  Defaulting to 1.')
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			av_int = 1.
		
		# number of size bins ---------------------------------------------------------------------------------
		
		try:
			#  counter's number of size bins (channels) per decade of particle size
			csbn = int((self.e230_m.toPlainText()))
			
		except: # give error message
			self.l203a.setText('Note - number of channels per decade of particle size should be a single number, defaulting to 128')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px magenta', 0., 0.)
				self.bd_pl = 1
				
			csbn = 128 # default value
		
		
		
		import plotter_counters
		importlib.reload(plotter_counters) # ensure latest version uploaded
		dir_path = self.l201.text() # name of folder with results
		
		plotter_counters.smps_plotter(0, dir_path, self, dryf, cdt, max_dt, sdt, max_size, uncert, delays, wfuncs, Hz, loss_func_str, losst, av_int, Q, tau, coi_maxD, csbn) # plot SMPS results
		
		return()

	# button to plot detection efficiency as a function of particle 
	# diameter for scanning mobility particle spectrometers
	@pyqtSlot()
	def on_click230_i(self):
		
		# reset message box
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		try:	
			# get particle diameter at 50 % counting efficiency and width 
			# factor for counting efficiency dependence on particle size
			sdt = ((self.e230_f.toPlainText()).split(','))
			sdt = [float(i) for i in sdt] 
			
		except: # give error message
			self.l203a.setText('Note - particle diameter (nm) at 50 % detection efficiency and factor for detection efficiency dependence on particle size should be two numbers separated by a comma, e.g.: 5, 1, defaulting to 5, 0.5')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			
			sdt = [5., 0.5]
						
		import plotter_counters
		dir_path = self.l201.text() # name of folder with results
		plotter_counters.count_eff_plot(0, dir_path, self, sdt) # plot
		
		return()
	

	# button to plot detection efficiency as a function of particle 
	# diameter for total particle number concentration counters
	@pyqtSlot()
	def on_click234_a(self):
		
		# reset message box
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		try:	
			# get particle diameter at 50 % counting efficiency and width 
			# factor for counting efficiency dependence on particle size
			sdt = ((self.e220_c.toPlainText()).split(','))
			sdt = [float(i) for i in sdt] 
			
		except: # give error message
			self.l203a.setText('Note - particle diameter (nm) at 50 % detection efficiency and factor for detection efficiency dependence on particle size should be two numbers separated by a comma, e.g.: 5, 1, defaulting to 5, 0.5')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			
			sdt = [5., 0.5]
						
		import plotter_counters
		dir_path = self.l201.text() # name of folder with results
		plotter_counters.count_eff_plot(0, dir_path, self, sdt) # plot
		
		return()
	
	# button to plot particle weighting by age due to instrument response time for SMPS
	@pyqtSlot()
	def on_click230_j(self):

		import plotter_counters
		
		# obtain the relevant inputs
		# reset message box
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		try:
			# get relevant inputs in list form
			ins = ((self.e230_i.toPlainText()).split(','))
			
			# obtain the delay times (s)
			delays = np.array((float(ins[0]), float(ins[1]), float(ins[3])))
			
			# obtain the weighting functions with delay time
			wfuncs = [str(ins[2]), str(ins[4])]
			
		except: # give error message
			self.l203a.setText('Note - there must be five inputs for weighting dependency on instrument response, all separated by commas.  The first, second and fourth must be numbers that represent the shortest response time at which particles counted, the response time at which greatest weighting is given to particle counts and the longest response time at which particles counted, respectively, with all these in seconds.  The third and fifth inputs should be the right hand side of the equations for weighting as a function of response time, with the third input providing this function between the shortest response time and the response time at maximum weighting and the fifth input providing this function between the response time at maximum weighting and the longest response time.  These functions should use t to represent response time, np for any numpy functions and python symbols for math functions such as multiply.  E.g. 0.1, 1.3, np.exp(t), 2.5, np.exp(np.flip(t-1.3)).  Defaulting to 1., 1., 1.*t, 1., 1.*t, which represents a 1 s response time with no mixing of particles of different ages.')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			delays = np.array((1., 1., 1.))
			wfuncs = ['1*t', '1*t']
						
		[w, t] = plotter_counters.resp_time_func(0, delays, wfuncs)

		return()
	
	# button to plot particle weighting by age due to instrument response time for CPC
	@pyqtSlot()
	def on_click222_b(self):

		import plotter_counters
		
		# obtain the relevant inputs
		# reset message box
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		try:
			# get relevant inputs in list form
			ins = ((self.e220_f.toPlainText()).split(','))
			
			# obtain the delay times (s)
			delays = np.array((float(ins[0]), float(ins[1]), float(ins[3])))
			
			# obtain the weighting functions with delay time
			wfuncs = [str(ins[2]), str(ins[4])]
			
		except: # give error message
			self.l203a.setText('Note - there must be five inputs for weighting dependency on instrument response, all separated by commas.  The first, second and fourth must be numbers that represent the shortest response time at which particles counted, the response time at which greatest weighting is given to particle counts and the longest response time at which particles counted, respectively, with all these in seconds.  The third and fifth inputs should be the right hand side of the equations for weighting as a function of response time, with the third input providing this function between the shortest response time and the response time at maximum weighting and the fifth input providing this function between the response time at maximum weighting and the longest response time.  These functions should use t to represent response time, np for any numpy functions and python symbols for math functions such as multiply.  E.g. 0.1, 1.3, np.exp(t), 2.5, np.exp(np.flip(t-1.3)).  Defaulting to 1., 1., 1.*t, 1., 1.*t, which represents a 1 s response time with no mixing of particles of different ages.')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			delays = np.array((1., 1., 1.))
			wfuncs = ['1*t', '1*t']
						
		[w, t] = plotter_counters.resp_time_func(0, delays, wfuncs)

		return()

	
	# button to plot loss rate of particles during passage through SMPS
	@pyqtSlot()
	def on_click230_k(self):

		import inlet_loss
		
		# obtain the relevant inputs
		# reset message box
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		# ------------------------------------------------------------------------------------------------
		
		try:	
			# get particle diameter at 50 % counting efficiency and width 
			# factor for counting efficiency dependence on particle size
			sdt = ((self.e230_f.toPlainText()).split(','))
			sdt = [float(i) for i in sdt] 
			
		except: # give error message
			self.l203a.setText('Note - particle diameter (nm) at 50 % detection efficiency and factor for detection efficiency dependence on particle size should be two numbers separated by a comma, e.g.: 5, 1, defaulting to 5, 0.5')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			
			sdt = [5., 0.5]
		
		# maximum particle size -------------------------------------------------------
		try:
			#  maximum particle size of counter
			max_size = float((self.e230_f.toPlainText()))
			
		except: # give error message
			self.l203a.setText('Note - maximum detectable particle diameter (nm) should be a single number, e.g. 3e5, defaulting to -1, which implies no maximum')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			max_size = -1
		# -------------------------------------------------------------------------------------
		
		# inlet loss function and time ------------------------------------------------------------
		try:
			# get relevant inputs in list form
			ins = ((self.e230_k.toPlainText()).split(';'))
			
			# loss rate (fraction/s) as a function os particle size (um)
			loss_func_str = str(ins[0])

			losst = float(ins[1]) # time loss rate applies over (s)
			
		except: # give error message
			self.l203a.setText('Note - Loss rate (fraction s<sup>-1</sup>) as a function of particle size (um) and time of passage through inlet (s), should be a string (using Dp for diameter (um), np for numpy functions and python math symbols for math functions) followed by a  number, with the two separated by a semicolon, e.g.: np.append(10.**(-5.5-0.5*np.log10(Dp[Dp &lt;= 1.e-1])), 10.**(-4.2+0.8*np.log10(Dp[Dp &gt; 1.e-1]))); 5..  Defaulting to 0., 0. which implies no particle losses in inlet.')
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			loss_func_str = str('0.')
			losst = 0.
		# -----------------------------------------------------------------------------------
		if (max_size == -1): # if defaulted
			max_size = 1.e3
		# radius range (um)
		xn = np.logspace(np.log10((sdt[0]*1.e-3)/2.), np.log10(max_size/2.), int(1e3))
		xn = xn.reshape(1, -1)
		
		sd_lrate = inlet_loss.inlet_loss(3, [], xn, [], loss_func_str, losst, 0)
		
		if (sd_lrate[0:5] == 'Error'):
			self.l203a.setText(sd_lrate)
	
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1

		return()


	# button to plot loss rate of particles during inlet passage
	@pyqtSlot()
	def on_click222_h(self):
		
		import inlet_loss
		
		# obtain the relevant inputs
		# reset message box
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		try:	
			# get particle diameter at 50 % counting efficiency and width 
			# factor for counting efficiency dependence on particle size
			sdt = ((self.e220_c.toPlainText()).split(','))
			sdt = [float(i) for i in sdt] 
			
		except: # give error message
			self.l203a.setText('Note - particle diameter (nm) at 50 % detection efficiency and factor for detection efficiency dependence on particle size should be two numbers separated by a comma, e.g.: 5, 1, defaulting to 5, 0.5')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			
			sdt = [5., 0.5]
		
		# maximum particle size -------------------------------------------------------
		try:
			#  maximum particle size of counter
			max_size = float((self.e220_d.toPlainText()))
			
		except: # give error message
			self.l203a.setText('Note - maximum detectable particle diameter (nm) should be a single number, e.g. 3e5, defaulting to -1, which implies no maximum')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			max_size = -1
		# -------------------------------------------------------------------------------------
		
		# inlet loss function and time ------------------------------------------------------------
		try:
			# get relevant inputs in list form
			ins = ((self.e220_h.toPlainText()).split(';'))
			
			# loss rate (fraction/s) as a function os particle size (um)
			loss_func_str = str(ins[0])
			losst = float(ins[1]) # time loss rate applies over (s)
			
		except: # give error message
			self.l203a.setText('Note - Loss rate (fraction s<sup>-1</sup>) as a function of particle size (um) and time of passage through inlet (s), should be a string (using Dp for diameter (um), np for numpy functions and python math symbols for math functions) followed by a  number, with the two separated by a semi-colon, e.g.: np.append(10.**(-5.5-0.5*np.log10(Dp[Dp &lt;= 1.e-1])), 10.**(-4.2+0.8*np.log10(Dp[Dp &gt; 1.e-1]))); 5..  Defaulting to 0., 0. which implies no particle losses in inlet.')
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			
			loss_func_str = str('0.')
			losst = 0.
		# -----------------------------------------------------------------------------------
		if (max_size == -1): # if defaulted
			max_size = 1.e3
		# radius range (um)
		xn = np.logspace(np.log10((sdt[0]*1.e-3)/2.), np.log10(max_size/2.), int(1e3))
		xn = xn.reshape(1, -1)

		sd_lrate = inlet_loss.inlet_loss(3, [], xn, [], loss_func_str, losst, 0)

		if (sd_lrate[0:5] == 'Error'):
			self.l203a.setText(sd_lrate)
	
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1

		return()

	# button to plot probability distribution function 
	# demonstrating mass:charge resolution
	@pyqtSlot()
	def on_click290_aa(self):
		# reset error message
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		dir_path = self.l201.text() # name of folder with results
		
		# get resolution inputs
		try:
			res_in = [] # empty list to contain values
			for i in ((self.e280.toPlainText().strip(' ').split(','))):
				res_in.append(float(i))
		except:
			self.l203a.setText('Note - failed to interpret input values for starting point and width of probability distribution function that represents mass:charge resolution, therefore defaulting to 1.0, 0.3')
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			res_in = [1.0, 0.3] # default

		y_MW = np.arange(0., 1000., res_in[1]/3.)
		
		import plotter_CIMS
		
		[_, _, _, _] = plotter_CIMS.write_mzres(3, res_in, y_MW)
		
	@pyqtSlot() # button to plot sensitivity to molar mass
	def on_click290_a(self):
		# reset error message
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		dir_path = self.l201.text() # name of folder with results
		
		y_MM = np.arange(1000.)
		
		# get sensitivity (Hz/ppt) dependence on molar mass
		try:
			sensit = str((self.e283.toPlainText()))
		except:
			sensit = 'np.ones(len(y_MW))' # default
		# means that the edit label text has not been 
		# changed from the description
		if (sensit[0:3] == 'Sen' or sensit == ''):
			sensit = 'np.ones(len(y_MM))'
		
		import plotter_CIMS
		
		blank = plotter_CIMS.write_sens2mm(3, sensit, y_MM)


	@pyqtSlot() # button to save results in CIMS format at all times through experiment
	def on_click290_aaa(self):
		# reset error message
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		self.dir_path = self.l201.text() # name of folder with results
		
		# get resolution inputs
		try:
			self.resol_in = [] # empty list to contain values
			for i in ((self.e280.toPlainText().strip(' ').split(','))):
				self.resol_in.append(float(i))
		except:
			self.l203a.setText('Note - failed to interpret input values for starting point and width of probability distribution function that represents mass:charge resolution, therefore defaulting to 1.0, 0.3')
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			self.resol_in = [1.0, 0.3] # default
		
		# get ioniser type
		try:
			self.iont = ((self.e282.toPlainText()).split(','))
			if (iont[0][0:3] == 'Ion' or len(iont) != 2):
				1./'a' # force exception
			if (int(iont[1]) != 0 and int(iont[1]) != 1):
				1./'a' # force exception
			if (iont[0] != 'I' and iont[0] != 'N'):
				1./'a' # force exception

		except: # give error message
			self.l203a.setText('Note - ionising agent, whether molar mass of ionising agent to be included in molar mass of components, should be a letter (I for iodide or N for nitrate) followed by a comma, followed by a number (1 to add agent molar mass to component mass or 0 not to).  But this information could not be correctly detected, so defaulting to I, 0')
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1

			self.iont = ['I', 0] # default
		
		# get sensitivity (Hz/ppt) dependence on molar mass
		try:
			self.sensit = str((self.e283.toPlainText()))
		except:
			self.sensit = 'np.ones(len(y_MW))' # default
		
		if (self.sensit[0:3] == 'Sen' or sensit == ''): # means that the edit label text has not been changed from the description
			self.sensit = 'np.ones(len(y_MW))'
		
		import plotter_CIMS
		
		plotter_CIMS.write_CIMS_output(self)
	
	# button to plot mass spectrum replication of 
	# chemical ionisation mass spectrometer
	@pyqtSlot()
	def on_click290(self):
	
		# reset error message
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		self.dir_path = self.l201.text() # name of folder with results
		
		# get resolution inputs
		try:
			res_in = [] # empty list to contain values
			for i in ((self.e280.toPlainText().strip(' ').split(','))):
				res_in.append(float(i))
		except:
			self.l203a.setText(str('Note - failed to interpret input ' +
			'values for starting point and width of probability ' +
			'distribution function that represents mass:charge ' +
			'resolution, therefore defaulting to 1.0, 0.3'))
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
			res_in = [1.0, 0.3] # default
			
		# get time to plot
		try:
			tn = float((self.e281.toPlainText()))
		except:
			try:
				tn = str((self.e281.toPlainText()))
				if tn != 'all times':
					tn = 0. # default
			except:
				tn = 0. # default
		
		# get ioniser type
		try:
			iont = ((self.e282.toPlainText()).split(','))
			if (iont[0][0:3] == 'Ion' or len(iont) != 2):
				1./'a' # force exception
			if (int(iont[1]) != 0 and int(iont[1]) != 1):
				1./'a' # force exception
			if (iont[0] != 'I' and iont[0] != 'N'):
				1./'a' # force exception

		except: # give error message
			self.l203a.setText(str('Note - ionising agent, whether ' +
			'molar mass of ionising agent to be included in molar mass ' +
			'of components, should be a letter (I for iodide or N for ' +
			'nitrate) followed by a comma, followed by a number (1 to ' +
			'add agent molar mass to component mass or 0 not to). ' +
			'But this information could not be correctly detected, ' +
			'so defaulting to I, 0'))
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1

			iont = ['I', 0] # default
			
		# get sensitivity (Hz/ppt) dependence on molar mass
		try:
			sensit = str((self.e283.toPlainText()))
		except:
			sensit = '1.' # default
		# means that the edit label text has not been changed 
		# from the description
		if (sensit[0:3] == 'Sen' or sensit == ''):
			sensit = '1.' # default
		
		# import required plotting function
		import plotter_CIMS
		plotter_CIMS.plotter_CIMS(self, res_in, tn, iont, sensit)

	@pyqtSlot() # button to check supplied values
	def on_clickb80(self):
		
		# get values to check from drop down button
		input_check_text = self.b80s.currentText() # drop down selection

		if (input_check_text == 'Photolysis Rates'):
			import plotter_simulate_tab
			plotter_simulate_tab.plotter_taf(self)

		if (input_check_text == 'Particle Number Size Distributions'):
			import plotter_nsd # for plotting supplied number size distributions
		
			# path for pickle file
			input_by_sim = str(self.PyCHAM_path + '/PyCHAM/pickle.pkl')

			# get the most recent model variables
			with open(input_by_sim, 'rb') as pk:
				[y0, Press,
				siz_stru, num_sb,
				 lowsize, uppsize, 
				std, Compt, injectt, Ct,
				seed_mw, seed_diss, seed_dens,
				dens_comp, dens, vol_comp, volP, 
				act_comp, act_user, 
				accom_comp, accom_val, uman_up, int_tol, 				new_partr, coag_on, inflectDp, pwl_xpre,
				 pwl_xpro, 
				inflectk, chamSA, Rader, p_char, 
				e_field, ser_H2O, 
				wat_hist, drh_str, erh_str, 
				z_prt_coeff, chamV] = pickle.load(pk)
				pk.close()
		
			# call on plotting script
			plotter_nsd.plotter_nsd(lowsize, num_sb, 
			uppsize, std, 
			0, self)

		if (input_check_text == 'Gas-phase Diffusion Coefficients'):
			import plotter_simulate_tab
			plotter_simulate_tab.plotter_gpdc(self)

		if (input_check_text == 'Gas-phase Mean Thermal Speeds'):
			import plotter_simulate_tab
			plotter_simulate_tab.plotter_gpmts(self)

		if (input_check_text == 'Molar Masses'):
			import plotter_simulate_tab
			plotter_simulate_tab.plotter_mm(self)
			
		if (input_check_text == 'Vapour Pressures'):
			import plotter_simulate_tab
			plotter_simulate_tab.plotter_vp(self)

		if (input_check_text == 'Nucleation Function'):
			import plotter_simulate_tab
			plotter_simulate_tab.plotter_nucfunc(self)

	@pyqtSlot() # button to retrieve and report component consumption
	def on_click224(self):
	
		self.dir_path = self.l201.text() # name of folder with results

		# get component name
		try:
			self.comp_names_to_plot = [comp_name for comp_name in (self.e224.text().split(','))]
			
		except: # give error message
			self.l203a.setText('Error - could not read chemical scheme name of component to estimate consumption of from box above')
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()
		try:
			# get starting time (hours)
			self.tmin = float((self.e224a.text()))
		except:
			self.tmin = 0.
		try:
			# get finishing time (hours)
			self.tmax = float((self.e224b.text()))
		except:
			self.tmax = 1.

		import consumption # function to estimate consumption
		consumption.cons(self, 0)

	@pyqtSlot() # button to retrieve and report yield
	def on_click225(self):
	
		self.dir_path = self.l201.text() # name of folder with results

		# get component name
		try:
			self.comp_names_to_plot = [comp_name for comp_name in (self.e224.text().split(','))]

		except: # give error message
			self.l203a.setText('Error - could not read chemical scheme name of component to estimate yield of from box above')
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()

		# get time range to estimate yield over
		try:
			# get starting time (hours)
			self.tmin = float((self.e224a.text()))
		except:
			self.tmin = 0.
		try:
			# get finishing time (hours)
			self.tmax = float((self.e224b.text()))
		except:
			self.tmax = 1.

		import consumption # function to estimate consumption
		consumption.cons(self, 1)
		return()

	@pyqtSlot()
	def on_click362(self): # for viewing properties of radicals
	
		# get number of top contributors
		# in case no user-defined value provided then default
		if (self.e360.text() == 'Provide the top number of components to plot'):
			self.rad_ord_num = 3
		else:
			self.rad_ord_num = int(self.e360.text())

		# get the name of the radical pool to plot from the drop down selection
		input_check_text = self.b361a.currentText()
		
		# if RO2 absolute concentrations
		if ('RO2 (alkyl peroxy radical) Concentration' in input_check_text):
			self.rad_mark = 0

		# if RO abolsute concentrations
		if ('RO (alkoxy radical) Concentration' in input_check_text):
			self.rad_mark = 1

		# if RO2 flux
		if ('RO2 (alkyl peroxy radical) Flux (molecules/cm3/s)' in input_check_text):
			self.rad_mark = 2

		# if RO flux
		if ('RO (alkoxy radical) Flux' in input_check_text):
			self.rad_mark = 3
		
		# get minimum carbon number of component to consider
		if str((self.e361.text()))[0:3] == 'Pro':
			# default to 1
			self.Cnum_thresh = 1
		else:
			self.Cnum_thresh = float((self.e361.text()))
		
		self.dir_path = self.l201.text() # name of folder with results


		import plotter_gp # required module
		if self.rad_mark == 0 or self.rad_mark == 1:
			# call on gas-phase plotter for radical pools
			plotter_gp.plotter_rad_pool(self) # plot results

		if (self.rad_mark == 2 or self.rad_mark == 3):
			# call on gas-phase plotter for radical flux
			plotter_gp.plotter_rad_flux(self) # plot results

		return()

	@pyqtSlot()
	def on_click400(self): # when observation file requires selection

		# button to get path to folder containing relevant files
		options = QFileDialog.Options()
		fil_nme = QFileDialog.getOpenFileName(self, str('Select file ' +
		'containing required observations'), './PyCHAM/output/')[0]
		self.l401.clear() # clear old label
		self.l401.setText(fil_nme) # set new label
		
		self.l203a.setText('')
		if (self.bd_pl == 1):
			self.l203a.setStyleSheet(0., '0px dashed magenta', 0., 0.)
			self.bd_pl = 2
		else:
			self.l203a.setStyleSheet(0., '0px solid magenta', 0., 0.)
			self.bd_pl = 1

		return()
	
	@pyqtSlot()
	def on_click401(self): # when observations and model results to be plotted

		# get path to observations
		self.xls_path = self.l401.text()
		# get path to model results
		self.mod_path = dir_path = self.l201.text()
	
		# get what to plot from drop-down button
		om_choice = self.b402.currentText()
		
		# flag for what is being plotted
		self.oandm = 0

		# if gas-phase concentrations
		if ('Gas-phase Concentrations' in om_choice):
			self.oandm = 1
			# get names of gas-phase components to plot
			self.gp_names = [str(i) for i in self.e205.text(). split(',')]
			# gas-phase concentration units
			self.gp_units = self.b206b.currentText()

		
		# if total particle concenrations
		if ('Particle Concentrations' in om_choice): 
			if ('Standard Results Plot' in om_choice):
				self.oandm = 1.1
			if ('Cumulative particle mass concentration without water' 
			in om_choice):
				self.oandm = 1.2
	
		# if Van Krevelen
		if ('Van Krevelen (Time Profile, Averaged Over All Hydrocarbons)' in om_choice):
			self.oandm = 2
		if ('Van Krevelen (All Individual Hydrocarbons at The Time Through Experiment (s) Provided Below)' in om_choice):
			self.oandm = 3
		if ('Van Krevelen (Time Profile, Averaged Over Non-methane Hydrocarbons)' in om_choice):
			self.oandm = 4
		if ('Van Krevelen (Non-methane Individual Hydrocarbons at The Time Through Experiment (s) Provided Below)' in om_choice):
			self.oandm = 5
		if ('Van Krevelen (Time Profile, Averaged Over _ Extension Hydrocarbons)' in om_choice):
			self.oandm = 6
		if ('Van Krevelen (_ Extension Individual Hydrocarbons at The Time Through Experiment (s) Provided Below)' in om_choice):
			self.oandm = 7

		# if mass defect
		if ('Mass Defect of All Components in Chemical Scheme' in om_choice):
			self.oandm = 8
		if ('Mass Defect of All Hydrocarbons Scaled to Concentrations at Time Through Experiment (s) Provided Below' in om_choice):
			self.oandm = 9

		# if CIMS mass spectrum
		if ('CIMS mass spectrum (as defined in Plot/Convolution/CIMS)' 
		in om_choice):
			self.oandm = 10

		# check on inputs and provide message for user if anything missing
		if (self.xls_path == ''):
			self.l203a.setText('Error - no path to observations selected')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()

		if (self.mod_path == ''):
			self.l203a.setText('Error - no path to model results selected')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()

		if (self.oandm == 1): # if on gas-phase temporal profiles
			if (self.gp_names == ['']):
				self.l203a.setText('Note, no components specified to be plotted for model results')
				if (self.bd_pl == 1):
					self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
					self.bd_pl = 2
				else:
					self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
					self.bd_pl = 1
		import plotter_xls
		
		if (self.oandm == 1): # if the gas-phase temporal profiles to be plotted
			plotter_xls.plotter_gp_mod_n_obs(self)
		# if the total particle concentrations temporal profiles 
		# to be plotted against the standard plot
		if (self.oandm == 1.1):
			plotter_xls.plotter_pp_mod_n_obs(self)
		if (self.oandm == 1.2):
			plotter_xls.plotter_pp_mod_n_obs(self)
		if (self.oandm >= 2 and self.oandm <= 7): # if the Van Krevelen to be plotted
			plotter_xls.plotter_VK_mod_n_obs(self)
		if (self.oandm == 8 or self.oandm == 9): # if mass defect plot
			plotter_xls.plotter_mass_defect(self)
		if (self.oandm == 10): # if CIMS mass spectrum
			self.on_click290() # call CIMS plotting button
	
		return()
		
	@pyqtSlot()
	def on_click500(self): # when experiment preparation outputs needed
		
		import plotter_xls
		
		# call function for estimating and saving outputs that are typically 
		# useful for designing experiments
		plotter_xls.plotter_exp_prep(self)
		
		self.l203a.setText('')
		if (self.bd_pl == 1):
			self.l203a.setStyleSheet(0., '0px dashed magenta', 0., 0.)
			self.bd_pl = 2
		else:
			self.l203a.setStyleSheet(0., '0px solid magenta', 0., 0.)
			self.bd_pl = 1

		return()

	def autorun(self): # function to automatically run PyCHAM
		# note that default model variables are already set 
		# before autorun function called

		# choose variable values
		param_range = self.param_const['param_ranges']		

		# filler for original starter simulation name
		starter_name0 = 'fill'
		
		if ((self.param_const['sim_num'] == 'set') and 
		(self.param_const['sim_type'] == 'finisher')):
	
			# tells mod_var_read that model variables not contained in a file
			self.inname = 'Not found'
		
			# loop through starter simulations
			for starteri in range(len(param_range['starter_paths'])):

				# get file name of starter	
				starter_name = param_range['starter_paths'][starteri]

				# set initial concentrations (ppb)
				# note that param_range['ys'] set in automated_setup_and_call.py
				# set starting concentration of components now (# molecules/cm3)
				self.param_const['ynow'] = param_range['ys'][starteri]

				# set temperature (K)
				self.param_const['temperature'] = param_range['temps'][starteri]
				# set relative humidity (0-1)
				self.param_const['rh'] = param_range['rhs'][starteri]
				# set pressure (Pa)
				self.param_const['p_init'] = param_range['pressures'][starteri]
				# set transmission factor for light (0-1)
				self.param_const['trans_fac'] = param_range['js'][starteri]

				if sys.platform == 'win32':
					save_path_start = 'C:\\Users\\Psymo\\OneDrive - The University of Manchester\\PyCHAM\\outputs\\interact\\'
				if sys.platform == 'darwin':
					save_path_start = '/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/PyCHAM/outputs/interact/'
				
				if ('loJ' in starter_name):
					self.param_const['res_file_name'] = str(save_path_start + starter_name + '_loJ')
				if ('meJ' in starter_name):
					self.param_const['res_file_name'] = str(save_path_start + starter_name + '_meJ')
				if ('hiJ' in starter_name):
					self.param_const['res_file_name'] = str(save_path_start + starter_name + '_hiJ')
		
				if ('loT' in starter_name):
					self.param_const['res_file_name'] = str(self.param_const['res_file_name'] + '_loT')
					if ('loT' not in starter_name0):
						# if temperature changes, need to renew property estimation 
						self.param_const['pars_skip'] = 0
				if ('meT' in starter_name):
					self.param_const['res_file_name'] = str(self.param_const['res_file_name'] + '_meT')
					if ('meT' not in starter_name0):
						# if temperature changes, need to renew property estimation 
						self.param_const['pars_skip'] = 0
				if ('hiT' in starter_name):
					self.param_const['res_file_name'] = str(self.param_const['res_file_name'] + '_hiT')
					if ('hiT' not in starter_name0):
						# if temperature changes, need to renew property estimation 
						self.param_const['pars_skip'] = 0

				# remember original starter simulation name
				starter_name0 = starter_name

				# remember original starting string in results save name
				res_file_name0 = self.param_const['res_file_name']

				# get names of influxing components in a useful list form (from string)
				const_influxers = self.param_const['const_infl'].split(',')

				# set influx rate for methane and carbon monoxide (ppb/s)
				Cinfl_0 = np.zeros((len(const_influxers)))

				if 'loCH4' in starter_name:
					# get index for methane
					CH4indx = const_influxers.index('CH4')
					Cinfl_0[CH4indx] = param_range['Cinfl'][CH4indx][0]	
					# get index for carbon monoxide
					COindx = const_influxers.index('CO')
					Cinfl_0[COindx] = param_range['Cinfl'][COindx][0]

				if 'meCH4' in starter_name:
					# get index for methane
					CH4indx = const_influxers.index('CH4')
					Cinfl_0[CH4indx] = param_range['Cinfl'][CH4indx][1]	
					# get index for carbon monoxide
					COindx = const_influxers.index('CO')
					Cinfl_0[COindx] = param_range['Cinfl'][COindx][1]

				if 'hiCH4' in starter_name:
					# get index for methane
					CH4indx = const_influxers.index('CH4')
					Cinfl_0[CH4indx] = param_range['Cinfl'][CH4indx][2]
					# get index for carbon monoxide
					COindx = const_influxers.index('CO')
					Cinfl_0[COindx] = param_range['Cinfl'][COindx][2]
		
				# get alpha-pinene index
				APindx = const_influxers.index('APINENE')

				# get benzene, NO and NO2 indices
				BZindx = const_influxers.index('BENZENE')
				NOindx = const_influxers.index('NO')
				NO2indx = const_influxers.index('NO2')
	
				# now loop through alpha-pinene influxes in combination with 
				# influxes of benzene, nitrogen oxide and nitrogen dioxide
				for iap in range(len(param_range['Cinfl'][APindx])):
					Cinfl_0[APindx] = param_range['Cinfl'][APindx][iap]
					if (iap == 0):
						self.param_const['res_file_name'] = str(self.param_const['res_file_name'] + '_loAP')
					if (iap == 1):
						self.param_const['res_file_name'] = str(self.param_const['res_file_name'] + '_meAP')
					if (iap == 2):
						self.param_const['res_file_name'] = str(self.param_const['res_file_name'] + '_hiAP')
					
					# loop through benzene influxes
					for ibz in range(len(param_range['Cinfl'][BZindx])):
						Cinfl_0[BZindx] = param_range['Cinfl'][BZindx][ibz]
						Cinfl_0[NOindx] = param_range['Cinfl'][NOindx][ibz]
						Cinfl_0[NO2indx] = param_range['Cinfl'][NO2indx][ibz]	
						if (ibz == 0):
							self.param_const['res_file_name'] = str(self.param_const['res_file_name'] + '_loBZ')
						if (ibz == 1):
							self.param_const['res_file_name'] = str(self.param_const['res_file_name'] + '_meBZ')
						if (ibz == 2):
							self.param_const['res_file_name'] = str(self.param_const['res_file_name'] + '_hiBZ')	
						# finally, either add (no interaction) or exclude (interaction) 
						# nint label from results save name
						if ('nint' in self.param_const['sch_name']):
							self.param_const['res_file_name'] = str(self.param_const['res_file_name'] + '_nint')
						
						# ensure correct format
						self.param_const['Cinfl'] = (str(Cinfl_0)).strip('[').strip(']')
						# loop through characters, white spaces with ; (components)
						cci = 1
						cc_all = -1
						for ccn in self.param_const['Cinfl']: # loop through characters
							cc_all += 1
							if (ccn == ' ' and cci % 2. > 0.):
								
								self.param_const['Cinfl'] = str(self.param_const['Cinfl'][0:cc_all]+';'+self.param_const['Cinfl'][cc_all+1::])
								cci += 1
								continue
							if (ccn == ' ' and cci % 2. < 1.):
								self.param_const['Cinfl'] = str(self.param_const['Cinfl'][0:cc_all]+';'+self.param_const['Cinfl'][cc_all+1::])
								cci += 1
								continue		
						print(self.param_const['res_file_name'])
						print(self.param_const['Cinfl'])
						print(self.param_const['const_infl'])
						print(self.param_const['temperature'])									
						print(self.param_const['rh'])
						print(self.param_const['p_init'])
						print(self.param_const['trans_fac'])

						if os.path.exists(self.param_const['res_file_name']):
							self.param_const['pars_skip'] = 0
							print('already exists!, so skipping')
							
						else:
							# establish parameters provided by user by calling mod_var_read
							import mod_var_read
							mod_var_read.mod_var_read(self)
							
							self.on_click2() # assign chemical scheme
							self.on_click3() # assign xml file
							self.on_click4() # provide model variables label
							
							self.on_click81sing() # run simulation
							# ensure that after the first call parsing and property estimation is skipped 
							self.param_const['pars_skip'] = 1

						# remove benzene tag, ready for new one
						if 'nint' in self.param_const['res_file_name']:
							self.param_const['res_file_name'] = self.param_const['res_file_name'][0:-10]
						else:
							self.param_const['res_file_name'] = self.param_const['res_file_name'][0:-5]
						
					# reset saving name to original
					self.param_const['res_file_name'] = res_file_name0
		
		if (self.param_const['sim_num'] == 'set') and (self.param_const['sim_type'] == 'starter'):

			# tells mod_var_read that model variables not in a model variables file 
			self.inname = 'Not found'

			# loop through light flux
			for Ji in range(len(param_range['trans_fac'])):
				
				# reset results save name
				res_name = '' 
				self.param_const['trans_fac'] = param_range['trans_fac'][Ji]

				if (Ji == 0):
					res_name = 'loJ' 
				if (Ji == 1):
					res_name = 'meJ'
				if (Ji == 2):
					res_name = 'hiJ' 

				# prepare for starting concentrations
				C0_0 = np.zeros((len(self.param_const['Comp0'].split(','))))
				
				# index for NO3 starting concentration in supplied 
				# values (constant between scenarios)
				#NO3i = (self.param_const['Comp0'].split(',').index('NO3'))
				# NO3 starting concentration (ppb)
				#C0_0[NO3i] = param_range['C0'][NO3i][0]

				# index for HNO3 starting concentration in supplied 
				# values (constant between scenarios)
				#HNO3i = (self.param_const['Comp0'].split(',').index('HNO3'))
				# HNO3 starting concentration (ppb)
				#C0_0[HNO3i] = param_range['C0'][HNO3i][0]

				# index for O3 starting concentration in supplied 
				# values (constant between scenarios)
				#O3i = (self.param_const['Comp0'].split(',').index('O3'))
				# O3 starting concentration (ppb)
				#C0_0[O3i] = param_range['C0'][O3i][0]

				# prepare for influxes
				Cinfl_0 = np.zeros((len(self.param_const['const_infl'].split(','))))
				# setup constant continuous influxes 
				# (not changing between scenarios)
				for const_infl_i in self.param_const['Cinfl_const_indx']:
					Cinfl_0[const_infl_i] = param_range['Cinfl'][const_infl_i][0]
				
				# loop through NOx emission
				for Ni in range(len(param_range['Cinfl'][0])):
				
					# reset results name to just the lighting part
					res_name = res_name[0:3]

					# index for NO starting concentration in supplied values
					NOi = self.param_const['Comp0'].split(',').index('NO')
					# NO starting concentration (ppb)
					C0_0[NOi] = param_range['C0'][NOi][Ni]

					# index for NO influx in supplied values
					NOi = self.param_const['const_infl'].split(',').index('NO')
					# NO emission rate (ppb/s)
					Cinfl_0[NOi] = param_range['Cinfl'][NOi][Ni]

					# index for NO2 starting concentration in supplied values
					NO2i = self.param_const['Comp0'].split(',').index('NO2')
					# NO2 starting concentration (ppb)
					C0_0[NO2i] = param_range['C0'][NO2i][Ni]

					# index for NO2 in supplied values
					NO2i = self.param_const['const_infl'].split(',').index('NO2')
					# NO2 emission rate (ppb/s)
					Cinfl_0[NO2i] = param_range['Cinfl'][NO2i][Ni]

					# index for O3 in supplied values
					O3i = self.param_const['const_infl'].split(',').index('O3')

					if (Ni == 0):
						res_name = str(res_name  + '_loNOx')
						# O3 transport rate (ppb/s)
						Cinfl_0[O3i] = param_range['Cinfl'][O3i][0]
					if (Ni == 1):
						res_name = str(res_name  + '_meNOx')
						# O3 transport rate (ppb/s)
						Cinfl_0[O3i] = param_range['Cinfl'][O3i][1]
					if (Ni == 2):
						res_name = str(res_name  + '_hiNOx')
						# O3 transport rate (ppb/s)
						Cinfl_0[O3i] = param_range['Cinfl'][O3i][2]

					# loop through CH4 and CO emission
					for Ci in range(len(param_range['Cinfl'][0])):
						
						# reset results name to the lighting and NOx part
						res_name = res_name[0:9]
						
						# index for CH4 starting concentration in supplied values
						CH4i = self.param_const['Comp0'].split(',').index('CH4')
						# CH4 starting concentration (ppb)
						C0_0[CH4i] = param_range['C0'][CH4i][Ci]

						# index for CH4 in supplied values
						CH4i = self.param_const['const_infl'].split(',').index('CH4')
						# CH4 emission rate (ppb/s)
						Cinfl_0[CH4i] = param_range['Cinfl'][CH4i][Ci]

						# index for CO starting concentration in supplied values
						COi = self.param_const['Comp0'].split(',').index('CO')
						# CO starting concentration (ppb)
						C0_0[COi] = param_range['C0'][COi][Ci]

						# index for CO in supplied values
						COi = self.param_const['const_infl'].split(',').index('CO')
						# CO emission rate (ppb/s)
						Cinfl_0[COi] = param_range['Cinfl'][COi][Ci]
										
						if (Ci == 0):
							res_name = str(res_name  + '_loCH4') 
						if (Ci == 1):
							res_name = str(res_name  + '_meCH4')
						if (Ci == 2):
							res_name = str(res_name  + '_hiCH4')
						
						for ti in range(len(param_range['temperature'])):

							# reset results name to the lighting, NOx and CH4 part
							res_name = res_name[0:15]

							self.param_const['temperature'] = param_range['temperature'][ti]
							if (ti == 0):
								res_name = str(res_name + '_loT') 
							if (ti == 1):
								res_name = str(res_name + '_meT')
							if (ti == 2):
								res_name = str(res_name + '_hiT')
							if (sys.platform == 'win32'):
								save_path_start = 'C:\\Users\\Psymo\\OneDrive - The University of Manchester\\PyCHAM\\outputs\\interact\\AP_BZ_MCM_PRAMAP_autoAPRAMBZ_scheme\\'
							if (sys.platform == 'darwin'):
								save_path_start = '/Users/user/Library/CloudStorage/OneDrive-TheUniversityofManchester/PyCHAM/outputs/interact/AP_BZ_MCM_PRAMAP_autoAPRAMBZ_scheme/'
							# check if this result already present
							print(str(save_path_start + res_name))
							
							# set results path
							self.param_const['res_file_name'] = str(save_path_start + res_name)

							if os.path.exists(self.param_const['res_file_name']):
								print('already exists!')
								continue
							
							# ensure correct format for mod_var_read
							self.param_const['C0'] = ''
						
							for i in C0_0:
								if i == '[' or i == ']':
									continue
								if i == ' ':
									continue
								self.param_const['C0'] = str(self.param_const['C0'] + str(i) + ', ')
							if self.param_const['C0'][-2::] == ', ':
								self.param_const['C0'] = self.param_const['C0'][0:-2] 	
							
							# ensure correct format
							self.param_const['Cinfl'] = (str(Cinfl_0)).strip('[').strip(']')
							# loop through characters, replacing odd white spaces with , 
							# (times) and even white spaces with ; (components)
							cci = 1
							cc_all = -1
							for ccn in self.param_const['Cinfl']: # loop through characters
								cc_all += 1
								if (ccn == ' ' and cci % 2. > 0.):
									
									self.param_const['Cinfl'] = str(self.param_const['Cinfl'][0:cc_all]+';'+self.param_const['Cinfl'][cc_all+1::])
									cci += 1
									continue
								if (ccn == ' ' and cci % 2. < 1.):
									self.param_const['Cinfl'] = str(self.param_const['Cinfl'][0:cc_all]+';'+self.param_const['Cinfl'][cc_all+1::])
									cci += 1
									continue
							
							# establish parameters provided by user by calling mod_var_read
							import mod_var_read
							mod_var_read.mod_var_read(self)
							
							self.on_click2() # assign chemical scheme
							self.on_click3() # assign xml file
							self.on_click4() # provide model variables label
							
							self.on_click81sing() # run simulation
						
							# ensure that after first simulation equation parsing is skipped
							self.param_const['pars_skip'] = 1
							
		if (type(self.param_const['sim_num']) == float or 
		type(self.param_const['sim_num']) == int and 
		self.param_const['sim_type'] != 'standard_call'):
	
			# prepare for randomness
			from numpy.random import default_rng
			rng = default_rng()

			# tells mod_var_read that model variables are 
			# not in a model variables file 
			self.inname = 'Not found'


			# loop through simulations
			for simi in range(self.param_const['sim_num']):
				
				# reset component concentrations
				self.param_const['Cinfl'] = ''

				# linear distribution in transmission factor of light (based on common sense, where 
				# 0=dark and 1=midday sunshine in Meditteranean in summer) 
				self.param_const['trans_fac'] = str('0_' + str((rng.integers(low=param_range['trans_fac'][0]*100., high=param_range['trans_fac'][1]*100., size=1))[0]/100.))
				
				# linear distribution in temperature (fig 11. of doi.org/10.1021/acsearthspacechem.1c00090) 
				self.param_const['temperature'] = (rng.integers(low=param_range['temperature'][0], high=param_range['temperature'][1], size=1))[0]
				# linear distribution in relative humidity (fig 11. of doi.org/10.1021/acsearthspacechem.1c00090) 
				self.param_const['rh'] = (rng.integers(low=param_range['rh'][0]*100., high=param_range['rh'][1]*100., size=1))[0]/100.
			
				# log-normal distribution for emission rate of gas phase: benzene, alpha-pinene, NOx, CO, SO2, CH4 (fig 10. of doi.org/10.1021/acsearthspacechem.1c00090)
				comp_cnt = 0 # count on components
				for compi in range(len(param_range['Cinfl'])):

					# create log-normal distribution for concentration range of this component
					minCi = param_range['Cinfl'][compi][0] # minimum concentration (ppb)
					maxCi = param_range['Cinfl'][compi][1] # maximum concentration (ppb)

					# linear distribution along log10 of range extremes
					lin_dis = np.linspace(np.log10(minCi), np.log10(maxCi), num=100)

					# randomly select concentration and raise to power 10
					conc_rand = 10**(lin_dis[(rng.integers(0, 99, size=1))[0]])

					if (comp_cnt == 0):
						self.param_const['Cinfl'] = str(self.param_const['Cinfl'] + str(conc_rand))
					else:
						self.param_const['Cinfl'] = str(self.param_const['Cinfl'] + ', ' + str(conc_rand))
			
					comp_cnt += 1 # count on components
				
				new_Cinfl = ''
				# now ensure continuous influx of components ends after 1 hour
				for ci in self.param_const['Cinfl'].split(','):

					new_Cinfl = str(new_Cinfl + str(ci) + ', 0. ;')
				# remove final ;			
				self.param_const['Cinfl'] = new_Cinfl[0:-1]	

				
				# log-normal distribution of seed particle concentration
				# create log-normal distribution for concentration range of this component
				minC = param_range['pconc'][0] # minimum concentration (# particles/cm3)
				maxC = param_range['pconc'][1] # maximum concentration (# particles/cm3)

				# linear distribution along log10 of range extremes
				lin_dis = np.linspace(np.log10(minC), np.log10(maxC), num=100)

				# randomly select total particle concentration influx rate and raise to power 10
				ptotal = 10**(lin_dis[(rng.integers(0, 99, size=1))[0]])

				# randomly select numbers between 0-1 to represent fraction in each size bin
				p0rnd = (rng.integers(0, 99, size=1))[0]
				p1rnd = (rng.integers(0, 99, size=1))[0]
				p2rnd = (rng.integers(0, 99, size=1))[0]
				ptotrnd = p0rnd+p1rnd+p2rnd
				self.param_const['pconc'] = str(str((p0rnd/ptotrnd)*ptotal) + ',' + str((p1rnd/ptotrnd)*ptotal) + ',' + str((p2rnd/ptotrnd)*ptotal) + '; 0,0,0')
				
				if (self.param_const['sim_type'] == 'finisher'):
					#import ast # for converting imported strings to list

					param_range['starter_paths'] = [] # prepare to list starter folders
					starter_names = [] # just the folder name

					# loop through starter simulations
					for starteri in range(len(param_range['starter_paths'])):

						self.param_const['res_file_name'] = str(starter_names[starteri] + '_' + str(simi))
						
						# withdraw concentrations (ppb in gas, # molecules/cm3 in particle and wall)
						#fname = str(param_range['starter_paths'][starteri] + '/concentrations_all_components_all_times_gas_particle_wall')
						#ystarter = (np.loadtxt(fname, delimiter=',', skiprows=1))[-1, :]

						#fname = str(param_range['starter_paths'][starteri] + '/model_and_component_constants')
						#const_in = open(fname)
						#for line in const_in.readlines():

							#if str(line.split(',')[0]) == 'factor_for_multiplying_ppb_to_get_molec/cm3_with_time':

								# find index of first [ and index of last ]
								#icnt = 0 # count on characters
								#for i in line:
								#	if i == '[':
								#		st_indx = icnt
								#		break
								#	icnt += 1 # count on characters
								#for cnt in range(10):
								#	if line[-cnt] == ']':
								#		fi_indx = -cnt+1
								#		break

								# conversion factor to change gas-phase concentrations from # molecules/cm3 
								# (air) into ppb
								#Cfactor = (ast.literal_eval(line[st_indx:fi_indx]))[-1]
						
							#for i in line.split(',')[1::]:
							#	if str(line.split(',')[0]) == 'number_of_components':
							#		num_comp = int(i)

						# convert ppb to # molecules/cm3
						#ystarter[0:num_comp] = ystarter[0:num_comp]*Cfactor
						
						# note that param_range['ys'] set in automated_setup_and_call.py
						# set starting concentration of components now (# molecules/cm3)
						self.param_const['ynow'] = param_range['ys'][starteri]

						# withdraw number-size distributions (# particles/cm3 (air))
						#fname = str(param_range['starter_paths'][starteri] +  '/particle_number_concentration_wet')
						#Nstarter = (np.loadtxt(fname, delimiter=',', skiprows=1))[-1, :]

						# note that param_range['Ns'] set in automated_setup_and_call.py
						# set starting concentration of particles now (# particles/cm3)
						self.param_const['Nnow'] = param_range['Ns'][starteri]

						# establish parameters provided by user by calling mod_var_read
						import mod_var_read
						mod_var_read.mod_var_read(self)
				
						self.on_click2() # assign chemical scheme
						self.on_click3() # assign xml file
						self.on_click4() # provide model variables label
				
						self.on_click81sing() # run simulation

			if (self.param_const['sim_type'] == 'starter'):

				# tells mod_var_read that model variables are 
				# not in a model variables file
				self.inname = 'Not found'

				# establish parameters provided by user by calling mod_var_read
				import mod_var_read
				mod_var_read.mod_var_read(self)
			
				param_const['sch_name']
				param_const['xml_name']

				self.on_click2() # assign chemical scheme
				self.on_click3() # assign xml file
				self.on_click4() # provide model variables label
			
				self.on_click81sing() # run simulation

		if (self.param_const['sim_type'] == 'standard_call'):

			self.inname = self.param_const['mod_var_name']
			# establish parameters provided by user by calling mod_var_read
			import mod_var_read
			mod_var_read.mod_var_read(self)
			
			# in case chemical scheme and xml 
			# file names contained in model variables file
			try:
				self.param_const['sch_name'] = self.sch_name
				self.param_const['xml_name'] = self.xml_name
			# in case chemical scheme and xml file names contained
			# in provided self.param_const variable
			except:
				self.sch_name = self.param_const['sch_name']
				self.xml_name = self.param_const['xml_name']
				self.inname = self.param_const['mod_var_name']
	
			self.on_click2() # assign chemical scheme
			self.on_click3() # assign xml file
			self.on_click4() # provide model variables label
			
			self.on_click81sing() # run simulation

		QWidget.close(self) # quit and close gui window

		return()

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
		
	def setStyleSheet(self, padn, brd, mrgn, spcn): # setting the style of the scroll area
		style_in = str('padding:  ' + str(padn) + ' px; border: ' + (brd) + ' px; margin: ' +str(mrgn) + ' px; spacing: ' + str(spcn))
		self.label.setStyleSheet(style_in)
	
	def setText(self, text): # the setText method
		self.label.setText(text) # setting text to the label
	
	def text(self): # interpreting text within the label method
		label_text = self.label.text()
		return(label_text)

	def clear(self): # the clear method
		self.label.clear() # clearing text

# function to allow opening of both file and directories
def getExistingFilesAndDirs(parent, caption, directory, 
                        filter='', initialFilter='', options=None):
	def updateText():
		# update the contents of the line edit widget with the selected files
		selected = []
		for index in view.selectionModel().selectedRows():
			selected.append('"{}"'.format(index.data()))
			lineEdit.setText(' '.join(selected))

	dialog = QFileDialog(parent, windowTitle=caption)
	dialog.setFileMode(dialog.ExistingFiles)
	if options:
		dialog.setOptions(options)
	dialog.setOption(dialog.DontUseNativeDialog, True)
	if directory:
		dialog.setDirectory(directory)
	if filter:
		dialog.setNameFilter(filter)
		if initialFilter:
			dialog.selectNameFilter(initialFilter)

	# by default, if a directory is opened in file listing mode, 
	# QFileDialog.accept() shows the contents of that directory, but we 
	# need to be able to "open" directories as we can do with files, so we 
	# just override accept() with the default QDialog implementation which 
	# will just return exec_()
	dialog.accept = lambda: QDialog.accept(dialog)

	# there are many item views in a non-native dialog, but the ones displaying 
	# the actual contents are created inside a QStackedWidget; they are a 
	# QTreeView and a QListView, and the tree is only used when the 
	# viewMode is set to QFileDialog.Details, which is not this case
	stackedWidget = dialog.findChild(QStackedWidget)
	view = stackedWidget.findChild(QListView)
	view.selectionModel().selectionChanged.connect(updateText)

	lineEdit = dialog.findChild(QLineEdit)
	# clear the line edit contents whenever the current directory changes
	dialog.directoryEntered.connect(lambda: lineEdit.setText(''))

	dialog.exec_()

	return(dialog.selectedFiles())
		
