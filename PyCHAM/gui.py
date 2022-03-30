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
'''The module that generates the Graphical User Interface (GUI) for PyCHAM, and connects that GUI with the core PyCHAM model'''
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
import re
import vol_contr_analys
import importlib
		
class PyCHAM(QWidget):

	def __init__(self):
		super().__init__()
		self.title = 'PyCHAM'
		self.left = 10
		self.top = 10
		self.width = 800
		self.height = 530
		
		self.initUI() # call on initialisation function to fill window
		
		return
		
	def initUI(self):
	
		# general window setup ---------------------------------------------------------------------
		self.setWindowTitle(self.title)
		self.setGeometry(self.left, self.top, self.width, self.height)
		
		# define grid layout
		grid = QGridLayout() 
		self.setLayout(grid)
		
		# title --------------------------------------------------------------------------------------------
		l0 = QLabel(self)
		l0.setOpenExternalLinks(True)
		l0.setText("Welcome to PyCHAM - please see <a href=\"http://www.github.com/simonom/PyCHAM\">README</a> for guidance")
		l0.setFont(QFont('Arial', 20))
		grid.addWidget(l0, 0, 0, 1, 1)
		
		# EUROCHAMP logo
		l0b= QLabel(self)
		pixmap = QPixmap('PyCHAM/logos_eurochamp2020-orange-noir.png')
		l0b.setPixmap(pixmap.scaled(100, 100, transformMode = Qt.SmoothTransformation))
		EUROlogo_hindx = 1 # horizontal position on grid
		grid.addWidget(l0b, 0, EUROlogo_hindx)
		# link to EUROCHAMP website
		bn1 = QPushButton('', self)
		bn1.setStyleSheet('background : transparent')
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
		bn1a.setStyleSheet('background : transparent')
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
		bn1b.setStyleSheet('background : transparent')
		bn1b.setToolTip('Visit the University of Manchester website')
		bn1b.clicked.connect(self.on_clickn1b)
		grid.addWidget(bn1b, 0, UoMlogo_hindx)
		
		# add tabs --------------------------------------------------------------------
		
		tabs = QTabWidget()
		tabs.addTab(self.NStab(), "Simulate")
		tabs.addTab(self.PLtab(), "Plot")
		grid.addWidget(tabs, 1, 0, 1, UoMlogo_hindx)
		
		# Quit pane ---------------------------------------------------------------------------------
		
		b89 = QPushButton('Quit', self)
		b89.setToolTip('Finish with PyCHAM and close this window')
		b89.clicked.connect(self.on_click89)
		grid.addWidget(b89, 1, UoMlogo_hindx)
		
		
		self.show()
		return
		
	def NStab(self): # New Simulation tab
	
		NSTab = QWidget()
		self.NSlayout = QGridLayout() 
		NSTab.setLayout(self.NSlayout)
		
		# folder and file starting column number
		ffscn = 0
		
		# folder selection for new simulation -------------------------------------------------------------------------
		b0 = QPushButton('Select Folder Containing Input Files', self)
		b0.setToolTip('Select the folder containing the required input files')
		b0.clicked.connect(self.on_click1)
		self.NSlayout.addWidget(b0, 0, ffscn, 1, 2)
		
		# default variables for all required input model variables -------------------------
		[sav_nam, chem_sch_mark, update_stp, tot_time, comp0, 
		y0, temp, tempt, RH, RHt, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, 
		pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, Compt, 
		injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedx,
		con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, 
		act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, 
		nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, 
		p_char, e_field, dil_fac, partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str, pcont, 
		Vwat_inc, seed_eq_wat, z_prt_coeff, chamV, 
		self] = def_mod_var.def_mod_var(0, self)
		
		# listing input files -----------------------------------------------------------
		l1 = QLabel(self)
		l1.setText("The following files have been found: ")
		self.NSlayout.addWidget(l1, 1, ffscn, 1, 2)
		
		l2 = QLabel(self)
		l2.setText("Chemical \nscheme: ")
		self.NSlayout.addWidget(l2, 2, ffscn, 1, 1)
			
		self.l3 = ScrollLabel(self)
		self.l3.setText(self.sch_name)
		self.NSlayout.addWidget(self.l3, 2, ffscn+1,  1, 1)
		
		b1 = QPushButton('Select Different Chemical scheme', self)
		b1.setToolTip('Select the file containing the required chemical scheme file')
		b1.clicked.connect(self.on_click2)
		self.NSlayout.addWidget(b1, 3, ffscn, 1, 2)
		
		l4 = QLabel(self)
		l4.setText("xml: ")
		self.NSlayout.addWidget(l4, 4, ffscn, 1, 1)
			
		self.l5 = ScrollLabel(self)
		self.l5.setText(self.xml_name)
		self.NSlayout.addWidget(self.l5, 4, ffscn+1, 1, 1)
			
		b2 = QPushButton('Select Different xml', self)
		b2.setToolTip('Select the file containing the required xml file')
		b2.clicked.connect(self.on_click3)
		self.NSlayout.addWidget(b2, 5, ffscn, 1, 2)
		
		l6 = QLabel(self)
		l6.setText("Model \nvariables: ")
		self.NSlayout.addWidget(l6, 6, ffscn, 1, 1)
			
		self.l7 = ScrollLabel(self)
		self.l7.setText(self.inname)
		self.NSlayout.addWidget(self.l7, 6, ffscn+1, 1, 1)
			
		b3 = QPushButton('Select Different Model variables', self)
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
		self.l9a.setText(sav_nam)
		self.varbox.addWidget(self.l9a, gen_row+1, 1)
		
		l10 = QLabel(self)
		l10.setText('Chemical scheme markers: ')
		self.varbox.addWidget(l10, gen_row+2, 0)
		self.l10a = QLabel(self)
		self.l10a.setText((str(chem_sch_mark)).replace('\'', '').replace(' ', '')[1:-1])
		self.varbox.addWidget(self.l10a, gen_row+2, 1)
		
		l11 = QLabel(self)
		l11.setText('Total experiment time (s): ')
		self.varbox.addWidget(l11, gen_row+3, 0)
		self.l11a = QLabel(self)
		self.l11a.setText((str(tot_time)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.l11a, gen_row+3, 1)
		
		l12 = QLabel(self)
		l12.setText('Update time interval (s): ')
		self.varbox.addWidget(l12, gen_row+4, 0)
		self.l12a = QLabel(self)
		self.l12a.setText((str(update_stp)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.l12a, gen_row+4, 1)
		
		l13 = QLabel(self)
		l13.setText('Recording time interval (s): ')
		self.varbox.addWidget(l13, gen_row+5, 0)
		self.l13a = QLabel(self)
		self.l13a.setText((str(save_step)).replace('\'', '').replace(' ', ''))
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
		self.l13_2a.setText((str(int_tol)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l13_2a, gen_row+7, 1)
		
		l13_2 = QLabel(self)
		l13_2.setText('Absolute and relative integration \ntolerances: ')
		self.varbox.addWidget(l13_2, gen_row+7, 0)
		self.l13_2a = QLabel(self)
		self.l13_2a.setText((str(int_tol)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
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
		self.l14a.setText((str(temp)).replace('\'', '').replace(' ', '')[1:-1])
		self.varbox.addWidget(self.l14a, env_row+1, 1)
		
		l15 = QLabel(self)
		l15.setText('Time(s) for chamber \ntemperature(s) (s): ')
		self.varbox.addWidget(l15, env_row+2, 0)
		self.l15a = QLabel(self)
		self.l15a.setText((str(tempt)).replace('\'', '').replace(' ', '')[1:-1])
		self.varbox.addWidget(self.l15a, env_row+2, 1)
		
		l16 = QLabel(self)
		l16.setText('Relative humidity (0-1): ')
		self.varbox.addWidget(l16, env_row+3, 0)
		self.l16a = QLabel(self)
		self.l16a.setText((str(RH)).replace('\'', '').replace(' ', ','))
		self.varbox.addWidget(self.l16a, env_row+3, 1)
		
		l16b = QLabel(self)
		l16b.setText('Relative humidity times (s): ')
		self.varbox.addWidget(l16b, env_row+4, 0)
		self.l16c = QLabel(self)
		self.l16c.setText((str(RHt)).replace('\'', '').replace(' ', ','))
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
		self.l17_1a.setText((str(dil_fac)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
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
		self.l20a.setText((str(pconc)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l20a, par_row+3, 1)
		
		l21 = QLabel(self)
		l21.setText('Particle concentration times (s): ')
		self.varbox.addWidget(l21, par_row+4, 0)
		self.l21a = QLabel(self)
		self.l21a.setText((str(pconct)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l21a, par_row+4, 1)
		
		l22 = QLabel(self)
		l22.setText('Molecular weight of seed particle \ncomponent (g/mol): ')
		self.varbox.addWidget(l22, par_row+5, 0)
		self.l22a = QLabel(self)
		self.l22a.setText((str(seed_mw)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l22a, par_row+5, 1)
		
		l23 = QLabel(self)
		l23.setText('Dissociation constant(s) of seed \ncomponent(s): ')
		self.varbox.addWidget(l23, par_row+6, 0)
		self.l23a = QLabel(self)
		self.l23a.setText((str(seed_diss)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l23a, par_row+6, 1)
		
		l24 = QLabel(self)
		l24.setText('Density of seed particles (g/cc): ')
		self.varbox.addWidget(l24, par_row+7, 0)
		self.l24a = QLabel(self)
		self.l24a.setText((str(seed_dens)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l24a, par_row+7, 1)
		
		l25 = QLabel(self)
		l25.setText('Name of seed component: ')
		self.varbox.addWidget(l25, par_row+8, 0)
		self.l25a = QLabel(self)
		self.l25a.setText((str(seed_name)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l25a, par_row+8, 1)
		
		l26 = QLabel(self)
		l26.setText('Mole fraction of non-water \ncomponents in \ndry seed particles: ')
		self.varbox.addWidget(l26, par_row+9, 0)
		self.l26a = QLabel(self)
		self.l26a.setText((str(seedx)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l26a, par_row+9, 1)
		
		l26b = QLabel(self)
		l26b.setText('Whether (1) or not (0) \nvolume of water included \nin seed particle \nnumber size distribution: ')
		self.varbox.addWidget(l26b, par_row+10, 0)
		self.l26b = QLabel(self)
		self.l26b.setText((str(Vwat_inc)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l26b, par_row+10, 1)
		
		l26c = QLabel(self)
		l26c.setText('Whether (1) or not (0) \nwater to be equilibrated \nwith seed particle \nprior to experiment start: ')
		self.varbox.addWidget(l26c, par_row+11, 0)
		self.l26c = QLabel(self)
		self.l26c.setText((str(Vwat_inc)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
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
		self.l29a.setText((str(space_mode)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
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
		self.l31a.setText((str(mean_rad)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l31a, par_row+16, 1)
		
		l32 = QLabel(self)
		l32.setText('Radius (cm) of newly nucleated \nparticles: ')
		self.varbox.addWidget(l32, par_row+17, 0)
		self.l32a = QLabel(self)
		self.l32a.setText((str(new_partr)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l32a, par_row+17, 1)
		
		l33 = QLabel(self)
		l33.setText('First nucleation \nparameterisation parameter: ')
		self.varbox.addWidget(l33, par_row+18, 0)
		self.l33a = QLabel(self)
		self.l33a.setText((str(nucv1)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l33a, par_row+18, 1)
		
		l34 = QLabel(self)
		l34.setText('Second nucleation \nparameterisation parameter: ')
		self.varbox.addWidget(l34, par_row+19, 0)
		self.l34a = QLabel(self)
		self.l34a.setText((str(nucv2)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l34a, par_row+19, 1)
		
		l35 = QLabel(self)
		l35.setText('Third nucleation \nparameterisation parameter: ')
		self.varbox.addWidget(l35, par_row+20, 0)
		self.l35a = QLabel(self)
		self.l35a.setText((str(nucv3)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l35a, par_row+20, 1)
		
		l36 = QLabel(self)
		l36.setText('Chemical scheme name of \nnucleating component: ')
		self.varbox.addWidget(l36, par_row+21, 0)
		self.l36a = QLabel(self)
		self.l36a.setText((str(nuc_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l36a, par_row+21, 1)
		
		l37 = QLabel(self)
		l37.setText('Whether (1) or not (0) to adapt \nintegration time interval \nto nucleation: ')
		self.varbox.addWidget(l37, par_row+22, 0)
		self.l37a = QLabel(self)
		self.l37a.setText((str(nuc_ad)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
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
		self.l38_5a.setText((str(pcont)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
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
		self.l40a.setText((str(comp0)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
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
		self.l42a.setText((str(con_infl_nam)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l42a, gas_row+3, 1)
		
		l43 = QLabel(self)
		l43.setText('Influx rate of components with \ncontinuous influx (ppb/s): ')
		self.varbox.addWidget(l43, gas_row+4, 0)
		self.l43a = QLabel(self)
		self.l43a.setText((str(con_infl_C)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l43a, gas_row+4, 1)
		
		l44 = QLabel(self)
		l44.setText('Times of component influx (s): ')
		self.varbox.addWidget(l44, gas_row+5, 0)
		self.l44a = QLabel(self)
		self.l44a.setText((str(con_infl_t)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
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
		self.l59a.setText((str(wall_on)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l59a, wall_row+1, 1)
		
		l60 = QLabel(self)
		l60.setText('Effective absorbing mass of \nwall (g/m3 (air)): ')
		self.varbox.addWidget(l60, wall_row+2, 0)
		self.l60a = QLabel(self)
		self.l60a.setText((str(Cw)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l60a, wall_row+2, 1)
		
		l61 = QLabel(self)
		l61.setText('Gas-wall mass transfer \ncoefficient (/s): ')
		self.varbox.addWidget(l61, wall_row+3, 0)
		self.l61a = QLabel(self)
		self.l61a.setText((str(kw)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
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
		l67.setText('Whether particle deposition to \nwall treated by Rader and \nMcMurry (1) or customised (0): ')
		self.varbox.addWidget(l67, wall_row+10, 0)
		self.l67a = QLabel(self)
		self.l67a.setText((str(Rader)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l67a, wall_row+10, 1)
		
		l68 = QLabel(self)
		l68.setText('Average number of charges per \nparticle (/particle): ')
		self.varbox.addWidget(l68, wall_row+11, 0)
		self.l68a = QLabel(self)
		self.l68a.setText((str(p_char)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l68a, wall_row+11, 1)
		
		l69 = QLabel(self)
		l69.setText('Average electric field inside chamber \n(g.m/A.s3): ')
		self.varbox.addWidget(l69, wall_row+12, 0)
		self.l69a = QLabel(self)
		self.l69a.setText((str(e_field)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l69a, wall_row+12, 1)
		
		# specific component properties ----------------------
		
		l70b = QLabel(self)
		l70b.setText('Specific Component Properties')
		l70b.setFont(QFont("Arial", 13, QFont.Bold))
		scp_row = wall_row+12
		self.varbox.addWidget(l70b, scp_row+0, 0)
		
		l70 = QLabel(self)
		l70.setText('Chemical scheme name of components \nto track process tendencies: ')
		self.varbox.addWidget(l70, scp_row+1, 0)
		self.l70a = QLabel(self)
		self.l70a.setText((str(dydt_trak)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l70a, scp_row+1, 1)
		
		l71 = QLabel(self)
		l71.setText('Chemical scheme name of components \nwith specified densities: ')
		self.varbox.addWidget(l71, scp_row+2, 0)
		self.l71a = QLabel(self)
		self.l71a.setText((str(dens_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l71a, scp_row+2, 1)
		
		l72 = QLabel(self)
		l72.setText('Specified densities of \ncomponents (g/cc): ')
		self.varbox.addWidget(l72, scp_row+3, 0)
		self.l72a = QLabel(self)
		self.l72a.setText((str(dens)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l72a, scp_row+3, 1)
		
		l73 = QLabel(self)
		l73.setText('Chemical scheme name of components \nwith specified vapour pressures: ')
		self.varbox.addWidget(l73, scp_row+4, 0)
		self.l73a = QLabel(self)
		self.l73a.setText((str(vol_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l73a, scp_row+4, 1)
		
		l74 = QLabel(self)
		l74.setText('Specified vapour pressures of \ncomponents (Pa): ')
		self.varbox.addWidget(l74, scp_row+5, 0)
		self.l74a = QLabel(self)
		self.l74a.setText((str(volP)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l74a, scp_row+5, 1)
		
		l75 = QLabel(self)
		l75.setText('Chemical scheme name of components \nwith specified activity coefficients: ')
		self.varbox.addWidget(l75, scp_row+6, 0)
		self.l75a = QLabel(self)
		self.l75a.setText((str(act_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l75a, scp_row+6, 1)
		
		l76 = QLabel(self)
		l76.setText('Specified activity coefficients of \ncomponents:  ')
		self.varbox.addWidget(l76, scp_row+7, 0)
		self.l76a = QLabel(self)
		self.l76a.setText((str(act_user)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l76a, scp_row+7, 1)

		l77 = QLabel(self)
		l77.setText('Chemical scheme name of components \nwith specified accommodation \ncoefficients: ')
		self.varbox.addWidget(l77, scp_row+8, 0)
		self.l77a = QLabel(self)
		self.l77a.setText((str(accom_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l77a, scp_row+8, 1)
		
		l78 = QLabel(self)
		l78.setText('Specified accommodation coefficients of \ncomponents:  ')
		self.varbox.addWidget(l78, scp_row+9, 0)
		self.l78a = QLabel(self)
		self.l78a.setText((str(accom_val)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l78a, scp_row+9, 1)
		
		
		# properties of model variables scroll area ----------------
		self.scrollwidget.setLayout(self.varbox)
		self.scroll.setWidget(self.scrollwidget)
		self.NSlayout.addWidget(self.scroll, 1, self.mvpn, 3, 3)
		
		# --------------------------------------------------------------------
		# label to let user know preparedness of simulation - displayed text updated in ui_check module below
		self.l80 = ScrollLabel(self)
		self.NSlayout.addWidget(self.l80, 4, self.mvpn, 1, 3)
		self.bd_st = 2 # border status

		# --------------------------------------------------------------------
		# drop down button to let users view a variety of variables 
		# determined by model variables
		self.b80s = QComboBox(self)
		self.b80s.addItem('Photolysis rates')
		self.b80s.addItem('Particle number size distributions')
		self.b80s.addItem('Gas-phase diffusion coefficients')
		self.b80s.addItem('Gas-phase mean thermal speeds')
		self.b80s.addItem('Molar masses')
		self.NSlayout.addWidget(self.b80s, 7, self.mvpn+0, 1, 2)
		
		# button to run checks on variables selected in drop down button
		self.b80 = QPushButton('Check Values', self)
		self.b80.setToolTip('See the values of the variables selected in the drop down button to the left')
		self.b80.clicked.connect(self.on_clickb80)
		self.NSlayout.addWidget(self.b80, 7, self.mvpn+2, 1, 1)

		# -------------------------------------------------------------------

		# label to let users know file combinations included in batch
		self.btch_str = 'File combinations included in batch: chemical scheme, xml, model variables\n'
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
		
		# running check on default model variables ---------------------------------------------------------------
		# let checking module know this is a first call
		self.chck_num = 1

		import ui_check # module for checking on model variables
		# check on inputs - note this loads the last saved pickle file and saves 
		# any change to this pickle file
		ui_check.ui_check(self)
		# finished check on model variables -------------------------------------------------------------------------
		
		# relative stretching (width-wise) of each column in Simuate tab
		self.NSlayout.setColumnStretch(0, 0.7)
		self.NSlayout.setColumnStretch(1, 1)
		self.NSlayout.setColumnStretch(2, 1)
		self.NSlayout.setColumnStretch(3, 1)
		self.NSlayout.setColumnStretch(4, 0.1)
		self.NSlayout.setColumnStretch(5, 1)

		return(NSTab)
	
	def PLtab(self): # Plotting tab
	
		PLTab = QWidget()
		self.PLlayout = QGridLayout() 
		PLTab.setLayout(self.PLlayout)
		
		# results folder and dialogue row --------------
		
		# label for name of results folder to plot:
		l200 = QLabel(self)
		l200.setText('Path to results folder to be used for plotting: ')
		l200.setWordWrap(True)
		#l200.setAlignment(Qt.AlignRight)
		self.PLlayout.addWidget(l200, 0, 0)
		
		self.l201 = ScrollLabel(self)
		cwd = os.getcwd() # current working directory
		path = str(cwd + '/PyCHAM/output/ex_chem_scheme/Example_Run_Output')
		self.l201.setText(path)
		self.PLlayout.addWidget(self.l201, 0, 1, 1, 1)
		
		b202 = QPushButton('Select new folder', self)
		b202.setToolTip('Select the folder containing the result files to plot')
		b202.clicked.connect(self.on_click202)
		self.PLlayout.addWidget(b202, 0, 2)
		
		# label for name of dialogue box:
		l202a = QLabel(self)
		l202a.setText('Messages from scripts: ')
		l202a.setWordWrap(True)
		self.PLlayout.addWidget(l202a, 0, 3)
		
		# status for border around plotting messages
		self.bd_pl = 3
		
		# plotting dialogue box:
		self.l203a = ScrollLabel(self)
		self.l203a.setText('No message currently')
		self.PLlayout.addWidget(self.l203a, 0, 4, 1, 1)
		
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
		PLtabs.addTab(self.RADtab(), "Radicals") # for radical chemicals
		PLtabs.addTab(self.INSTRtab(), "Convolution")
		PLtabs.addTab(self.OBStab(), "Observed")
		self.PLlayout.addWidget(PLtabs, 1, 0, 1, 5)
		
		# relative stretching (width-wise) of each column in Plot tab
		self.PLlayout.setColumnStretch(0, 1)
		self.PLlayout.setColumnStretch(1, 2)
		self.PLlayout.setColumnStretch(2, 0.8)
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

		# particle-phase concentrations temporal profiles -------------

		# button to plot temporal profile of total particle-phase concentration of supplied components
		self.b209 = QPushButton(str('Total particle-phase concentrations \n('+u'\u03BC'+'g/m'+u'\u00B3'+')'), self)
		self.b209.setToolTip('Plot particle-phase concentration temporal profile of these components')
		self.b209.clicked.connect(self.on_click209)
		self.PRIMlayout.addWidget(self.b209, 4, 0)

		# button to plot temporal profile of total particle-phase concentration excluding seed and water 
		self.b209a = QPushButton(str('Total particle-phase concentrations \n('+u'\u03BC'+'g/m'+u'\u00B3'+') excluding seed and water'), self)
		self.b209a.setToolTip('Plot total particle-phase concentration of all components except for seed and water')
		self.b209a.clicked.connect(self.on_click209a)
		#self.b209a.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.PRIMlayout.addWidget(self.b209a, 4, 1)
		
		# wall (from gas-wall partitioning) concentrations temporal profiles -------------
		
		# button to plot temporal profile of total particle-phase concentrations
		self.b212 = QPushButton('Wall concentrations (from gas-wall partitioning) ('+u'\u03BC'+'g/m'+u'\u00B3'+')', self)
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
		self.b215_a = QPushButton('Chamber Conditions (T, P, RH)', self)
		self.b215_a.setToolTip('Plot the temporal profile of chamber variables (temperature, pressure, relative humidity)')
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
		self.SEClayout.addWidget(self.separatorLine4, 3, 1, 7, 1)
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

		# volatility basis set ------------------
		
		# label for names of components to plot tracked change tendencies
		l219 = QLabel(self)
		l219.setText('Buttons below plot mass fractions of components grouped by vapour pressure at 298.15 K')
		l219.setWordWrap(True)
		self.VOLlayout.addWidget(l219, 4, 0, 2, 1)
		
		# button to plot temporal profile of volatility basis set mass fractions with water
		self.b220 = QPushButton('Plot Volatility Basis Set With Water', self)
		self.b220.setToolTip('Plot the temporal profile of volatility basis set mass fractions')
		self.b220.clicked.connect(self.on_click220)
		#self.b220.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.VOLlayout.addWidget(self.b220, 6, 0)
		
		# button to plot temporal profile of volatility basis set mass fractions without water
		self.b221 = QPushButton('Plot Volatility Basis Set Without Water', self)
		self.b221.setToolTip('Plot the temporal profile of volatility basis set mass fractions')
		self.b221.clicked.connect(self.on_click221)
		#self.b221.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.VOLlayout.addWidget(self.b221, 7, 0)
		
		# two-dimensional volatility basis set ------------------
		
		# input bar for time through experiment to plot the 2D VBS for
		self.e222 = QLineEdit(self)
		self.e222.setText('Provide the time (seconds) through experiment at which to plot the two-dimensional volatility basis set - the closest recorded time to this will be used')
		self.e222.setStyleSheet('qproperty-cursorPosition : 0')
		self.VOLlayout.addWidget(self.e222, 8, 0)
		
		# button to plot 2D VBS
		self.b223 = QPushButton('Plot 2D Volatility Basis Set', self)
		self.b223.setToolTip('Plot the two-dimensional volatility basis set (O:C ratio and vapour pressures) at the time through experiment specified above')
		self.b223.clicked.connect(self.on_click223)
		#self.b223.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.VOLlayout.addWidget(self.b223, 9, 0)

		return(VOLTab)

	def PARtab(self): # more detailed particle-phase plotting tab definition

		PARTab = QWidget()
		self.PARlayout = QGridLayout() 
		PARTab.setLayout(self.PARlayout)
		
		# input bar for number of components contributing
		# to particle-phase concentration
		self.e300 = QLineEdit(self)
		self.e300.setText('Provide the top number of components contributing to particle-phase')
		self.e300.setStyleSheet('qproperty-cursorPosition : 0')
		self.PARlayout.addWidget(self.e300, 0, 0, 1, 1)
	
		# button to plot particle-phase contributions
		self.b300 = QPushButton(str('Particle-phase contributions (%)'))
		self.b300.setToolTip('Show the contribution to particle-phase by the top contributors (number given in box above)')
		self.b300.clicked.connect(self.on_click300)
		#self.b300.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.PARlayout.addWidget(self.b300, 1, 0)

		# drop down button for type of contributors
		self.b300a = QComboBox(self)
		self.b300a.addItem('All components')
		self.b300a.addItem('Excluding Seed and Water')
		self.PARlayout.addWidget(self.b300a, 1, 1, 1, 1)

		# button to plot particle-phase surface concentration
		self.b301 = QPushButton(str('Particle surface area'))
		self.b301.setToolTip('Graph the temporal profile of particle surface area')
		self.b301.clicked.connect(self.on_click301)
		#self.b301.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.PARlayout.addWidget(self.b301, 2, 0)

		# drop down button for type of surface area
		self.b301a = QComboBox(self)
		self.b301a.addItem('All components')
		self.b301a.addItem('Seed Only')
		self.PARlayout.addWidget(self.b301a, 2, 1, 1, 1)

		# button to plot particle-phase mass contribution by different 
		# generations of oxidised organic molecules
		self.b302 = QPushButton(str('Organic molecule \ncontribution by generation'))
		self.b302.setToolTip('See the particle-phase mass contribution by different generations of oxidised organic molecules')
		self.b302.clicked.connect(self.on_click302)
		self.b302.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.PARlayout.addWidget(self.b302, 0, 2, 1, 1)

		# horizontal separator line -------------------------------
		self.separatorLine5 = QFrame()
		self.separatorLine5.setFrameShape(QFrame.HLine)
		self.separatorLine5.setFrameShadow(QFrame.Raised)
		self.PARlayout.addWidget(self.separatorLine5, 3, 0, 1, 2)
		self.separatorLine5.show()
	
		# section for consumption and yield calculations ------------------------------
		# input bar for component to estimate consumption for
		self.e224 = QLineEdit(self)
		self.e224.setText('Provide the chemical scheme name of the component to view consumption/yield for (result displayed in message box above)')
		self.e224.setStyleSheet('qproperty-cursorPosition : 0')
		self.PARlayout.addWidget(self.e224, 4, 0, 1, 2)

		# input bar for starting time to estimate consumption for
		self.e224a = QLineEdit(self)
		self.e224a.setText('Provide the starting time to calculate consumption for (hours)')
		self.e224a.setStyleSheet('qproperty-cursorPosition : 0')
		self.PARlayout.addWidget(self.e224a, 5, 0, 1, 2)

		# input bar for finshing time to estimate consumption for
		self.e224b = QLineEdit(self)
		self.e224b.setText('Provide the finishing time to calculate consumption for (hours)')
		self.e224b.setStyleSheet('qproperty-cursorPosition : 0')
		self.PARlayout.addWidget(self.e224b, 6, 0, 1, 2)

		# button to estimate consumption
		self.b224 = QPushButton(str('Consumption (' + u'\u03BC' + 'g/m' + u'\u00B3' +')'))
		self.b224.setToolTip('For the component specified above show the mass concentration consumed throughout the whole simulation in the message box above')
		self.b224.clicked.connect(self.on_click224)
		#self.b224.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.PARlayout.addWidget(self.b224, 7, 0, 1, 1)

		# button to estimate SOA yield
		self.b225 = QPushButton(str('SOA yield (fraction 0-1)'))
		self.b225.setToolTip('In the message box above show the SOA yield from the component given in the box above between the times stated in the boxes above.')
		self.b225.clicked.connect(self.on_click225)
		#self.b225.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.PARlayout.addWidget(self.b225, 7, 1, 1, 1)
		
		return(PARTab)

	def RADtab(self): # more detailed plotting tab for radical chemicals

		RADTab = QWidget()
		self.RADlayout = QGridLayout() 
		RADTab.setLayout(self.RADlayout)
		
		# input bar for number of components contributing
		# to radical population
		self.e360 = QLineEdit(self)
		self.e360.setText('Provide the top number of components contributing to a radical pool')
		self.e360.setStyleSheet('qproperty-cursorPosition : 0')
		self.RADlayout.addWidget(self.e360, 0, 0, 1, 2)

		# drop down button for type of radical
		self.b361a = QComboBox(self)
		self.b361a.addItem('RO2 (alkyl peroxy radical) Number Contribution (%)')
		self.b361a.addItem('RO (alkoxy radical) Number Contribution (%)')
		self.b361a.addItem('RO (alkoxy radical) Flux (molecules/cm3/s)')
		self.RADlayout.addWidget(self.b361a, 1, 0, 1, 1)
	
		# button to plot radical output
		self.b362 = QPushButton(str('Plot'))
		self.b362.setToolTip('Show the selected output by the top number of components (number given in box above)')
		self.b362.clicked.connect(self.on_click362)
		#self.b362.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.RADlayout.addWidget(self.b362, 1, 1)

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
		self.e280.setText('Mass:charge resolution determinants: interval of probability peaks (m:z), width of distribution per peak.  Defaults to 1., 0.3 ')
		self.CIMSscrolllayout.addWidget(self.e280, 0, 0)
		
		# input for time to show mass spectrum for
		self.e281 = QTextEdit(self)
		self.e281.setText('Time through experiment to show mass spectrum for (s)')
		self.CIMSscrolllayout.addWidget(self.e281, 1, 0)
		
		# type of ionisation source
		self.e282 = QTextEdit(self)
		self.e282.setText('Ionisation source (I for iodide, N for nitrate, defaults to I), whether to add molar mass to mass of ionised components (1 for yes), 0 (for no), defaults to 0 for no')
		self.CIMSscrolllayout.addWidget(self.e282, 2, 0)
		
		# sensitivity dependence on molar mass
		self.e283 = QTextEdit(self)
		self.e283.setText('Sensitivity (Hz/ppt) dependence on molar mass (g/mol), use y_MW to denote molar mass (g/mol) of components.  Defaults to 1.0 which implies no dependency on molar mass.')
		self.CIMSscrolllayout.addWidget(self.e283, 3, 0)
		
		# button to plot probability distribution function demonstrating mass:charge resolution
		self.b290_aa = QPushButton('PDF graph', self)
		self.b290_aa.setToolTip('Plot the probability distribution function causing mass:charge resolution')
		self.b290_aa.clicked.connect(self.on_click290_aa)
		self.CIMSscrolllayout.addWidget(self.b290_aa, 0, 1)

		# button to plot sensitivity dependence on molar mass
		self.b290_a = QPushButton('Sensitivity', self)
		self.b290_a.setToolTip('Plot the sensitivity to molar mass')
		self.b290_a.clicked.connect(self.on_click290_a)
		self.CIMSscrolllayout.addWidget(self.b290_a, 3, 1)
	
		# button to plot mass spectrum
		self.b290 = QPushButton('CIMS observations', self)
		self.b290.setToolTip('Plot the mass spectrum as observed by a chemical ionisation mass spectrometer')
		self.b290.clicked.connect(self.on_click290)
		self.CIMSscrolllayout.addWidget(self.b290, 4, 0)
		
		# properties of CIMS scroll area ----------------
		self.scrollwidget.setLayout(self.CIMSscrolllayout)
		self.scroll.setWidget(self.scrollwidget)
		self.CIMSlayout.addWidget(self.scroll, 0, 0, 1, 1)
	
		return(CIMSTab)
	
	def OBStab(self): # comparing with observations

		OBSTab = QWidget()
		self.OBSlayout = QGridLayout() 
		OBSTab.setLayout(self.OBSlayout)

		# show path to observations file --------------------------------------------------------
		self.l401 = ScrollLabel(self)
		cwd = os.getcwd() # current working directory
		path = str(cwd + '/PyCHAM/output/26_10.xlsx')
		self.l401.setText(path)
		self.OBSlayout.addWidget(self.l401, 0, 0, 1, 1)

		# select observations file --------------------------------------------------------------
		b400 = QPushButton('Select File Containing Observations', self)
		b400.setToolTip('Select the file containing the required observations')
		b400.clicked.connect(self.on_click400)
		self.OBSlayout.addWidget(b400, 0, 1, 1, 1)

		# plot observations and model results
		# select observations file --------------------------------------------------------------
		b401 = QPushButton('Observed && Modelled', self)
		b401.setToolTip('Plot Observations and Model Results as Described in Quick tab')
		b401.clicked.connect(self.on_click401)
		self.OBSlayout.addWidget(b401, 1, 0, 1, 1)
		
		return(OBSTab)

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
			self.fab = 0 # remember that single simulation widgets not showing
		if (self.atb == 1): # if showing then remove add to batch button
			self.b82.deleteLater()
			self.atb = 0 # remember that add to batch button not showing
		# remove any old message from previous run
		if (len(self.l81b.text()) > 0):
			if (self.l81b.text()[0:17] != 'File combinations'):
				self.l81b.setText('')
				self.l81b.setStyleSheet(0., '0px', 0., 0.)
		
		# prepare by enforcing default variables
		# default variables for all required input model variables -------------------------
		[sav_nam, chem_sch_mark, update_stp, tot_time, comp0, 
		y0, temp, tempt, RH, RHt, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, 
		pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, Compt, 
		injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedx, 
		con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, 
		act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, 
		nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, 
		p_char, e_field, dil_fac, partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str, pcont, 
		Vwat_inc, seed_eq_wat, z_prt_coeff, chamV, self] = def_mod_var.def_mod_var(0, self)
		
		# then open default variables, ready for modification
		input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
		
		with open(input_by_sim, 'rb') as pk:
			[sav_nam, chem_sch_mark, update_stp, 
			tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on,
			Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
			save_step, Compt, injectt, Ct, seed_name,
			seed_mw, seed_diss, seed_dens, seedx,
			con_infl_nam, con_infl_t, con_infl_C, 
			dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
			accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
			nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
			inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, 
			wat_hist, drh_str, erh_str, pcont, Vwat_inc, seed_eq_wat, z_prt_coeff, 
			chamV] = pickle.load(pk)
			pk.close()
		
		# button to get path to folder containing relevant files
		options = QFileDialog.Options()
		fol_nme = QFileDialog.getExistingDirectory(self, "Select Folder Containing Required Input Files", "./PyCHAM/input/")
		
		if (fol_nme == ''): # if no folder selected (e.g. because selection cancelled)
			return()
		
		# list of files (and any folders) here
		dir_con = os.listdir(fol_nme)
		
		# unknown names
		self.sch_name = self.xml_name = self.inname = 'Not found'
		
		# look for corresponding files here
		for i in dir_con:
			if ('chem' in i): # chemical scheme file
				self.sch_name = str(fol_nme+'/'+i)
			if ('xml' in i): # xml file
				self.xml_name = str(fol_nme+'/'+i)
			if ('var' in i): # model variables file
				self.inname = str(fol_nme+'/'+i)
			
		# read in model variables of this model variables file and store to pickle
		import mod_var_read
		mod_var_read.mod_var_read(self)

		# end this function if an error thrown by reading of model variables
		if (self.bd_st == 1 or self.bd_st == 2):
			return()
			
		# updating scroll labels showing path to files
		self.l3.setText(self.sch_name)
		self.l5.setText(self.xml_name)
		self.l7.setText(self.inname)
			
		self.show()
		
		import mod_var_up # update displayed model variables to those of the selected model variables file
		mod_var_up.mod_var_up(self)
		
		# running check on model variables ---------------------------------------------------------------
		import ui_check # module for checking on model variables
		# check on inputs - note this loads the last saved pickle file and saves any change to this pickle file
		ui_check.ui_check(self)
		# finished check on model variables --------------------------------------------------------------		

		return()
		
	@pyqtSlot()
	def on_click2(self): # when different chemical scheme file requires selection
	
		self.sch_name, _ = QFileDialog.getOpenFileName(self, "Select Chemical Scheme File", "./PyCHAM/input/") # get path of file
		self.l3.clear() # clear old label
		self.l3.setText(str(self.sch_name))
		self.l3.show()
	
		# running check on model variables ---------------------------------------------------------------
		import ui_check # module for checking on model variables
		# check on inputs - note this loads the last saved pickle file and saves any change to this pickle file
		ui_check.ui_check(self)
		# finished check on model variables --------------------------------------------------------------
	
	@pyqtSlot()
	def on_click3(self): # when different xml file requires selection
		
		self.xml_name, _ = QFileDialog.getOpenFileName(self, "Select xml File", "./PyCHAM/input/") # get path of file
		self.l5.clear() # clear old label
		self.l5.setText(str(self.xml_name))
		self.l5.show()
			
	@pyqtSlot()
	def on_click4(self): # when different model variables file requires selection
	
		if (self.fab == 1): # if showing, remove single simulation widgets
			self.b81.deleteLater()
			self.l81.deleteLater()
			self.fab = 0 # remember that single simulation widgets not showing
		if (self.atb == 1): # if showing then remove add to batch button
			self.b82.deleteLater()
			self.atb = 0 # remember that add to batch button not showing
		# remove any old 'Simulation complete' message from previous run
		if ((self.l81b.text() == 'Simulation complete') or (self.l81b.text() == 'Simulations complete')):
			self.l81b.setText('')

		# user chooses path of file to model variables file
		self.inname, _ = QFileDialog.getOpenFileName(self, "Select Model Variables File", "./PyCHAM/input/")
		
		if (self.inname == ''): # if no file selected, e.g. because selection was cancelled
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
	
		import mod_var_up # update displayed model variables to those of the selected model variables file
		mod_var_up.mod_var_up(self)
		
		# running check on model variables ---------------------------------------------------------------
		import ui_check # module for checking on model variables
		# check on inputs - note this loads the last saved pickle file and saves any change to this pickle file
		ui_check.ui_check(self)
		# finished check on model variables --------------------------------------------------------------

		return()
			
	
	@pyqtSlot()	
	def on_click81b(self): # when button to run simulation pressed
		
		for sim_num in range(self.btch_no):
			
			# --------------------------------------------
			if ((self.btch_no > 1) and (sim_num == self.btch_no-1)): # reach end of list in batch mode
				# once all simulations done tidy up
				# clear the list of file combinations for simulations in batch
				self.output_list = [] # reset list of output paths
				# return to single simulation mode
				self.btch_no = 1
				self.btch_str = 'File combinations included in batch: chemical scheme, xml, model variables\n'
				
				# remove old progress message
				self.l81b.setText('')
				# tell user that simulations finished
				self.l81b.setText(str('Simulations complete'))
				return(err_mess)
			
			# --------------------------------------------
		
			# reset to default variables to allow any new variables to arise
			# from the current model variables file only
			[sav_nam, chem_sch_mark, update_stp, 
			tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on, Cw, 
			kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, 
			space_mode, std, mean_rad, save_step, Compt, injectt, 
			Ct, seed_name, seed_mw, seed_diss, seed_dens, seedx, 
			con_infl_nam, con_infl_t, con_infl_C, dydt_trak, 
			dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, 
			accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, 
			nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
			inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, 
			ser_H2O, wat_hist, drh_str, erh_str, pcont, Vwat_inc, 
			seed_eq_wat, z_prt_coeff, chamV, self] = def_mod_var.def_mod_var(0, self)
			
			# get text from batch list label
			btch_list = self.btch_str
			if ((sim_num+1) < (self.btch_no-1)): # if prior to final item in list
				# get location of text relevant to this simulation
				txt_st = btch_list.find(str(str(sim_num+1)+'.'+'\n'))+len(str(str(sim_num+1)+'.'+'\n'))
				txt_fi = btch_list.find(str(str(sim_num+2)+'.'+'\n'))
				txtn = btch_list[txt_st:txt_fi]
			else: # if final item in list
				# get location of text relevant to this simulation
				txt_st = btch_list.find(str(str(sim_num+1)+'.'+'\n'))+len(str(str(sim_num+1)+'.'+'\n'))
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
			input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
			with open(input_by_sim, 'rb') as pk:
				[sav_nam, chem_sch_mark, update_stp, 
				tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on,
				Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
				save_step, Compt, injectt, Ct, seed_name,
				seed_mw, seed_diss, seed_dens, seedx,
				con_infl_nam, con_infl_t, con_infl_C, 
				dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
				accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
				nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
				inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, 
				wat_hist, drh_str, erh_str, pcont, Vwat_inc, seed_eq_wat, z_prt_coeff, 
				chamV] = pickle.load(pk)
				pk.close() # close pickle file
			
			# let check know this is a second call
			self.chck_num = 2

			# run another check on inputs - means any changes made by default are set
			import ui_check; ui_check.ui_check(self)

			# reset check number to a first call
			self.chck_num = 1

			# saving path - copied from saving module
			dir_path = os.getcwd() # current working directory
			output_root = 'PyCHAM/output'
			filename = os.path.basename(self.sch_name)
			filename = os.path.splitext(filename)[0]
			output_by_sim = os.path.join(dir_path, output_root, filename, sav_nam)
			
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
			err_log = str(os.getcwd() + '/PyCHAM/err_log.txt')
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
			
			err_mess = self.act_81(output_by_sim, sim_num) # call on function to simulate

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
	def act_81(self, output_by_sim, sim_num): # start the simulation
	
		from middle import middle # prepare to communicate with main program
		
		note_messf = 0 # cancel note message flag
		
		for prog in middle(self): # call on modules to solve problem
			
		
			if (isinstance(prog, str)): # check if it's a message
				mess = prog
				if (mess[0:5] == 'Error'): # if it's an error message
					# remove the progress bar
					self.progress.deleteLater()
					return(mess)
				else:
					self.l81b.setText(mess)
					note_messf = 1 # flag that a note message has been generated 
			
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
			
				self.progress.setValue(prog) # get progress
				
			QApplication.processEvents() # allow progress bar/message panel to update
		
		# remove the progress bar after each simulation
		self.progress.deleteLater()
		
		# set the path to folder to plot results to the latest simulation results
		self.l201.setText(output_by_sim)
		# if this point reached then no error message generated
		mess = ''
		return(mess)
	
	@pyqtSlot()
	def on_click81sing(self): # when single simulation button pressed
	
		cs_file = self.l3.text()
		xml_file = self.l5.text()
		mv_file = self.l7.text()
		self.btch_no = 1 # ensure just one simulation when single simulation button pressed
		self.btch_str = (str(str(self.btch_no)+'.\n' + cs_file + '\n' + xml_file + '\n' + mv_file + '\n'))
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
	def on_click202(self): # when different model results folder requires selection

		# button to get path to folder containing relevant files
		options = QFileDialog.Options()
		fol_nme = QFileDialog.getExistingDirectory(self, "Select Folder Containing Required Input Files", "./PyCHAM/output/")
		
		self.l201.clear() # clear old label
		self.l201.setText(fol_nme)
		
		# check whether required files present here
		try:
			# name of file where experiment constants saved
			fname = str(fol_nme + '/model_and_component_constants')
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
			self.b220.setEnabled(True)
			self.b221.setEnabled(True)
			self.b223.setEnabled(True)
			
			# clear error message
			self.l203a.setText('')
			self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
			self.bd_pl = 3
			
		except:
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
			self.b220.setEnabled(False)
			self.b221.setEnabled(False)
			self.b223.setEnabled(False)
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

		plotter.plotter(0, dir_path, uc, self) # plot results
	
	@pyqtSlot() # button to plot gas-phase concentration temporal profile
	def on_click206(self):	
		
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
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
	
	
	# button to plot total particle-phase concentration for individual components temporal profile
	@pyqtSlot()
	def on_click209(self):	
		
		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		# get names of components to plot
		comp_names = [str(i) for i in self.e205.text(). split(',')]
		
		import plotter_pp
		dir_path = self.l201.text() # name of folder with results
		plotter_pp.plotter(0, dir_path, comp_names, self) # plot results

	# button to plot temporal profile of total particle-phase concentration
	# excluding seed and water
	@pyqtSlot()
	def on_click209a(self):	
		
		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		# get names of components to plot
		comp_names = []
		
		import plotter_pp
		dir_path = self.l201.text() # name of folder with results
		plotter_pp.plotter(3, dir_path, comp_names, self) # plot results

	# button to plot temporal profile of particle-phase contribution 
	# by top contributors to particle-phase concentration
	@pyqtSlot()
	def on_click300(self):	
		
		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
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

		caller = 8
		
		# names of components (filler)
		comp_names = []

		import plotter_pp
		dir_path = self.l201.text() # name of folder with results
		plotter_pp.plotter(caller, dir_path, comp_names, self) # plot results
	
	@pyqtSlot() # button to plot temporal profile of total concentration of 
	# components that have gas-wall partitioned to wall
	def on_click212(self):	
		
		# clear dialogue
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		# get names of components to plot
		comp_names = [str(i) for i in self.e205.text(). split(',')]
		
		import plotter_wp
		dir_path = self.l201.text() # name of folder with results
		plotter_wp.plotter(0, dir_path, comp_names, self) # plot results
	
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
		top_num = [int(i) for i in self.e217a.text().split(',')]

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
		plotter_ct.plotter_ind(0, dir_path, comp_names, top_num, uc, self) # plot results
	
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
		
	# button to plot volatility basis set mass fractions with water
	@pyqtSlot()
	def on_click220(self):

		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		import vol_contr_analys	
		dir_path = self.l201.text() # name of folder with results
		vol_contr_analys.plotter_wiw(0, dir_path, self, 0) # plot results
	
	# button to plot volatility basis set mass fractions without water
	@pyqtSlot()
	def on_click221(self):

		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
			
		import vol_contr_analys
		dir_path = self.l201.text() # name of folder with results
		vol_contr_analys.plotter_wiw(0, dir_path, self, 1) # plot results
	
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

	@pyqtSlot() # button to plot probability distribution function demonstrating mass:charge resolution
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
		
		y_MW = np.arange(1000.)
		
		# get sensitivity (Hz/ppt) dependence on molar mass
		try:
			sensit = str((self.e283.toPlainText()))
		except:
			sensit = 'np.ones(len(y_MW))' # default
		if (sensit[0:3] == 'Sen' or sensit == ''): # means that the edit label text has not been changed from the description
			sensit = 'np.ones(len(y_MW))'
		
		import plotter_CIMS
		
		blank = plotter_CIMS.write_sens2mm(3, sensit, y_MW)
		
	@pyqtSlot() # button to plot mass spectrum replication of chemical ionisation mass spectrometer
	def on_click290(self):
	
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
			
		# get time to plot
		try:
			tn = float((self.e281.toPlainText()))
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
			self.l203a.setText('Note - ionising agent, whether molar mass of ionising agent to be included in molar mass of components, should be a letter (I for iodide or N for nitrate) followed by a comma, followed by a number (1 to add agent molar mass to component mass or 0 not to).  But this information could not be correctly detected, so defaulting to I, 0')
			
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
		if (sensit[0:3] == 'Sen' or sensit == ''): # means that the edit label text has not been changed from the description
			sensit = '1.' # default
		
		# import required plotting function
		import plotter_CIMS
		plotter_CIMS.plotter_CIMS(dir_path, res_in, tn, iont, sensit)

	@pyqtSlot() # button to check supplied values
	def on_clickb80(self):
		
		# get values to check from drop down button
		input_check_text = self.b80s.currentText() # drop down selection

		if (input_check_text == 'Photolysis rates'):
			import plotter_simulate_tab
			plotter_simulate_tab.plotter_taf(self)

		if (input_check_text == 'Particle number size distributions'):
			import plotter_nsd # for plotting supplied number size distributions
		
			# path for pickle file
			input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')

			# get the most recent model variables
			with open(input_by_sim, 'rb') as pk:
				[sav_nam, chem_sch_mark, update_stp, 
				tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on,
				Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, 
				std, mean_rad, save_step, Compt, injectt, Ct, seed_name,
				seed_mw, seed_diss, seed_dens, seedx,
				con_infl_nam, con_infl_t, con_infl_C, 
				dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
				accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
				nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
				inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, 
				wat_hist, drh_str, erh_str, pcont, Vwat_inc, seed_eq_wat, z_prt_coeff, 
				chamV] = pickle.load(pk)
				pk.close()
		
			# call on plotting script
			plotter_nsd.plotter_nsd(lowsize, num_sb, uppsize, mean_rad, std, pmode, pconc, 
			space_mode, 0, pconct)

		if (input_check_text == 'Gas-phase diffusion coefficients'):
			import plotter_simulate_tab
			plotter_simulate_tab.plotter_gpdc(self)

		if (input_check_text == 'Gas-phase mean thermal speeds'):
			import plotter_simulate_tab
			plotter_simulate_tab.plotter_gpmts(self)

		if (input_check_text == 'Molar masses'):
			import plotter_simulate_tab
			plotter_simulate_tab.plotter_mm(self)

	@pyqtSlot() # button to retrieve and report component consumption
	def on_click224(self):
	
		dir_path = self.l201.text() # name of folder with results

		# get component name
		try:
			comp_chem_schem_name = str((self.e224.text()))

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
		consumption.cons(comp_chem_schem_name, dir_path, self, 0)

	@pyqtSlot() # button to retrieve and report yield
	def on_click225(self):
	
		dir_path = self.l201.text() # name of folder with results

		# get component name
		try:
			comp_chem_schem_name = str((self.e224.text()))

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
		consumption.cons(comp_chem_schem_name, dir_path, self, 1)
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
		
		# if RO2
		if ('RO2 (alkyl peroxy radical)' in input_check_text):
			self.rad_mark = 0

		# if RO
		if ('RO (alkoxy radical) Number' in input_check_text):
			self.rad_mark = 1

		if ('RO (alkoxy radical) Flux' in input_check_text):
			self.rad_mark = 2
		
		self.dir_path = self.l201.text() # name of folder with results


		import plotter_gp # required module
		if self.rad_mark == 0 or self.rad_mark == 1:
			# call on gas-phase plotter for radical pools
			plotter_gp.plotter_rad_pool(self) # plot results

		if (self.rad_mark == 2):
			# call on gas-phase plotter for radical flux
			plotter_gp.plotter_rad_flux(self) # plot results

		return()

	@pyqtSlot()
	def on_click400(self): # when observation file requires selection

		# button to get path to folder containing relevant files
		options = QFileDialog.Options()
		fil_nme = QFileDialog.getOpenFileName(self, "Select File Containing Required Observations", "./PyCHAM/output/")[0]
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
		# get names of gas-phase components to plot
		self.gp_names = [str(i) for i in self.e205.text(). split(',')]
		# gas-phase concentration units
		self.gp_units = self.b206b.currentText()
		
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
		
		if (self.gp_names == ['']):
			self.l203a.setText('Note, no components specified to be plotted for model results')
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
				self.bd_pl = 1
		import plotter_xls
		plotter_xls.plotter_gp_mod_n_obs(self)		
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