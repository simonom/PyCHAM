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
		#bn1.setToolTip('Visit the EUROCHAMP website')
		#bn1.clicked.connect(self.on_clickn1)
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
		[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str] = def_mod_var.def_mod_var(0)
		
		# listing input files -----------------------------------------------------------
		l1 = QLabel(self)
		l1.setText("The following files have been found: ")
		self.NSlayout.addWidget(l1, 1, ffscn, 1, 2)
		
		l2 = QLabel(self)
		l2.setText("Chemical \nscheme: ")
		self.NSlayout.addWidget(l2, 2, ffscn, 1, 1)
			
		self.l3 = ScrollLabel(self)
		self.l3.setText(sch_name)
		self.NSlayout.addWidget(self.l3, 2, ffscn+1,  1, 1)
		
		b1 = QPushButton('Select Different Chemical scheme', self)
		b1.setToolTip('Select the file containing the required chemical scheme file')
		b1.clicked.connect(self.on_click2)
		self.NSlayout.addWidget(b1, 3, ffscn, 1, 2)
		
		l4 = QLabel(self)
		l4.setText("xml: ")
		self.NSlayout.addWidget(l4, 4, ffscn, 1, 1)
			
		self.l5 = ScrollLabel(self)
		self.l5.setText(xml_name)
		self.NSlayout.addWidget(self.l5, 4, ffscn+1, 1, 1)
			
		b2 = QPushButton('Select Different xml', self)
		b2.setToolTip('Select the file containing the required xml file')
		b2.clicked.connect(self.on_click3)
		self.NSlayout.addWidget(b2, 5, ffscn, 1, 2)
		
		l6 = QLabel(self)
		l6.setText("Model \nvariables: ")
		self.NSlayout.addWidget(l6, 6, ffscn, 1, 1)
			
		self.l7 = ScrollLabel(self)
		self.l7.setText(inname)
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
		l13_2.setText('The absolute and relative integration \ntolerances for use: ')
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
		l26.setText('Volume ratio of components in \nseed particles: ')
		self.varbox.addWidget(l26, par_row+9, 0)
		self.l26a = QLabel(self)
		self.l26a.setText((str(seedVr)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l26a, par_row+9, 1)
		
		l27 = QLabel(self)
		l27.setText('Smallest radius size bin \nboundary (um): ')
		self.varbox.addWidget(l27, par_row+10, 0)
		self.l27a = QLabel(self)
		self.l27a.setText((str(lowsize)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l27a, par_row+10, 1)
		
		l28 = QLabel(self)
		l28.setText('Largest radius size bin \nboundary (um): ')
		self.varbox.addWidget(l28, par_row+11, 0)
		self.l28a = QLabel(self)
		self.l28a.setText((str(uppsize)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l28a, par_row+11, 1)
		
		l29 = QLabel(self)
		l29.setText('Method for spacing size bins (lin) \nor (log): ')
		self.varbox.addWidget(l29, par_row+12, 0)
		self.l29a = QLabel(self)
		self.l29a.setText((str(space_mode)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l29a, par_row+12, 1)
		
		l30 = QLabel(self)
		l30.setText('Standard deviation for particle \nnumber size distribution: ')
		self.varbox.addWidget(l30, par_row+13, 0)
		self.l30a = QLabel(self)
		self.l30a.setText((str(std)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l30a, par_row+13, 1)
		
		l31 = QLabel(self)
		l31.setText('Mean radius (um) for particle \nnumber size distribution: ')
		self.varbox.addWidget(l31, par_row+14, 0)
		self.l31a = QLabel(self)
		self.l31a.setText((str(mean_rad)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l31a, par_row+14, 1)
		
		l32 = QLabel(self)
		l32.setText('Radius (um) of newly nucleated \nparticles: ')
		self.varbox.addWidget(l32, par_row+15, 0)
		self.l32a = QLabel(self)
		self.l32a.setText((str(new_partr)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l32a, par_row+15, 1)
		
		l33 = QLabel(self)
		l33.setText('First nucleation \nparameterisation parameter: ')
		self.varbox.addWidget(l33, par_row+16, 0)
		self.l33a = QLabel(self)
		self.l33a.setText((str(nucv1)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l33a, par_row+16, 1)
		
		l34 = QLabel(self)
		l34.setText('Second nucleation \nparameterisation parameter: ')
		self.varbox.addWidget(l34, par_row+17, 0)
		self.l34a = QLabel(self)
		self.l34a.setText((str(nucv2)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l34a, par_row+17, 1)
		
		l35 = QLabel(self)
		l35.setText('Third nucleation \nparameterisation parameter: ')
		self.varbox.addWidget(l35, par_row+18, 0)
		self.l35a = QLabel(self)
		self.l35a.setText((str(nucv3)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l35a, par_row+18, 1)
		
		l36 = QLabel(self)
		l36.setText('Chemical scheme name of \nnucleating component: ')
		self.varbox.addWidget(l36, par_row+19, 0)
		self.l36a = QLabel(self)
		self.l36a.setText((str(nuc_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l36a, par_row+19, 1)
		
		l37 = QLabel(self)
		l37.setText('Whether (1) or not (0) to adapt \nintegration time interval \nto nucleation: ')
		self.varbox.addWidget(l37, par_row+20, 0)
		self.l37a = QLabel(self)
		self.l37a.setText((str(nuc_ad)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l37a, par_row+20, 1)
		
		l38 = QLabel(self)
		l38.setText('Whether (1) or not (0) to serialise \ngas-particle partitioning of water: ')
		self.varbox.addWidget(l38, par_row+21, 0)
		self.l38a = QLabel(self)
		self.l38a.setText((str(ser_H2O)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l38a, par_row+21, 1)
		
		l38_1 = QLabel(self)
		l38_1.setText('Whether (1) or not (0) to model \nparticle coagulation: ')
		self.varbox.addWidget(l38_1, par_row+22, 0)
		self.l38_1a = QLabel(self)
		self.l38_1a.setText((str(coag_on)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l38_1a, par_row+22, 1)
		
		l38_2 = QLabel(self)
		l38_2.setText('Particle-phase history with \nrespect to water partitioning (0 or 1): ')
		self.varbox.addWidget(l38_2, par_row+23, 0)
		self.l38_2a = QLabel(self)
		self.l38_2a.setText((str(wat_hist)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l38_2a, par_row+23, 1)
		
		l38_3 = QLabel(self)
		l38_3.setText('Deliquescence relative humidity \nas a function of temperature: ')
		self.varbox.addWidget(l38_3, par_row+24, 0)
		self.l38_3a = QLabel(self)
		self.l38_3a.setText((str(drh_str)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l38_3a, par_row+24, 1)
		
		l38_4 = QLabel(self)
		l38_4.setText('Efflorescence relative humidity \nas a function of temperature: ')
		self.varbox.addWidget(l38_4, par_row+25, 0)
		self.l38_4a = QLabel(self)
		self.l38_4a.setText((str(erh_str)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l38_4a, par_row+25, 1)
		
		# gas inputs ----------------
		
		l39a = QLabel(self)
		l39a.setText('Gas')
		l39a.setFont(QFont("Arial", 13, QFont.Bold))
		gas_row = par_row+26
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
		self.l45a.setText((str(const_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
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
		self.l49a.setText((str(light_stat)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l49a, light_row+1, 1)
		
		l50 = QLabel(self)
		l50.setText('Time through simulation that \nlight status attained (s): ')
		self.varbox.addWidget(l50, light_row+2, 0)
		self.l50a = QLabel(self)
		self.l50a.setText((str(light_time)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l50a, light_row+2, 1)
		
		l51 = QLabel(self)
		l51.setText('Time of day experiment starts \n(s since 00:00): ')
		self.varbox.addWidget(l51, light_row+3, 0)
		self.l51a = QLabel(self)
		self.l51a.setText((str(daytime)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l51a, light_row+3, 1)
		
		l52 = QLabel(self)
		l52.setText('Latitude of experiment (degrees): ')
		self.varbox.addWidget(l52, light_row+4, 0)
		self.l52a = QLabel(self)
		self.l52a.setText((str(lat)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l52a, light_row+4, 1)
		
		l53 = QLabel(self)
		l53.setText('Longitude of experiment (degrees): ')
		self.varbox.addWidget(l53, light_row+5, 0)
		self.l53a = QLabel(self)
		self.l53a.setText((str(lat)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l53a, light_row+5, 1)
		
		l54 = QLabel(self)
		l54.setText('Path to customised (non-MCM) \nactinic flux file: ')
		self.varbox.addWidget(l54, light_row+6, 0)
		self.l54a = QLabel(self)
		self.l54a.setText((str(af_path)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l54a, light_row+6, 1)
		
		l55 = QLabel(self)
		l55.setText('Path to file containing absorption \ncross-sections and quantum yields: ')
		self.varbox.addWidget(l55, light_row+7, 0)
		self.l55a = QLabel(self)
		self.l55a.setText((str(photo_path)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l55a, light_row+7, 1)
		
		l56 = QLabel(self)
		l56.setText('Day number of experiment \n(number of days since the \npreceding 31st December): ')
		self.varbox.addWidget(l56, light_row+8, 0)
		self.l56a = QLabel(self)
		self.l56a.setText((str(dayOfYear)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l56a, light_row+8, 1)
		
		l57 = QLabel(self)
		l57.setText('Transmission factor for natural \nsunlight (0-1 fraction): ')
		self.varbox.addWidget(l57, light_row+9, 0)
		self.l57a = QLabel(self)
		self.l57a.setText((str(tf)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l57a, light_row+9, 1)
		
		l58 = QLabel(self)
		l58.setText('Whether to (1) or not to (0) adapt \nintegration time interval and initial \ncondition update to changing \nnatural light intensity: ')
		self.varbox.addWidget(l58, light_row+10, 0)
		self.l58a = QLabel(self)
		self.l58a.setText((str(light_ad)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l58a, light_row+10, 1)
		
		# wall inputs ---------------------------
		
		l59b = QLabel(self)
		l59b.setText('Walls')
		l59b.setFont(QFont("Arial", 13, QFont.Bold))
		wall_row = light_row+11
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
		
		l67 = QLabel(self)
		l67.setText('Whether particle deposition to \nwall treated by Rader and \nMcMurry (1) or customised (0): ')
		self.varbox.addWidget(l67, wall_row+9, 0)
		self.l67a = QLabel(self)
		self.l67a.setText((str(Rader)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l67a, wall_row+9, 1)
		
		l68 = QLabel(self)
		l68.setText('Average number of charges per \nparticle (/particle): ')
		self.varbox.addWidget(l68, wall_row+10, 0)
		self.l68a = QLabel(self)
		self.l68a.setText((str(p_char)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l68a, wall_row+10, 1)
		
		l69 = QLabel(self)
		l69.setText('Average electric field inside chamber \n(g.m/A.s3): ')
		self.varbox.addWidget(l69, wall_row+11, 0)
		self.l69a = QLabel(self)
		self.l69a.setText((str(e_field)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.l69a, wall_row+11, 1)
		
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
		import ui_check # module for checking on model variables
		# check on inputs - note this loads the last saved pickle file and saves any change to this pickle file
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
		l202a.setText('Messages from plotting scripts: ')
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
		PLtabs.addTab(self.PRIMtab(), "Basic")
		PLtabs.addTab(self.SECtab(), "Detailed")
		PLtabs.addTab(self.INSTRtab(), "Under Development - not fully functional")
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


		# input bar for names of components to plot temporal profiles of
		self.e205 = QLineEdit(self)
		self.e205.setText('Chemical scheme names of components to be plotted using the buttons below (use H2O for water)')
		# show left most point first
		self.e205.setStyleSheet('qproperty-cursorPosition : 0')
		self.PRIMlayout.addWidget(self.e205, 1, 0)

		# gas-phase concentrations temporal profiles -------------
		
		# button to plot temporal profile of gas-phase concentrations
		self.b206 = QPushButton('Gas-phase concentrations', self)
		self.b206.setToolTip('Plot gas-phase concentration temporal profile for the specified components')
		self.b206.clicked.connect(self.on_click206)
		self.PRIMlayout.addWidget(self.b206, 2, 0)
		
		 # particle-phase concentrations temporal profiles -------------

		
		# button to plot temporal profile of total particle-phase concentrations
		self.b209 = QPushButton('Total particle-phase concentrations', self)
		self.b209.setToolTip('Plot particle-phase concentration temporal profile of these components')
		self.b209.clicked.connect(self.on_click209)
		self.PRIMlayout.addWidget(self.b209, 3, 0)
		
		# wall (from gas-wall partitioning) concentrations temporal profiles -------------
		
		# button to plot temporal profile of total particle-phase concentrations
		self.b212 = QPushButton('Wall concentrations (from gas-wall partitioning)', self)
		self.b212.setToolTip('Plot the temporal profile of wall concentration (from gas-wall partitioning) for the specified components')
		self.b212.clicked.connect(self.on_click212)
		self.PRIMlayout.addWidget(self.b212, 4, 0)
		
		# wall (from particle deposition to wall) concentrations temporal profiles -------------
		
		# button to plot temporal profile of wall concentrations (from particle deposition to wall)
		self.b215 = QPushButton('Wall concentrations (from particle deposition to wall)', self)
		self.b215.setToolTip('Plot the temporal profile of the wall concentration (from particle deposition to wall) for the specified components')
		self.b215.clicked.connect(self.on_click215)
		self.PRIMlayout.addWidget(self.b215, 5, 0)
		
		# chamber conditions temporal profiles -------------
		
		# button to plot temporal profile of chamber environmental variables
		self.b215_a = QPushButton('Chamber Conditions (T, P, RH)', self)
		self.b215_a.setToolTip('Plot the temporal profile of chamber variables (temperature, pressure, relative humidity)')
		self.b215_a.clicked.connect(self.on_click215_a)
		self.PRIMlayout.addWidget(self.b215_a, 0, 1)
		
		return(PRIMTab)
	
	def SECtab(self): # more detailed plotting tab definition

		SECTab = QWidget()
		self.SEClayout = QGridLayout() 
		SECTab.setLayout(self.SEClayout)
		
		# input bar for names of components to plot change tendencies
		self.e217 = QLineEdit(self)
		self.e217.setText('Provide the chemical scheme names of components for plotting change tendencies')
		self.e217.setStyleSheet('qproperty-cursorPosition : 0')
		self.SEClayout.addWidget(self.e217, 0, 0)
		
		# button to plot temporal profile of change tendencies
		self.b218 = QPushButton('Plot change tendencies', self)
		self.b218.setToolTip('Plot the rate of change due to relevant processes')
		self.b218.clicked.connect(self.on_click218)
		self.SEClayout.addWidget(self.b218, 1, 0)
		
		# volatility basis set ------------------
		
		# label for names of components to plot tracked change tendencies
		l219 = QLabel(self)
		l219.setText('Buttons below plot mass fractions of components grouped by vapour pressure at 298.15 K')
		l219.setWordWrap(True)
		self.SEClayout.addWidget(l219, 2, 0, 1, 1)
		
		# button to plot temporal profile of volatility basis set mass fractions with water
		self.b220 = QPushButton('Plot Volatility Basis Set With Water', self)
		self.b220.setToolTip('Plot the temporal profile of volatility basis set mass fractions')
		self.b220.clicked.connect(self.on_click220)
		#self.b220.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.SEClayout.addWidget(self.b220, 3, 0)
		
		# button to plot temporal profile of volatility basis set mass fractions without water
		self.b221 = QPushButton('Plot Volatility Basis Set Without Water', self)
		self.b221.setToolTip('Plot the temporal profile of volatility basis set mass fractions')
		self.b221.clicked.connect(self.on_click221)
		#self.b221.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.SEClayout.addWidget(self.b221, 4, 0)
		
		# two-dimensional volatility basis set ------------------
		
		# input bar for time through experiment to plot the 2D VBS for
		self.e222 = QLineEdit(self)
		self.e222.setText('Provide the time (seconds) through experiment at which to plot the two-dimensional volatility basis set - the closest recorded time to this will be used')
		self.e222.setStyleSheet('qproperty-cursorPosition : 0')
		self.SEClayout.addWidget(self.e222, 5, 0)
		
		# button to plot 2D VBS
		self.b223 = QPushButton('Plot 2D Volatility Basis Set', self)
		self.b223.setToolTip('Plot the two-dimensional volatility basis set (O:C ratio and vapour pressures) at the time through experiment specified above')
		self.b223.clicked.connect(self.on_click223)
		#self.b223.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.SEClayout.addWidget(self.b223, 6, 0)
		
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
	
		CPCTab = QWidget()
		self.CPClayout = QGridLayout() 
		CPCTab.setLayout(self.CPClayout)
	
		# input for equilibrium humidity on reaching the condensing 
		# section of the condensation particle counter
		self.e220 = QTextEdit(self)
		self.e220.setText('Relative humidity (fraction (0-1)) on reaching CPC condensing unit (particles assumed to equilibrate to this) (defaults to 0.65)')
		self.CPClayout.addWidget(self.e220, 0, 0)
		
		# input for minimum particle concentration detectable
		self.e220_a = QTextEdit(self)
		self.e220_a.setText('False background counts (used as the minimum detectable particle concentration) (# particles cm<sup>-3</sup>, defaults to 1.e-2)')
		self.CPClayout.addWidget(self.e220_a, 0, 1)
		
		# input for maximum particle concentration detectable
		self.e220_b = QTextEdit(self)
		self.e220_b.setText('Maximum detectable particle concentration (# particles cm<sup>-3</sup> (e.g. 3.e5), defaults to -1, which implies no maximum)')
		self.CPClayout.addWidget(self.e220_b, 0, 2)
		
		# input for coincidence correction
		self.e220_j = QTextEdit(self)
		self.e220_j.setText('Coincidence inputs: volumetric flow rate through counting unit (cm<sup>3</sup> s<sup>-1</sup>), instrument dead time (s) and upper limit of actual particle concentration (# particles cm<sup>-3</sup>) this can be applied to; should be three numbers separated by a comma (e.g. 5., 2.e-6, 3.e5).  Defaults to -1, -1, -1, which indicates no convolution needed for coincidence (e.g. because coincidence already corrected for by instrument)')
		self.CPClayout.addWidget(self.e220_j, 0, 3)
		
		# input for detection efficiency curve of counter
		self.e220_c = QTextEdit(self)
		self.e220_c.setText('Particle diameter (nm) at 50 % detection efficiency (>0 nm), factor for detection efficiency dependency on particle size (>0) (determines the range of particle sizes affected by reduced detection efficiency and assumes a sigmoid function.  Defaults to 5, 0.5)')
		self.CPClayout.addWidget(self.e220_c, 0, 4)
		
		# input for maximum detectable size of particle
		self.e220_d = QTextEdit(self)
		self.e220_d.setText('Maximum detectable particle diameter (nm), e.g. 3.e5.  Defaults to -1 which implies no maximum')
		self.CPClayout.addWidget(self.e220_d, 0, 5)
		
		# input for uncertainty in total particle number concentration (%)
		self.e220_e = QTextEdit(self)
		self.e220_e.setText('Uncertainty (%) around total particle number concentration, defaults to 10')
		self.CPClayout.addWidget(self.e220_e, 0, 6)
		
		# input for response time function
		self.e220_f = QTextEdit(self)
		self.e220_f.setText('Inputs for accounting for instrument response time and any mixing of particles from different times entering inlet (five inputs in total, all separated by a comma: i) shortest delay (s) in particles reaching counting unit, ii) delay at which weighting of particles at a maximum (s), iii) function of weighting against delay time (s) (use t for time and np for numpy) for particles between the shortest delay (i) and the delay at which weighting at maximum (ii), iv) longest delay (s) in particles reaching counting unit, v) function of weighting against delay time (s)  (use t for time and np for numpy) for particles between the delay at which weighting at maximum (ii) and the longest delay (iv).  For example: 0.1, 1.3, np.exp(t), 2.5, np.exp(np.flip(t-1.3)).  Note that weighting is normalised by the integral so that the final integral is one.  Defaults to: 1., 1., 1.*t, 1., 1.*t, which represents a response time of 1 s with no mixing of particles of different ages.')
		self.CPClayout.addWidget(self.e220_f, 0, 7)
		
		# input for frequency of instrument output
		self.e220_g = QTextEdit(self)
		self.e220_g.setText('Frequency of instrument output (Hz), defaults to 1.')
		self.CPClayout.addWidget(self.e220_g, 0, 8)
		
		# input for particle loss during inlet passage
		self.e220_h = QTextEdit(self)
		self.e220_h.setText('Loss rate (fraction s<sup>-1</sup>) as a function of particle size (um) (using Dp for diameter (um), np for numpy functions and python math symbols for math functions); time of passage through inlet (s).  E.g.: np.append(10.**(-5.5-0.5*np.log10(Dp[Dp<=1.e-1])), 10.**(-4.2+0.8*np.log10(Dp[Dp>1.e-1]))); 5..  Defaults to 0., 0., which implies no particle losses in inlet')
		self.CPClayout.addWidget(self.e220_h, 0, 9)
		
		# input for averaging interval (s)
		self.e220_i = QTextEdit(self)
		self.e220_i.setText('Averaging interval (s).  Defaults to 1.')
		self.CPClayout.addWidget(self.e220_i, 0, 10)
		
		# button to plot counting efficiency dependence on particle size 
		self.b220_0 = QPushButton('Counting \nefficiency \ndependence \non size', self)
		self.b220_0.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.b220_0.setToolTip('Counting efficiency dependence on particle size')
		self.b220_0.clicked.connect(self.on_click234_a)
		self.CPClayout.addWidget(self.b220_0, 1, 4)
		
		# button to plot weighting as a function of response time
		self.b220_f = QPushButton('Weighting \ndependency \non response \ntime', self)
		self.b220_f.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.b220_f.setToolTip('Plot the weighting of particles by age due to the response time function')
		self.b220_f.clicked.connect(self.on_click222_b)
		self.CPClayout.addWidget(self.b220_f, 1, 7)
		
		# button to plot inlet loss rate as a function of particle diameter
		self.b220_h = QPushButton('Inlet loss \nrate with \nparticle \ndiameter', self)
		self.b220_h.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.b220_h.setToolTip('Plot the inlet loss rate as a function of particle diameter')
		self.b220_h.clicked.connect(self.on_click222_h)
		self.CPClayout.addWidget(self.b220_h, 1, 9)
		
		# button to plot temporal profile of total particle number concentration
		self.b220_a = QPushButton('CPC \nobservations', self)
		self.b220_a.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		self.b220_a.setToolTip('Plot the total particle number concentration as observed by a condensation particle counter')
		self.b220_a.clicked.connect(self.on_click222)
		self.CPClayout.addWidget(self.b220_a, 1, 10)
		
		
	
		return(CPCTab)
	
	def SMPStab(self): # instrument comparison plotting tab definition
	
		SMPSTab = QWidget()
		self.SMPSlayout = QGridLayout() 
		SMPSTab.setLayout(self.SMPSlayout)
	
		# input for whether to use wet or dried particles
		self.e230 = QTextEdit(self)
		self.e230.setText('0 for dried particles or 1 for not dried particles')
		self.SMPSlayout.addWidget(self.e230, 0, 0)
		
		# input for minimum particle concentration detectable
		self.e232 = QTextEdit(self)
		self.e232.setText('Minimum detectable particle concentration (particles cm<sup>-3</sup>)')
		self.SMPSlayout.addWidget(self.e232, 0, 1)
		
		# input for counting efficiency curve of counter
		self.e233 = QTextEdit(self)
		self.e233.setText('Particle diameter (nm) at 50 % counting efficiency (>0 nm), factor for counting efficiency dependency on particle size (>0) (determines the range of particle sizes affected by reduced counting efficiency and assumes a sigmoid function)')
		self.SMPSlayout.addWidget(self.e233, 0, 2)
		
		# input for minimum size particle size range of counter
		self.e234 = QTextEdit(self)
		self.e234.setText('Minimum detectable particle diameter (nm)')
		self.SMPSlayout.addWidget(self.e234, 0, 3)
		
		# input for maximum size particle size range of counter
		self.e235 = QTextEdit(self)
		self.e235.setText('Maximum detectable particle diameter (nm)')
		self.SMPSlayout.addWidget(self.e235, 0, 4)
		
		# input for number of size bins of counter
		self.e236 = QTextEdit(self)
		self.e236.setText('Number of size bins within the detectable particle diameter range (assumed to be logarithmically spaced)')
		self.SMPSlayout.addWidget(self.e236, 0, 5)
		
		# input for assumed density of particles
		self.e237 = QTextEdit(self)
		self.e237.setText('Assumed density of particles (g cm<sup>-3</sup>)')
		self.SMPSlayout.addWidget(self.e237, 0, 6)
		
		# button to plot counting efficiency dependence on particle size 
		self.b234 = QPushButton('Counting efficiency curve', self)
		self.b234.setToolTip('Counting efficiency dependence on particle size')
		self.b234.clicked.connect(self.on_click234)
		self.SMPSlayout.addWidget(self.b234, 1, 2)
		
		# button to plot temporal profile of number size distribution
		self.b233 = QPushButton('SMPS observations', self)
		self.b233.setToolTip('Plot the number size distribution as observed by a particle counter')
		self.b233.clicked.connect(self.on_click233)
		self.SMPSlayout.addWidget(self.b233, 1, 6)
		
		return(SMPSTab)
	
	def CIMStab(self): # instrument comparison plotting tab definition
	
		CIMSTab = QWidget()
		self.CIMSlayout = QGridLayout() 
		CIMSTab.setLayout(self.CIMSlayout)
	
		# input for minimum particle concentration detectable
		self.e280 = QTextEdit(self)
		self.e280.setText('Resolution of mass:charge ratio')
		self.CIMSlayout.addWidget(self.e280, 0, 0)
	
		return(CIMSTab)
	
	#@pyqtSlot() # eurochamp website under development 10/02/2021
	#def on_clickn1(self): # EUROCHAMP website
	#	import webbrowser
	#	webbrowser.open('https://www.eurochamp.org/Eurochamp2020.aspx')
	#	return()
		
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
		if (self.atb == 1): # if showing remove add to batch button
			self.b82.deleteLater()
			self.atb = 0 # remember that add to batch button not showing
		# remove any old message from previous run
		if (len(self.l81b.text()) > 0):
			if (self.l81b.text()[0:17] != 'File combinations'):
				self.l81b.setText('')
				self.l81b.setStyleSheet(0., '0px', 0., 0.)
		
		# prepare by enforcing default variables
		# default variables for all required input model variables -------------------------
		[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str] = def_mod_var.def_mod_var(0)
		
		# then open default variables, ready for modification
		input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
		
		with open(input_by_sim, 'rb') as pk:
			[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, 
			tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on,
			Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
			save_step, const_comp, Compt, injectt, Ct, seed_name,
			seed_mw, seed_diss, seed_dens, seedVr,
			light_stat, light_time, daytime, lat, lon, af_path, 
			dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, 
			dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
			accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
			nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
			inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str] = pickle.load(pk)
			pk.close()
		
		# button to get path to folder containing relevant files
		options = QFileDialog.Options()
		fol_nme = QFileDialog.getExistingDirectory(self, "Select Folder Containing Required Input Files", "./PyCHAM/input/")
		
		if (fol_nme == ''): # if no folder selected (e.g. because selection cancelled)
			return()
		
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
		list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str]
		with open(input_by_sim, 'wb') as pk:
			pickle.dump(list_vars, pk) # pickle
			pk.close() # close
			
		# read in model variables of this model variables file and store to pickle
		import mod_var_read
		mod_var_read.mod_var_read()
			
		# updating scroll labels showing path to files
		self.l3.setText(sch_name)
		self.l5.setText(xml_name)
		self.l7.setText(inname)
			
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
	
		# prepare by opening default inputs, ready for modification
		input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
		
		with open(input_by_sim, 'rb') as pk:
			[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, 
			tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on,
			Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
			save_step, const_comp, Compt, injectt, Ct, seed_name,
			seed_mw, seed_diss, seed_dens, seedVr,
			light_stat, light_time, daytime, lat, lon, af_path, 
			dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, 
			dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
			accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
			nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
			inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str] = pickle.load(pk)
			pk.close()
	
		sch_name, _ = QFileDialog.getOpenFileName(self, "Select Chemical Scheme File", "./PyCHAM/input/") # get path of file
		self.l3.clear() # clear old label
		self.l3.setText(str(sch_name))
		self.l3.show()

		# pickle with new chemical scheme name	
		list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str]
		with open(input_by_sim, 'wb') as pk:
			pickle.dump(list_vars, pk) # pickle
			pk.close() # close
	
	@pyqtSlot()
	def on_click3(self): # when different xml file requires selection
	
		# prepare by opening default inputs, ready for modification
		input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
		
		with open(input_by_sim, 'rb') as pk:
			[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, 
			tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on,
			Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
			save_step, const_comp, Compt, injectt, Ct, seed_name,
			seed_mw, seed_diss, seed_dens, seedVr,
			light_stat, light_time, daytime, lat, lon, af_path, 
			dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, 
			dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
			accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
			nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
			inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str] = pickle.load(pk)
			pk.close()
	
		xml_name, _ = QFileDialog.getOpenFileName(self, "Select xml File", "./PyCHAM/input/") # get path of file
		self.l5.clear() # clear old label
		self.l5.setText(str(xml_name))
		self.l5.show()
		
		# pickle with new chemical scheme name	
		list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str]
		with open(input_by_sim, 'wb') as pk:
			pickle.dump(list_vars, pk) # pickle
			pk.close() # close
			
	@pyqtSlot()
	def on_click4(self): # when different model variables file requires selection
	

		if (self.fab == 1): # if showing, remove single simulation widgets
			self.b81.deleteLater()
			self.l81.deleteLater()
			self.fab = 0 # remember that single simulation widgets not showing
		if (self.atb == 1): # if showing remove add to batch button
			self.b82.deleteLater()
			self.atb = 0 # remember that add to batch button not showing
		# remove any old 'Simulation complete' message from previous run
		if ((self.l81b.text() == 'Simulation complete') or (self.l81b.text() == 'Simulations complete')):
			self.l81b.setText('')

		# prepare by opening existing inputs, ready for modification
		input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
		with open(input_by_sim, 'rb') as pk:
			[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, 
			tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on,
			Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
			save_step, const_comp, Compt, injectt, Ct, seed_name,
			seed_mw, seed_diss, seed_dens, seedVr,
			light_stat, light_time, daytime, lat, lon, af_path, 
			dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, 
			dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
			accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
			nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
			inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str] = pickle.load(pk)
			pk.close() # close pickle file

		# user chooses path of file to model variables file
		inname, _ = QFileDialog.getOpenFileName(self, "Select Model Variables File", "./PyCHAM/input/")
		
		if (inname == ''): # if no file selected, e.g. because selection was cancelled
			return()
		
		self.l7.clear() # clear old label
		self.l7.setText(inname)
		self.l7.show()

		# pickle with new file name	
		list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str]
		with open(input_by_sim, 'wb') as pk:
			pickle.dump(list_vars, pk) # pickle
			pk.close() # close
	
		# read in model variables of this model variables file and store to pickle
		import mod_var_read
		mod_var_read.mod_var_read()
	
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
			[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str] = def_mod_var.def_mod_var(0)
		
			
			
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
			sch_name = txtn[0]
			xml_name = txtn[1]
			inname = txtn[2]
				
			# pickle with new file names	
			list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str]
				
			# path to pickle file
			input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
			with open(input_by_sim, 'wb') as pk:
				pickle.dump(list_vars, pk) # pickle
				pk.close() # close
			
			# read in model variables of this model variables file 
			# (as identified by the inname variable of the pickle file) and store to pickle
			import mod_var_read
			mod_var_read.mod_var_read()
			
			# get the save path name variables
			input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
			with open(input_by_sim, 'rb') as pk:
				[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, 
				tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on,
				Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, 
				save_step, const_comp, Compt, injectt, Ct, seed_name,
				seed_mw, seed_diss, seed_dens, seedVr,
				light_stat, light_time, daytime, lat, lon, af_path, 
				dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, 
				dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
				accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, 
				nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, 
				inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, wat_hist, drh_str, erh_str] = pickle.load(pk)
				pk.close() # close pickle file
			
			# saving path - copied from saving module
			dir_path = os.getcwd() # current working directory
			output_root = 'PyCHAM/output'
			filename = os.path.basename(sch_name)
			filename = os.path.splitext(filename)[0]
			output_by_sim = os.path.join(dir_path, output_root, filename, sav_nam)
			
			# display the progress label, including the name of the saving path
			if (self.btch_no == 1): # single run mode
				self.l81b.setText('')
				self.l81b.setText(str('Progress through simulation saving to: \n' + str(output_by_sim) + '\n' + str(sim_num+1) + ' of ' + str(self.btch_no)))
			if (self.btch_no > 1): # batch run mode
				self.l81b.setText('')
				self.l81b.setText(str('Progress through simulation saving to: \n' + str(output_by_sim) + '\n' + str(sim_num+1) + ' of ' + str(self.btch_no-1)))
			
			# disable then remove start simulation button
			self.b81.setEnabled(False)
			self.b81.deleteLater()
			
			# if showing, disable then remove add to batch button
			if (self.atb == 1): 
				self.b82.setEnabled(False)
				self.b82.deleteLater()
				# remove or label
				self.l81.deleteLater()
				self.atb = 0
			
			# path to error log
			err_log = str(os.getcwd() + '/PyCHAM/err_log.txt')
			if (sim_num == 0): # # delete any existing error log and create new log
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
	def act_81(self, output_by_sim, sim_num): # action the simulation
	
		from middle import middle # prepare to communicate with main programme
		
		note_messf = 0 # cancel note message flag
		
		for prog in middle(): # call on modules to solve problem
			
		
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
						self.l81b.setText(str('Progress through simulation saving to: \n' + str(output_by_sim) + '\n' + str(sim_num+1) + ' of ' + str(self.btch_no)))
					if (self.btch_no > 1): # batch run mode
						self.l81b.setText('')
						self.l81b.setText(str('Progress through simulation saving to: \n' + str(output_by_sim) + '\n' + str(sim_num+1) + ' of ' + str(self.btch_no-1)))
			
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
			self.l81.deleteLater()
			self.fab = 0 # remember that single simulation widgets not showing
		if (self.atb == 1): # if showing, disable and remove add to batch button
			self.b82.setEnabled(False)
			self.b82.deleteLater()
			self.atb = 0 # remember that add to batch button not showing
		if (self.ssb == 1):# if showing, remove series of simulations button
			self.b83.deleteLater()
			self.ssb = 0 # remember this button removed
		
		self.l81b.setText('') # clear progress message
		
		# action to simulate
		err_mess = self.on_click81b()
		
		# once all simulations done in single simulation mode, tidy up
		# clear the list of file combinations for simulations in batch
		self.output_list = [] # reset list of output paths
		# return to single simulation mode
		self.btch_no = 1
		self.btch_str = 'File combinations included in batch: chemical scheme, xml, model variables\n'

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
		
		# action to simulate
		err_mess = self.on_click81b()

	
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
		sys.exit() # end programme and release all memory
		
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
		
		dir_path = self.l201.text() # name folder with results
		plotter.plotter(0, dir_path, self) # plot results
	
	@pyqtSlot() # button to plot gas-phase concentration temporal profile
	def on_click206(self):	
		
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		# get names of components to plot
		comp_names = [str(i) for i in self.e205.text(). split(',')]
		
		import plotter_gp
		dir_path = self.l201.text() # name folder with results
		plotter_gp.plotter(0, dir_path, comp_names, self) # plot results
	
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
		comp_names = [str(i) for i in self.e217.text(). split(',')]
		
		import plotter_ct
		dir_path = self.l201.text() # name of folder with results
		plotter_ct.plotter(0, dir_path, comp_names, self) # plot results
	
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
				self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
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
			self.l203a.setText('Note - Loss rate (fraction s<sup>-1</sup>) as a function of particle size (um) and time of passage through inlet (s), should be a string (using Dp for diameter (um), np for numpy functions and python math symbols for math functions) followed by a  number, with the two separated by a semi-colon, e.g.: np.append(10.**(-5.5-0.5*np.log10(Dp[Dp<=1.e-1])), 10.**(-4.2+0.8*np.log10(Dp[Dp>1.e-1]))); 5..  Defaulting to 0., 0. which implies no particle losses in inlet.')
			
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
	
	@pyqtSlot() # button to plot number size distribution replication of particle sizer and counter
	def on_click233(self):
	
		# reset error message
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		try: 
			dryf = int(self.e230.toPlainText()) # whether dry or not
			
		except: # give error message
			self.l203a.setText('Error - whether particles are dried (0) or not (1) should be a single integer')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1

			return()
		
		# false background counts ----------------------------------------------------------------
		try:
			cdt = float(self.e232.toPlainText()) # concentration detection limit (particles/cm3)
			
		except: # give error message
			self.l203a.setText('Error - particle number concentration detection limit of counter should be a single number')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1

			return()
		# ----------------------------------------------------------------------------------------------
		
		# diameter at 50 % detection efficiency and width factor for sigmoidal 
		# curve of detection efficiency against size ---------------------------------
		try:
			# get particle diameter at 50 % counting efficiency and width 
			# factor for counting efficiency dependence on particle size
			sdt = ((self.e233.toPlainText()).split(','))
			sdt = [float(i) for i in sdt] 
			
		except: # give error message
			self.l203a.setText('Error - particle diameter (nm) at 50 % detection efficiency and factor for detection efficiency dependence on particle size should be two numbers separated by a comma, e.g.: 5, 1')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()
		
		if (len(sdt) != 2): # if not two numbers show error
			self.l203a.setText('Error - particle diameter (nm) at 50 % detection efficiency and width factor for detection efficiency dependence on particle size should be two numbers separated by a comma, e.g.: 5, 1')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()
		
		# -----------------------------------------------------------------------------------------
		
		try:
			#  minimum particle size of counter
			min_size = float((self.e234.toPlainText()))
			
		except: # give error message
			self.l203a.setText('Error - minimum detectable particle diameter (nm) should be a single number, e.g. 5')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()
		
		# maximum particle size -------------------------------------------------------
		try:
			#  maximum particle size of counter
			max_size = float((self.e235.toPlainText()))
			
		except: # give error message
			self.l203a.setText('Error - maximum size of particle size range (nm) should be a single number, e.g. 500')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()
		# -------------------------------------------------------------------------------------
		
		try:
			#  number of size bins of counter
			csbn = int((self.e236.toPlainText()))
			
		except: # give error message
			self.l203a.setText('Error - number of size bins within detectable particle size range (nm) should be a single number, e.g. 128')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()
		
		try:
			#  assumed density of particles (g/cm3)
			p_rho = float((self.e237.toPlainText()))
			
		except: # give error message
			self.l203a.setText('Error - assumed density of particles (g/cm3) should be a single number, e.g. 1.3')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()
		
		import plotter_counters
		importlib.reload(plotter_counters) # ensure latest version uploaded
		dir_path = self.l201.text() # name of folder with results
		plotter_counters.plotter(0, dir_path, self, dryf, cdt, sdt, min_size, max_size, csbn, p_rho) # plot results
		
		return()

	# button to plot detection efficiency as a function of particle 
	# diameter for number size distribution counters
	@pyqtSlot()
	def on_click234(self):
		
		# reset message box
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		try:	
			# get particle diameter at 50 % counting efficiency and width 
			# factor for counting efficiency dependence on particle size
			sdt = ((self.e233.toPlainText()).split(','))
			sdt = [float(i) for i in sdt] 
			
		except: # give error message
			self.l203a.setText('Error - particle diameter (nm) at 50 % detection efficiency and width factor for detection efficiency dependence on particle size should be two numbers separated by a comma, e.g.: 5, 1')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()
		
		if (len(sdt) != 2): # if not two numbers show error
			self.l203a.setText('Error - particle diameter (nm) at 50 % detection efficiency and width factor for detection efficiency dependence on particle size should be two numbers separated by a comma, e.g.: 5, 1')
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()
						
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
	
	# button to plot particle weighting by age due to instrument response time
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

	# button to plot loss rate of particles during inlet passage
	@pyqtSlot()
	def on_click222_h(self):

		import inlet_loss
		
		# obtain the relevant inputs
		# reset message box
		self.l203a.setStyleSheet(0., '0px dashed red', 0., 0.)
		self.l203a.setText('')
		
		# particle diameter at 50 % detection efficiency -----------------------------------
		
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
		
		# ------------------------------------------------------------------------------------------------
		
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
			self.l203a.setText('Note - Loss rate (fraction s<sup>-1</sup>) as a function of particle size (um) and time of passage through inlet (s), should be a string (using Dp for diameter (um), np for numpy functions and python math symbols for math functions) followed by a  number, with the two separated by a semi-colon, e.g.: np.append(10.**(-5.5-0.5*np.log10(Dp[Dp<=1.e-1])), 10.**(-4.2+0.8*np.log10(Dp[Dp>1.e-1]))); 5..  Defaulting to 0., 0. which implies no particle losses in inlet.')
			
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
		
		inlet_loss.inlet_loss(3, [], xn, [], loss_func_str, [], 0)

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