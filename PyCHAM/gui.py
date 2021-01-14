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
		tabs.addTab(self.NStab(), "Simulation Setup")
		tabs.addTab(self.PLtab(), "Plotting")
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
		[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O] = def_mod_var.def_mod_var(0)
		
		# listing input files -----------------------------------------------------------
		l1 = QLabel(self)
		l1.setText("The following files have been found: ")
		self.NSlayout.addWidget(l1, 1, ffscn, 1, 2)
		
		l2 = QLabel(self)
		l2.setText("Chemical scheme: ")
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
		l6.setText("Model variables: ")
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
		self.NSlayout.addWidget(l8, 0, self.mvpn, 1, 2)
		
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
		
		
		# Chamber Environment ----------
		
		l14a = QLabel(self)
		l14a.setText('Chamber Environment')
		l14a.setFont(QFont("Arial", 13, QFont.Bold))
		env_row = gen_row+6
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
		self.l16a.setText((str(RH)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.l16a, env_row+3, 1)
		
		l17 = QLabel(self)
		l17.setText('Chamber air pressure (Pa): ')
		self.varbox.addWidget(l17, env_row+4, 0)
		self.l17a = QLabel(self)
		self.l17a.setText((str(Press)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.l17a, env_row+4, 1)
		
		# particle properties ----------------
		
		l18a = QLabel(self)
		l18a.setText('Particle Properties')
		l18a.setFont(QFont("Arial", 13, QFont.Bold))
		par_row = env_row+5
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
		
		# gas inputs ----------------
		
		l39a = QLabel(self)
		l39a.setText('Gas')
		l39a.setFont(QFont("Arial", 13, QFont.Bold))
		gas_row = par_row+22
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

		
		
		# properties of model variables scroll area ----------------
		self.scrollwidget.setLayout(self.varbox)
		self.scroll.setWidget(self.scrollwidget)
		self.NSlayout.addWidget(self.scroll, 1, self.mvpn, 3, 2)
		
		# label to let user know preparedness of simulation - displayed text updated in ui_check module below
		self.l80 = ScrollLabel(self)
		self.NSlayout.addWidget(self.l80, 4, self.mvpn, 1, 2)
		self.bd_st = 2 # border status
		
		# label to let users know file combinations included in batch
		self.btch_str = 'File combinations included in batch: chemical scheme, xml, model variables\n'
		# begin count on number of simulations in batch
		self.btch_no = 1
		
		# running check on default model variables ---------------------------------------------------------------
		import ui_check # module for checking on model variables
		# check on inputs - note this loads the last saved pickle file and saves any change to this pickle file
		ui_check.ui_check(self)
		# finished check on model variables -------------------------------------------------------------------------
		
		# relative stretching (width-wise) of each column in Simulation Setup tab
		self.NSlayout.setColumnStretch(0, 0.7)
		self.NSlayout.setColumnStretch(1, 1)
		self.NSlayout.setColumnStretch(2, 1)
		self.NSlayout.setColumnStretch(3, 1)
		self.NSlayout.setColumnStretch(4, 1)

		return(NSTab)
	
	def PLtab(self): # Plotting tab
		PLTab = QWidget()
		layout = QGridLayout() 
		PLTab.setLayout(layout)
		
		# label for name of results folder to plot:
		l81a = QLabel(self)
		l81a.setText('Path to results folder to be used for plotting: ')
		l81a.setWordWrap(True)
		layout.addWidget(l81a, 0, 0)

		return(PLTab)
	
	@pyqtSlot()
	def on_clickn1(self): # EUROCHAMP website
		import webbrowser
		webbrowser.open('https://www.eurochamp.org/Eurochamp2020.aspx')
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
		
		# prepare by opening existing variables, ready for modification
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
		list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O]
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
		self.l3.setText(str(sch_name))
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
		self.l5.clear() # clear old label
		self.l5.setText(str(xml_name))
		self.l5.show()
		
		# pickle with new chemical scheme name	
		list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O]
		with open(input_by_sim, 'wb') as pk:
			pickle.dump(list_vars, pk) # pickle
			pk.close() # close
			
	@pyqtSlot()
	def on_click4(self): # when different model variables file requires selection
	

		# prepare by opening existing inputs, ready for modification
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

		# user chooses path of file to model variables file
		inname, _ = QFileDialog.getOpenFileName(self, "Select Model Variables File", "./PyCHAM/input/")
		
		if (inname == ''): # if no file selected, e.g. because selection was cancelled
			return()
		
		self.l7.clear() # clear old label
		self.l7.setText(inname)
		self.l7.show()

		# pickle with new file name	
		list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O]
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
	def on_click81(self): # when button to run simulation pressed
	
		for sim_num in range(self.btch_no):
		
			# if in batch mode, prepare the relevant inputs
			if (self.btch_no > 1):
			
				if (sim_num == self.btch_no-1): # reach end of list
					return()
					
				# reset to default variables to allow any new variables to arise
				# from the current model variables file only
				[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O] = def_mod_var.def_mod_var(0)
			
			
				# get text from batch list label
				btch_list = self.l82.text()
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
				list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O]
				# path to pickle file
				input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
				with open(input_by_sim, 'wb') as pk:
					pickle.dump(list_vars, pk) # pickle
					pk.close() # close
	
				# read in model variables of this model variables file and store to pickle
				import mod_var_read
				mod_var_read.mod_var_read()
			
							
			import middle # prepare to communicate with main programme
			middle.middle() # call on modules to solve problem
		
		# running check on model variables - note this useful here as will stop the user running
		# exactly the same simulation again -------------------------------------------------------------
		import ui_check # module for checking on model variables
		# check on inputs - note this loads the last saved pickle file and saves any change to this pickle file
		ui_check.ui_check(self)
		# finished check on model variables --------------------------------------------------------------
	
		
	@pyqtSlot()
	def on_click82(self): # when button to add to batch pressed
	
		# replace 'Start Simulation' button with 'Start Batch Processing' button
		self.b81.deleteLater()
		self.b81 = QPushButton('Start Batch Processing', self)
		self.b81.setToolTip('Start the series of simulations')
		self.b81.clicked.connect(self.on_click81)
		self.NSlayout.addWidget(self.b81, 7, self.mvpn, 1, 2)
	
		# label to let user know simulation setups included in the batch
		self.l82 = ScrollLabel(self)
		self.NSlayout.addWidget(self.l82, 6, self.mvpn, 1, 2)
		
		cs_file = self.l3.text()
		xml_file = self.l5.text()
		mv_file = self.l7.text()
		
		self.btch_str = (str(self.btch_str + str(self.btch_no)+'.\n' + cs_file + '\n' + xml_file + '\n' + mv_file + '\n'))
		self.l82.setText(self.btch_str)
		
		self.btch_no += 1 # increase number in batch
		
		

	@pyqtSlot() # button to quit software
	def on_click89(self):
		QWidget.close(self)

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