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
	
	def text(self): # interpreting text within the label method
		label_text = self.label.text()
		return(label_text)

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
		grid.addWidget(l0, 0, 0, 1, 3)
		
		# EUROCHAMP logo
		l0b= QLabel(self)
		pixmap = QPixmap('PyCHAM/logos_eurochamp2020-orange-noir.png')
		l0b.setPixmap(pixmap.scaled(100, 100, transformMode = Qt.SmoothTransformation))
		EUROlogo_hindx = 4 # horizontal position on grid
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
		NCASlogo_hindx = 6 # horizontal position on grid
		grid.addWidget(l0c, 0, NCASlogo_hindx)
		# link to EUROCHAMP website
		bn1a = QPushButton('', self)
		bn1a.setStyleSheet('background : transparent')
		bn1a.setToolTip('Visit the NCAS website')
		bn1a.clicked.connect(self.on_clickn1a)
		grid.addWidget(bn1a, 0, NCASlogo_hindx)
		
		
		# horizontal line
		self.separatorLine = QFrame()
		self.separatorLine.setFrameShape(QFrame.HLine)
		self.separatorLine.setFrameShadow(QFrame.Raised)
		grid.addWidget(self.separatorLine, 1, 0, 1, 7)
		self.separatorLine.show()
		
		# vertical lines
		self.separatorLine = QFrame()
		self.separatorLine.setFrameShape(QFrame.VLine)
		self.separatorLine.setFrameShadow(QFrame.Raised)
		vline_hindx = 1 # starting grid index for vertical lines
		grid.addWidget(self.separatorLine, 1, vline_hindx, 12, 1)
		self.separatorLine.show()

		self.separatorLine = QFrame()
		self.separatorLine.setFrameShape(QFrame.VLine)
		self.separatorLine.setFrameShadow(QFrame.Raised)
		grid.addWidget(self.separatorLine, 1, vline_hindx+2, 12, 1)
		self.separatorLine.show()
		
		self.separatorLine = QFrame()
		self.separatorLine.setFrameShape(QFrame.VLine)
		self.separatorLine.setFrameShadow(QFrame.Raised)
		grid.addWidget(self.separatorLine, 1, vline_hindx+4, 12, 1)
		self.separatorLine.show()
		
		# default variables for all required input model variables -------------------------
		# stored to pickle file and output here
		[sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O] = def_mod_var.def_mod_var(0)
		
		
		# file selection pane -------------------------------------------------------------------------
		b0 = QPushButton('Folder Containing Input Files', self)
		b0.setToolTip('Select the folder containing the required input files')
		b0.clicked.connect(self.on_click1)
		grid.addWidget(b0, 2, 0)
		
		l1 = QLabel(self)
		l1.setText("The following files have been found: ")
		grid.addWidget(l1, 3, 0)
		
		l2 = QLabel(self)
		l2.setText("Chemical scheme file: ")
		grid.addWidget(l2, 4, 0)
			
		self.l3 = ScrollLabel(self)
		self.l3.setText(sch_name)
		grid.addWidget(self.l3, 5, 0)
			
		b1 = QPushButton('Select Different File', self)
		b1.setToolTip('Select the file containing the required chemical scheme file')
		b1.clicked.connect(self.on_click2)
		grid.addWidget(b1, 6, 0)
		
		l4 = QLabel(self)
		l4.setText("xml file: ")
		grid.addWidget(l4, 7, 0)
			
		self.l5 = ScrollLabel(self)
		self.l5.setText(xml_name)
		grid.addWidget(self.l5, 8, 0)
			
		b2 = QPushButton('Select new file', self)
		b2.setToolTip('Select the file containing the required xml file')
		b2.clicked.connect(self.on_click3)
		grid.addWidget(b2, 9, 0)
		
		l6 = QLabel(self)
		l6.setText("Model variables file: ")
		grid.addWidget(l6, 10, 0)
			
		self.l7 = ScrollLabel(self)
		self.l7.setText(inname)
		grid.addWidget(self.l7, 11, 0)
			
		b3 = QPushButton('Select new file', self)
		b3.setToolTip('Select the file containing the required model variables file')
		b3.clicked.connect(self.on_click4)
		grid.addWidget(b3, 12, 0)
		
		# model variables pane -------------------------------------------------------------------------
		mvpn = 2 # column number for pyqt grid
		
		
		# button to load model variables from selected Model Variables file
		b4 = QPushButton(r'Load Model Variables From Selected File', self)
		b4.setToolTip('Match model variables to those in the selected model variables file')
		b4.clicked.connect(self.on_click5)
		grid.addWidget(b4, 2, mvpn)
		
		l8 = QLabel(self)
		l8.setText("The following model variables have been found: ")
		grid.addWidget(l8, 3, mvpn)
		
		# begin a scroll area to contain model variables
		self.scroll = QScrollArea()
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
		self.e9 = QLineEdit(self)
		self.e9.setText(sav_nam)
		self.varbox.addWidget(self.e9, gen_row+2, 0)
		
		l10 = QLabel(self)
		l10.setText('Chemical scheme markers: ')
		self.varbox.addWidget(l10, gen_row+3, 0)
		self.e10 = QLineEdit(self)
		self.e10.setText((str(chem_sch_mark)).replace('\'', '').replace(' ', '')[1:-1])
		self.varbox.addWidget(self.e10, gen_row+4, 0)
		
		l11 = QLabel(self)
		l11.setText('Total experiment time (s): ')
		self.varbox.addWidget(l11, gen_row+5, 0)
		self.e11 = QLineEdit(self)
		self.e11.setText((str(tot_time)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.e11, gen_row+6, 0)
		
		l12 = QLabel(self)
		l12.setText('Update time interval (s): ')
		self.varbox.addWidget(l12, gen_row+7, 0)
		self.e12 = QLineEdit(self)
		self.e12.setText((str(update_stp)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.e12, gen_row+8, 0)
		
		l13 = QLabel(self)
		l13.setText('Recording time interval (s): ')
		self.varbox.addWidget(l13, gen_row+9, 0)
		self.e13 = QLineEdit(self)
		self.e13.setText((str(save_step)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.e13, gen_row+10, 0)
		
		# Chamber Environment ----------
		
		l14a = QLabel(self)
		l14a.setText('Chamber Environment')
		l14a.setFont(QFont("Arial", 13, QFont.Bold))
		env_row = gen_row+11
		self.varbox.addWidget(l14a, env_row+0, 0)
		
		
		l14 = QLabel(self)
		l14.setText('Chamber temperature(s) (K): ')
		self.varbox.addWidget(l14, env_row+1, 0)
		self.e14 = QLineEdit(self)
		self.e14.setText((str(temp)).replace('\'', '').replace(' ', '')[1:-1])
		self.varbox.addWidget(self.e14, env_row+2, 0)
		
		l15 = QLabel(self)
		l15.setText('Time(s) for chamber temperature(s) (s): ')
		self.varbox.addWidget(l15, env_row+3, 0)
		self.e15 = QLineEdit(self)
		self.e15.setText((str(tempt)).replace('\'', '').replace(' ', '')[1:-1])
		self.varbox.addWidget(self.e15, env_row+4, 0)
		
		l16 = QLabel(self)
		l16.setText('Relative humidity (0-1): ')
		self.varbox.addWidget(l16, env_row+5, 0)
		self.e16 = QLineEdit(self)
		self.e16.setText((str(RH)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.e16, env_row+6, 0)
		
		l17 = QLabel(self)
		l17.setText('Chamber air pressure (Pa): ')
		self.varbox.addWidget(l17, env_row+7, 0)
		self.e17 = QLineEdit(self)
		self.e17.setText((str(Press)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.e17, env_row+8, 0)
		
		# Particle properties ----------------
		
		l18a = QLabel(self)
		l18a.setText('Particle Properties')
		l18a.setFont(QFont("Arial", 13, QFont.Bold))
		par_row = env_row+9
		self.varbox.addWidget(l18a, par_row+0, 0)
		
		
		l18 = QLabel(self)
		l18.setText('Size structure (MC=0, FM=1): ')
		self.varbox.addWidget(l18, par_row+1, 0)
		self.e18 = QLineEdit(self)
		self.e18.setText((str(siz_stru)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.e18, par_row+2, 0)
		
		l19 = QLabel(self)
		l19.setText('Number of particle size bins: ')
		self.varbox.addWidget(l19, par_row+3, 0)
		self.e19 = QLineEdit(self)
		self.e19.setText((str(num_sb)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.e19, par_row+4, 0)
		
		l20 = QLabel(self)
		l20.setText('Particle concentrations (#/cc): ')
		self.varbox.addWidget(l20, par_row+5, 0)
		self.e20 = QLineEdit(self)
		self.e20.setText((str(pconc)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e20, par_row+6, 0)
		
		l21 = QLabel(self)
		l21.setText('Particle concentration times (s): ')
		self.varbox.addWidget(l21, par_row+7, 0)
		self.e21 = QLineEdit(self)
		self.e21.setText((str(pconct)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e21, par_row+8, 0)
		
		l22 = QLabel(self)
		l22.setText('Molecular weight of seed particle \ncomponent (g/mol): ')
		self.varbox.addWidget(l22, par_row+9, 0)
		self.e22 = QLineEdit(self)
		self.e22.setText((str(seed_mw)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e22, par_row+10, 0)
		
		l23 = QLabel(self)
		l23.setText('Dissociation constant(s) of seed \ncomponent(s): ')
		self.varbox.addWidget(l23, par_row+11, 0)
		self.e23 = QLineEdit(self)
		self.e23.setText((str(seed_diss)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e23, par_row+12, 0)
		
		l24 = QLabel(self)
		l24.setText('Density of seed particles (g/cc): ')
		self.varbox.addWidget(l24, par_row+13, 0)
		self.e24 = QLineEdit(self)
		self.e24.setText((str(seed_dens)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e24, par_row+14, 0)
		
		l25 = QLabel(self)
		l25.setText('Name of seed component: ')
		self.varbox.addWidget(l25, par_row+15, 0)
		self.e25 = QLineEdit(self)
		self.e25.setText((str(seed_name)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e25, par_row+16, 0)
		
		l26 = QLabel(self)
		l26.setText('Volume ratio of components in seed \nparticles: ')
		self.varbox.addWidget(l26, par_row+17, 0)
		self.e26 = QLineEdit(self)
		self.e26.setText((str(seedVr)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e26, par_row+18, 0)
		
		l27 = QLabel(self)
		l27.setText('Smallest radius size bin \nboundary (um): ')
		self.varbox.addWidget(l27, par_row+19, 0)
		self.e27 = QLineEdit(self)
		self.e27.setText((str(lowsize)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e27, par_row+20, 0)
		
		l28 = QLabel(self)
		l28.setText('Largest radius size bin \nboundary (um): ')
		self.varbox.addWidget(l28, par_row+21, 0)
		self.e28 = QLineEdit(self)
		self.e28.setText((str(uppsize)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e28, par_row+22, 0)
		
		l29 = QLabel(self)
		l29.setText('Method for spacing size bins (lin) \nor (log): ')
		self.varbox.addWidget(l29, par_row+23, 0)
		self.e29 = QLineEdit(self)
		self.e29.setText((str(space_mode)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e29, par_row+24, 0)
		
		l30 = QLabel(self)
		l30.setText('Standard deviation for \nparticle number \nsize distribution: ')
		self.varbox.addWidget(l30, par_row+25, 0)
		self.e30 = QLineEdit(self)
		self.e30.setText((str(std)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e30, par_row+26, 0)
		
		l31 = QLabel(self)
		l31.setText('Mean radius (um) for \nparticle number \nsize distribution: ')
		self.varbox.addWidget(l31, par_row+27, 0)
		self.e31 = QLineEdit(self)
		self.e31.setText((str(mean_rad)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e31, par_row+28, 0)
		
		l32 = QLabel(self)
		l32.setText('Radius (um) of newly nucleated \nparticles: ')
		self.varbox.addWidget(l32, par_row+29, 0)
		self.e32 = QLineEdit(self)
		self.e32.setText((str(new_partr)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e32, par_row+30, 0)
		
		l33 = QLabel(self)
		l33.setText('First nucleation \nparameterisation parameter: ')
		self.varbox.addWidget(l33, par_row+31, 0)
		self.e33 = QLineEdit(self)
		self.e33.setText((str(nucv1)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e33, par_row+32, 0)
		
		l34 = QLabel(self)
		l34.setText('Second nucleation \nparameterisation parameter: ')
		self.varbox.addWidget(l34, par_row+33, 0)
		self.e34 = QLineEdit(self)
		self.e34.setText((str(nucv2)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e34, par_row+34, 0)
	
		l35 = QLabel(self)
		l35.setText('Third nucleation \nparameterisation parameter: ')
		self.varbox.addWidget(l35, par_row+35, 0)
		self.e35 = QLineEdit(self)
		self.e35.setText((str(nucv3)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e35, par_row+36, 0)
		
		l36 = QLabel(self)
		l36.setText('Chemical scheme name of \nnucleating component: ')
		self.varbox.addWidget(l36, par_row+37, 0)
		self.e36 = QLineEdit(self)
		self.e36.setText((str(nuc_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e36, par_row+38, 0)
		
		l37 = QLabel(self)
		l37.setText('Whether (1) or not (0) to adapt \nintegration time interval to \nnucleation: ')
		self.varbox.addWidget(l37, par_row+39, 0)
		self.e37 = QLineEdit(self)
		self.e37.setText((str(nuc_ad)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e37, par_row+, 0)
		
		l38 = QLabel(self)
		l38.setText('Whether (1) or not (0) to serialise \ngas-particle partitioning of water: ')
		self.varbox.addWidget(l38, par_row+41, 0)
		self.e38 = QLineEdit(self)
		self.e38.setText((str(ser_H2O)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e38, par_row+42, 0)

		
		# properties of model variables scroll area
		self.scrollwidget.setLayout(self.varbox)
		self.scroll.setWidget(self.scrollwidget)
		grid.addWidget(self.scroll, 4, mvpn, 6, 1)
		
		# button to update model variables to those shown in GUI boxes
		b79 = QPushButton('Update Model Variables \n To Those Edited In Text Boxes', self)
		b79.setToolTip('In case model variables in the text boxes above have been edited, update here')
		b79.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		b79.clicked.connect(self.on_click79)
		grid.addWidget(b79, 10, mvpn)
		
		# button to check model variables from selected Model Variables file
		b80 = QPushButton(r'Check Model Variables', self)
		b80.setToolTip('Run a check to ensure model variables are valid')
		b80.clicked.connect(self.on_click80)
		grid.addWidget(b80, 11, mvpn)
		
		
		b81 = QPushButton(r'Run Model', self)
		b81.setToolTip('Start the simulation')
		b81.clicked.connect(self.on_click81)
		grid.addWidget(b81, 12, mvpn)
		
		# Action pane ---------------------------------------------------------------------
		apn = 4 # column number for qt grid
		
		# label for name of results folder to plot:
		l81a = QLabel(self)
		l81a.setText('Path to results folder to be used for plotting: ')
		l81a.setWordWrap(True)
		grid.addWidget(l81a, 2, apn)
		
		self.l81b = ScrollLabel(self)
		cwd = os.getcwd() # current working directory
		path = str(cwd + '/PyCHAM/output/ex_chem_scheme/Example_Run_output')
		self.l81b.setText(path)
		grid.addWidget(self.l81b, 3, apn)
			
		b81c = QPushButton('Select new folder', self)
		b81c.setToolTip('Select the folder containing the result files to plot')
		b81c.clicked.connect(self.on_click81c)
		grid.addWidget(b81c, 4, apn)

		b82 = QPushButton('Create Standard Results Plot', self)
		b82.setToolTip('Create standard results plot ')
		b82.clicked.connect(self.on_click82)
		grid.addWidget(b82, 5, apn)
		
		# gas-phase concentrations temporal profiles -------------
		
		# label for names of components to plot temporal profile of gas-phase for
		l83 = QLabel(self)
		l83.setText('Chemical scheme name of component(s) for plotting gas-phase concentration temporal profile: ')
		l83.setWordWrap(True)
		grid.addWidget(l83, 6, apn)
		
		# input bar for names of components to plot temporal profiles of gas-phase concentration
		self.e83 = QLineEdit(self)
		grid.addWidget(self.e83, 7, apn)
		
		# button to plot temporal profile of gas-phase concentrations
		b83 = QPushButton('Plot temporal profile of \ngas-phase concentrations', self)
		b83.setToolTip('Plot gas-phase concentration temporal profile')
		b83.clicked.connect(self.on_click83)
		b83.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		grid.addWidget(b83, 8, apn)
		
		# particle-phase concentrations temporal profiles -------------
		
		# label for names of components to plot temporal profile of total particle-phase concentration for
		l84 = QLabel(self)
		l84.setText('Chemical scheme name of component(s) for plotting total particle-phase concentration temporal profile: ')
		l84.setWordWrap(True)
		grid.addWidget(l84, 9, apn)
		
		# input bar for names of components to plot temporal profiles of gas-phase concentration
		self.e84 = QLineEdit(self)
		grid.addWidget(self.e84, 10, apn)
		
		# button to plot temporal profile of total particle-phase concentrations
		b84 = QPushButton('Plot temporal profile of total \nparticle-phase concentrations', self)
		b84.setToolTip('Plot particle-phase concentration temporal profile')
		b84.clicked.connect(self.on_click84)
		b84.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		grid.addWidget(b84, 11, apn)
		
		# wall (from gas-wall partitioning) concentrations temporal profiles -------------
		
		# label for names of components to plot temporal profile of total particle-phase concentration for
		l85 = QLabel(self)
		l85.setText('Chemical scheme name of component(s) \n for plotting wall concentration \n (from gas-wall partitioning) \n temporal profile: ')
		#l85.setWordWrap(True)
		grid.addWidget(l85, 2, apn+2)
		
		# input bar for names of components to plot temporal profiles of gas-phase concentration
		self.e85 = QLineEdit(self)
		grid.addWidget(self.e85, 3, apn+2)
		
		# button to plot temporal profile of total particle-phase concentrations
		b85 = QPushButton('Plot temporal profile of \nwall concentrations \n(from gas-wall partitioning)', self)
		b85.setToolTip('Plot wall concentration temporal profile')
		b85.clicked.connect(self.on_click85)
		b85.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		grid.addWidget(b85, 4, apn+2)
		
		# wall (from particle deposition to wall) concentrations temporal profiles -------------
		
		# label for names of components to plot temporal profile for
		l86 = QLabel(self)
		l86.setText('Chemical scheme name of component(s) \n for plotting wall concentration \n (from particle deposition to wall) \n temporal profile: ')
		grid.addWidget(l86, 5, apn+2)
		
		# input bar for names of components to plot temporal profiles of gas-phase concentration
		self.e86 = QLineEdit(self)
		grid.addWidget(self.e86, 6, apn+2)
		
		# button to plot temporal profile of total particle-phase concentrations
		b86 = QPushButton('Plot temporal profile of \nwall concentrations \n(from particle deposition to wall)', self)
		b86.setToolTip('Plot wall concentration temporal profile')
		b86.clicked.connect(self.on_click86)
		b86.setStyleSheet('background-color : white; border-width : 1px; border-radius : 7px; border-color: silver; padding: 2px; border-style : solid')
		grid.addWidget(b86, 7, apn+2)
		
		# Quit pane ---------------------------------------------------------------------------------
		qpn = apn+2 # pane number
		
		b89 = QPushButton('Quit', self)
		b89.setToolTip('Finish with PyCHAM and close this window')
		b89.clicked.connect(self.on_click89)
		grid.addWidget(b89, 12, qpn)
		

		self.show()
		return
	
	@pyqtSlot()
	def on_clickn1a(self):
		
		import webbrowser
		webbrowser.open('https://ncas.ac.uk')
		return()
	
	@pyqtSlot()
	def on_clickn1(self):
		
		import webbrowser
		webbrowser.open('https://www.eurochamp.org/Eurochamp2020.aspx')
		return()
	
	@pyqtSlot()
	def on_click1(self): # selecting folder containing input files
	
		@pyqtSlot()
		def click1_up(self): # update labels following selection of folder
		
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
		
		# button to get path to folder containing relevant files
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

		inname, _ = QFileDialog.getOpenFileName(self, "Select Model Variables File", "./PyCHAM/input/") # get path of file
		
		self.l7.clear() # clear old label
		self.l7.setText(inname)
		self.l7.show()
		
		return()
		
		
	@pyqtSlot()
	def on_click5(self): # interrogating the model variables file to check inputs
	
		@pyqtSlot()
		def click5_up(self): # update the display of model variables following click of check button since checking may have altered some
			
			# contents of model variables scroll area
			self.e9 = QLineEdit(self)
			self.e9.setText(sav_nam)
			self.varbox.addWidget(self.e9, 0, 1)
		
			self.e10 = QLineEdit(self)
			self.e10.setText((str(chem_sch_mark)).replace('\'', '').replace(' ', '')[1:-1])
			self.varbox.addWidget(self.e10, 1, 1)
		
			self.e11 = QLineEdit(self)
			self.e11.setText((str(tot_time)).replace('\'', '').replace(' ', ''))
			self.varbox.addWidget(self.e11, 2, 1)
		
			self.e12 = QLineEdit(self)
			self.e12.setText((str(update_stp)).replace('\'', '').replace(' ', ''))
			self.varbox.addWidget(self.e12, 3, 1)
		
			self.e13 = QLineEdit(self)
			self.e13.setText((str(save_step)).replace('\'', '').replace(' ', ''))
			self.varbox.addWidget(self.e13, 4, 1)
			
			self.e14 = QLineEdit(self)
			self.e14.setText((str(temp)).replace('\'', '').replace(' ', '')[1:-1])
			self.varbox.addWidget(self.e14, 5, 1)
			
			self.e15 = QLineEdit(self)
			self.e15.setText((str(tempt)).replace('\'', '').replace(' ', '')[1:-1])
			self.varbox.addWidget(self.e15, 6, 1)
			
			self.e16 = QLineEdit(self)
			self.e16.setText((str(RH)).replace('\'', '').replace(' ', ''))
			self.varbox.addWidget(self.e16, 7, 1)
			
			self.e17 = QLineEdit(self)
			self.e17.setText((str(Press)).replace('\'', '').replace(' ', ''))
			self.varbox.addWidget(self.e17, 8, 1)
			
			self.e18 = QLineEdit(self)
			self.e18.setText((str(siz_stru)).replace('\'', '').replace(' ', ''))
			self.varbox.addWidget(self.e18, 9, 1)
			
			self.e19 = QLineEdit(self)
			self.e19.setText((str(num_sb)).replace('\'', '').replace(' ', ''))
			self.varbox.addWidget(self.e19, 10, 1)
			
			# process pconc textbox presentation
			pconc_str = ''
			for ic in range(pconc.shape[1]): # loop through columns (times)
			
				if (pmode == 0): # if in modal form
					pconc_str = str(pconc_str+str(pconc[:, ic]).replace(' ', ':').replace('[', '').replace(']', ''))
				if (pmode == 1): # if in explicit form
					pconc_str = str(pconc_str+str(pconc[:, ic]).replace(' ', ',').replace('[', '').replace(']', ''))
			
				# if more times, separate with a semi-colon
				if ((ic+1) < pconc.shape[1]):
					pconc_str = str(pconc_str+';')
			
			self.e20 = QLineEdit(self)
			self.e20.setText(pconc_str)
			self.varbox.addWidget(self.e20, 11, 1)
			
			self.e21 = QLineEdit(self)
			self.e21.setText((str(pconct)).replace('\'', '').replace(' ', '').replace(',', ';').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e21, 12, 1)
			
			self.e22 = QLineEdit(self)
			self.e22.setText((str(seed_mw)).replace('\'', '').replace(' ', '').replace(',', ';').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e22, 13, 1)
			
			self.e23 = QLineEdit(self)
			self.e23.setText((str(seed_diss)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e23, 14, 1)
			
			self.e24 = QLineEdit(self)
			self.e24.setText((str(seed_dens)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e24, 15, 1)
			
			self.e25 = QLineEdit(self)
			self.e25.setText((str(seed_name)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e25, 16, 1)
			
			self.e26 = QLineEdit(self)
			self.e26.setText((str(seedVr)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e26, 17, 1)
			
			self.e27 = QLineEdit(self)
			self.e27.setText((str(lowsize)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e27, 18, 1)
			
			self.e28 = QLineEdit(self)
			self.e28.setText((str(uppsize)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e28, 19, 1)
			
			self.e29 = QLineEdit(self)
			self.e29.setText((str(space_mode)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e29, 20, 1)
			
			self.e30 = QLineEdit(self)
			self.e30.setText((str(std)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e30, 21, 1)
			
			self.e31 = QLineEdit(self)
			self.e31.setText((str(mean_rad)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e31, 22, 1)
			
			self.e32 = QLineEdit(self)
			self.e32.setText((str(new_partr)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e32, 23, 1)
			
			self.e33 = QLineEdit(self)
			self.e33.setText((str(nucv1)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e33, 24, 1)
			
			self.e34 = QLineEdit(self)
			self.e34.setText((str(nucv2)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e34, 25, 1)
			
			self.e35 = QLineEdit(self)
			self.e35.setText((str(nucv3)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e35, 26, 1)
			
			self.e36 = QLineEdit(self)
			self.e36.setText((str(nuc_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e36, 27, 1)
			
			self.e37 = QLineEdit(self)
			self.e37.setText((str(nuc_ad)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e37, 28, 1)
			
			self.e38 = QLineEdit(self)
			self.e38.setText((str(ser_H2O)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
			self.varbox.addWidget(self.e38, 29, 1)
			
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
					if i == ';':
						time_cnt += 1 # increase time count
					if (time_cnt == 1 and i == ','):
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
			
		# display inputs
		click5_up(self)
		
	@pyqtSlot()	
	def on_click79(self): # button to update model variables updated in text fields
	
		sav_nam = self.e9.text()
		self.e9 = QLineEdit(self)
		self.e9.setText(sav_nam)
		self.varbox.addWidget(self.e9, 0, 1)
		
		chem_sch_mark = self.e10.text().split(',')
		self.e10 = QLineEdit(self)
		self.e10.setText((str(chem_sch_mark)).replace('\'', '').replace(' ', '')[1:-1])
		self.varbox.addWidget(self.e10, 1, 1)
		
		tot_time = float(self.e11.text())
		self.e11 = QLineEdit(self)
		self.e11.setText((str(tot_time)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.e11, 2, 1)
		
		update_stp = float(self.e12.text())
		self.e12 = QLineEdit(self)
		self.e12.setText((str(update_stp)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.e12, 3, 1)
		
		save_step = float(self.e13.text())
		self.e13 = QLineEdit(self)
		self.e13.setText((str(save_step)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.e13, 4, 1)
				
		temp = [float(i) for i in ((self.e14.text().strip()).split(','))]
		self.e14 = QLineEdit(self)
		self.e14.setText((str(temp)).replace('\'', '').replace(' ', '')[1:-1])
		self.varbox.addWidget(self.e14, 5, 1)
		
		tempt = [float(i) for i in ((self.e15.text().strip()).split(','))]
		self.e15 = QLineEdit(self)
		self.e15.setText((str(tempt)).replace('\'', '').replace(' ', '')[1:-1])
		self.varbox.addWidget(self.e15, 6, 1)
	
		RH = float(self.e16.text())
		self.e16 = QLineEdit(self)
		self.e16.setText((str(RH)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.e16, 7, 1)
		
		Press = float(self.e17.text())
		self.e17 = QLineEdit(self)
		self.e17.setText((str(Press)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.e17, 8, 1)
		
		siz_stru = int(self.e18.text())
		self.e18 = QLineEdit(self)
		self.e18.setText((str(siz_stru)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.e18, 9, 1)
	
		num_sb = int(self.e19.text())
		self.e19 = QLineEdit(self)
		self.e19.setText((str(num_sb)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.e19, 10, 1)
	
		value = str(self.e20.text()) # get the string version of pconc from text box
		
		# update the form of particle number size distribution (modal or explicit)
		time_cnt = 1 # track number of times
		sb_cnt = 1 # track number of size bins
		mode_cnt = 1 # track number of modes
		for i in value:
			if i == ';':
				time_cnt += 1 # increase time count
			if (time_cnt == 1 and i == ','):
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
		
		
		# process pconc textbox presentation
		pconc_str = ''
		for ic in range(pconc.shape[1]): # loop through columns (times)
			
			if (pmode == 0): # if in modal form
				pconc_str = str(pconc_str+str(pconc[:, ic]).replace(' ', ':').replace('[', '').replace(']', ''))
			if (pmode == 1): # if in explicit form
				pconc_str = str(pconc_str+str(pconc[:, ic]).replace(' ', ',').replace('[', '').replace(']', ''))
			
			# if more times, separate with a semi-colon
			if ((ic+1) < pconc.shape[1]):
				pconc_str = str(pconc_str+';')
			
		self.e20 = QLineEdit(self)
		self.e20.setText(pconc_str)
		self.varbox.addWidget(self.e20, 11, 1)
		
		
		value = str(self.e21.text()) # get the string version of pconc from text box
		# process times of particle number size distribution appearance
		time_cnt = 1 # track number of times
		for i in value:
			if (i == ';'):
				time_cnt += 1 # increase time count
					
		# times in columns
		pconct = np.zeros((1, time_cnt))
		pconct[0, :] = [float(i) for i in ((value.strip()).split(';'))]
	
		self.e21 = QLineEdit(self)
		self.e21.setText((str(pconct)).replace('\'', '').replace(' ', ';').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e21, 12, 1)
		
		seed_mw = float(self.e22.text())
		self.e22 = QLineEdit(self)
		self.e22.setText((str(seed_mw)).replace('\'', '').replace(' ', ''))
		self.varbox.addWidget(self.e22, 13, 1)
		
		seed_diss = [float(i) for i in (self.e23.text().split(','))]
		self.e23 = QLineEdit(self)
		self.e23.setText((str(seed_diss)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e23, 14, 1)
		
		seed_dens = float(self.e24.text())
		self.e24 = QLineEdit(self)
		self.e24.setText((str(seed_dens)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e24, 15, 1)
		
		seed_name = [str(i) for i in (self.e25.text().split(','))]
		self.e25 = QLineEdit(self)
		self.e25.setText((str(seed_name)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e25, 16, 1)
		
		seedVr = [int(i) for i in (self.e26.text().split(','))]
		self.e26 = QLineEdit(self)
		self.e26.setText((str(seedVr)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e26, 17, 1)
		
		lowsize =float(self.e27.text())
		self.e27 = QLineEdit(self)
		self.e27.setText((str(lowsize)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e27, 18, 1)
		
		uppsize =float(self.e28.text())
		self.e28 = QLineEdit(self)
		self.e28.setText((str(uppsize)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e28, 19, 1)
		
		space_mode =str(self.e29.text())
		self.e29 = QLineEdit(self)
		self.e29.setText((str(space_mode)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e29, 20, 1)
		
		std = [float(i) for i in (self.e30.text().split(','))]
		self.e30 = QLineEdit(self)
		self.e30.setText((str(std)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e30, 21, 1)
		
		mean_rad = [float(i) for i in (self.e31.text().split(','))]
		self.e31 = QLineEdit(self)
		self.e31.setText((str(mean_rad)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e31, 22, 1)
		
		new_partr = float(self.e32.text())
		self.e32 = QLineEdit(self)
		self.e32.setText((str(new_partr)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e32, 23, 1)
		
		nucv1 = float(self.e33.text())
		self.e33 = QLineEdit(self)
		self.e33.setText((str(nucv1)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e33, 24, 1)
		
		nucv2 = float(self.e34.text())
		self.e34 = QLineEdit(self)
		self.e34.setText((str(nucv2)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e34, 25, 1)
		
		nucv3 = float(self.e35.text())
		self.e35 = QLineEdit(self)
		self.e35.setText((str(nucv3)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e35, 26, 1)
		
		nuc_comp = [str(i) for i in (self.e36.text().split(','))]
		self.e36 = QLineEdit(self)
		self.e36.setText((str(nuc_comp)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e36, 27, 1)
		
		nuc_ad = int(self.e37.text())
		self.e37 = QLineEdit(self)
		self.e37.setText((str(nuc_ad)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e37, 28, 1)
		
		ser_H2O = int(self.e38.text())
		self.e38 = QLineEdit(self)
		self.e38.setText((str(ser_H2O)).replace('\'', '').replace(' ', '').replace('[', '').replace(']', ''))
		self.varbox.addWidget(self.e38, 29, 1)
		
		# prepare for pickling
		list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, tot_time, comp0, y0, temp, tempt, RH, Press, wall_on, Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, uppsize, space_mode, std, mean_rad, save_step, const_comp, Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, seedVr, light_stat, light_time, daytime, lat, lon, af_path, dayOfYear, photo_path, tf, light_ad, con_infl_nam, con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, accom_comp, accom_val, uman_up, int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O]

		input_by_sim = str(os.getcwd() + '/PyCHAM/pickle.pkl')
		with open(input_by_sim, 'wb') as pk: # the file to be used for pickling
			pickle.dump(list_vars, pk) # pickle
			pk.close() # close
	
	@pyqtSlot()	
	def on_click80(self): # button to check validity of model variables
		
		@pyqtSlot()
		def click80_up(self): # update the display of model variables following click of check button since checking may have altered some
			
			# contents of model variables scroll area
			self.e9 = QLineEdit(self)
			self.e9.setText(sav_nam)
			self.varbox.addWidget(self.e9, 0, 1)
		
			self.e10 = QLineEdit(self)
			self.e10.setText((str(chem_sch_mark)).replace('\'', '').replace(' ', '')[1:-1])
			self.varbox.addWidget(self.e10, 1, 1)
		
			self.e11 = QLineEdit(self)
			self.e11.setText((str(tot_time)).replace('\'', '').replace(' ', ''))
			self.varbox.addWidget(self.e11, 2, 1)
		
			self.e12 = QLineEdit(self)
			self.e12.setText((str(update_stp)).replace('\'', '').replace(' ', ''))
			self.varbox.addWidget(self.e12, 3, 1)
		
			self.e13 = QLineEdit(self)
			self.e13.setText((str(save_step)).replace('\'', '').replace(' ', ''))
			self.varbox.addWidget(self.e13, 4, 1)
			
			self.show()
			return()
		
		import ui_check # check on inputs module
		# check on inputs - note this loads the last saved pickle file and saves any change to this pickle file
		ui_check.ui_check()
		
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
		
		# update displayed model variables
		click80_up(self)
		
	@pyqtSlot()	
	def on_click81(self): # button to run simulation
		import middle # prepare to communicate with main programme
		middle.middle() # call on modules to solve problem
		
	# plotting options ----------------------------------------------------------------------------------------
	
	@pyqtSlot()
	def on_click81c(self): # when different model results folder requires selection

		# button to get path to folder containing relevant files
		options = QFileDialog.Options()
		fol_nme = QFileDialog.getExistingDirectory(self, "Select Folder Containing Required Input Files", "./PyCHAM/output/")
		
		self.l81b.clear() # clear old label
		self.l81b.setText(fol_nme)
		print(self.l81b.text())
		self.l81b.show()
		
		return()
	
	@pyqtSlot() # button to plot standard results graphically
	def on_click82(self):	
		import plotter
		dir_path = self.l81b.text() # name folder with results
		plotter.plotter(0, dir_path) # plot results
		
	@pyqtSlot() # button to plot gas-phase concentration temporal profile
	def on_click83(self):	
		
		# get names of components to plot
		comp_names = [str(i) for i in self.e83.text(). split(',')]
		
		import plotter_gp
		dir_path = self.l81b.text() # name folder with results
		plotter_gp.plotter(0, dir_path, comp_names) # plot results
	
	@pyqtSlot() # button to plot total particle-phase concentration temporal profile
	def on_click84(self):	
		
		# get names of components to plot
		comp_names = [str(i) for i in self.e84.text(). split(',')]
		
		import plotter_pp
		dir_path = self.l81b.text() # name of folder with results
		plotter_pp.plotter(0, dir_path, comp_names) # plot results
	
	@pyqtSlot() # button to plot temporal profile of total concentration of 
	# components that have gas-wall partitioned to wall
	def on_click85(self):	
		
		# get names of components to plot
		comp_names = [str(i) for i in self.e85.text(). split(',')]
		
		import plotter_wp
		dir_path = self.l81b.text() # name of folder with results
		plotter_wp.plotter(0, dir_path, comp_names) # plot results
	
	@pyqtSlot() # button to plot temporal profile of total concentration of 
	# components on wall due to particle deposition to wall
	def on_click86(self):	
		
		# get names of components to plot
		comp_names = [str(i) for i in self.e86.text(). split(',')]
		
		import plotter_wp_part
		dir_path = self.l81b.text() # name of folder with results
		plotter_wp_part.plotter(0, dir_path, comp_names) # plot results
		
	# quiting option ------------------------------------------------------------------------------------------
	
	@pyqtSlot() # button to quit software
	def on_click89(self):
		QWidget.close(self)
		
	def openFileNameDialog(self): # allows opening of system's directory navigator
		options = QFileDialog.Options()
		fname, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
		return(fname)
