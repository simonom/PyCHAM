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
'''estimating consumption of a component over simulation'''
# a module to calculate the consumed mass concentration of a 
# component introduced artificially to chamber, so not produced
# by means other than through injection

import scipy.constants as si # scientific constants
import retr_out # retrieving information
import numpy as np # for arithmetic

def cons(dir_path, self, caller):

	# inputs: -------------------------------
	# dir_path - path to results
	# self - reference to GUI
	# caller - flag for calling function
	# ---------------------------------------

	# retrieve results
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_mw, Nwet, comp_names, y_MV, _, wall_on, space_mode, indx_plot, 
		comp0, _, PsatPa, OC, H2Oi, seedi, _, _, _, tot_in_res, 
		_) = retr_out.retr_out(dir_path)
	

	# loop through components to plot to check they are available
	for comp_name in (self.comp_names_to_plot):
		
		fname = str(dir_path+ '/' +comp_name +'_rate_of_change')
		try: # try to open
			dydt = np.loadtxt(fname, delimiter = ',', skiprows = 1) # skiprows = 1 omits header	
		except:
			mess = str('Please note, a change tendency record for the component ' + str(comp_name) + ' was not found, was it specified in the tracked_comp input of the model variables file?  Please see README for more information.')
			self.l203a.setText(mess)
			
			# set border around error message
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			
			return()

	# if all files are available, then proceed without error message
	mess = str('')
	self.l203a.setText(mess)
			
	if (self.bd_pl < 3):
		self.l203a.setStyleSheet(0., '0px solid red', 0., 0.)
		self.bd_pl == 3

	# note that penultimate column in dydt is gas-particle 
	# partitioning and final column is gas-wall partitioning, whilst
	# the first row contains chemical reaction numbers
	# extract the change tendency due to gas-particle partitioning
	gpp = dydt[1::, -2]
	# extract the change tendency due to gas-wall partitioning
	gwp = dydt[1::, -1]
	# sum chemical reaction gains
	crg = np.zeros((dydt.shape[0]-1, 1))
	# sum chemical reaction losses
	crl = np.zeros((dydt.shape[0]-1, 1))
	for ti in range(dydt.shape[0]-1): # loop through times
		indx = dydt[ti+1, 0:-2] < 0 # indices of reactions that lose component
		crl[ti] = dydt[ti+1, 0:-2][indx].sum()
	
	# convert change tendencies from molecules/cm3/s to ug/m3/s
	crl = ((crl/si.N_A)*y_mw[comp_names.index(self.comp_names_to_plot[0])])*1.e12
	indxt = (timehr>=self.tmin)*(timehr<=self.tmax)

	# integrate chemical losses over period of interest (ug/m3) for total consumed
	cons = np.abs((np.sum(crl[indxt[0:-1]])/len(indxt[0:-1]))*((timehr[indxt][-1]-timehr[indxt][0])*3.6e3))

	cons = yrec[:, comp_names.index(self.comp_names_to_plot[0])]
	cons = (((cons[0]-cons[-1])*Cfac[0])/si.N_A)*y_mw[comp_names.index(self.comp_names_to_plot[0])]*1.e12

	if (caller == 0): # call from the consumption button
	
		self.l203a.setText(str('Consumption of ' + self.comp_names_to_plot[0] + ': ' + str(cons) + ' ' + u'\u03BC' + 'g/m' + u'\u00B3'))
		# set border around message
		if (self.bd_pl == 1):
			self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
			self.bd_pl = 2
		else:
			self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
			self.bd_pl = 1
		return() # return now

	if (caller == 1): # call from the yield button

		# concentrations of components in particle phase at end of simulation (# molecules/cm3)
		SOA = yrec[-2, num_comp:num_comp*(num_sb-wall_on+1)]
		# remove seed and water in all size bins
		SOA[seedi[0]::num_comp] = 0.
		SOA[H2Oi::num_comp] = 0.
		
		# convert from # molecules/cm3 to ug/m3
		SOA = ((SOA/si.N_A)*np.tile(y_mw, (num_sb-wall_on)))*1.e12

		# sum for total (ug/m3)
		SOA = np.sum(SOA)

		yld = SOA/cons

	
		self.l203a.setText(str('Yield of ' + self.comp_names_to_plot[0] + ': ' + str(yld)))
		# set border around message
		if (self.bd_pl == 1):
			self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
			self.bd_pl = 2
		else:
			self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
			self.bd_pl = 1
		return() # return now

	
	
	return()