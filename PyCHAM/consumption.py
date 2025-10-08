########################################################################
#                                                                      #
# Copyright (C) 2018-2025                                              #
# Simon O'Meara : simon.omeara@manchester.ac.uk                        #
#                                                                      #
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
'''estimating consumption of a component over simulation'''
# a module to calculate the consumed mass concentration of a 
# component introduced artificially to chamber, so not produced
# by means other than through injection

import scipy.constants as si # scientific constants
import numpy as np # for arithmetic

def cons(self, caller):

	# inputs: -------------------------------
	# self.dir_path - path to results
	# self - reference to GUI
	# caller - flag for calling function
	# ---------------------------------------

	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = np.zeros((self.ro_obj.yrec.shape[0], self.ro_obj.yrec.shape[1]))
	yrec[:, :] = self.ro_obj.yrec[:, :]
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = np.zeros((self.ro_obj.Nrec_wet.shape[0], self.ro_obj.Nrec_wet.shape[1]))
	Nwet[:, :] = self.ro_obj.Nrec_wet[:, :]
	Ndry = np.zeros((self.ro_obj.Nrec_dry.shape[0], self.ro_obj.Nrec_dry.shape[1]))
	Ndry[:, :] = self.ro_obj.Nrec_dry[:, :]
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = np.squeeze((np.array((self.ro_obj.comp_MW))).reshape(-1, 1))
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec = np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	space_mode = self.ro_obj.spacing
	Cfac = self.ro_obj.cfac
	Cfac = (np.array(Cfac)).reshape(-1, 1) # convert to numpy array from list
	inname = self.ro_obj.vp
	
	dydt_list = [] # prepare to store change tendencies
	# loop through components to plot to check they are available
	for comp_name in (self.comp_names_to_plot):
		comp_name = comp_name.strip() # remove any white space
		
		fname = str(self.dir_path + '/' + comp_name +'_rate_of_change')
		try: # try to open
			dydt_list.append(np.loadtxt(fname, delimiter = ',', skiprows = 1)) # skiprows = 1 omits header	
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

	compi = 0 # count on components

	# get times of interest
	indxt = (timehr>=self.tmin)*(timehr<=self.tmax)

	# prepare to sum chemical reaction losses
	cons = np.zeros((sum(indxt)-1))

	# loop through components to plot
	for comp_name in (self.comp_names_to_plot):
		comp_name = comp_name.strip() # remove any white space
		
		# prepare to store chemical reaction losses in # molecules/cm^3/s
		# note that length of time axis in dydt_list[compi] has -1
		# because of header in dydt_list[compi]
		crl = np.zeros((dydt_list[compi].shape[0]-1, 1))

		for ti in range(dydt_list[compi].shape[0]-1): # loop through times
			indx = dydt_list[compi][ti+1, 0:-2] < 0 # indices of reactions that lose component
			crl[ti] = dydt_list[compi][ti+1, 0:-2][indx].sum()
		
		# convert change tendencies from # molecules/cm^3/s to ug/m^3/s
		crl = ((crl/si.N_A)*y_MW[comp_names.index(self.comp_names_to_plot[compi].strip())])*1.e12
		
		# integrate chemical losses over time intervals within 
		# the period of interest (ug/m^3) for total consumed,
		# averaging (arithmetic mean) the change tendency over 
		# each time step
		cons[:] += ((crl[indxt, 0][0:-1]+crl[indxt, 0][1::])/2.)*(np.diff(timehr[indxt])*3.6e3)
		
		compi += 1 # count on components
	

	if (caller == 0): # call from the consumption button
		
		self.l203a.setText(str('Total consumption of ' + str(self.comp_names_to_plot) + 
						 ': ' + str(sum(-1*cons)) + ' ' + u'\u03BC' + 'g/m' + u'\u00B3'))
		# set border around message
		if (self.bd_pl == 1):
			self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
			self.bd_pl = 2
		else:
			self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
			self.bd_pl = 1
		return() # return now

	if (caller == 1): # call from the yield button
		# concentrations of components in particle-phase at end of user-defined time interval (# molecules/cm3)
		SOA = yrec[timehr==max(timehr[indxt]), num_comp:num_comp*(num_sb-wall_on+1)]

		# remove seed and water in all size bins
		SOA[0, seedi[0]::num_comp] = 0.
		SOA[0, H2Oi::num_comp] = 0.
		
		# convert from # molecules/cm^3 to ug/m^3
		SOA = ((SOA/si.N_A)*np.tile(y_MW, (num_sb-wall_on)))*1.e12

		# sum for total (ug/m^3)
		SOAfi = np.sum(SOA)

		# concentrations of components in particle phase at start of time interval (# molecules/cm^3)
		SOA = yrec[timehr==min(timehr[indxt]), num_comp:num_comp*(num_sb-wall_on+1)]
		# remove seed and water in all size bins
		SOA[0, seedi[0]::num_comp] = 0.
		SOA[0, H2Oi::num_comp] = 0.
		
		# convert from # molecules/cm^3 to ug/m^3
		SOA = ((SOA/si.N_A)*np.tile(y_MW, (num_sb-wall_on)))*1.e12

		# sum for total (ug/m^3)
		SOAst = np.sum(SOA)

		# get model variables file name
		try:
			inname = inname[::-1][0:inname[::-1].index('/')][::-1]
		except:
			inname = inname[::-1][0:inname[::-1].index('\\')][::-1]
		
		inname = str(self.dir_path + '/inputs/' + inname)	
		
		inputs = open(inname, mode= 'r' ) # open model variables file
		in_list = inputs.readlines() # read file and store everything into a list
		inputs.close() # close file

		# loop through supplied model variables to interpret
		for i in range(len(in_list)):
			
			# ----------------------------------------------------
			# if commented out continue to next line
			if (in_list[i][0] == '#'):
				continue
			try:
				key, value = in_list[i].split('=') # split values from keys
			except:
				continue
			if key.strip() == 'dil_fac': # dilution factor rate
				dil_fac = np.array(([float(i) for i in (((value.strip()).split(',')))]))
				
			if key.strip() == 'dil_fact': # dilution factor rate times through experiment (s)
				dil_fact = np.array(([float(i) for i in (((value.strip()).split(',')))]))

		# in case there is dilution during this interval
		try: # will work if dilt given
			# get final dilution factor in this time inetrval
			dil_facfi = dil_fac[sum(dil_fact<=max(timehr[indxt]))]
			# get first dilution factor in this time inetrval
			dil_facst = dil_fac[sum(dil_fact<=min(timehr[indxt]))]
			dil_fac = np.mean(dil_fac[dil_facst:dil_facfi+1])
		except:
			try:
				dil_fac = dil_fac			
			except:
				dil_fac = 0.


		# integrate loss of SOA due to dilution over time interval
		tint = (max(timehr[indxt])-min(timehr[indxt]))*3600. # time interval (s)
		# note calculation of mean SOA in this time interval
		SOA_loss_by_dil = ((SOAfi+SOAst)/2.)*tint*dil_fac

		yld = (SOA_loss_by_dil+(SOAfi-SOAst))/sum(-1*cons)

	
		self.l203a.setText(str('Mass SOA yield of ' + str(self.comp_names_to_plot) + ' (%) : ' + str(yld*100.)  + ' (%)'))
		# set border around message
		if (self.bd_pl == 1):
			self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
			self.bd_pl = 2
		else:
			self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
			self.bd_pl = 1
		return() # return now

		
	return()
