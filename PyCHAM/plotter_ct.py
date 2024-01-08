##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2024 Simon O'Meara : simon.omeara@manchester.ac.uk                #
#                                                                                       											 #
#    All Rights Reserved.                                                                									 #
#    This file is part of PyCHAM                                                         									 #
#                                                                                        											 #
#    PyCHAM is free software: you can redistribute it and/or modify it under              #
#    the terms of the GNU General Public License as published by the Free Software        #
#    Foundation, either version 3 of the License, or (at your option) any later           #
#    version.                                                                            										 #
#                                                                                        											 #
#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT                #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS        #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more               #
#    details.                                                                            										 #
#                                                                                        											 #
#    You should have received a copy of the GNU General Public License along with         #
#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                 							 #
#                                                                                        											 #
##########################################################################################
'''plots results for the change tendency temporal profiles of specified components'''
# simulation results are represented graphically

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import os
import numpy as np
import scipy.constants as si

def plotter(caller, dir_path, comp_names_to_plot, self):
	
	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0) or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# comp_names_to_plot - chemical scheme names of components to plot
	# self - reference to PyCHAM
	# --------------------------------------------------------------------------

	# chamber condition ---------------------------------------------------------
		
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
	y_mw = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	y_MV = (np.array((self.ro_obj.comp_MV))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec= np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	x = self.ro_obj.cen_size
	space_mode = self.ro_obj.spacing
	Cfac = self.ro_obj.cfac
	wall_on = self.ro_obj.wf
	PsatPa = self.ro_obj.vpPa
	OC = self.ro_obj.O_to_C
	
	# no record of change tendency for final experiment time point
	timehr = timehr[0:-1]

	# loop through components to plot to check they are available
	for comp_name in (comp_names_to_plot):
		
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
			
			plt.ioff() # turn off interactive mode
			plt.close() # close figure window
			
			return()
	
	
	# if all files are available, then proceed without error message
	mess = str('')
	self.l203a.setText(mess)
			
	if (self.bd_pl < 3):
		self.l203a.setStyleSheet(0., '0px solid red', 0., 0.)
		self.bd_pl == 3
	
	# prepare figure
	plt.ion() # display figure in interactive mode
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
	
	for comp_name in (comp_names_to_plot): # loop through components to plot
		
		if (comp_name != 'HOMRO2' and comp_name != 'RO2'):
			ci = comp_names.index(comp_name) # get index of this component
			
		if (comp_name == 'HOMRO2'):
			# get mean molecular weight of HOMRO2 (g/mol)
			counter = 0 # count on HOMRO2 components
			mw_extra = 0. # count on HOMRO2 molecular weights (g/mol)
			for cni in range(len(comp_names)): # loop through all components
				if 'API_' in comp_names[cni] or 'api_' in comp_names[cni]:
					if 'RO2' in comp_names[cni]:
						mw_extra += y_mw[cni]
						counter += 1
			y_mw = np.concatenate((y_mw, (np.array(mw_extra/counter)).reshape(1)))
			ci = len(y_mw)-1
			
		if (comp_name == 'RO2'):
			# get indices of RO2
			RO2indx = (np.array((ro_obj.gi['RO2i'])))
			# get mean molecular weight of RO2 (g/mol)
			mw_extra = np.sum((np.array((y_mw)))[RO2indx])/len(RO2indx)
						
			y_mw = np.concatenate((y_mw, (np.array(mw_extra)).reshape(1)))
			ci = len(y_mw)-1
		
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
			indx = dydt[ti+1, 0:-2] > 0 # indices of reactions that produce component
			crg[ti] = dydt[ti+1, 0:-2][indx].sum()
			indx = dydt[ti+1, 0:-2] < 0 # indices of reactions that lose component
			crl[ti] = dydt[ti+1, 0:-2][indx].sum()
			 
		# convert change tendencies from molecules/cm3/s to ug/m3/s
		gpp = ((gpp/si.N_A)*y_mw[0, ci])*1.e12
		gwp = ((gwp/si.N_A)*y_mw[0, ci])*1.e12
		crg = ((crg/si.N_A)*y_mw[0, ci])*1.e12
		crl = ((crl/si.N_A)*y_mw[0, ci])*1.e12
			 
		# plot temporal profiles of change tendencies due to chemical 
		# reaction production and loss, gas-particle partitioning and gas-wall partitioning
		ax0.plot(timehr, gpp, label = str('gas-particle partitioning '+ comp_name))
		ax0.plot(timehr, gwp, label = str('gas-wall partitioning '+ comp_name))
		ax0.plot(timehr, crg, label = str('chemical reaction gain '+ comp_name))
		ax0.plot(timehr, crl, label = str('chemical reaction loss '+ comp_name))
		ax0.yaxis.set_tick_params(direction = 'in')
		
		ax0.set_title('Change tendencies, where a tendency to decrease \ngas-phase concentrations is treated as negative')
		ax0.set_xlabel('Time through experiment (hours)')
		ax0.set_ylabel('Change tendency ($\mathrm{\{mu} g\, m^{-3}\, s^{-1}}$)')
		
		ax0.yaxis.set_tick_params(direction = 'in')
		ax0.xaxis.set_tick_params(direction = 'in')
		
		ax0.legend()
			 

	return()
	
# plot change tendency due to individual chemical reactions
def plotter_ind(caller, dir_path, comp_names_to_plot, top_num, uc, self):
	
	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0) or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# comp_names_to_plot - chemical scheme names of components to plot
	# top_num - top number of chemical reactions to plot
	# uc - units to use for change tendency
	# self - reference to GUI
	# --------------------------------------------------------------------------

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
	y_mw = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	y_MV = (np.array((self.ro_obj.comp_MV))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec= np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	x = self.ro_obj.cen_size
	space_mode = self.ro_obj.spacing
	Cfac = self.ro_obj.cfac
	wall_on = self.ro_obj.wf
	PsatPa = self.ro_obj.vpPa
	OC = self.ro_obj.O_to_C
	
	# loop through components to plot to check they are available
	for comp_name in (comp_names_to_plot):
		
		fname = str(dir_path+ '/' +comp_name +'_rate_of_change')
		try: # try to open
			dydt = np.loadtxt(fname, delimiter = ',', skiprows = 0) # skiprows = 0 skips first header	
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
			
			plt.ioff() # turn off interactive mode
			plt.close() # close figure window
			
			return()
	
	# if all files are available, then proceed without error message
	mess = str('')
	self.l203a.setText(mess)
			
	if (self.bd_pl < 3):
		self.l203a.setStyleSheet(0., '0px solid red', 0., 0.)
		self.bd_pl == 3
	
	# prepare figure
	plt.ion() # display figure in interactive mode
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
	for comp_name in (comp_names_to_plot): # loop through components to plot
		
		if (comp_name != 'HOMRO2' and comp_name != 'RO2'):
			ci = comp_names.index(comp_name) # get index of this component
			
		if (comp_name == 'HOMRO2'):
			# get mean molecular weight of HOMRO2 (g/mol)
			counter = 0 # count on HOMRO2 components
			mw_extra = 0. # count on HOMRO2 molecular weights (g/mol)
			for cni in range(len(comp_names)): # loop through all components
				if 'API_' in comp_names[cni] or 'api_' in comp_names[cni]:
					if 'RO2' in comp_names[cni]:
						mw_extra += y_mw[cni]
						counter += 1
			y_mw = np.concatenate((y_mw, (np.array(mw_extra/counter)).reshape(1)))
			ci = len(y_mw)-1
			
		if (comp_name == 'RO2'):
			# get indices of RO2
			RO2indx = (np.array((ro_obj.gi['RO2i'])))
			# get mean molecular weight of RO2 (g/mol)
			mw_extra = np.sum((np.array((y_mw)))[RO2indx])/len(RO2indx)
						
			y_mw = np.concatenate((y_mw, (np.array(mw_extra)).reshape(1)))
			ci = len(y_mw)-1
		
		# note that penultimate column in dydt is gas-particle 
		# partitioning and final column is gas-wall partitioning, whilst
		# the first row contains chemical reaction numbers
		
		# prepare to store results of change tendency due to chemical reactions
		res = np.zeros((dydt.shape[0], dydt.shape[1]-2))
		res[:, :] = dydt[:, 0:-2] # get chemical reaction numbers and change tendencies
		
		import scipy.constants as si

		if (uc == 0):
			Cfaca = (np.array(Cfac)).reshape(-1, 1) # convert to numpy array from list
			# convert change tendencies from # molecules/cm3/s to ppb
			res[1::, :] = (res[1::, :]/Cfaca[1::])
			ct_units = str('(ppb/s)')
		if (uc == 1):			
			# convert change tendencies from # molecules/cm3/s to ug/m3/s
			res[1::, :] = ((res[1::, :]/si.N_A)*y_mw[ci])*1.e12
			ct_units = str('(' + u'\u03BC' + 'g/m' +u'\u00B3' + '/s)')
		if (uc == 2):
			# keep change tendencies as # molecules/cm3/s
			res[1::, :] = res[1::, :]
			ct_units = str('\n(' + u'\u0023' + ' molecules/cm' +u'\u00B3' + '/s)')
		
		# identify most active chemical reactions
		# first sum total change tendency over time (ug/m3/s)
		res_sum = np.abs(np.sum(res[1::, :], axis=0))
		 
		# sort in ascending order
		res_sort = np.sort(res_sum)
		
		if (len(res_sort) < top_num[0]):

			# if less reactions are present than the number requested inform user
			mess = str('Please note that ' + str(len(res_sort)) + ' relevant reactions were found, although a maximum of ' + str(top_num[0]) + ' were requested by the user.')
			self.l203a.setText(mess)
		
		# get all reactions out of the used chemical scheme --------------------------------------------------------------------
		import sch_interr # for interpeting chemical scheme
		import re # for parsing chemical scheme
		import scipy.constants as si
		
		try:
			sch_name = self.ro_obj.sp
			inname = self.ro_obj.vp
			f_open_eqn = open(sch_name, mode= 'r' ) # open model variables file
			inputs = open(inname, mode= 'r' ) # open model variables file
		except:
			# get index of slashes in path to scheme
			try:
				slashi = (sch_name[::-1]).index('/') # in case saved on UNIX
			except:
				slashi = (sch_name[::-1]).index('\\') # in case saved on windows
			# get scheme name
			sch_name = sch_name[-slashi::]
			sch_name = str(dir_path + '/inputs/' + sch_name)
			
			f_open_eqn = open(sch_name, mode= 'r' ) # open chemical scheme file

			# get index of slashes in path to model variables
			try:
				slashi = (inname[::-1]).index('/') # in case saved on UNIX
			except:
				slashi = (inname[::-1]).index('\\') # in case saved on windows
			# get scheme name
			inname = inname[-slashi::]
			inname = str(dir_path + '/inputs/' + inname)
			
			inputs = open(inname, mode= 'r' ) # open model variables file
		inname = self.ro_obj.vp
		
		
		# read the file and store everything into a list
		total_list_eqn = f_open_eqn.readlines()
		f_open_eqn.close() # close file
		
		in_list = inputs.readlines() # read file and store everything into a list
		inputs.close() # close file

		# start by assuming default chemical scheme markers
		self.chem_sch_mrk = ['{', 'RO2', '+', 'C(ind_', ')','' , '&', '' , '', ':', '}', ';', '']

		for i in range(len(in_list)): # loop through supplied model variables to interpret

			# ----------------------------------------------------
			# if commented out continue to next line
			if (in_list[i][0] == '#' or in_list[i][0:12] == 'param_ranges'):
				continue
			try:
				key, value = in_list[i].split('=') # split values from keys
				# model variable name - a string with bounding white space removed
				key = key.strip()
				# ----------------------------------------------------

				if key == 'chem_scheme_markers' and (value.strip()): # formatting for chemical scheme
					self.chem_sch_mrk = [str(i).strip() for i in (value.split(','))]
			except:
				continue
		# interrogate scheme to list equations
		[rrc, rrc_name, RO2_names, self] = sch_interr.sch_interr(total_list_eqn, self)	
	
		for cnum in range(np.min([top_num[0], len(res_sort)])): # loop through chemical reactions
			
			# identify this chemical reaction
			cindx = np.where((res_sort[-(cnum+1)] ==  res_sum) == 1)[0]
			
			for indx_two in (cindx):
				
				reac_txt = str(self.eqn_list[int(res[0, indx_two])]) # get equation text
				# plot, note the +1 in the label to bring label into MCM index
				ax0.plot(timehr[0:-1], res[1::, indx_two], label = str(' Eq. # ' + str(int(res[0, indx_two])+1) + ':  ' + reac_txt))
		
		ax0.yaxis.set_tick_params(direction = 'in')
		
		ax0.set_title('Change tendencies, where a tendency to decrease \ngas-phase concentrations is negative')
		ax0.set_xlabel('Time through experiment (hours)')
		ax0.set_ylabel(str('Change tendency ' + ct_units))
		
		ax0.yaxis.set_tick_params(direction = 'in')
		ax0.xaxis.set_tick_params(direction = 'in')
		
		ax0.legend()
			 

	return()

# output total production of component
def plotter_prod(caller, dir_path, comp_names_to_plot, tp, uc, self):
	
	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0) or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# comp_names_to_plot - chemical scheme names of components to plot
	# tp - times to calculate between (hours)
	# uc - units to use for change tendency
	# self - reference to GUI
	# --------------------------------------------------------------------------

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
	y_mw = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	y_MV = (np.array((self.ro_obj.comp_MV))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec= np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	x = self.ro_obj.cen_size
	space_mode = self.ro_obj.spacing
	Cfac = self.ro_obj.cfac
	wall_on = self.ro_obj.wf
	PsatPa = self.ro_obj.vpPa
	OC = self.ro_obj.O_to_C

	# loop through components due to be plotted, to check they are available
	for comp_name in (comp_names_to_plot):
		
		fname = str(dir_path+ '/' +comp_name +'_rate_of_change')
		try: # try to open
			dydt = np.loadtxt(fname, delimiter = ',', skiprows = 0) # skiprows = 0 skips first header	
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
			
			plt.ioff() # turn off interactive mode
			plt.close() # close figure window
			
			return()

	# if all files are available, then proceed without error message
	mess = str('')
	self.l203a.setText(mess)
			
	if (self.bd_pl < 3):
		self.l203a.setStyleSheet(0., '0px solid red', 0., 0.)
		self.bd_pl == 3

	for comp_name in (comp_names_to_plot): # loop through components to plot
		
		if (comp_name != 'HOMRO2' and comp_name != 'RO2'):
			ci = comp_names.index(comp_name) # get index of this component
			
		if (comp_name == 'HOMRO2'):
			# get mean molecular weight of HOMRO2 (g/mol)
			counter = 0 # count on HOMRO2 components
			mw_extra = 0. # count on HOMRO2 molecular weights (g/mol)
			for cni in range(len(comp_names)): # loop through all components
				if 'API_' in comp_names[cni] or 'api_' in comp_names[cni]:
					if 'RO2' in comp_names[cni]:
						mw_extra += y_mw[cni]
						counter += 1
			y_mw = np.concatenate((y_mw, (np.array(mw_extra/counter)).reshape(1)))
			ci = len(y_mw)-1
			
		if (comp_name == 'RO2'):
			# get indices of RO2
			RO2indx = (np.array((ro_obj.gi['RO2i'])))
			# get mean molecular weight of RO2 (g/mol)
			mw_extra = np.sum((np.array((y_mw)))[RO2indx])/len(RO2indx)
						
			y_mw = np.concatenate((y_mw, (np.array(mw_extra)).reshape(1)))
			ci = len(y_mw)-1
		
		# note that penultimate column in dydt is gas-particle 
		# partitioning and final column is gas-wall partitioning, whilst
		# the first row contains chemical reaction numbers
		
		# prepare to store results of change tendency due to chemical reactions
		res = np.zeros((dydt.shape[0], dydt.shape[1]-2))
		res[:, :] = dydt[:, 0:-2] # get chemical reaction numbers and change tendencies
		
		if (uc == 0):
			Cfaca = (np.array(Cfac)).reshape(-1, 1) # convert to numpy array from list
			# convert change tendencies from # molecules/cm3/s to ppb/s
			res[1::, :] = (res[1::, :]/Cfaca[1::])
			ct_units = str('ppb')
		if (uc == 1):			
			# convert change tendencies from # molecules/cm3/s to ug/m3/s
			res[1::, :] = ((res[1::, :]/si.N_A)*y_mw[ci])*1.e12
			ct_units = str(u'\u03BC' + 'g/m' +u'\u00B3')
		if (uc == 2):
			# keep change tendencies as # molecules/cm3/s
			res[1::, :] = res[1::, :]
			ct_units = str(u'\u0023' + ' molecules/cm' +u'\u00B3')
		
		# remove reactions that destroy component
		res[1::, :][res[1::, :] < 0] = 0.

		# sum production rates over production reactions
		res_sum = np.sum(res[1::, :], axis=1)
		
		# retain only the required time period
		tindx = (timehr >= tp[0])*(timehr < tp[1])
		res_sum = res_sum[tindx[0:-1]]

		# time intervals over this period (s)
		tindx = (timehr >= tp[0])*(timehr <= tp[1])
		tint = (timehr[tindx][1::]-timehr[tindx][0:-1])*3600.
		# if one short then concatenate assuming same interval
		if len(tint) < (len(res_sum)):
			tint = np.concatenate((tint, np.reshape(np.array(tint[-1]), 1)))
		
		# integrate production rate to get total production
		res_sum = np.sum(res_sum*tint)

		# display to message board
		mess = str(mess + 'Total production of ' + str(comp_name) + ': ' + str(res_sum) + ' ' + ct_units)
		self.l203a.setText(mess)
			
		# set border around message
		if (self.bd_pl == 1):
			self.l203a.setStyleSheet(0., '2px dashed magenta', 0., 0.)
			self.bd_pl = 2
		else:
			self.l203a.setStyleSheet(0., '2px solid magenta', 0., 0.)
			self.bd_pl = 1
		return() # return now

	return()

def plotter_reac_ratios(self):
	
	# inputs: ------------------------------------------------------------------
	# self - reference to PyCHAM
	# --------------------------------------------------------------------------

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
	y_mw = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	y_MV = (np.array((self.ro_obj.comp_MV))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec= np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	x = self.ro_obj.cen_size
	space_mode = self.ro_obj.spacing
	Cfac = self.ro_obj.cfac
	wall_on = self.ro_obj.wf
	PsatPa = self.ro_obj.vpPa
	OC = self.ro_obj.O_to_C


	# prepare figure
	plt.ion() # display figure in interactive mode
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))	
	
	# empty results array for numerator and denominator
	num = np.zeros((len(timehr), 1))
	den = np.zeros((len(timehr), 1))

	# loop through file names in results path
	for fnamei in os.listdir(self.dir_path):
		if (fnamei[-15::] == '_rate_of_change' and np.sum(self.num_reac_num != -1.e6)+np.sum(self.den_reac_num != -1.e6)>0 ):
			
			dydt = np.loadtxt(str(self.dir_path + '/' + fnamei), delimiter = ',', skiprows = 1) # skiprows = 1 omits header
			res = dydt[:, 0:-2] # get chemical reaction numbers and change tendencies
			# transform recation numbers from starting at 0 to (pythonic) to starting at 1 (user-supplied)
			res[0, :] += 1

			rn_cnt = 0 # count on numerator reaction numbers
			for reac_num in self.num_reac_num:
				
				if sum(res[0, :] == reac_num) == 1:
					# add to numerator
					num[:, :] += np.abs(res[:, res[0, :] == reac_num])
					# remove this reaction from list
					self.num_reac_num[rn_cnt] = -1.e6
				rn_cnt += 1 # count on reaction number

			rn_cnt = 0 # count on denominator reaction numbers
			for reac_num in self.den_reac_num:
				if sum(res[0, :] == reac_num) == 1:
					# add to denominator
					den[:, :] += np.abs(res[:, res[0, :] == reac_num])
					# remove this reaction from list
					self.den_reac_num[rn_cnt] = -1.e6
				rn_cnt += 1 # count on reaction number

	# plot
	ax0.plot(timehr[den[:,0]>0.], num[den>0.]/den[den>0.])
	return()

# plot resevoir sizes of carbon
def plotter_carb_res(self):
	
	# inputs: ------------------------------------------------------------------
	# self - reference to PyCHAM
	# --------------------------------------------------------------------------

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
	y_mw = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	y_MV = (np.array((self.ro_obj.comp_MV))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec= np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	x = self.ro_obj.cen_size
	space_mode = self.ro_obj.spacing
	Cfac = self.ro_obj.cfac
	wall_on = self.ro_obj.wf
	PsatPa = self.ro_obj.vpPa
	OC = self.ro_obj.O_to_C
	mv_path = self.ro_obj.vp
	yrec_p2w = self.ro_obj.part_to_wall

	# prepare figure
	plt.ion() # display figure in interactive mode
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))

	init_flag = 0 # flag for whether initial concentrations are present


	# estimate total carbon influxed by user settings ------------
	inputs = open(mv_path, mode= 'r' ) # open model variables file
	in_list = inputs.readlines() # read file and store everything into a list
	inputs.close() # close file
	
	# get model variable value for influxed components
	for i in range(len(in_list)): # loop through supplied model variables to interpret

		try:
			key, value = in_list[i].split('=') # split values from keys
		except:
			continue

		# model variable name - a string with bounding white space removed
		key = key.strip()

		if key == 'Comp0' and (value.strip()): # names of components present at experiment start
			comp0 = [str(i).strip() for i in (value.split(','))]			

		if key == 'C0' and (value.strip()): # initial concentrations of components present at experiment start (ppb)
			y0 = [float(i) for i in (value.split(','))]
			init_flag = 1
		if key == 'const_infl' and (value.strip()): # names of components with continuous influx
			# check if this is a path to a file containing continuous influxes, if not treat as a list of component names
			if '/' in value or '\\' in value: # treat as path to file containing continuous influxes
				self.const_infl_path = str(value.strip())
				self = const_infl_open(self)
				
			else: # treat as list of components 
				self.con_infl_nam = [str(i).strip() for i in (value.split(','))]

		if key == 'const_infl_t' and (value.strip()): # times of continuous influxes (s)
			self.con_infl_t = [float(i.strip()) for i in (value.split(','))]
			self.con_infl_t = np.array((self.con_infl_t))

		if key == 'Cinfl' and (value.strip()): # influx rate of components with continuous influx (ppb/s)
			comp_count = 1 # count number of components
			time_count = 1 # track number of times
			for i in value:
				if (i==';'):
					comp_count += 1 # record number of components
				if (i==',' and comp_count == 1):
					time_count += 1 # record number of times
			self.con_infl_C = np.zeros((comp_count, time_count))
				
			try:	
				for i in range(comp_count): # loop through components
					for ii in range(time_count): # loop through times
						self.con_infl_C[i, ii] = float((((value.split(';')[i]).split(',')))[ii].strip())
							
			# in case semicolons and commas messed up on input, note this will invoke an 
			# error message from the user input check module
			except:
				self.con_infl_C = np.empty(0)
		
		if key == 'dil_fac' and (value.strip()): # dilution factor rate
			self.dil_fac = np.array(([float(i) for i in (((value.strip()).split(',')))]))
				
		if key == 'dil_fact' and (value.strip()): # dilution factor rate times through experiment (s)
			self.dil_fact = np.array(([float(i) for i in (((value.strip()).split(',')))]))

	# get carbon reservoir at start of simulation
	if init_flag == 1:
		# molar mass of carbon in each component present at start (g/mol)
		mm_C = np.zeros((len(comp0), 1))

		for i in range(len(comp0)): # loop through starting components
			compi = comp_names.index(comp0[i])
			mm_C[i, 0] = (rel_SMILES[compi].count('C')+rel_SMILES[compi].count('c'))*12.0107
		# now convert influxes to # molecules/cm3 from ppb
		y0 = y0*Cfac[0]
		# now divide by Avogadro's constant to convert # molecules/cm3 to mol/cm3
		y0 = y0/si.N_A
		# now multiply by molar mass of carbon and g/cm3 to ug/m3 conversion factor to get ug/m3
		y0 = y0*mm_C*1.e12

		# sum over components (ug/m3)
		y0 = (np.sum(y0)).reshape(1, 1)
	else:
		y0 = np.zeros((1, 1))


	# now convert influx of components (ppb/s) into influx of carbon (ug/m3/s) -------
	# molar mass of carbon in each influxed component (g/mol)
	mm_C = np.zeros((len(self.con_infl_nam), 1))
	for i in range(len(self.con_infl_nam)): # loop through influxed components
		compi = comp_names.index(self.con_infl_nam[i])
		mm_C[i, 0] = (rel_SMILES[compi].count('C')+rel_SMILES[compi].count('c'))*12.0107
	
	# now convert influxes to # molecules/cm3/s from ppb/s
	self.con_infl_C = self.con_infl_C*Cfac[0]
	# now divide by Avogadro's constant to convert # molecules/cm3/s to mol/cm3/s
	self.con_infl_C = self.con_infl_C/si.N_A
	# now multiply by molar mass of carbon to get g/cm3/s
	self.con_infl_C = self.con_infl_C*mm_C
	# convert g/cm3/s to ug/m3/s
	self.con_infl_C = self.con_infl_C*1.e12
	
	# sum over components
	self.con_infl_C = (np.sum(self.con_infl_C, axis = 0)).reshape(1, -1)

	# only keep the influxes within the simulated time
	self.con_infl_C = self.con_infl_C[0, self.con_infl_t<timehr[-1]*3600.]

	# only keep the times within the simulated time
	self.con_infl_t = self.con_infl_t[self.con_infl_t<timehr[-1]*3600.]
	
	# now, get time intervals (s)
	# append end simulation time to times
	self.con_infl_t = np.concatenate((self.con_infl_t, (np.array((timehr[-1]*3600.)).reshape(1))))
	# time intervals (s)
	t_int = self.con_infl_t[1::]-self.con_infl_t[0:-1]
	# integrate over time intervals (ug/m3)
	self.con_infl_C = self.con_infl_C*(t_int.reshape(1, -1))
	
	# consider reservoir present at start of simulation
	self.con_infl_C = np.concatenate((y0, self.con_infl_C), axis=1)

	# sum with time
	user_influx_C = np.cumsum(self.con_infl_C, axis=1)
	
	# ---------------------------------------
	# get carbon present in gas-phase at all times in order to calculate 
	# cumulative gas-phase loss through ventilation
	# convert gas-phase concentration from ppb to # molecules/cm3
	yrecg = yrec[:, 0:num_comp]*Cfac[0]
	# convert to mol/cm3
	yrecg = yrecg/si.N_A 

	# get carbon number of every component
	Cnum = np.zeros((len(rel_SMILES)))
	for compi in range(len(rel_SMILES)):
		Cnum[compi] = rel_SMILES[compi].count('C')+rel_SMILES[compi].count('c')

	# get molar mass of carbon per component (g of C/mol)
	Cnum = (Cnum*12.0107).reshape(1, -1)

	# get mass concentration of components (ug/m3 of C)
	yrecg = (yrecg*Cnum)*1.e12
	
	# sum over components (ug/m3 of C)
	yrecg = np.sum(yrecg, axis=1)
	
	# take average over a time interval to best represent what the 
	# ode solver works on
	yrecg[1::] = (yrecg[1::]+yrecg[0:-1])/2.

	# prepare for aligning dilution factors with recorded times
	dil_fac_align = np.zeros((len(yrecg)))

	# loop through record times to assign correct dilution factor, note
	# that  used < rather than <= because we want dilution factor up to that 
	# point, not at it.  By the same, principle, the first point 
	# (at time=0), must have a dilution factor of 0 because no time has 
	# passed to allow air exchange
	for ti in range(1, len(timehr)):
		
		dil_fac_align[ti] = (self.dil_fac[self.dil_fact<timehr[ti]*3.6e3])[-1]

	# every point represents the resorvoir as it was
	# up to that point, so we should go all the way up to the final recorded 
	# time point.  Note that in the first point 0 time had passed, so 0 carbon
	# can have been lost to air excahnge.  Therefore, include zero time passed in
	# time interval
	t_int = np.zeros((len(timehr)))
	t_int[1::] = np.diff(timehr*3.6e3)
	dil_fac_align = dil_fac_align*t_int

	# multiply by carbon concentration to get mass concentration removed 
	# by air exchange (ug/m3)
	ax_removed = np.cumsum(dil_fac_align*yrecg)

	# calculation of carbon lost through particle loss during air exchange --------
	if (num_sb-wall_on > 0):
		# concentration of all components in the particle phase (# molecules/cm3)
		yrecp = yrec[:, num_comp:num_comp*(num_sb-wall_on+1)]
		# convert # molecules/cm3 to mol/cm3
		yrecp = yrecp/si.N_A

		# sum over size bins
		for i in range(1, num_sb-wall_on):
			yrecp[:, 0:num_comp] += yrecp[:, num_comp*i:num_comp*(i+1)]
		
		# get mass concentration of components (note conversion from g/cm3 to ug/m3 of C)
		yrecp = (yrecp[:, 0:num_comp]*Cnum)*1.e12
	
		# sum over components (ug/m3 of C)
		yrecp = np.sum(yrecp, axis=1)
	
		# take average over a time step to best represent what the 
		# model works on for carbon lost through particle air exchange
		yrecp[1::] = (yrecp[1::]+yrecp[0:-1])/2.

		# multiply by dilution factor to get mass concentration removed 
		# by air exchange (ug/m3)
		ax_part_removed = np.cumsum(dil_fac_align*yrecp)
	else:
		ax_part_removed = np.zeros((len(timehr)))

	# particles on wall ---------------------------------------------
	if (num_sb-wall_on > 0 and wall_on > 0):

		# convert concentrations of components on wall due to particle deposition to wall
		# from # molecules/cm3 to mol/cm3
		yrec_p2w = yrec_p2w/si.N_A

		# sum over particle size bins
		for i in range(1, num_sb-wall_on):
			yrec_p2w[:, 0:num_comp] += yrec_p2w[:, num_comp*i:num_comp*(i+1)]
		
		# get mass concentration of carbon in components (note the conversion 
		# from g/cm3 to ug/m3 of C)
		yrec_p2w = (yrec_p2w[:, 0:num_comp]*Cnum)*1.e12
	
		# sum over components (ug/m3 of C)
		yrec_p2w = np.sum(yrec_p2w, axis=1)
	else:
		yrec_p2w = np.zeros((len(timehr)))

	# vapours on wall ---------------------------------------------
	if (wall_on > 0):
		# wall concentrations (# molecules/cm3)
		yrec_w = yrec[:, -num_comp*wall_on::]

		# convert concentrations of components on wall due to particle deposition to wall
		# from # molecules/cm3 to mol/cm3
		yrec_w = yrec_w/si.N_A

		# sum over wall bins
		for i in range(1, wall_on):
			yrec_w[:, 0:num_comp] += yrec_w[:, num_comp*i:num_comp*(i+1)]
		
		# get mass concentration of components (note conversion of g/cm3 
		# to ug/m3 of C)
		yrec_w = (yrec_w[:, 0:num_comp]*Cnum)*1.e12
	
		# sum over components (ug/m3 of C)
		yrec_w = np.sum(yrec_w[:, 0:num_comp], axis=1)
	else: # if no wall
		yrec_w = np.zeros((len(timehr)))

	# amount in gas phase at any one time ----------------------
	if (num_comp > 0): # if components are present
		# gas-phase concentrations (note conversion from ppb to # molecules/cm3)
		yrec_g = yrec[:, 0:num_comp]*Cfac[0]

		# convert concentrations of components in gas
		# from # molecules/cm3 to mol/cm3
		yrec_g = yrec_g/si.N_A

		# get mass concentration of components (ug/m3 of C)
		yrec_g = (yrec_g[:, 0:num_comp]*Cnum)*1.e12
	
		# sum over components (ug/m3 of C)
		yrec_g = np.sum(yrec_g, axis=1)

	else: # if no components in gas phase
		yrec_g = np.zeros((len(timehr)))

	# amount in particle phase at any one time ----------------------
	if (num_sb-wall_on > 0): # if components are present
		# particle-phase concentrations (# molecules/cm3)
		yrec_p = yrec[:, num_comp:num_comp*(num_sb-wall_on+1)]

		# convert concentrations of components
		# from # molecules/cm3 to mol/cm3
		yrec_p = yrec_p/si.N_A

		# sum over particle size bins
		for i in range(1, num_sb-wall_on-1):
			yrec_p[:, 0:num_comp] += yrec_p[:, num_comp*i:num_comp*(i+1)]

		# get mass concentration of components (ug/m3 of C)
		yrec_p = (yrec_p[:, 0:num_comp]*Cnum)*1.e12
	
		# sum over components (ug/m3 of C)
		yrec_p = np.sum(yrec_p[:, 0:num_comp], axis=1)

	else: # if no components in gas phase
		yrec_p = np.zeros((len(timehr)))


	ax0.plot((self.con_infl_t/3.6e3), user_influx_C[0, :], 'k', label = 'total in')
	ax0.stackplot(timehr, ax_removed, ax_part_removed, yrec_p2w, yrec_w, yrec_g, yrec_p, labels = ['exchange of gas', 'exchange of particle', 'particle on wall', 'vapour on wall', 'gas phase', 'particle phase'])
	ax0.set_xlabel('Time through experiment (hours)')
	ax0.set_ylabel(str('Cumulative Concentration (' + u'\u03BC' + 'g/m' +u'\u00B3' + ')'))
		
	ax0.yaxis.set_tick_params(direction = 'in')
	ax0.xaxis.set_tick_params(direction = 'in')
		
	ax0.legend()
	
	return()	
	

def const_infl_open(self): # define function to read in values relevant to constant influxes

	import openpyxl
	import os

	wb = openpyxl.load_workbook(filename = self.const_infl_path)
	sheet = wb['const_infl']
	# component names are in first column, times are in headers of first row		
	ic = 0 # count on row iteration
	
	# prepare to store component names
	self.con_infl_nam = []
	
	for i in sheet.iter_rows(values_only=True): # loop through rows
			if (ic == 0): # get times of influx (s through experiment)
				self.con_infl_t = np.array((i[1::]))
				# prepare to store emission rates
				self.con_infl_C = np.zeros((1, len(self.con_infl_t)))
				
			# get names of components (matching chemical scheme names) 
			# and their emission rates (ppb/s)
			else:
				# append component name
				self.con_infl_nam.append(i[0])
				
				# emission rates
				if (ic>self.con_infl_C.shape[0]): # if we need to concatenate
					self.con_infl_C = np.concatenate((self.con_infl_C, np.array((i[1::])).reshape(1, -1)), axis=0)
				else:
					self.con_infl_C[ic-1] = i[1::]
				
				
			ic += 1 # count on row iteration
		
	wb.close() # close excel file


	return(self)

# display properties of individual components
def plotter_individ_prop(self):
	
	# inputs: ------------------------------------------------------------------
	# self - reference to PyCHAM
	# --------------------------------------------------------------------------

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
	y_mw = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	y_MV = (np.array((self.ro_obj.comp_MV))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec= np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	x = self.ro_obj.cen_size
	space_mode = self.ro_obj.spacing
	Cfac = self.ro_obj.cfac
	wall_on = self.ro_obj.wf
	PsatPa = self.ro_obj.vpPa
	PsatPa0 = self.ro_obj.vpPa0
	OC = self.ro_obj.O_to_C
	mv_path = self.ro_obj.vp
	yrec_p2w = self.ro_obj.part_to_wall

	if ('Molar Mass' in self.single_comp_prop):

		# get molar mass of component of interest
		mm_interest = y_mw[0, comp_names.index(self.mm_comp_name)]
	
		# display in text area of GUI
		self.l203a.setText(str('Molar mass of ' + str(self.mm_comp_name) + ': ' + str(mm_interest) + ' g/mol'))

	if ('saturation vapour pressure at starting temperature' in self.single_comp_prop):

		PsatPa0 = np.squeeze(PsatPa0) # ensure minimum number of dimensions

		# get saturation vapour pressure at starting temperature of simulation
		vp0_interest = PsatPa0[comp_names.index(self.mm_comp_name)]
	
		# display in text area of GUI
		self.l203a.setText(str('Vapour pressure at starting temperature for ' + str(self.mm_comp_name) + ': ' + str(vp0_interest) + ' Pa'))

	if ('saturation vapour pressure at 298.15 K' in self.single_comp_prop):

		# get saturation vapour pressure at 298.15 K
		vp_interest = PsatPa[comp_names.index(self.mm_comp_name)]
	
		# display in text area of GUI
		self.l203a.setText(str('Vapour pressure at 298.15 K for ' + str(self.mm_comp_name) + ': ' + str(vp_interest) + ' Pa'))

	return(self)
