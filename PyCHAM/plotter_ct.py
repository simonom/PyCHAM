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
'''plots results for the change tendency temporal profiles of specified components'''
# simulation results are represented graphically

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import os
import retr_out
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
	# retrieve results
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_mw, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, _, _, _, _, _, _, _) = retr_out.retr_out(dir_path)
	
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
			 
		# convert change tendencies from molecules/cc/s to ug/m3/s
		gpp = ((gpp/si.N_A)*y_mw[ci])*1.e12
		gwp = ((gwp/si.N_A)*y_mw[ci])*1.e12
		crg = ((crg/si.N_A)*y_mw[ci])*1.e12
		crl = ((crl/si.N_A)*y_mw[ci])*1.e12
			 
		# plot temporal profiles of change tendencies due to chemical 
		# reaction production and loss, gas-particle partitioning and gas-wall partitioning
		ax0.plot(timehr, gpp, label = str('gas-particle partitioning '+ comp_name))
		ax0.plot(timehr, gwp, label = str('gas-wall partitioning '+ comp_name))
		ax0.plot(timehr, crg, label = str('chemical reaction gain '+ comp_name))
		ax0.plot(timehr, crl, label = str('chemical reaction loss '+ comp_name))
		ax0.yaxis.set_tick_params(direction = 'in')
		
		ax0.set_title('Change tendencies, where a tendency to decrease \ngas-phase concentrations is treated as negative')
		ax0.set_xlabel('Time through experiment (hours)')
		ax0.set_ylabel('Change tendency ($\mathrm{\mu g\, m^{-3}\, s^{-1}}$)')
		
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

	# chamber condition ---------------------------------------------------------
	# retrieve results
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_mw, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, _, _, _, _, _, _, ro_obj) = retr_out.retr_out(dir_path)
	
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

		sch_name = ro_obj.sp
		inname = ro_obj.vp
		
		f_open_eqn = open(sch_name, mode= 'r' ) # open model variables file
		# read the file and store everything into a list
		total_list_eqn = f_open_eqn.readlines()
		f_open_eqn.close() # close file
		
		inputs = open(inname, mode= 'r' ) # open model variables file
		in_list = inputs.readlines() # read file and store everything into a list
		inputs.close() # close file
		for i in range(len(in_list)): # loop through supplied model variables to interpret

			# ----------------------------------------------------
			# if commented out continue to next line
			if (in_list[i][0] == '#'):
				continue
			key, value = in_list[i].split('=') # split values from keys
			# model variable name - a string with bounding white space removed
			key = key.strip()
			# ----------------------------------------------------

			if key == 'chem_scheme_markers' and (value.strip()): # formatting for chemical scheme
				chem_sch_mrk = [str(i).strip() for i in (value.split(','))]

		# interrogate scheme to list equations
		[eqn_list, aqeqn_list, eqn_num, rrc, rrc_name, 
			RO2_names] = sch_interr.sch_interr(total_list_eqn, chem_sch_mrk)	
	
		for cnum in range(np.min([top_num[0], len(res_sort)])): # loop through chemical reactions
			
			# identify this chemical reaction
			cindx = np.where((res_sort[-(cnum+1)] ==  res_sum) == 1)[0]
			
			for indx_two in (cindx):
				
				reac_txt = str(eqn_list[int(res[0, indx_two])]) # get equation text
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

	# chamber condition ---------------------------------------------------------
	# retrieve results
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_mw, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, _, _, _, _, _, _, ro_obj) = retr_out.retr_out(dir_path)

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

	# chamber condition ---------------------------------------------------------
	# retrieve results
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_mw, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, _, _, _, _, _, _, _) = retr_out.retr_out(self.dir_path)


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