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
'''opening and reading any user-provided observation file'''
# called from eqn_pars.py and (if in simulation preparation mode) 
# ui_check.py

# define function
def obs_file_open(self):

	# inputs: --
	# self - PyCHAM object
	# ----------

	import openpyxl
	import os
	import numpy as np

	self.obs_file = str(self.obs_file)
	try: # in case full path provided
		wb = openpyxl.load_workbook(filename = self.obs_file)
	except: # in case in same folder as model variables file
		pd_indx = self.inname[::-1].index('/')
		pd = self.inname[0:-pd_indx]
		self.obs_file = str(pd +self.obs_file)
		wb = openpyxl.load_workbook(filename = self.obs_file)
	sheet = wb['PyCHAMobs']
	# time (seconds) is in first column, other 
	# variables in later columns
	ic = 0 # count on row iteration
	# column indices for temperature and 
	# relative humidity
	temper_indx = 0
	rh_indx = 0
	# index for columns of components
	comp_indx = []
	# loop through rows
	for i in sheet.iter_rows(values_only = True):
		if (ic == 0):
			# get chemical scheme names of variables
			names_xlsx = i[1::]
			# prepare for variable names
			self.obs_comp = []
			# keep count on column index
			col_num = 0
			for oc in names_xlsx:
				# keep count on column index
				col_num += 1
				if (oc == None):
					continue
				if (oc == 'Temperature (K)'):
					temper_indx = col_num
					self.TEMP = []
					self.tempt = []
					continue
				if (oc == 'RH (0-1)'):
					rh_indx = col_num
					self.RH = []
					self.RHt = []
					continue
				# if none of the above, then 
				# assume it's a chemical components
				self.obs_comp.append(oc)
				comp_indx.append(col_num)
			
			# in case of simulation preparation mode, set
			# these components as ones to track change 
			# tendencies of
			if (self.sim_ci_file != []):
				self.dydt_trak = self.obs_comp
					
			# get number of components to consider
			nc_obs = len(self.obs_comp)
			# prepare to store observations, not 
			# forgetting a column for time
			self.obs = np.zeros((1, nc_obs+1))
			comp_indx = np.array((comp_indx)).astype('int')
			rh_indx = np.array((rh_indx)).astype('int')
			temper_indx = np.array((
				temper_indx)).astype('int')
		else:
			if (ic > 1): # add row onto data array
				self.obs = np.concatenate((self.obs, 
					(np.zeros((1, nc_obs+1)))), 
					axis=0)
			
			self.obs[ic-1, 1::] = np.array((i))[
				0:col_num+1][comp_indx]
			# get times (s)
			self.obs[ic-1, 0] = np.array((i))[0]

			# temperature and relative humidity observations
			self.TEMP.append(np.array((i))[temper_indx])
			self.tempt.append(i[0])
			self.RH.append(np.array((i))[rh_indx])
			self.RHt.append(i[0])
			
		ic += 1 # count on row iteration

	# ensure numpy array for rh arrays
	self.RH = np.array((self.RH))
	self.RHt = np.array((self.RHt))
	wb.close() # close excel file
	
	
	if hasattr(self, 'sim_ci_file'):
		self.ci_array = np.zeros((len(self.obs_comp)+1, 
			1)).astype('str')
		self.ci_array[0, 0] = 'molec/cm3/s' # units
		self.ci_array[1::, 0] = self.obs_comp # component names
	
	# get indices of components with concentrations to fit 
	# observations
	self.obs_comp_i = np.zeros((len(self.obs_comp))).astype('int')
	
	# when called from ui_check, comp_namelist may not yet have 
	# been established
	if hasattr(self, 'comp_namelist'):
		for i in range(len(self.obs_comp)):
			self.obs_comp_i[i] = self.comp_namelist.index(
				self.obs_comp[i])

	return(self)
