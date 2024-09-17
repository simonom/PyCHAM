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
'''code for alterantive to eqn_pars when eqn_pars being skipped'''

def eqn_pars_skipper(self): # define function

	import numpy as np # matrix functions
	import scipy.constants as si
	from water_calc import water_calc
	# imports: ---------------------------------
	# self - reference to PyCHAM class
	# ------------------------------------------
	# if observation file provided for constraint
	if hasattr(self, 'obs_file') and self.obs_file != []:

		from obs_file_open import obs_file_open
		self = obs_file_open(self)

	# get index of components with continuous influx/concentration -----------
	# empty array for storing index of components with continuous influx
	self.con_infl_indx = np.zeros((len(self.con_infl_nam)))
	self.con_C_indx = np.zeros((len(self.const_comp))).astype('int')
	delete_row_list = [] # prepare for removing rows of unrecognised components

	icon = 0 # count on constant influxes

	for i in range (len(self.con_infl_nam)):
		
		# water not included explicitly in chemical schemes but accounted 
		# for later in init_conc
		if (self.con_infl_nam[icon] == 'H2O'):

			# if water not included explicitly in chemical schemes 
			# (note it is accounted for later in init_conc)
			if (self.H2O_in_cs == 2):
				self.con_infl_indx[icon] = int(self.comp_num)
				icon += 1 # count on constant influxes
				continue
			else: # if water in chemical scheme
				self.con_infl_indx[icon] = self.comp_namelist.index(
				self.con_infl_nam[icon])
				icon += 1 # count on constant influxes
				continue

		# if we want to remove continuous influxes 
		# not present in the chemical scheme
		if (self.remove_influx_not_in_scheme == 1):

			try:
				# index of where components with continuous 
				# influx occur in list of components
				self.con_infl_indx[icon] = self.comp_namelist.index(
				self.con_infl_nam[icon])
			except:
				# remove names of unrecognised components
				self.con_infl_nam = np.delete(self.con_infl_nam, (icon), axis=0)
				# remove emissions of unrecognised components
				self.con_infl_C = np.delete(self.con_infl_C, (icon), axis = 0)
				# remove empty indices of unrecognised components
				self.con_infl_indx = np.delete(
				self.con_infl_indx, (icon), axis = 0)
				icon -= 1 # count on constant influxes
		else:

			try:
				# index of where components with continuous 
				# influx occur in list of components
				self.con_infl_indx[i] = self.comp_namelist.index(
				self.con_infl_nam[i])
			except:
				erf = 1 # raise error

				err_mess = str('Error: continuous influx component with name ' +str(self.con_infl_nam[i]) + ' has not been identified in the chemical scheme, please check it is present and the chemical scheme markers are correct')
	
		icon += 1 # count on continuous influxes

	if len(self.con_infl_nam) > 0:

		# get ascending indices of components with continuous influx
		sindx = self.con_infl_indx.argsort()
		# order continuous influx indices, influx rates and 
		# names ascending by component index
		self.con_infl_indx = (self.con_infl_indx[sindx]).astype('int')
		self.con_infl_C = (self.con_infl_C[sindx])
		self.con_infl_nam = (self.con_infl_nam[sindx])

	# check if water in continuous influx components 
	if ('H2O' in self.con_infl_nam):

		self.H2Oin = 1 # flag for water influx

		# index of water in continuous influx
		wat_indx = self.con_infl_nam.tolist().index('H2O')

		# get influx rate of water
		self.con_infl_H2O = self.con_infl_C[wat_indx, :].reshape(1, -1)

		# do not allow continuous influx of water in the standard ode 
		# solver, instead deal with it inside the water ode-solver
		self.con_infl_C = np.delete(self.con_infl_C, wat_indx, axis=0)
		self.con_infl_indx = np.delete(self.con_infl_indx, wat_indx, axis=0)
		self.con_infl_nam = np.delete(self.con_infl_nam, wat_indx, axis=0)

	else:
		self.H2Oin = 0 # flag for no water influx

	# components with constant concentration
	for i in range (len(self.const_comp)):
		try:
			# index of where constant concentration components occur in list 
			# of components
			self.con_C_indx[i] = self.comp_namelist.index(self.const_comp[i])
		except:
			# if water then we know it will be the next component to 
			# be appended to the component list
			if (self.const_comp[i] == 'H2O'):
				self.con_C_indx[i] = len(self.comp_namelist)
			else: # if not water
				erf = 1 # raise error
				err_mess = str('Error: constant concentration component with name ' + str(self.const_comp[i]) + ' has not been identified in the chemical scheme, please check it is present and the chemical scheme markers are correct')

	# use eqn_pars output from previous simulation
	rowvals = self.rowvals; colptrs = self.colptrs 
	comp_num = self.comp_num; Jlen = self.Jlen  
	erf = 0.; err_mess = ''

	# if not using vapour pressures saved to file
	# get vapour pressure at first temperature (# molecules/cm3)
	self.Psat = self.Psat_rec0

	# get saturation vapour pressure of water (log10(atm))
	[_, Psat_water, _] = water_calc(self.TEMP[0], self.RH[0], si.N_A)	

	# convert to Pa
	Psat_water = (10**Psat_water)*101325.
	# convert to # molecules/cm3
	self.Psat[0, self.comp_namelist.index('H2O')] = Psat_water*(si.N_A/((si.R*1.e6)*self.TEMP[0]))

	return(rowvals, colptrs, comp_num, 
		Jlen, erf, err_mess)
