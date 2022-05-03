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
'''function to estimate total inputs of components'''
# for estimating the total input of components to the simulation
# based on the user-supplied model variables

import scipy.constants as si # for scientific constants
import numpy as np # for arithmetic

def tot_in(init_conc, Cfac, comp0, comp_namelist, y_mw,
		const_infl_t, tot_time, Compt, self): # define function

	# inputs: ----------------------------------------------
	# init_conc - initial concentrations (ppb)
	# Cfac - factor to convert ppb to # molecules/cm3
	# comp0 - chemical scheme names of components present 
	#	at simulation start
	# comp_namelist - list of all chemical scheme names
	# y_mw - molar mass of components (g/mol)
	# self.con_infl_nam - names of components with continuous influx
	# const_infl_t - times at which injections occur (s)
	# tot_time - total experiment time (s)
	# Compt - name of components injected instantaneously 
	#	after start of experiment
	# self - reference to PyCHAM
	# ------------------------------------------------------

	# indices of components injected
	tot_in_res_indx = []
	# total concentration of components injected
	tot_in_res_con = []

	ccnt = 0 # count on components

	for cnam in comp0: # loop through components present initially

		
		ci = comp_namelist.index(cnam) # index within all components
		tot_in_res_indx.append(ci) # remember component index

		# initial input ug/m3, note *1e12 converts g/cm3 to ug/m3
		Czero = ((init_conc[ccnt]*Cfac)/si.N_A)*y_mw[ci]*1.e12
		tot_in_res_con.append(Czero[0])
		
		ccnt += 1 # count on components

	# instantaneous injection -------------------------------------------
	Compti = [] # indices in record for components injected instantaneously
	# index for this record of components injected instantaneously after experiment start
	# loop through chemical scheme names of components injected instantaneously
	for cnam in Compt:
		if comp_namelist.index(cnam) in tot_in_res_indx:
			Compti.append(tot_in_res_indx.index(comp_namelist.index(cnam)))
			continue
		else:
			Compti.append(len(tot_in_res_indx))
			tot_in_res_indx.append(comp_namelist.index(cnam))
			tot_in_res_con.append(0.)

	# continuous injection ---------------------------------------------
	self.cont_inf_reci = [] # index inside record
	# loop through components with continuous influx
	for cnam in self.con_infl_nam: 
		if (cnam != 'H2O'): # omit H2O as this dealt with separately
			if comp_namelist.index(cnam) in tot_in_res_indx:
				self.cont_inf_reci.append(tot_in_res_indx.index(comp_namelist.index(cnam)))
				continue
			else:
				self.cont_inf_reci.append(len(tot_in_res_indx))
				tot_in_res_indx.append(comp_namelist.index(cnam))
				tot_in_res_con.append(0.)
	
	# array to hold all information (component indices in first 
	# column, influxes in second column)
	tot_in_res = np.array((tot_in_res_con))
	
	return(tot_in_res, Compti, tot_in_res_indx)