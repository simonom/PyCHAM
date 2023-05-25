##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2023 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
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
'''module to prepare PyCHAM for partitioning variable calculation (particle and wall)'''

# module responsible for preparing inputs to the calculation of the gas-particle 
# partitioning coefficient in kimt_calc

import numpy as np
import scipy.constants as si
import diff_vol_est

def prep(y_mw, TEMP, num_comp, act_comp, act_user, acc_comp, 
	accom_coeff_user, num_sb, num_asb, Pnow, 
	Pybel_object, name_SMILE, self):
	
	# ------------------------------------------------------------------
	# inputs:
	# y_mw - molecular weight of components (g/mol) (num_comp,1)
	# TEMP - temperature of chamber at start of experiment (K)
	# num_comp - number of components
	# self.testf - flag for whether in normal mode (0) or testing mode (1) or 
	#	plotting of gas-phase diffusion coefficients mode (2)
	# self.Cw - effective absorbing mass of wall (g/m3 (air))
	# act_comp - names of components (corresponding to chemical scheme name) with 
	# 			activity coefficient stated in act_user
	# act_user - user-specified activity coefficients of components with names given in
	#			act_comp
	# accom_comp - names of components with accommodation coefficient set by the user
	# accom_coeff_user - accommodation coefficient set by the user
	# self.comp_namelist - names of components as stated in the chemical scheme
	# num_sb - number of size bins (including wall)
	# num_asb - number of actual size bins excluding wall
	# Pnow - air pressure inside chamber (Pa)
	# Pybel_object - Pybel objects for components
	# name_SMILE - SMILE strings of components
	# self - reference to program
	# self.kw - rate of transfer of components to wall (/s)
	# -----------------------------------------------------------------
	
	# start by assuming no error message
	err_mess = ''

	if (self.testf == 1): # if in testing mode (for test_front.py)
		return(0,0,0,0,0,0,0, err_mess) # return dummies
	
	surfT = 72. # assume surface tension of water (g/s2==mN/m==dyn/cm) for all particles
	
	# dynamic viscosity of air (kg/m.s), eq. 4.54 of Jacobson 2005
	dyn_visc = 1.8325e-5*((416.16/(TEMP+120.))*(TEMP/296.16)**1.5)
	
	ma = 28.966e-3 # molecular weight of air (kg/mol) (eq. 16.17 Jacobson 2005)
	
	# air density (kg/m3 (air)), ideal gas law
	rho_a =  (Pnow*ma)/((si.R)*TEMP)
	
	# mean thermal speed of each component (m/s) (eq. 2.3 Jacobson 2005)
	# note that we need the weight of one molecule, which is why y_mw is divided by
	# Avogadro's constant, and we need it in kg, which is why we multiply by 1e-3
	therm_sp = ((8.*si.k*TEMP)/(np.pi*(y_mw/si.N_A)*1.e-3))**0.5
	
	# mean free path (m) for each component (15.24 of Jacobson 2005)
	mfp = (2.*dyn_visc/(rho_a*therm_sp)).reshape(-1, 1)

	nv = (Pnow/(si.R*TEMP))*si.N_A # concentration of molecules (# molecules/m3)
	
	# get diffusion volumes
	diff_vol = diff_vol_est.diff_vol_est(Pybel_object)
	
	# append water and core (water from Table 4.1 of the Taylor (1993) textbook 
	# Multicomponent Mass Transfer, ISBN: 0-471-57417-1)
	diff_vol = (np.append(diff_vol, np.array((13.1, 1.)))).reshape(-1, 1)
	
	# diffusion coefficient (m2/s) of components in gas phase (air), eq 4.1.4 of
	# the Taylor (1993) textbook 
	# Multicomponent Mass Transfer, ISBN: 0-471-57417-1, note diffusion volume for air 
	# (19.7) taken from Table 4.1 of Taylor (1993) and mw of air converted to g/mol from 
	# kg/mol.  This is a replication of the original method from Fuller et al. (1969): 
	# doi.org/10.1021/j100845a020
	Dstar_org = 1.013e-2*TEMP**1.75*(((y_mw+ma*1.e3)/(y_mw*ma*1.e3))**0.5)/(Pnow*(diff_vol**(1./3.)+19.7**(1./3.))**2.)
	
	# convert to cm2/s
	Dstar_org = Dstar_org*1.e4

	if (self.testf == 2):
		import matplotlib.pyplot as plt
		from matplotlib.colors import BoundaryNorm
		from matplotlib.ticker import MaxNLocator
		from matplotlib.colors import LinearSegmentedColormap # for customised colormap
		import matplotlib.ticker as ticker # set colormap tick labels to standard notation
		plt.ion() # allow plotting
		# prepare plot
		fig, (ax0) = plt.subplots(1, 1, figsize = (14, 7))
		# plot gas-phase diffusion coefficients (cm2/s)
		ax0.plot(np.arange(len(self.comp_namelist)), Dstar_org, '+')
		ax0.set_ylabel(r'Gas-phase diffusion coeffiecient (cm$\rm{^{2}}\,$s$\rm{^{-1}}$)', fontsize = 14)
		ax0.set_xlabel(r'Component name', fontsize = 14)
		# set location of x ticks
		ax0.set_xticks(np.arange(len(self.comp_namelist)))
		ax0.set_xticklabels(self.comp_namelist, rotation = 90)
		ax0.set_title(str('Gas-phase diffusion coeffiecients at ' + str(TEMP) + ' K and ' + str(Pnow) + ' Pa'), fontsize = 14)
		err_mess = 'Stop'

	if (self.testf == 3):
		import matplotlib.pyplot as plt
		from matplotlib.colors import BoundaryNorm
		from matplotlib.ticker import MaxNLocator
		from matplotlib.colors import LinearSegmentedColormap # for customised colormap
		import matplotlib.ticker as ticker # set colormap tick labels to standard notation
		plt.ion() # allow plotting
		# prepare plot
		fig, (ax0) = plt.subplots(1, 1, figsize = (14, 7))
		# plot gas-phase diffusion coefficients (cm2/s)
		ax0.plot(np.arange(len(self.comp_namelist)), therm_sp, '+')
		ax0.set_ylabel(r'Gas-phase mean thermal speed (m$\,$s$\rm{^{-1}}$)', fontsize = 14)
		ax0.set_xlabel(r'Component name', fontsize = 14)
		# set location of x ticks
		ax0.set_xticks(np.arange(len(self.comp_namelist)))
		ax0.set_xticklabels(self.comp_namelist, rotation = 90)
		ax0.set_title(str('Gas-phase mean thermal speeds at ' + str(TEMP) + ' K'), fontsize = 14)
		err_mess = 'Stop'
	
	# accommodation coefficient of components in each size bin
	accom_coeff = np.ones((num_comp, num_sb))*1.e0
	
	# list containing accommodation coefficients that are functions
	accom_coeff_func = []
	
	ac_indx = []
	for i in range(len(acc_comp)): # user-defined accommodation coefficients
		# get index of component stated
		ac_indx.append(self.comp_namelist.index(acc_comp[i].strip()))

	# check for any accommodation coefficients set by user
	if len(ac_indx)>0:
		for i in range(len(ac_indx)):

			# if it is a constant (not a function, which would be a string)
			if isinstance(accom_coeff_user[i], str) == False:
				accom_coeff[ac_indx] = accom_coeff_user[i]
			# if it is a function, it will be a string and needs making available to the 
			# kimt_calc module
			if isinstance(accom_coeff_user[i], str) == True:
				accom_coeff_func.append(str('accom_coeff[' + str(ac_indx[i]) + ', :]' + ' = ' + accom_coeff_user[i]))
	
	# generate module that contains any accommodation coefficient functions, note, do 
	# this even if no functions supplied so that the accomm_coeff_calc is updated and
	# accurate for this simulation
	f = open('PyCHAM/accom_coeff_calc.py', mode='w')
	f.write('##########################################################################################\n')
	f.write('#                                                                                        											 #\n')
	f.write('#    Copyright (C) 2018-2023 Simon O\'Meara : simon.omeara@manchester.ac.uk                  				 #\n')
	f.write('#                                                                                       											 #\n')
	f.write('#    All Rights Reserved.                                                                									 #\n')
	f.write('#    This file is part of PyCHAM                                                         									 #\n')
	f.write('#                                                                                        											 #\n')
	f.write('#    PyCHAM is free software: you can redistribute it and/or modify it under              						 #\n')
	f.write('#    the terms of the GNU General Public License as published by the Free Software       					 #\n')
	f.write('#    Foundation, either version 3 of the License, or (at your option) any later          						 #\n')
	f.write('#    version.                                                                            										 #\n')
	f.write('#                                                                                        											 #\n')
	f.write('#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT                						 #\n')
	f.write('#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       			 #\n')
	f.write('#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              				 #\n')
	f.write('#    details.                                                                            										 #\n')
	f.write('#                                                                                        											 #\n')
	f.write('#    You should have received a copy of the GNU General Public License along with        					 #\n')
	f.write('#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                 							 #\n')
	f.write('#                                                                                        											 #\n')
	f.write('##########################################################################################\n')
	f.write('\'\'\'module for calculating accommodation coefficients, automatically generated by kimt_prep\'\'\'\n')
	f.write('\n')
	f.write('# code that expresses and performs the functions for accommodation \n')
	f.write('# coefficients that are given by the user in the model variables file \n')
	f.write('# and that are intended to be calculated real-time via the kimt_calc \n')
	f.write('# module \n')
	f.write('\n')

	# following part is the function (there should be an indent at the start of each line
	# following the def line - suggest using 1 tab)
	f.write('def accom_coeff_func(accom_coeff, radius):\n')
	f.write('\n')
	f.write('	# ------------------------------------------------------------------ \n')
	f.write('	# inputs:\n')
	f.write('	# accom_coeff - array containing accommodation coefficients for all \n')
	f.write('	# components\n')
	f.write('	# radius - radii of size bins (m)\n')
	f.write('	# ------------------------------------------------------------------ \n')
	f.write('\n')
	f.write('	# calculate accommodation coefficients \n')
	# code to calculate accommodation coefficients as given by user 
	for line in accom_coeff_func:
		f.write('	%s \n' %line)
	f.write('\n')
	f.write('	return(accom_coeff)\n')
	f.close()
	
	# activity coefficient of components - affects the particle- and wall-phase
	act_coeff = np.ones((1, num_comp))
	for i in range(len(act_comp)): # user-defined activity coefficients
		# get index of component stated
		ac_indx = self.comp_namelist.index(act_comp[i].strip())
		act_coeff[0, ac_indx] = act_user[i].strip()
	
	# in preparation for use in ode solver, tile activity coefficients over
	# particle and wall bins
	act_coeff = np.tile(act_coeff, (num_sb, 1))
	# convert Cw (effective absorbing mass of wall) from g/m3 (air) to 
	# # molecules/cm3 (air), assuming a molecular weight of 200 g/mol (*1.e-6 to convert from
	# /m3 (air) to /cm3 (air))
	self.Cw = ((self.Cw*1.e-6)/200.)*si.N_A
	# ensure effective absorbing mass of walls represents walls in rows and components in columns
	self.Cw = np.tile((self.Cw.reshape(-1, 1)), (1, num_comp))
	# in case user has not given a Cw value for every wall
	if (self.Cw.shape[0] == 1 and self.wall_on > 1):
		self.Cw = np.tile((self.Cw), (self.wall_on, 1))

	# start of mass transfer rate coefficient to wall -----------------
	
	# empty array ready to hold useful mass transfer rate coefficient to wall matrix
	kwn = np.zeros((self.wall_on, num_comp))
	pre_count = 0 # count on prescribed components
	mask = np.zeros((self.wall_on, num_comp)) # masking matrix
	mask = mask == 0.
	
	if (self.wall_on > 0): # if wall present
		# loop through user-defined inputs
		for wi in range(self.kw.shape[0]): # loop through walls
			for ci in range(self.kw.shape[1]): # loop through components	
				if self.kw[wi,ci] == -1.e-6: # flag for empty space in array
					continue
				
				if self.kw[wi,ci] == -1.e-7: # coefficient specified for a component

					# becuase water may not be in chemical scheme
					if self.wmtc_names[pre_count] == 'H2O':
						kwn[wi, self.H2Oi] = self.wmtc[pre_count]
						mask[wi, self.H2Oi] = False
						pre_count += 1 # count on prescribed components
						continue

					if (self.remove_influx_not_in_scheme == 1):
						
						try:
							# coefficient of this component on this wall
							kwn[wi, self.comp_namelist.index(self.wmtc_names[pre_count])] = self.wmtc[pre_count]
							mask[wi, self.comp_namelist.index(self.wmtc_names[pre_count])] = False
						except:
							# just ignore if told to by user-defined 
							# self.remove_influx_not_in_scheme, and leave as zero
							continue
					if (self.remove_influx_not_in_scheme == 0):
						try:
							# coefficient of this component on this wall
							kwn[wi, self.comp_namelist.index(self.wmtc_names[pre_count])] = self.wmtc[pre_count]
							mask[wi, self.comp_namelist.index(self.wmtc_names[pre_count])] = False
						except:
							# give error message
							err_mess = str('Error: component ' + str(self.wmtc_names[pre_count]) + ' has a gas-wall mass transfer coefficient but has not been identified in the chemical scheme')
							break

					pre_count += 1 # count on prescribed components

				else: # coefficient specified for all components
					kwn[wi, mask[wi, :]] = self.kw[wi, ci]
	
	# finally convert kw to kwn
	self.kw = np.zeros((((self.wall_on, num_comp))))
	self.kw[:, :] = kwn[:, :]

	# end of mass transfer rate coefficient to wall -------------------

	R_gas = si.R # ideal gas constant (kg.m2.s-2.K-1.mol-1)
	NA = si.Avogadro # Avogadro's constant (molecules/mol)

	return(mfp, accom_coeff, therm_sp, surfT, act_coeff, 
			R_gas, NA, diff_vol, Dstar_org, err_mess, self)
