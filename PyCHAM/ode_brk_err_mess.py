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
'''edits the output file if the ODE solver fails to solve within the minimum integration time interval'''
# if the ODE solver can't produce positive results within the minimum integration 
# time interval then this module writes a file that contains the fluxes through 
# all reactions of the chemical scheme and the partitioning fluxes of components
# with negative concentrations

import numpy as np # for arithmetic
import datetime # for dealing with times
	
def ode_brk_err_mess(y0, neg_names, rindx, y_arr, y_rind, rstoi, 
			pstoi, rrc, nreac, nprod, num_comp, 
			num_asb, Cw, Psat, act_coeff, kw, neg_comp_indx,
			N_perbin, core_diss, kelv_fac, kimt, 
			eqn_num, rindx_aq, y_rind_aq, y_arr_aq, 
			rstoi_aq, pstoi_aq, nreac_aq, nprod_aq, call, 
			H2Oi, y, self):
	
	# inputs: ------------------------------------------------------------------
	# y0 - concentrations prior to attempted integration (# molecules/cm3)
	# neg_names - chemical scheme names of components with negative concentrations
	# rindx - index of reactants per equation
	# y_arr - index for matrix used to arrange concentrations of gas-phase reactants, 
	#	enabling calculation of reaction rate coefficients 	
	# y_rind - index of y relating to reactants for reaction rate 
	# 	coefficient equation
	# rstoi - stoichiometry of reactants
	# pstoi - stoichiometry of products
	# rrc - reaction rate coefficient (/s)
	# nreac - number of reactants per equation
	# nprod - number of products per equation
	# self.wall_on - whether wall considerations turned on or off
	# num_comp - number of components
	# num_asb - number of actual size bins (excluding wall)
	# Cw - effective absorbing mass concentration of wall 
	#	(molecules/cc (air))
	# Psat - pure component saturation vapour pressures (molecules/cc)
	# act_coeff - activity coefficient of components
	# kw - mass transfer coefficient to wall (/s)
	# neg_comp_indx - indices of components with negative concentrations
	# N_perbin - number concentration of particles per size bin (#/cc)
	# self.seedi - index of seed components
	# core_diss - dissociation of seed components
	# kelv_fac - kelvin factor for particle
	# kimt - mass transfer coefficient for gas-particle partitioning (s)
	# eqn_num - number of gas- and particle-phase reactions
	# rindx_aq - index of aqueous-phase reactants
	# y_rind_aq - index of y relating to particle-phase reactants for reaction rate 
	# 	coefficient equation
	# y_arr_aq - index for matrix used to arrange concentrations of particle-phase 
	#	reactants, enabling calculation of reaction rate coefficients
	# rstoi_aq - stoichiometry of particle-phase reactants
	# pstoi_aq - stoichiometry of particle-phase products
	# nreac_aq - number of reactants per equation for aqueous-phase reactions
	# nprod_aq - number of products per equation for aqueous-phase reactions
	# call - flag for whether module called after ode_solv_wat or ode_solv
	# H2Oi - index for water
	# y - concentrations (molecules/cm3) after solver
	# self - reference to program
	# --------------------------------------------------------------------------

	# create new  file to store fluxes
	f = open('PyCHAM/ODE_solver_break_relevant_fluxes.txt', mode='w')
	f.write('##########################################################################################\n')
	f.write('#                                                                                        											 #\n')
	f.write('#    Copyright (C) 2018-2022 Simon O\'Meara : simon.omeara@manchester.ac.uk                  				 #\n')
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
	f.write('# File created at %s by ode_brk_err_mess\n' %(datetime.datetime.now()))

	if (call == 0): # if called after ode_solv_wat
		# record gas and particle-phase concentrations of water (molecules/cm3/s)
		f.write('\n')
		f.write('Gas and particle-phase concentrations of water prior to solver (molecules/cm3) with size bin numbers starting at 1: \n')
		for sbi in range(num_asb+1):
			if (sbi == 0):
				f.write(str('Gas: ' + str(y0[H2Oi]) + '\n'))
			if (sbi > 0):
				f.write(str('Size bin ' + str(sbi) + ': ' + str(y0[H2Oi+sbi*num_comp]) + '\n'))
		f.write('\n')
		f.write('Gas and particle-phase concentrations of water following solver (molecules/cm3) with size bin numbers starting at 1: \n')
		for sbi in range(num_asb+1):
			if (sbi == 0):
				f.write(str('Gas: ' + str(y[H2Oi]) + '\n'))
			if (sbi > 0):
				f.write(str('Size bin ' + str(sbi) + ': ' + str(y[H2Oi+sbi*num_comp]) + '\n'))

	if (eqn_num[0] > 0):# if gas-phase reactions present
		# gas-phase reactions -------------------------
		# empty array to hold relevant concentrations for
		# reaction rate coefficient calculation
		rrc_y = np.ones((rindx.shape[0]*rindx.shape[1]))
		rrc_y[y_arr] = y0[y_rind]
		rrc_y = rrc_y.reshape(rindx.shape[0], rindx.shape[1], order = 'C')
		# reaction rate (molecules/cm3/s) 
		rr = rrc[0:rindx.shape[0]]*((rrc_y**rstoi).prod(axis=1))
		rr = rr.reshape(-1, 1) # allow multiplication across multiple columns
		# loss of reactants
		reac_loss_rate = rr*rstoi # loss values (# molecules/cm3/s)
		# gain of products
		prod_gain_rate = rr*pstoi # gain values (# molecules/cm3/s)

		f.write('\n')
		f.write('Gas-phase reaction fluxes with equation numbers starting at 1 (# molecules/cm3/s):\n')

		for i in range(rindx.shape[0]):
			f.write(str('Eq. ' + str(i+1) + ' reac: ' + str(reac_loss_rate[i, 0:nreac[i]]) + '\n'))
			f.write(str('Eq. ' + str(i+1) + ' prod: ' + str(prod_gain_rate[i, 0:nprod[i]]) + '\n'))

	if (eqn_num[1] > 0):# if particle-phase reactions present

		# particle-phase reactions -------------------------
		# tile aqueous-phase reaction rate coefficients
		rr_aq = np.tile(rrc[rindx.shape[0]::], num_asb)
		# prepare for aqueous-phase concentrations
		rrc_y = np.ones((rindx_aq.shape[0]*rindx_aq.shape[1]))
		rrc_y[y_arr_aq] = y0[y_rind_aq]
		rrc_y = rrc_y.reshape(rindx_aq.shape[0], rindx_aq.shape[1], order = 'C')
		# reaction rate (molecules/cc/s) 
		rr = rr_aq*((rrc_y**rstoi_aq).prod(axis=1))
		rr = rr.reshape(-1, 1) # allow multiplication across multiple columns
		# loss of reactants
		reac_loss_rate = rr*rstoi_aq # prepare loss values
		# gain of products
		prod_gain_rate = rr*pstoi_aq # gain values (# molecules/cm3/s)
		
		f.write('\n')
		f.write('Particle-phase reaction fluxes with both size bin numbers and equation numbers starting at 1 (# molecules/cm3/s):\n')
		for sbi in range(num_asb): # loop through size bins
			for i in range(eqn_num[1]): # loop through equations
				f.write(str('size bin ' + str(sbi+1) + ', eq. ' + str(i+1) + ', reac: ' + str(reac_loss_rate[sbi*eqn_num[1]+i, 0:nreac_aq[i]]) + '\n'))
				f.write(str('size bin ' + str(sbi+1) + ', eq. ' + str(i+1)  + ', prod: ' + str(prod_gain_rate[sbi*eqn_num[1]+i, 0:nprod_aq[i]]) + '\n'))
	
	if (self.wall_on == 1): # include fluxes of trouble components to wall if wall is considered
		
		f.write('\n')
		f.write('Fluxes (molecules/cm3/s) of components with negative values output by ODE solver to (-) or from (+) wall\n')
		
		# concentration on wall (molecules/cc (air))
		Csit = y0[num_comp*(num_asb+1):num_comp*(num_asb+2)]
		# saturation vapour pressure on wall (molecules/cc (air))
		# note, just using the top rows of Psat and act_coeff
		# as do not need the repetitions over size bins
		if (Cw > 0.):
			Csit = Psat[0, neg_comp_indx]*(Csit[neg_comp_indx]/Cw)*act_coeff[0, neg_comp_indx]
			# rate of transfer (molecules/cc/s)
			dd_trouble = -1.*kw*(y0[neg_comp_indx]-Csit)
			
		if (Cw == 0.):
			dd_trouble = np.zeros((len(neg_comp_indx))) # rate of transfer (molecules/cc/s)

		for i in range(len(neg_comp_indx)): # loop through trouble components
			f.write(str(str(neg_names[i]) + ' : ' + str(dd_trouble[i]) + '\n'))	
		
	else:
		f.write('\n')
		f.write('Wall not turned on so no gas-wall partitioning fluxes reported\n')
		
	if (num_asb > 0): # if particle size bins present

		f.write('\n')
		f.write('Gas-particle partitioning fluxes (molecules/cm3/s) for each component with negative concentrations following call to ODE solver, where a negative flux represents loss from the gas-phase and positive represents gain to the gas-phase.  For each component, flux to the smallest size bin is the first value, and flux to the largest size bin is the final value.\n')
		
		# get particle-phase concentrations (molecules/cm3/s)
		ymat = (y0[num_comp:num_comp*(num_asb+1)]).reshape(num_asb, num_comp)
		# force all components in size bins with no particle to zero
		ymat[N_perbin[:, 0] == 0, :] = 0.
		# total particle-phase concentration per size bin (molecules/cc (air))		
		csum = ((ymat.sum(axis=1)-ymat[:, self.seedi].sum(axis=1))+((ymat[:, self.seedi]*core_diss).sum(axis=1)).reshape(-1)).reshape(-1, 1)
		# tile over components
		csum = np.tile(csum, [1, len(neg_comp_indx)])
		ymat = ymat[:, neg_comp_indx] # keep just the components with negative values output by ODE solver
		# container for gas-phase concentrations at particle surface (molecules/cm3)
		Csit = np.zeros((num_asb, len(neg_comp_indx)))
		# index of particles containing components
		isb = csum[:,0] > 0
		# mole fractions at particle surface
		Csit[isb, :] = (ymat[isb, :]/csum[isb, :])
		# filter just the components with negative concentrations following call to ODE solver
		Psat = Psat[:, neg_comp_indx]
		act_coeff = act_coeff[:, neg_comp_indx]
		kimt = kimt[:, neg_comp_indx]
		# gas-phase concentration of components at particle surface (molecules/cm3)
		Csit[isb, :] = Csit[isb, :]*Psat[isb, :]*kelv_fac[isb]*act_coeff[isb, :]
		# gas-particle partitioning rate (molecules/cm3/s)
		dd_trouble = -1.*kimt*(y0[neg_comp_indx].reshape(1, -1)-Csit)
		
		for i in range(len(neg_comp_indx)): # loop through trouble components
			f.write(str(str(neg_names[i]) + ': ' + str(dd_trouble[:, i]) + '\n'))		

	else:
		f.write('\n')
		f.write('Particle size bins not present so no gas-particle partitioning fluxes reported\n')


	f.close() # close file	

	return()
