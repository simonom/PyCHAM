##########################################################################################
#                                                                                        #
#    Copyright (C) 2018-2024 Simon O'Meara : simon.omeara@manchester.ac.uk               #
#                                                                                        #
#    All Rights Reserved.                                                                #
#    This file is part of PyCHAM                                                         #
#                                                                                        #
#    PyCHAM is free software: you can redistribute it and/or modify it under             #
#    the terms of the GNU General Public License as published by the Free Software       #
#    Foundation, either version 3 of the License, or (at your option) any later          #
#    version.                                                                            #
#                                                                                        #
#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT               #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              #
#    details.                                                                            #
#                                                                                        #
#    You should have received a copy of the GNU General Public License along with        #
#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                #
#                                                                                        #
##########################################################################################
'''generate the module to solve ODEs'''
# writes a module based on the supplied chemical scheme and user inputs,
# since the solved ODEs can represent gas-phase photochemistry,
# gas-particle partitioning and gas-wall partitioning, note use %f when
# writing floats from inputs

import datetime

# function to generate the ordinary differential equation (ODE)
# solver file
def ode_gen(int_tol, rowvals, num_comp, num_asb, testf, self):
	
	# inputs: ------------------------------------------------
	# self.con_infl_indx - indices of components with continuous influx
	# int_tol - integration tolerances
	# rowvals - indices of rows for Jacobian
	# self.wall_on - marker for whether to consider wall 
	# 	partitioning
	# num_comp - number of components in chemical scheme, plus 
	# any additional components needing consideration
	# num_asb - number of actual size bins (excluding wall)
	# testf - marker for whether in test mode or not
	# self.eqn_num - number of gas- and particle-phase reactions
	# self.dil_fac - fraction of chamber air extracted/s
	# self - reference to PyCHAM
	# -------------------------------------------------------
	
	# create new  file to store solver module
	f = open(self.PyCHAM_path + '/PyCHAM/ode_solv.py', mode='w')
	f.write('##########################################################################################\n')
	f.write('#                                                                                        #\n')
	f.write('#    Copyright (C) 2018-2024 Simon O\'Meara : simon.omeara@manchester.ac.uk               #\n')
	f.write('#                                                                                        #\n')
	f.write('#    All Rights Reserved.                                                                #\n')
	f.write('#    This file is part of PyCHAM                                                         #\n')
	f.write('#                                                                                        #\n')
	f.write('#    PyCHAM is free software: you can redistribute it and/or modify it under             #\n')
	f.write('#    the terms of the GNU General Public License as published by the Free Software       #\n')
	f.write('#    Foundation, either version 3 of the License, or (at your option) any later          #\n')
	f.write('#    version.                                                                            #\n')
	f.write('#                                                                                        #\n')
	f.write('#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT               #\n')
	f.write('#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       #\n')
	f.write('#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              #\n')
	f.write('#    details.                                                                            #\n')
	f.write('#                                                                                        #\n')
	f.write('#    You should have received a copy of the GNU General Public License along with        #\n')
	f.write('#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                #\n')
	f.write('#                                                                                        #\n')
	f.write('##########################################################################################\n')
	f.write('\'\'\'solution of ODEs, generated by eqn_pars.py\'\'\'\n')
	f.write('# module to solve system of ordinary differential equations (ODEs) using solve_ivp of Scipy \n')
	f.write('# File Created at %s\n' %(datetime.datetime.now()))	
	f.write('\n')
	f.write('import numpy as np\n')
	f.write('import scipy.sparse as SP\n')
	f.write('from scipy.integrate import solve_ivp\n')
	f.write('\n')	
	f.write('# define function\n')
	f.write('def ode_solv(y, integ_step, rrc, \n')
	f.write('	Cinfl_now, \n')
	f.write('	rowvals, colptrs, num_comp, num_sb,\n')
	f.write('	act_coeff, core_diss, kelv_fac, kimt, num_asb,\n')
	f.write('	jac_mod_len, jac_part_hmf_indx, rw_indx, N_perbin, jac_part_H2O_indx,\n')
	f.write('	H2Oi, self):\n')
	f.write('\n')
	f.write('	# inputs: -------------------------------------\n')
	f.write('	# y - initial concentrations (# molecules/cm3)\n')
	f.write('	# integ_step - the maximum integration time step (s)\n')
	f.write('	# self.rindx_g - index of reactants per equation\n')
	f.write('	# self.pindx_g - index of products per equation\n')
	f.write('	# self.rstoi_g - stoichiometry of reactants\n')
	f.write('	# self.pstoi_g - stoichiometry of products\n')
	f.write('	# self.nreac_g - number of reactants per equation\n')
	f.write('	# self.nprod_g - number of products per equation\n')
	f.write('	# rrc - reaction rate coefficient\n')
	f.write('	# self.jac_stoi_g - stoichiometries relevant to Jacobian\n')
	f.write('	# self.njac_g - number of Jacobian elements affected per equation\n')
	f.write('	# self.jac_den_indx_g - index of component denominators for Jacobian\n')
	f.write('	# self.jac_indx_g - index of Jacobian to place elements per equation (rows)\n')
	f.write('	# Cinfl_now - influx of components with continuous influx \n')
	f.write('	#		(# molecules/cm3/s)\n')
	f.write('	# self.y_arr_g - index for matrix used to arrange concentrations of gas-phase reactants, \n')
	f.write('	#	enabling calculation of reaction rate coefficients \n')
	f.write('	# self.y_rind_g - index of y relating to reactants for reaction rate \n')
	f.write('	# 	coefficient equation\n')
	f.write('	# self.uni_y_rind_g - unique index of reactants \n')
	f.write('	# self.y_pind_g - index of y relating to products\n')
	f.write('	# self.uni_y_pind_g - unique index of products \n')
	f.write('	# self.reac_col_g - column indices for sparse matrix of reaction losses\n')
	f.write('	# self.prod_col_g - column indices for sparse matrix of production gains\n')
	f.write('	# self.rstoi_flat_g - 1D array of reactant stoichiometries per equation\n')
	f.write('	# self.pstoi_flat_g - 1D array of product stoichiometries per equation\n')
	f.write('	# rr_arr - index for reaction rates to allow reactant loss\n')
	f.write('	# 	calculation\n')
	f.write('	# rr_arr_p - index for reaction rates to allow reactant loss\n')
	f.write('	# 	calculation\n')
	f.write('	# rowvals - row indices of Jacobian elements\n')
	f.write('	# colptrs - indices of  rowvals corresponding to each column of the\n') 
	f.write('	# 	Jacobian\n')
	f.write('	# num_comp - number of components\n')
	f.write('	# num_sb - number of size bins\n')
	f.write('	# self.wall_on - flag saying whether to include wall partitioning\n')
	f.write('	# self.Psat - pure component saturation vapour pressures (# molecules/cm3)\n')
	f.write('	# self.Cw - effective absorbing mass concentration of wall (# molecules/cm3) \n')
	f.write('	# act_coeff - activity coefficient of components\n')
	f.write('	# self.jac_wall_indxn - index of inputs to Jacobian by wall partitioning\n')
	f.write('	# self.seedi - index of seed material\n')
	f.write('	# core_diss - dissociation constant of seed material\n')
	f.write('	# kelv_fac - kelvin factor for particles\n')
	f.write('	# kimt - mass transfer coefficients for gas-particle partitioning (s) and gas-wall partitioning (/s)\n')
	f.write('	# num_asb - number of actual size bins (excluding wall)\n')
	f.write('	# self.jac_part_indxn - index for sparse Jacobian for particle influence \n')
	f.write('	# self.jac_extr_indx - index for sparse Jacobian for air extraction influence \n')
	f.write('	# self.rindx_aq - index of aqueous-phase reactants \n')
	f.write('	# self.eqn_num - number of gas- and aqueous-phase reactions \n')
	f.write('	# jac_mod_len - modification length due to high fraction of component(s)\n')
	f.write('	# 	in particle phase\n')
	f.write('	# jac_part_hmf_indx - index of Jacobian affected by water\n')
	f.write('	#	 in the particle phase\n')
	f.write('	# rw_indx - indices of rows affected by water in particle phase\n')
	f.write('	# N_perbin - number concentration of particles per size bin (#/cm3)\n')
	f.write('	# jac_part_H2O_indx - sparse Jacobian indices for the effect of\n')
	f.write('	#	particle-phase water on all other components\n')
	f.write('	# H2Oi - index for water\n')
	f.write('	# self.dil_fac - dilution factor for chamber (fraction of chamber air removed/s)\n')
	f.write('	# self.RO2_indx - index of organic peroxy radicals\n')
	f.write('	# self.comp_namelist - chemical scheme names of components\n')
	f.write('	# self.Psat_Pa - saturation vapour pressure of components (Pa) at starting\n')
	f.write('	#	temperature of chamber\n')
	f.write('	# self - reference to program \n')
	f.write('	# ---------------------------------------------\n')
	f.write('\n')
	
	# the module if needed for testing
	if (testf > 0):
		f.write('	\n')
		f.write('	# gas-particle partitioning-----------------\n')
		f.write('	# transform particle phase concentrations into\n')
		f.write('	# size bins in rows, components in columns\n')
		f.write('	ymat = (y[num_comp:num_comp*(num_asb+1), 0]).reshape(num_asb, num_comp)\n')
		f.write('	# total particle-phase concentration per size bin (molecules/cm3 (air))\n')
		f.write('	csum = ((ymat.sum(axis=1)-ymat[:, self.seedi])+((ymat[:, self.seedi]*core_diss))).reshape(-1, 1)\n')
		f.write('	# size bins with contents \n')
		f.write('	isb = (csum[:, 0]>0.)\n')
		f.write('	\n')
		f.write('	# container for gas-phase concentrations at particle surface\n')
		f.write('	Csit = np.zeros((num_asb, num_comp))\n')
		f.write('	# mole fraction of components at particle surface\n')
		f.write('	Csit[isb, :] = (ymat[isb, :]/csum[isb, :])\n')
		f.write('	\n')
		f.write('	return(Csit)\n')
		f.close() # close file
		return()


	# testing with 16 size bins and the MCM alpha-pinene chemical scheme
	# showed that using the vectorised Python code gave just 1 %
	# increase in wall clock time compared to using numba, and won't
	# give the gradual slow down in computation that arises with numba
	# when many integration time steps are set, furthermore it means that 
	# for fast integration systems, time isn't wasted on compilation, 
	# therefore we use the
	# vectorised form by default and the numba version is commented out 
	# below
	
	f.write('	def dydt(t, y): # define the ODE(s)\n')
	f.write('		\n')
	f.write('		# inputs: ----------------\n')
	f.write('		# y - concentrations (# molecules/cm3), note when using\n')
	f.write('		#	scipy integrator solve_ivp, this should have shape\n')
	f.write('		#	(number of elements, 1)\n')
	f.write('		# t - time interval to integrate over (s)\n')
	f.write('		# ---------------------------------------------\n')
	f.write('		\n')
	f.write('		# ensure y is correct shape\n')
	f.write('		if (y.shape[1] > 1):\n')
	f.write('			y = y[:, 0].reshape(-1, 1)\n')
	f.write('		# empty array to hold rate of change per component (this is the returned value from dydt)\n')
	f.write('		dd = np.zeros((y.shape[0], 1))\n')
	f.write('		\n')
	
	
	if (self.eqn_num[0] > 0): # if gas-phase reactions present
		f.write('		# gas-phase reactions -------------------------\n')
		f.write('		# empty array to hold relevant concentrations for\n')
		f.write('		# reaction rate coefficient calculation\n')
		f.write('		rrc_y = np.ones((self.rindx_g.shape[0]*self.rindx_g.shape[1]))\n')
		f.write('		rrc_y[self.y_arr_g] = y[self.y_rind_g, 0]\n')
		f.write('		rrc_y = rrc_y.reshape(self.rindx_g.shape[0], self.rindx_g.shape[1], order = \'C\')\n')
		f.write('		# reaction rate (molecules/cm3/s) \n')
		f.write('		rr = rrc[0:self.rindx_g.shape[0]]*((rrc_y**self.rstoi_g).prod(axis=1))\n')
		f.write('		# loss of reactants\n')
		f.write('		data = rr[self.rr_arr_g]*self.rstoi_flat_g # prepare loss values\n')
		f.write('		# convert to sparse matrix\n')
		f.write('		loss = SP.csc_matrix((data, self.y_rind_g, self.reac_col_g))\n')
		f.write('		# register loss of reactants\n')
		f.write('		dd[self.uni_y_rind_g, 0] -= np.array((loss.sum(axis = 1))[self.uni_y_rind_g])[:, 0]\n')
		f.write('		# gain of products\n')
		f.write('		data = rr[self.rr_arr_p_g]*self.pstoi_flat_g # prepare loss values\n')
		f.write('		# convert to sparse matrix\n')
		f.write('		loss = SP.csc_matrix((data, self.y_pind_g, self.prod_col_g))\n')
		f.write('		# register gain of products\n')
		f.write('		dd[self.uni_y_pind_g, 0] += np.array((loss.sum(axis = 1))[self.uni_y_pind_g])[:, 0]\n')
		f.write('		\n')

	if (self.eqn_num[1] > 0): # if particle-phase reactions present
		f.write('		# particle-phase reactions -------------------------\n')
		f.write('		\n')
		f.write('		# empty array to hold relevant concentrations for\n')
		f.write('		# reaction rate coefficient calculation\n')
		f.write('		# tile aqueous-phase reaction rate coefficients\n')
		f.write('		rr_aq = np.tile(rrc[self.rindx_g.shape[0]::], num_asb)\n')
		f.write('		# prepare aqueous-phase concentrations\n')
		f.write('		rrc_y = np.ones((self.rindx_aq.shape[0]*self.rindx_aq.shape[1]))\n')
		f.write('		rrc_y[self.y_arr_aq] = y[self.y_rind_aq, 0]\n')
		f.write('		rrc_y = rrc_y.reshape(self.rindx_aq.shape[0], self.rindx_aq.shape[1], order = \'C\')\n')
		f.write('		# reaction rate (# molecules/cm3/s) \n')
		f.write('		rr = rr_aq*((rrc_y**self.rstoi_aq).prod(axis=1))\n')
		f.write('		# loss of reactants\n')
		f.write('		data = rr[self.rr_arr_aq]*self.rstoi_flat_aq # prepare loss values\n')
		f.write('		# convert to sparse matrix\n')
		f.write('		loss = SP.csc_matrix((data, self.y_rind_aq, self.reac_col_aq))\n')
		f.write('		# register loss of reactants\n')
		f.write('		dd[self.uni_y_rind_aq, 0] -= np.array((loss.sum(axis = 1))[self.uni_y_rind_aq])[:, 0]\n')
		f.write('		# gain of products\n')
		f.write('		data = rr[self.rr_arr_p_aq]*self.pstoi_flat_aq # prepare loss values\n')
		f.write('		# convert to sparse matrix\n')
		f.write('		loss = SP.csc_matrix((data, self.y_pind_aq, self.prod_col_aq))\n')
		f.write('		# register gain of products\n')
		f.write('		dd[self.uni_y_pind_aq, 0] += np.array((loss.sum(axis = 1))[self.uni_y_pind_aq])[:, 0]\n')
		f.write('		\n')

	if (self.eqn_num[2] > 0): # if surface reactions present
		f.write('		# surface-phase reactions -------------------------\n')
		f.write('		\n')	
		f.write('		# empty array to hold relevant concentrations for\n')
		f.write('		# reaction rate coefficient calculation\n')
		f.write('		# tile surface-phase reaction rate coefficients\n')
		f.write('		rr_su = np.tile(rrc[-self.rindx_su.shape[0]::], self.wall_on)\n')
		f.write('		# prepare surface-phase concentrations\n')
		f.write('		rrc_y = np.ones((self.rindx_su.shape[0]*self.rindx_su.shape[1]))\n')
		f.write('		rrc_y[self.y_arr_su] = y[self.y_rind_su, 0]\n')
		f.write('		rrc_y = rrc_y.reshape(self.rindx_su.shape[0], self.rindx_su.shape[1], order = \'C\')\n')
		f.write('		# reaction rate (# molecules/cm3/s) \n')
		f.write('		rr = rr_su*((rrc_y**self.rstoi_su).prod(axis=1))\n')
		f.write('		# loss of reactants\n')
		f.write('		data = rr[self.rr_arr_su]*self.rstoi_flat_su # prepare loss values\n')
		f.write('		# flatten data in column-major order (rows are reaction) (columns are surfaces)\n')
		f.write('		data = data.flatten(order = \'F\')\n')
		f.write('		# convert to sparse matrix\n')
		f.write('		loss = SP.csc_matrix((data, self.y_rind_su, self.reac_col_su))\n')
		f.write('		# register loss of reactants\n')
		f.write('		dd[self.uni_y_rind_su, 0] -= np.array((loss.sum(axis = 1))[self.uni_y_rind_su])[:, 0]\n')
		f.write('		# gain of products\n')
		f.write('		data = rr[self.rr_arr_p_su]*self.pstoi_flat_su # prepare loss values\n')
		f.write('		# flatten data in column-major order (rows are reaction) (columns are surfaces)\n')
		f.write('		data = data.flatten(order = \'F\')\n')
		f.write('		# convert to sparse matrix\n')
		f.write('		loss = SP.csc_matrix((data, self.y_pind_su, self.prod_col_su))\n')
		f.write('		# register gain of products\n')
		f.write('		dd[self.uni_y_pind_su, 0] += np.array((loss.sum(axis = 1))[self.uni_y_pind_su])[:, 0]\n')
		f.write('		\n')
	
	if (len(self.con_infl_indx) > 0): # if a component has a continuous gas-phase influx
		
		f.write('		# account for components with continuous gas-phase influx\n')	
		f.write('		dd[[self.con_infl_indx], 0] += Cinfl_now[:, 0]\n')
		if ('H2O' in self.con_infl_nam): 
			f.write('		if self.odsw_flag == 0: # if water not solved separately\n')
			f.write('			dd[H2Oi, 0] += self.Cinfl_H2O_now\n')

	if (any(self.obs_comp_i)): # if a component has a fixed concentration
		
		f.write('		# account for components with fixed gas-phase concentration\n')	
		f.write('		dd[[self.obs_comp_i], 0] = 0.\n')
	
	if (any(self.dil_fac > 0.)): # if chamber air being extracted
		f.write('		# account for continuous extraction of chamber air\n')
		f.write('		# index for estimating dilution factors \n')
		f.write('		df_indx = np.ones((dd.shape[0])).astype(\'int\') \n')
		f.write('		if (self.odsw_flag == 0): # if water solver not used \n')
		f.write('			dd[H2Oi::num_comp, 0] -= y[H2Oi::num_comp, 0]*self.dil_fac_H2O_now\n')
		f.write('		# water diluted either in water solver or above \n')
		f.write('		df_indx[H2Oi::num_comp] = 0 \n')
		f.write('		# cannot dilute what is on wall \n')
		f.write('		df_indx[num_comp*(num_sb-self.wall_on+1)::] = 0  \n')
		f.write('		# in case particle-phase components should not be diluted, \n')
		f.write('		# e.g. when observations already account for dilution, as in obs_file_open \n')
		f.write('		if (self.pp_dil == 0):\n')
		f.write('			df_indx[num_comp:num_comp*(num_sb-self.wall_on+1)] = 0  \n')
		f.write('		df_indx = (df_indx == 1) # transform to Boolean array \n')
		f.write('		dd[df_indx, 0] -= y[df_indx, 0]*self.dil_fac_now\n')
		f.write('		\n')
		
	# note the following needs two indents (as for the reaction section), so that it
	# sits within the dydt function
	if (num_asb > 0 or self.wall_on > 0): # include gas-particle partitioning in ode solver
		if (num_asb > 0 and self.wall_on > 0): 
			f.write('		# gas-particle and gas-wall partitioning-----------------\n')
		if (num_asb > 0 and self.wall_on == 0):
			f.write('	\n')
			f.write('		# gas-particle partitioning-----------------\n')
		if (num_asb == 0 and self.wall_on > 0): 
			f.write('		# gas-wall partitioning-----------------\n')
			
		f.write('		# transform component concentrations in particles and walls\n')
		f.write('		# into size bins in rows, components in columns\n')
		f.write('		ymat = (y[num_comp::, 0]).reshape(num_sb, num_comp)\n')
		f.write('		\n')
		if (num_asb > 0):
			f.write('		# for particles, force all components in bins with no particle to zero\n')
			f.write('		ymat[0:num_asb, :][N_perbin[:, 0] == 0, :] = 0\n')
			f.write('		\n')
			f.write('		# for particles, calculate total particle-phase concentration per size bin (# molecules/cm3 (air))\n')
			f.write('		csum = ((ymat[0:num_asb, :].sum(axis=1)-ymat[0:num_asb, self.seedi].sum(axis=1))+((ymat[0:num_asb, self.seedi]*core_diss).sum(axis=1)).reshape(-1)).reshape(-1, 1)\n')
			f.write('		# tile total particle-phase concentration over components (# molecules/cm3 (air))\n')
			f.write('		csum = np.tile(csum, [1, num_comp])\n')
			if (self.wall_on > 0):
				f.write('		# concatenate wall bin total concentrations to total particle-phase concentration (# molecules/cm3)\n')
				f.write('		csum = np.concatenate((csum, self.Cw), axis=0)\n')
				f.write('		\n')
		if (num_asb == 0 and self.wall_on > 0):
			f.write('		# rename wall bin total concentrations (# molecules/cm3)\n')
			f.write('		csum = self.Cw\n')
			f.write('		\n')
		if (num_asb > 0):
			f.write('		# size bins with contents\n')
			f.write('		isb = (csum[0:num_asb, 0] > 0.)\n')
			f.write('		\n')
		if (self.wall_on > 0):
			f.write('		# wall bins with contents\n')
			f.write('		wsb = (self.Cw[:, 0] > 0.)\n')
			f.write('		\n')
		if (num_asb > 0 and self.wall_on > 0):
			f.write('		# container for gas-phase concentrations at particle surface and at wall surface\n')
		if (num_asb > 0 and self.wall_on == 0):
			f.write('		# container for gas-phase concentrations at particle surface\n')
		if (num_asb == 0 and self.wall_on > 0):
			f.write('		# container for gas-phase concentrations at wall surface\n')
		f.write('		Csit = np.zeros((num_sb, num_comp))\n')
		f.write('		\n')
		if (num_asb > 0):
			f.write('		# mole fractions of components at particle surface\n')
			f.write('		Csit[0:num_asb, :][isb, :] = (ymat[0:num_asb, :][isb, :]/csum[0:num_asb, :][isb, :])\n')
		if (self.wall_on > 0):
			f.write('		# mole fraction of components on walls, note that Cw included in csum above\n')
			f.write('		Csit[num_asb::, :][wsb, :] = (ymat[num_asb::, :][wsb, :]/csum[num_asb::, :][wsb, :])\n')
		f.write('		\n')
		if (num_asb > 0):
			#f.write('		if any(isb):\n')
			f.write('		# gas-phase concentration of components at\n')
			f.write('		# particle surface (# molecules/cm3 (air))\n')
			f.write('		Csit[0:num_asb, :][isb, :] = Csit[0:num_asb, :][isb, :]*self.Psat[0:num_asb, :][isb, :]*kelv_fac[isb]*act_coeff[0:num_asb, :][isb, :]\n')	
			f.write('		# partitioning rate (# molecules/cm3/s)\n')
			f.write('		dd_all = kimt[0:num_asb, :]*(y[0:num_comp, 0].reshape(1, -1)-Csit[0:num_asb, :])\n')
			f.write('		# gas-phase change\n')
			f.write('		dd[0:num_comp, 0] -= dd_all.sum(axis=0)\n')
			f.write('		# particle change\n')
			f.write('		dd[num_comp:num_comp*(num_asb+1), 0] += (dd_all.flatten())\n')
			f.write('		\n')
		if (self.wall_on > 0):
			f.write('		if any(wsb):\n')
			f.write('			# gas-phase concentration of components at\n')
			f.write('			# wall surface (# molecules/cm3 (air))\n')
			f.write('			Csit[num_asb::, :][wsb, :] = Csit[num_asb::, :][wsb, :]*self.Psat[num_asb::, :][wsb, :]*act_coeff[num_asb::, :][wsb, :]\n')	
			f.write('			# partitioning rate (# molecules/cm3/s)\n')
			f.write('			dd_all = kimt[num_asb::, :]*(y[0:num_comp, 0].reshape(1, -1)-Csit[num_asb::, :])\n')
			f.write('			# gas-phase change (summed over all wall bins)\n')
			f.write('			dd[0:num_comp, 0] -= dd_all.sum(axis=0)\n')
			f.write('			# wall change\n')
			f.write('			dd[num_comp*(num_asb+1)::, 0] += (dd_all.flatten())\n')
			f.write('		\n')
		
	
	f.write('		dd = (dd[:, 0]).reshape(num_sb+1, num_comp)\n')
	f.write('		# force all components in size bins with no particle to zero\n')
	f.write('		if (num_asb > 0):\n')
	f.write('			dd[1:num_asb+1, :][N_perbin[:, 0] == 0, :] = 0\n')
	f.write('		# return to array, note that consistent with the solve_ivp manual, this ensures dd is\n')
	f.write('		# a vector rather than matrix, since y0 is a vector\n')
	f.write('		dd = dd.flatten()\n')
	f.write('		nzindx = dd != 0.\n')
	f.write('		return (dd)\n')
	f.write('\n')
	
	# set the Jacobian
	f.write('	def jac(t, y): # define the Jacobian\n')
	f.write('		\n')
	f.write('		# inputs: ----------------\n')
	f.write('		# y - concentrations (# molecules/cm3), note when using scipy integrator solve_ivp, this should have shape (number of elements, 1)\n')
	f.write('		# t - time interval to integrate over (s)\n')
	f.write('		# ---------------------------------------------\n')
	f.write('		\n')
	f.write('		# ensure y is correct shape\n')
	f.write('		if (y.ndim == 2):\n')
	f.write('			if (y.shape[1] > 1):\n')
	f.write('				y = y[:, 0].reshape(-1, 1)\n')
	f.write('		if (y.ndim <= 1):\n')
	f.write('			y = y.reshape(-1, 1)\n')
	f.write('		\n')
	
	f.write('		# elements of sparse Jacobian matrix\n')
	if (num_asb > 0): # include any particle-phase modifiers
		f.write('		data = np.zeros((%s+jac_mod_len))\n' %len(rowvals))
	else: # don't include any particle-phase modifiers
		f.write('		data = np.zeros((%s))\n' %len(rowvals))
	f.write('		\n')
	
	if (self.eqn_num[0] > 0): # if gas-phase reactions present
		f.write('		for i in range(self.rindx_g.shape[0]): # gas-phase reaction loop\n')
		f.write('			# reaction rate (# molecules/cm3/s)\n')
		f.write('			rr = rrc[i]*(y[self.rindx_g[i, 0:self.nreac_g[i]], 0].prod())\n')
		f.write('			# prepare Jacobian inputs\n')
		f.write('			jac_coeff = np.zeros((self.njac_g[i, 0]))\n')
		f.write('			# only fill Jacobian if reaction rate sufficient\n')
		f.write('			if (rr != 0.):\n')
		f.write('				jac_coeff = (rr*(self.jac_stoi_g[i, 0:self.njac_g[i, 0]])/\n')
		f.write('				(y[self.jac_den_indx_g[i, 0:self.njac_g[i, 0]], 0]))\n')
		f.write('			data[self.jac_indx_g[i, 0:self.njac_g[i, 0]]] += jac_coeff\n')
		f.write('		\n')
	
	if (self.eqn_num[1] > 0): # if particle-phase reactions present
		f.write('		n_aqr = self.nreac_aq.shape[0] # number of aqueous-phase reactions \n')
		f.write('		aqi = 0 # aqueous-phase reaction counter\n')
		f.write('		\n')
		f.write('		for i in range(self.rindx_g.shape[0], rrc.shape[0]): # aqueous-phase reaction loop\n')
		f.write('			# reaction rate (molecules/cm3/s)\n')
		f.write('			rr = rrc[i]*(y[self.rindx_aq[aqi::n_aqr, 0:self.nreac_aq[aqi]], 0].prod(axis=1))\n')
		f.write('			# spread along affected components\n')
		f.write('			rr = rr.reshape(-1, 1)\n')
		f.write('			rr = (np.tile(rr, int(self.njac_aq[aqi, 0]/(num_sb-self.wall_on)))).flatten(order=\'C\')\n')
		f.write('			# prepare Jacobian inputs\n')
		f.write('			jac_coeff = np.zeros((self.njac_aq[aqi, 0]))\n')
		f.write('			nzi = (rr != 0)\n')
		f.write('			jac_coeff[nzi] = (rr[nzi]*((self.jac_stoi_aq[aqi, 0:self.njac_aq[aqi, 0]])[nzi])/\n')
		f.write('				((y[self.jac_den_indx_aq[aqi, 0:self.njac_aq[aqi, 0]], 0])[nzi]))\n')
		f.write('			# stack size bins\n')
		f.write('			jac_coeff = jac_coeff.reshape(int(num_sb-self.wall_on), int(self.njac_aq[aqi, 0]/(num_sb-self.wall_on)), order=\'C\')\n')
		f.write('			data[self.jac_indx_aq[aqi::n_aqr, 0:(int(self.njac_aq[aqi, 0]/(num_sb-self.wall_on)))]] += jac_coeff\n')
		f.write('			\n')
		f.write('			aqi += 1 # keep count on aqueous-phase reactions \n')
	f.write('		\n')
	
	if (num_asb > 0): # include gas-particle partitioning in Jacobian
		f.write('		# gas-particle partitioning\n')
		f.write('		part_eff = np.zeros((%s))\n' %((num_comp)*(num_asb+1)+((num_comp)*(num_asb*2))))
		f.write('		if (sum(N_perbin[:, 0]) > 0.): # if any particles present \n')
		f.write('			part_eff[0:%s:%s] = -kimt[0:num_asb, :].sum(axis=0) # effect of gas on gas\n' %(num_comp*(num_asb+1), (num_asb+1)))
		f.write('		\n')
		f.write('		# empty array for any particle-on-gas and particle-on-particle effects on water in the particle-phase for rows of Jacobian\n')
		f.write('		part_eff_rw = np.zeros((len(jac_part_hmf_indx)))\n')
		f.write('		# empty array for any particle-on-gas and particle-on-particle effects of water in the particle-phase on non-water components in the particle-phase for columns of Jacobian\n')
		f.write('		part_eff_cl = np.zeros((len(jac_part_H2O_indx)))\n')
		f.write('		# starting index for jacobian row inputs for effect on water\n')
		f.write('		sti_rw = 0 \n')
		f.write('		\n')
		f.write('		# transform particle phase concentrations into\n')
		f.write('		# size bins in rows, components in columns\n')
		f.write('		ymat = (y[num_comp:num_comp*(num_asb+1), 0]).reshape(num_asb, num_comp)\n')
		f.write('		ymat[N_perbin[:, 0] == 0, :] = 0 # ensure zero components where zero particles\n')
		f.write('		# total particle-phase concentration per size bin (molecules/cm3 (air))\n')
		f.write('		csum = ymat.sum(axis=1)-ymat[:, self.seedi].sum(axis=1)+(ymat[:, self.seedi]*core_diss).sum(axis=1)\n')
		f.write('		\n')
		f.write('		# effect of particle on gas\n')
		f.write('		for isb in range(int(num_asb)): # size bin loop\n')
		f.write('			if (csum[isb] > 0): # if components present in this size bin\n')
		f.write('				# effect of gas on particle\n')
		f.write('				part_eff[1+isb:num_comp*(num_asb+1):num_asb+1] = +kimt[isb, :]\n')
		f.write('				# start index\n')
		f.write('				sti = int((num_asb+1)*num_comp+isb*(num_comp*2))\n')
		f.write('				# diagonal index\n')
		f.write('				diag_indxg = sti+np.arange(0, num_comp*2, 2).astype(\'int\')\n')
		f.write('				diag_indxp = sti+np.arange(1, num_comp*2, 2).astype(\'int\')\n')
		f.write('				# prepare for diagonal (component effect on itself)\n')
		f.write('				diag = kimt[isb, :]*self.Psat[0, :]*act_coeff[0, :]*kelv_fac[isb, 0]*(-(csum[isb]-ymat[isb, :])/(csum[isb]**2.)) \n')
		f.write('				# implement to part_eff\n')
		f.write('				part_eff[diag_indxg] -= diag\n')
		f.write('				part_eff[diag_indxp] += diag\n')
		f.write('				\n')
		f.write('				if (rw_indx[isb] > -1): # if water in this size bin \n') 
		f.write('					# prepare for row(s) (particle-phase non-water component effects on water in particle phase)\n')
		f.write('					rw = kimt[isb, rw_indx[isb]]*self.Psat[0, rw_indx[isb]]*act_coeff[0, rw_indx[isb]]*kelv_fac[isb, 0]*(-(-ymat[isb, rw_indx[isb]])/(csum[isb]**2.)) \n')
		f.write('					# indices\n')
		f.write('					indxg = sti_rw+np.arange(0, ((num_comp-1)*2), 2).astype(\'int\')\n')
		f.write('					indxp = sti_rw+np.arange(1, ((num_comp-1)*2), 2).astype(\'int\')\n')
		f.write('					# implement to part_eff_rw\n')
		f.write('					part_eff_rw[indxg] -= rw\n')
		f.write('					part_eff_rw[indxp] += rw\n')
		f.write('					\n')
		f.write('					# prepare for column(s) (particle-phase water effect on non-water in particle phase)\n')
		f.write('					#cl = kimt[isb, :]*self.Psat[0, :]*act_coeff[0, :]*kelv_fac[isb, 0]*(-(-ymat[isb, :])/(csum[isb]**2.))\n')
		f.write('					#cl = np.zeros((num_comp))\n')
		f.write('					# remove water\n')
		f.write('					#cl = np.concatenate((cl[0:H2Oi], cl[H2Oi+1::]))\n')
		f.write('					#indxg = sti_rw+np.arange(0, (num_comp-1)).astype(\'int\')\n')
		f.write('					#indxp = sti_rw+np.arange((num_comp-1), (num_comp-1)*2).astype(\'int\')\n')
		f.write('					# implement to part_eff_cl\n')
		f.write('					#part_eff_cl[indxg] -= cl\n')
		f.write('					#part_eff_cl[indxp] += cl\n')
		f.write('					\n')
		f.write('					# starting index update\n')
		f.write('					sti_rw += (num_comp-1)*2\n')
		f.write('		\n')
		f.write('		data[self.jac_part_indxn] += part_eff # diagonal\n')
		f.write('		data[jac_part_hmf_indx] += part_eff_rw # rows\n')
		f.write('		#data[jac_part_H2O_indx] += part_eff_cl # columns\n')
		f.write('		\n')
		
	if (self.wall_on > 0): # include gas-wall partitioning in Jacobian
		f.write('		wsb = 0 # count on wall bins\n')
		# holder for wall effect - note that 1st term is gas-on-gas, 2nd term is gas-on-wall, 3rd term is wall-on-gas and 4th term is wall-on-wall
		f.write('		# holder for wall effect\n')
		f.write('		wall_eff = np.zeros((%s))\n' %((num_comp)+(num_comp*self.wall_on)+(num_comp*self.wall_on)+(num_comp*self.wall_on)))
		f.write('		# effect of gas on gas \n')
		f.write('		wall_eff[0:%s:%s] = -np.sum(kimt[num_asb::, :], axis=0) \n' %(num_comp*(self.wall_on+1), self.wall_on+1))
		f.write('		for wsb in range(int(self.wall_on)): # wall bin loop\n')
		f.write('			if (self.Cw[wsb, 0] > 0.):\n')
		f.write('				# effect of gas on wall \n')
		f.write('				wall_eff[wsb+1:%s:%s] = +kimt[num_asb+wsb, :] \n' %(num_comp*(self.wall_on+1), self.wall_on+1))
		f.write('				# effect of wall on gas\n')
		f.write('				wall_eff[wsb*2*num_comp+num_comp*(self.wall_on+1):num_comp*(self.wall_on+1)+(wsb+1)*2*num_comp:2] = +kimt[num_asb::, :][wsb, :]*(self.Psat[num_asb::, :][wsb, :]*act_coeff[num_asb::, :][wsb, :]/self.Cw[wsb, :]) \n')
		f.write('				# effect of wall on wall\n')
		f.write('				wall_eff[wsb*2*num_comp+num_comp*(self.wall_on+1)+1:num_comp*(self.wall_on+1)+(wsb+1)*2*num_comp:2] = -kimt[num_asb::, :][wsb, :]*(self.Psat[num_asb::, :][wsb, :]*act_coeff[num_asb::, :][wsb, :]/self.Cw[wsb, :]) \n')
		f.write('		data[self.jac_wall_indxn] += wall_eff\n')
		f.write('		\n')
	if (any(self.dil_fac > 0)): # include extraction of chamber air in ode Jacobian
		f.write('		data[self.jac_extr_indx] -= 1.*self.dil_fac_now\n')
		f.write('		\n')
	# note that continuous influx should only be included in Jacobian if the influx is dependent
	# on the concentration of a component(s), otherwise it differentiates to zero when the 
	#ordinary differential equation is expressed
	#if (len(self.con_infl_indx) > 0.): # include continuous influx of gases
		#f.write('		Cinfl_gr_zero = Cinfl_now[:, 0] > 0. # influxes over zero\n')
		#f.write('		data[self.jac_cont_infl_indx][Cinfl_gr_zero] += Cinfl_now[Cinfl_gr_zero, 0]/(y[self.con_infl_indx, 0][Cinfl_gr_zero])\n')	
	f.write('		# create Jacobian\n')
	f.write('		j = SP.csc_matrix((data, rowvals, colptrs))\n')
	f.write('		\n')
	f.write('		return(j)\n')
	f.write('	\n')
	f.write('	# set ODE solver tolerances\n')
	f.write('	atol = %s\n'%int_tol[0])
	f.write('	rtol = %s\n'%int_tol[1])
	f.write('	self.ode_cnt = 0\n')
	f.write('	\n')
	f.write('	# check for underflow issues\n')
	f.write('	# reaction rate coefficient calculation\n')
	f.write('	#rrc_y = np.ones((self.rindx_g.shape[0]*self.rindx_g.shape[1]))\n')
	f.write('	#rrc_y[self.y_arr_g] = y[self.rindx_g]\n')
	f.write('	#rrc_y = rrc_y.reshape(self.rindx_g.shape[0], self.rindx_g.shape[1], order = \'C\')\n')
	f.write('	# reaction rate coefficient zeroed wherever product of reactant concentrations is zero (including where underflow causes zero, thereby preventing underflows breaking the solver which appears to be an issue on less powerful machines such as HP Spectre Folio) (/s) \n')
	f.write('	#rrc[((rrc_y**self.rstoi_g).prod(axis=1)) == 0.0] = 0.\n')
	f.write('	\n')
	f.write('	# call on the ODE solver, note y contains the initial condition(s) (molecules/cm3 (air)) and must be 1D even though y in dydt and jac has shape (number of elements, 1)\n')
	
	f.write('	sol = solve_ivp(dydt, [0, integ_step], y, atol = atol, rtol = rtol, method = \'BDF\', t_eval = [integ_step], vectorized = True, jac = jac)\n')
	f.write('	\n')
	f.write('	if (sol.status == -1): # if integration step failed, then we want to reduce the time step and try again \n')
	f.write('		y[0] = -1.e6\n')
	f.write('	else:\n')	
	f.write('		# force all components in size bins with no particle to zero\n')
	
	f.write('		y = np.squeeze(sol.y)\n')
	
	f.write('		y = y.reshape(num_sb+1, num_comp)\n')
	f.write('		if (num_asb > 0):\n')
	f.write('			y[1:num_asb+1, :][N_perbin[:, 0] == 0, :] = 0\n')
	f.write('		# return to array\n')
	f.write('		y = y.flatten()\n')
	f.write('		\n')
	f.write('	# return concentration(s) and time(s) following integration\n')
	f.write('	return(y, sol.t)\n')
	f.close() # close file
