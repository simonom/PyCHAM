'''solution of ODEs for water gas-particle partitioning'''
# module to solve system of ordinary differential equations (ODEs) using solve_ivp of Scipy 

import numpy as np
import scipy.sparse as SP
from scipy.integrate import solve_ivp

# define function
def ode_solv(y, integ_step, rindx, pindx, rstoi, pstoi, 
	nreac, nprod, rrc, jac_stoi, njac, jac_den_indx, 
	jac_indx, Cinfl_now, y_arr, y_rind, uni_y_rind, 
	y_pind, uni_y_pind, reac_col, prod_col, 
	rstoi_flat, pstoi_flat, rr_arr, rr_arr_p,
	rowvals, colptrs, num_comp, num_sb,
	wall_on, Psat, Cw, act_coeff, kw, jac_wall_indx,
	seedi, core_diss, kelv_fac, kimt, num_asb,
	jac_part_indx,
	rindx_aq, pindx_aq, rstoi_aq, pstoi_aq,
	nreac_aq, nprod_aq, jac_stoi_aq, njac_aq, jac_den_indx_aq, jac_indx_aq,
	y_arr_aq, y_rind_aq, uni_y_rind_aq, y_pind_aq, uni_y_pind_aq,
	reac_col_aq, prod_col_aq, rstoi_flat_aq,
	pstoi_flat_aq, rr_arr_aq, rr_arr_p_aq, eqn_num, jac_mod_len,
	jac_part_hmf_indx, rw_indx, N_perbin, jac_part_H2O_indx, H2Oi):

	# inputs: -------------------------------------
	# y - initial concentrations (moleucles/cm3)
	# integ_step - the maximum integration time step (s)
	# rindx - index of reactants per equation
	# pindx - index of products per equation
	# rstoi - stoichiometry of reactants
	# pstoi - stoichiometry of products
	# nreac - number of reactants per equation
	# nprod - number of products per equation
	# rrc - reaction rate coefficient
	# jac_stoi - stoichiometries relevant to Jacobian
	# njac - number of Jacobian elements affected per equation
	# jac_den_indx - index of component denominators for Jacobian
	# jac_indx - index of Jacobian to place elements per equation (rows)
	# Cinfl_now - influx of components with constant influx 
	#		(molecules/cc/s)
	# y_arr - index for matrix used to arrange concentrations, 
	#	enabling calculation of reaction rate coefficients 
	# y_rind - index of y relating to reactants for reaction rate 
	# 	coefficient equation
	# uni_y_rind - unique index of reactants 
	# y_pind - index of y relating to products
	# uni_y_pind - unique index of products 
	# reac_col - column indices for sparse matrix of reaction losses
	# prod_col - column indices for sparse matrix of production gains
	# rstoi_flat - 1D array of reactant stoichiometries per equation
	# pstoi_flat - 1D array of product stoichiometries per equation
	# rr_arr - index for reaction rates to allow reactant loss
	# 	calculation
	# rr_arr_p - index for reaction rates to allow reactant loss
	# 	calculation
	# rowvals - row indices of Jacobian elements
	# colptrs - indices of  rowvals corresponding to each column of the
	# 	Jacobian
	# num_comp - number of components
	# num_sb - number of size bins
	# wall_on - flag saying whether to include wall partitioning
	# Psat - pure component saturation vapour pressures (molecules/cc)
	# Cw - effective absorbing mass concentration of wall (molecules/cc) 
	# act_coeff - activity coefficient of components
	# kw - mass transfer coefficient to wall (/s)
	# jac_wall_indx - index of inputs to Jacobian by wall partitioning
	# seedi - index of seed material
	# core_diss - dissociation constant of seed material
	# kelv_fac - kelvin factor for particles
	# kimt - mass transfer coefficient for gas-particle partitioning (s)
	# num_asb - number of actual size bins (excluding wall)
	# jac_part_indx - index for sparse Jacobian for particle influence 
	# eqn_num - number of gas- and aqueous-phase reactions 
	# jac_mod_len - modification length due to high fraction of component(s)
	# 	in particle phase
	# jac_part_hmf_indx - index of Jacobian affected by water
	#	 in the particle phase
	# rw_indx - indices of rows affected by water in particle phase
	# N_perbin - number concentration of particles per size bin (#/cc)
	# jac_part_H2O_indx - sparse Jacobian indices for the effect of
	#	particle-phase water on all other components
	# H2Oi - index for water
	# ---------------------------------------------

	def dydt(t, y): # define the ODE(s)
		
		# inputs: ----------------
		# y - water concentrations (molecules/cc), note when using scipy integrator solve_ivp, 
		#	this should have shape (number of elements, 1)
		# t - time interval to integrate over (s)
		# ---------------------------------------------
		
		# ensure y is correct shape
		if (y.shape[1] > 1):
			y = y[:, 0].reshape(-1, 1)
		# empty array to hold rate of change per component
		dd = np.zeros((y.shape[0], 1))
		
		# gas-particle partitioning-----------------
		
		# update the water particle-phase concentrations in the matrix of particle-phase 
		# concentrations for all components
		ymat[:, H2Oi] = y[1::, 0]
		
		
		# total particle-phase concentration per size bin (molecules/cc (air))
		csum = ((ymat.sum(axis=1)-ymat[:, seedi].sum(axis=1))+((ymat[:, seedi]*core_diss).sum(axis=1)).reshape(-1)).reshape(-1, 1)
		
		# size bins with contents 
		isb = (csum[:, 0] > 0.)
		
		if (any(isb)): # if particle-phase components present
			# mole fraction of water at particle surface
			Csit = (y[1::, 0][isb]/csum[isb, 0])
			# gas-phase concentration of water at particle surface (molecules/cc (air))
			Csit = Csit*Psat[isb, H2Oi]*kelv_fac[isb, 0]*act_coeff[isb, H2Oi]
			# partitioning rate (molecules/cc/s)
			dd_all = (kimt[isb, H2Oi]*(y[0, 0]-Csit)).reshape(-1, 1)
			dd[0, 0] -= sum(dd_all) # gas-phase change
			
			dd[1::, 0][isb] += (dd_all.flatten()) # particle change
			
		
		# force all components in size bins with no particle to zero
		if (num_asb > 0):
			dd[1:num_asb+1, 0][N_perbin[:, 0] == 0] = 0.
		# return to array, note that consistent with the solve_ivp manual, this ensures dd is
		# a vector rather than matrix, since y0 is a vector
		dd = dd.flatten()
		return (dd)

	def jac(t, y): # define the Jacobian
		
		# inputs: ----------------
		# y - concentrations (molecules/cc), note when using scipy integrator solve_ivp, this should have shape (number of elements, 1)
		# t - time interval to integrate over (s)
		# ---------------------------------------------
		
		# ensure y is correct shape
		if (y.ndim == 2):
			if (y.shape[1] > 1):
				y = y[:, 0].reshape(-1, 1)
		if (y.ndim <= 1):
			y = y.reshape(-1, 1)
		
		# elements of sparse Jacobian matrix
		data = np.zeros(((num_asb+1)**2-num_asb*(num_asb-1)))
		
		# the row and column indices for the sparse Jacobian
		rowvals = np.arange(num_asb+1)
		rowvals_app = np.zeros((2, num_asb))
		rowvals_app[1, :] = np.arange(1, num_asb+1)
		rowvals_app = rowvals_app.flatten(order='F')
		rowvals = np.concatenate((rowvals, rowvals_app)).astype('int')
		colptrs = np.array((0, num_asb+1))
		colptrs_app = np.arange(1, num_asb+1)*2
		colptrs = np.concatenate((colptrs, colptrs[-1]+colptrs_app))
		
		# gas-particle partitioning
		if (sum(N_perbin[:, 0]) > 0.): # if any particles present 
			data[0] = -kimt[:, H2Oi].sum(axis=0) # effect of gas on gas
		
		y[1::][N_perbin[:, 0] == 0] = 0 # ensure zero water where zero particles
		
		# update the particle-phase concentrations of water (molecules/cc (air))
		ymat[:, H2Oi] = y[1::, 0]
		# total particle-phase concentration per size bin (molecules/cc (air))
		csum = ymat.sum(axis=1)-ymat[:, seedi].sum(axis=1)+(ymat[:, seedi]*core_diss).sum(axis=1)
		
		# effect of particle-on-gas and particle-on-particle
		for isb in range(int(num_asb)): # size bin loop
			if (csum[isb] > 0): # if components present in this size bin
				# effect of gas on particle
				data[1+isb] += kimt[isb, H2Oi]
				# prepare for diagonal (component effect on itself)
				diag = kimt[isb, H2Oi]*Psat[0, H2Oi]*act_coeff[0, H2Oi]*kelv_fac[isb, 0]*(-(csum[isb]-y[1+isb, :])/(csum[isb]**2.)) 
				# implement to part_eff
				data[(num_asb+1)+isb*2] -= diag
				data[(num_asb+1)+isb*2+1] += diag

		# create Jacobian
		j = SP.csc_matrix((data, rowvals, colptrs))
		
		return(j)
	
	# set ODE solver tolerances
	atol = 0.001
	rtol = 0.0001
	
	y0 = y # remember initial component concentrations (molecules/cc (air))
	
	# isolate just the water concentrations
	y_w = y[H2Oi:num_comp*(num_asb+1):num_comp]
	
	# transform particle phase concentrations into size bins in rows and components in columns
	ymat = (y[num_comp:num_comp*(num_asb+1)]).reshape(num_asb, num_comp)
	# force all components in size bins with no particle to zero
	ymat[N_perbin[:, 0] == 0, :] = 0.
	
	# call on the ODE solver, note y contains the initial condition(s) (molecules/cc (air)) 
	# and must be 1D even though y in dydt and jac has shape (number of elements, 1)
	sol = solve_ivp(dydt, [0, integ_step], y_w, atol = atol, rtol = rtol, method = 'Radau', t_eval = [integ_step], vectorized = True, jac = jac)
	
	# force all components in size bins with no particle to zero
	y_w = np.squeeze(sol.y)
	y_w = y_w.reshape(num_asb+1, 1)
	if (num_asb > 0):
		y_w[1:num_asb+1, 0][N_perbin[:, 0] == 0] = 0.
	# return to array
	y_w = y_w.flatten()
	
	# implement new water gas- and particle-phase concentrations (molecules/cc (air))
	y0[H2Oi:num_comp*(num_asb+1):num_comp] = y_w
	# updated concentrations (molecules/cc (air))
	y[:] = y0[:]
	
	# return concentration(s) and time(s) following integration
	return(y, sol.t)
