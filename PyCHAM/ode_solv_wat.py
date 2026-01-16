########################################################################
#                                                                      #
# Copyright (C) 2018-2026                                              #
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
'''solution of ODEs for water gas-particle partitioning'''
# module to solve system of ordinary differential equations (ODEs) 
# using solve_ivp of Scipy 

import numpy as np
import scipy.sparse as SP
from scipy.integrate import solve_ivp

# define function
def ode_solv(y, integ_step, Cinfl_now,
	rowvals, colptrs, num_comp, num_sb,
	act_coeff, kelv_fac, kimt, num_asb,
	jac_mod_len, jac_part_hmf_indx, rw_indx, 
	N_perbin, jac_part_H2O_indx, H2Oi, self):

	# inputs: -------------------------------------
	# y - initial concentrations (molecules/cm^3)
	# integ_step - the maximum integration time step (s)
	# rrc - reaction rate coefficient
	# Cinfl_now - influx of components with constant influx 
	#		(# molecules/cm^3/s)
	# rowvals - row indices of Jacobian elements
	# colptrs - indices of  rowvals corresponding to each column of the
	# 	Jacobian
	# num_comp - number of components
	# num_sb - number of size bins
	# self.wall_on - flag saying whether to include wall partitioning
	# self.Psat - pure component saturation vapour pressures (# molecules/cm3)
	# act_coeff - activity coefficient of components
	# self.seedi - index of seed material
	# self.diss_wrtw - dissociation constant of components with respect to water
	# kelv_fac - kelvin factor for particles
	# kimt - mass transfer coefficient for gas-particle partitioning (s)
	# num_asb - number of actual size bins (excluding wall)
	# self.eqn_num - number of gas- and aqueous-phase reactions 
	# jac_mod_len - modification length due to high fraction of component(s)
	# 	in particle phase
	# jac_part_hmf_indx - index of Jacobian affected by water
	#	 in the particle phase
	# rw_indx - indices of rows affected by water in particle phase
	# N_perbin - number concentration of particles per size bin (#/cm3)
	# jac_part_H2O_indx - sparse Jacobian indices for the effect of
	#	particle-phase water on all other components
	# H2Oi - index for water
	# self - reference to program
	# ---------------------------------------------

	def dydt(t, y): # define the ODE(s)
		
		# inputs: ----------------
		# y - water concentrations (molecules/cm3), note when 
		# 	using scipy integrator solve_ivp, 
		#	this should have shape (number of elements, 1)
		# t - time interval to integrate over (s)
		# ---------------------------------------------
		
		# ensure y is correct shape
		if (y.shape[1] > 1):
			y = y[:, 0].reshape(-1, 1)
		
		# empty array to hold rate of change per component
		dd = np.zeros((y.shape[0], 1))
		
		# check for continuous gas-phase inputs for water
		if (self.H2Oin == 1):
			dd[0, 0] += self.Cinfl_H2O_now
		# check for continuous dilution
		if (self.dil_fac_H2Og_now > 0):
			# gas-phase water
			dd[0, 0] -= y[0, 0]*self.dil_fac_H2Og_now

		if (self.dil_fac_now > 0):
			# particle-phase water
			dd[1:num_sb+1, 0] -= y[1:num_sb+1, 0]*self.dil_fac_now
				
		# gas-particle partitioning-----------------
		
		# update the water particle-phase concentrations in the matrix 
		# of particle-phase concentrations for all components
		ymat[:, H2Oi] = y[1::, 0]
		
		# total particle-phase concentration per size bin (# molecules/cm^3 (air))
		csum = ((ymat*self.diss_wrtw).sum(axis=1)).reshape(-1, 1)
		
		isb = (csum[:, 0] != 0.) # indices of size bins with contents 
		
		if (any(isb)): # if particle-phase components present

			# mole fraction of water at particle surface
			Csit = (y[1::, 0][isb]/csum[isb, 0])
			# gas-phase concentration of water at particle surface
			# (# molecules/cm^3 (air))
			Csit = Csit*self.Psat[0:num_sb-self.wall_on, :][
				isb, H2Oi]*kelv_fac[isb, 0]*act_coeff[
				0:num_sb-self.wall_on, :][isb, H2Oi]

			# partitioning rate (# molecules/cm^3/s)
			dd_all = (kimt[0:num_asb, :][isb, H2Oi]*(y[0, 0]-Csit)).reshape(-1, 1)

			# only let gas-phase water change if not set to be constant concentration
			if (H2Oi not in self.conCindxn):
				dd[0, 0] -= sum(dd_all) # gas-phase change
			
			dd[1::, 0][isb] += (dd_all.flatten()) # particle change
			
		# force all components in size bins with no particle to zero
		if (num_asb > 0):
			dd[1:num_asb+1, 0][N_perbin[:, 0] == 0] = 0.
		# return to array, note that consistent with the solve_ivp manual, 
		# this ensures dd is
		# a vector rather than matrix, since y00 is a vector
		dd = dd.flatten()

		return (dd)

	def jac(t, y): # define the Jacobian
		
		# inputs: ----------------
		# y - concentrations (molecules/cm^3), note when using scipy integrator 
		#	solve_ivp, this should have shape (number of elements, 1)
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
		
		# update the particle-phase concentrations of water (molecules/cm^3 (air))
		ymat[:, H2Oi] = y[1::, 0]

		# total particle-phase concentration per size bin (# molecules/cm^3 (air))
		csum = ((ymat*self.diss_wrtw).sum(axis=1))
		
		
		# effect of particle-on-gas and particle-on-particle
		for isb in range(int(num_asb)): # size bin loop
			if (csum[isb] > 0): # if components present in this size bin
				# effect of gas on particle
				data[1+isb] += kimt[isb, H2Oi]
				# prepare for diagonal (component effect on itself)
				diag = kimt[isb, H2Oi]*self.Psat[0, H2Oi]*act_coeff[
					0, H2Oi]*kelv_fac[isb, 0]*(-(csum[isb]-y[
					1+isb, :])/(csum[isb]**2.)) 
				# implement to part_eff
				data[(num_asb+1)+isb*2] -= diag
				data[(num_asb+1)+isb*2+1] += diag

		# create Jacobian
		j = SP.csc_matrix((data, rowvals, colptrs))
		
		return(j)
	
	
	# set ODE solver (integration) tolerances
	atol = 1.e-4
	rtol = 1.e-5
	
	# isolate just the water concentrations
	y_w = y[H2Oi:num_comp*(num_asb+1):num_comp]

	# transform particle phase concentrations into size bins in 
	# rows and components in columns
	ymat = (y[num_comp:num_comp*(num_asb+1)]).reshape(num_asb, num_comp)
	# force all components in size bins with no particle to zero
	ymat[N_perbin[:, 0] == 0, :] = 0.
	
	# call on the ODE solver, note y contains the initial condition(s) (molecules/cm^3 (air)) 
	# and must be 1D even though y in dydt and jac has shape (number of elements, 1)
	# as stated in GMD paper, BDF method used as this known to deal with stiff problems well
	sol = solve_ivp(dydt, [0, integ_step], y_w, atol = atol, rtol = rtol, 
	method = 'Radau', t_eval = [integ_step], vectorized = True, jac = jac)
	
	# force all components in size bins with no particle to zero
	y_w = np.squeeze(sol.y)
	y_w = y_w.reshape(num_asb+1, 1)

	if (num_asb > 0):
		y_w[1:num_asb+1, 0][N_perbin[:, 0] == 0] = 0.
	# return to array
	y_w = y_w.flatten()

	# incorporate new water concentrations (molecules/cm^3)
	# implement new water gas- and particle-phase concentrations (molecules/cm^3 (air))
	y[H2Oi:num_comp*((num_sb-self.wall_on)+1):num_comp] = y_w

	# return concentration(s) and time(s) following integration
	return(y, sol.t)
