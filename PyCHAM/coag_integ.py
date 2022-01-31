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
'''module to solve coagulation through integration'''
# using the differential equation for coagulation given in Eq. 15.2
# of Jacobson (2005), the new number size distribution is solved 

import numpy as np
import scipy as si
from scipy.integrate import solve_ivp

def integ_coag(N_perbin, T, x, xb, nsb, integ_step):
	
	# inputs: -------------------------
	# N_perbin - initial concentration of particles per size bin (# particles/cm3)
	# T - chamber temperature (K)
	# x - radius of particles per size bin (um)
	# xb - radius of particle size bin bounds (um)
	# nsb - number of size bins
	# integ_step - integration step (s)
	# -----------------------------------
	
	def dNdt(t, N_perbin): # define ordinary differential equation function for number-conserving solution
		
		N_perbinr = np.tile(N_perbin.reshape(1, -1), [nsb, 1])
		N_perbinc = np.tile(N_perbin.reshape(-1, 1), [1, nsb])
	
		# empty array for rate of change (# particles/cm3.s)
		dd = np.zeros((nsb))
		# size bin loop
		for isb in range(nsb):
	
			# identify pairs that coagulate to give particle in this size bin
			pindx = (Vcoag >= Vb[isb])*(Vcoag < Vb[isb+1])
			# rate of change (# particles/cm3.s)
			dd[isb] = 0.5*(np.sum(np.sum(Beta[pindx]*N_perbinr[pindx]*N_perbinc[pindx])))-(N_perbin[isb]*np.sum(Beta[isb, :]*N_perbin))
		
		return(dd)
		
		
	Vb = (4./3.)*np.pi*(xb**3.) # volume bound per size bin (um3)
	# volume (um3) per size bin
	V = (4./3.)*np.pi*(x**3.)
	# empty array for volume of combined pairs (um3)
	Vcoag = np.tile(V.reshape(1, -1), [nsb, 1])+np.tile(V.reshape(-1, 1), [1, nsb])
	
	
	# suggested kernel calculation from Eq. 15.16 Jacobson 2005 ------
	# dynamic viscosity of air (g/m.s) (Eq. 4.54 Jacobson 2005)
	# note the conversion of the universal gas constant from 
	# kg.m2/s2.mol.K to g.m2/s2.mol.K
	na = 5./(16.*si.constants.N_A*3.673e-10**2.)*((28.966*(si.constants.R*1.e3)*T/np.pi))**0.5
	# convert to kg/m.s
	na = na*1.e-3 
	Beta = (8.*si.constants.k*T)/(3.*na) # coagulation kernel (m3/particle.s)
	# convert to cm3/particle.s
	Beta = Beta*1.e6
	Beta = np.ones((nsb, nsb))*Beta # spread over all size bins
	
	# set ODE solver (integration) tolerances
	atol = 1.e-4
	rtol = 1.e-5
	
	# make 1 dimensional
	N_perbin = np.squeeze(N_perbin)

	# call on the ODE solver
	sol = solve_ivp(dNdt, [0, integ_step], N_perbin, atol = atol, rtol = rtol, method = 'RK45', t_eval = [integ_step])

	N_perbin = (sol.y).reshape(1, -1) # get results (# particles/cm3)
	
	
	return(N_perbin)

# number of size bins
#nsb = 20
# initial concentration of particles per size bin (# particles/cm3)
#N_perbin = np.zeros((1, nsb))
#N_perbin[0, 0] = 2.e7 # assume monodisperse population
#T = 298. # chamber temperature (K)
#x = np.logspace(np.log10(3.e-3), np.log10(2.e-2), num = nsb) # radius of particles per size bin (um)
# size bin radius bounds (um)
#xb = np.zeros((nsb+1))
#xb[1:-1] = x[0:-1]+(x[1::]-x[0:-1])/2.
#xb[-1] = x[-1]+(x[-1]-xb[-2])
#xb[-1] = x[-1]*1.e6

#integ_step = 12.*3600. # integration step (s)
#sumt = 0. # total time (s)
#tot_time = 12.*3600.
#num_steps = (tot_time/integ_step)
#num_steps = int(num_steps)

# empty array for record of number concentration (# particles/cm3)
#N_perbint = np.zeros((num_steps+2, nsb))
#N_perbint[0, :] = N_perbin
# count on steps
#it = 0

# loop through times (s)
#while (sumt <= tot_time ):
#	N_perbint[it+1, :] = integ_coag(N_perbint[it, :], T, x, xb, nsb, integ_step)
#	sumt += integ_step # total time (s)
#	it += 1 # count on steps
	

#xb[-1] = xb[-1]*1.e-6 # reverse amplification of final size bin bound
#xb[0] = x[0]/2. # attain a lowermost size bin bound greater than 0
#dNdlogDp = N_perbint/((np.diff(np.log10(xb*2))).reshape(1, -1))

#import matplotlib.pyplot as plt

#plt.loglog(x*2., dNdlogDp[0, :], '-x')
#plt.loglog(x*2., dNdlogDp[-1, :], '-+')
#plt.show()