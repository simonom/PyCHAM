'''module to solve coagulation through integration'''
# Using the differential equation for coagulation given in Eq. 15.2
# of Jacobson (2005), the new number size distribution is solved.  
# The volume-conserving result is found using Eq. 15.2 of Jacobson (2005) to
# estimate the number fraction per size bin contributing to/being lost from size
# bins to estimate the molecular concentration being exchanged
# Because the coagulation kernel is a constant for integration but is size-dependent
# and coagulation changes particle sizes, this is one source of imprecision that is
# reduced by having a high ratio of number of size bins to total particle size range.
# Another source of imprecision is the assumption that coagulation pairs have
# the same volume throughout the integration.  This flaw is highlighted when
# considering an initial number size distribution of all particles in one size bin.  If the
# size bin is sufficiently wide self-coagulation will produce larger particles in the 
# same size bin.  Eventually self-coagulation produces particles in a larger size bin
# but this is not recognised as the initial volumes of particles inform the size
# of coagulated particles.  This problem is minimised through increasing ratio of
# size bin number to total size range and through reducing time interval.  This problem
# is entirely removed if a monomer size distribution (Fig. 15.1 Jacobson 2005) is 
# present.
# The final source of imprecision for the volume-conserving method is that the number
# concentration of particles is not conserved, therefore number concentration is
# estimated assuming the particles in each size bin have the same volume per 
# particle as the volume per particle before the integration.
# All the above sources of imprecision are removed when a monomer size distribution
# of particles is present.  However, if this requires interpolation/extrapolation, then
# imprecision may be introduced through smoothing of particle composition (component
# particle-phase concentrations).  The monomer size distribution is not included here 
# because on a desktop computer (as at 2021), the computation time is prohibitively
# high when the size range of particles is within that expected for use with PyCHAM.
# For example, a radius range of three orders of magnitude would require a trillion
# size bins with a monomer volume size distribution.


import numpy as np
import scipy as si
from scipy.integrate import solve_ivp

def integ_coag(N_perbin, T, x, xb, nsb, integ_step, y, y_MV, num_comp, msdf):
	
	# inputs: -------------------------
	# N_perbin - initial concentration of particles per size bin (# particles/cm3)
	# T - chamber temperature (K)
	# x - radius of particles per size bin (um)
	# xb - radius of particle size bin bounds (um)
	# nsb - number of size bins
	# integ_step - integration step (s)
	# y - component concentrations (# molecules/cm3)
	# y_MV - molar volume of water (um3/mol)
	# num_comp - number of components
	# msdf - flag for monomer size distribution (1=on, 0 = off)
	# -----------------------------------
	
	# volume-conserving part -----------------------------------------------------------
	
	def dydt(t, y): # define ordinary differential equation function for volume-conserving solution
		
		# inputs: ------------------------
		# t - time to integrate over (s)
		# y - concentration of components (# molecules/cm3)
		# ----------------------------------
		# reshape y into size bins in rows and components in columns
		y = y.reshape(nsb, num_comp)
		
		# estimate number particle concentration based 
		# on new concentrations of components
		# here, concentrations are converted to moles/cm3, then um3/cm3 
		# and then divided by the original
		# single particle volume to get the new number of particles
		Npb = np.sum((y/si.constants.N_A)*y_MV, axis = 1)/V0
		
		Npbr = np.tile(Npb.reshape(1, -1), [nsb, 1])
		Npbc = np.tile(Npb.reshape(-1, 1), [1, nsb])
	
		# empty array for rate of change (# particles/cm3.s)
		dd = np.zeros((nsb, num_comp))
		
		# fraction tracking from size bins
		frac_trak = np.zeros((nsb, 2))
		
		# size bin loop
		for isb in range(nsb):
			
			# empty array to hold indices of contributing size bins tiled over rows
			pindxr = np.zeros((nsb, nsb))
			pindxr = (pindxr == 1) # convert to Boolean
			# empty array to hold indices of contributing size bins tiled over columns
			pindxc = np.zeros((nsb, nsb))
			pindxc = (pindxc == 1) # convert to Boolean
			
			# identify pairs that coagulate to give particle in this size bin
			pindx = (Vcoag0 >= Vb[isb])*(Vcoag0 < Vb[isb+1])
			
			# ensure coagulation with itself can not produce a new volume in this size bin
			pindx[isb, isb] = 0
			
			# prepare the production result matrix
			prod = np.zeros((nsb, nsb))
			# prepare the fraction result matrices
			rfrac = np.zeros((nsb, nsb))
			cfrac = np.zeros((nsb, nsb))
			
			# number production term from each constituent of size bin pairs (# particles/cm3.s)
			prod[pindx] = 0.5*Beta[pindx]*Npbr[pindx]*Npbc[pindx]
			
			pindx = (prod > 0.) # the non-zero production indices
			
			if ((pindx > 0.).any()): # if anything produced
			
				# fraction of size bins contributing to production	
				# fraction from size bins tiled over rows
				pindxr[:, :] = pindx[:, :]
				pindxr[:, isb] = 0 # ensure no self-counting
				rfrac[pindxr] = (prod[pindxr]/Npbr[pindxr])
				
				# fraction from size bins tiled over columns
				pindxc[:, :] = pindx[:, :]
				pindxc[isb, :] = 0 # ensure no self-counting
				cfrac[pindxc] = (prod[pindxc]/Npbc[pindxc])
				
				# fraction gain rate from each size bin (fraction/s)
				tot_frac = np.sum(rfrac, axis=0).reshape(-1, 1)+np.sum(cfrac, axis = 1).reshape(-1, 1)
				
				frac_trak[:, 0] += tot_frac[:, 0]
				
				
				# rate of change due to production (# molecules/cm3.s)
				dd[isb, :] += np.sum(y[:, :]*np.tile(tot_frac, [1, num_comp]), axis = 0)
				

			if (Npb[isb] > 0.): # if anything to lose
					
				# loss term (# particles/cm3.s)
				loss = Npb[isb]*(Beta[isb, :]*Npb)
				
				# ignore coagulation that produces particles in this size bin
				loss[(Vcoag0[isb, :] < Vb[isb+1])] = 0.
				
				# sum losses across all size bins
				loss = sum(loss)
				
				# fraction of this size bin lost (fraction/s)
				tot_frac = (loss/Npb[isb])
				
				frac_trak[isb, 1] += tot_frac
				
				# rate of change due to loss (# molecules/cm3.s)
				dd[isb, :] -= (y[isb, :]*tot_frac)
			
		# flatten y and dd so that components change then size bins
		y = y.flatten().reshape(-1, 1)
		y = y[:, 0]
		dd = dd.flatten()
		
		return(dd)
	
	# preparation for integration part --------------------------------------------------------------------------
	nsb0 = nsb # remember original number of size bins
	V0 = (4./3.)*np.pi*(x**3.) # single particle volume (um3) per size bin
	
	# if change to monomer size distribution wanted, then make appropriate changes ---------
	if (msdf == 1):
	
		# normalise component concentrations by size bin width
		# first turn into matrix with size bins in rows, components in columns
		ymat = y.reshape(nsb, num_comp)
		
		# shrink uppermost bin bound (um)
		xb[-1] = x[-1]*2.1
		
		# normalise by bin width (# molecules/cm3/um)
		ymat = ymat/(np.tile(np.diff(xb).reshape(-1, 1), [1, num_comp]))
		
		# smallest size bin single particle volume (um3)
		ssbv = V0[0]
		# largest size bin single particle volume (um3)
		lsbv = V0[-1]
		
		# number of size bins for a monomer size distribution
		nsb_md = int(np.ceil(lsbv/ssbv)) 
		
		# single particle volumes (um3) in monomer size distribution
		V0 = np.arange(ssbv, (nsb_md*ssbv+ssbv/10.), ssbv)
		# single particle radius (um) in monomer size distribution
		xn = ((3./(4.*np.pi))*V0)**(1./3.)
		# new number of size bins
		nsb = len(xn)
		# size bin bounds in monomer size distribution (um)
		xbn = np.zeros((nsb+1))
		xbn[1:-1] = xn[0:-1]+np.diff(xn)/2.
		xbn[0] = xb[0]
		xbn[-1] = xn[-1]*2.1
		
		# empty array for new concentrations normalised by bin width (# molecules/cm3/um)
		ymatn = np.zeros((nsb, num_comp))
		
		# interpolate to new size bins 
		for ic in range(num_comp): # loop through components
			# new concentrations normalised by size bin width (# molecules/cm3/um)
			ymatn[:, ic] = np.interp(xn, x, ymat[:, ic])
		
		# interpolated concentrations (# molecules/cm3)
		ymatn = ymatn*(np.tile(np.diff(xbn).reshape(-1, 1), [1, num_comp]))
	
		# flatten so that components change then size bins
		y = ymatn.flatten()
		
		Vb = (4./3.)*np.pi*(xbn**3.) # volume bound per size bin (um3)
		# empty array for volume of combined pairs (um3)
		Vcoag0 = np.tile(V0.reshape(1, -1), [nsb, 1])+np.tile(V0.reshape(-1, 1), [1, nsb])
		
		
	else: # if no change to monomer size distribution
		
		Vb = (4./3.)*np.pi*(xb**3.) # volume bound per size bin (um3)
		# empty array for volume of combined pairs (um3)
		Vcoag0 = np.tile(V0.reshape(1, -1), [nsb, 1])+np.tile(V0.reshape(-1, 1), [1, nsb])
	
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
	
	# ensure 1 dimensional
	if (y.ndim > 1):
		y = y[:, 0]
	
	# call on the ODE solver
	sol = solve_ivp(dydt, [0, integ_step], y, atol = atol, rtol = rtol, method = 'RK45', t_eval = [integ_step])

	y = (sol.y) # get results (# molecules/cm3)

	# reshape y into size bins in rows and components in columns
	y = y.reshape(nsb, num_comp)
	
	# estimate number particle concentration based 
	# on new concentrations of components
	# here, concentrations are converted to moles/cm3, then um3/cm3 
	# and then divided by the original
	# single particle volume to get the new number of particles (# particles/cm3)
	N_perbin = np.sum((y/si.constants.N_A)*y_MV, axis = 1)/V0
	
	# if monomer size distribution chosen the revert to original size distribution
	if (msdf == 1):
		# first normalise component concentrations by size bin width (# molecules/cm3/um)
		y = y/(np.tile(np.diff(xbn).reshape(-1, 1), [1, num_comp]))
		# normalise particle number concentration by size bin width (# particles/cm3/um)
		N_perbin = N_perbin/np.diff(xbn)
		
		# empty interpolation results array (# molecules/cm3/um)
		ymatn = np.zeros((nsb0, num_comp))
		
		for ic in range(num_comp): # loop through components
			# interpolate back to original size bin centres
			ymatn[:, ic] = np.interp(x, xn, y[:, ic])
		
		# (# particles/cm/um)
		
		N_perbin = np.interp(x, xn, N_perbin)
		
		# return to absolute concentrations (not normalised)
		# (# molecules/cm3)
		ymatn = ymatn*(np.tile(np.diff(xb).reshape(-1, 1), [1, num_comp]))
		N_perbin = N_perbin.reshape(-1, 1)*(np.diff(xb).reshape(-1, 1))
		y = ymatn # rename variable
	
	# flatten y so that components change then size bins
	y = y.flatten()
	
	print(N_perbin.reshape(1, -1))
	#import ipdb; ipdb.set_trace()
	
	return(N_perbin.reshape(1, -1), y, x)