'''module to test that coagulation functions correctly'''
# module for testing coagulation, assumes calling from the PyCHAM home directory

import numpy as np
import sys
import os
# ensure modules can be seen 
# (assumes calling function is in the home folder)
sys.path.append(str(os.getcwd() + '/PyCHAM'))
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import scipy.constants as si

import coag as coag


def test_coag(): # define function


	# start of coagulation kernel test -------------------------------------------------
	print('first part of the coagulation test is to reproduce Fig. 15.7 of Jacobson (2005)')
	RH = 0.50 # relative humidity (fraction (0-1))
	T = 298.0 # temperature (K)
	sbr = np.logspace(-8.0, -4.6, num=50, endpoint=True, base=10.0) # size bin radius (m)
	sbVi = (4.0/3.0)*np.pi*sbr**3.0 # single particle volume for i sizes (m3)
	M = 200. # molecular weight of components (g/mol) 
	rho = 1. # components densities (g/cm3) in a 1D array
	rint = np.array((1.0e-8, 1.0e-5)) # size of interest (m)
	# concentration of particles per size bin (sbr) (particle/cc (air))
	num_part = np.ones((1, len(sbr)))
	# molecular concentration (molecules/cc (air)), arranged by component in rows and size 
	# bins (sbr) in columns
	num_molec = ((rho*1.0e6)*((1.0/M)*si.N_A))*sbVi*num_part
	# concentration of particles per size bin (sbr) (# particle/cm3 (air))
	num_part_rint = np.ones((1, len(rint)))
	# molecular concentration (# molecules/cm3 (air)), arranged by component in rows and size 
	# bins (sbr) in columns
	num_molec_rint = ((rho*1.0e6)*((1.0/M)*si.N_A))*((4.0/3.0)*np.pi*rint**3.0)*num_part_rint
	tint = 1. # time interval coagulation occurs over (s)
	sbbound = sbVi[0:-1]+(sbVi[1::]-sbVi[0:-1])/2.0 # size bin volume boundaries (m3)
	sbbound = np.append(0.0, sbbound)
	sbbound = np.append(sbbound, sbVi[-1]*1.e6)
	# radius bounds (m)
	rbou = ((3./(4.*np.pi))*sbbound)**(1./3.)
	num_comp = 1 # number of components
	vdWon = 0 # saying whether the van der Waals kernel should be calculated or ignored
	rad0 = sbr # original radius at size bin centre (um)
	PInit = 1.e5 # pressure inside chamber (Pa)
	testf = 1 # unit testing flag on
	sbVj = (4./3.)*np.pi*rint**3. # single particle volume for j sizes (m3)
	coag_on = 1 # flag to ensure coagulation considered
	siz_str = 0 # size structure flag
	wall_on = 0 # whether the wall is considered

	# call on coagulation module to complete kernel testing
	coag.coag(RH, T, sbr, sbVi, M, rint, num_molec, num_part, tint, sbbound, rbou,
		num_comp, vdWon, rho, sbVi, rad0, PInit, testf, num_molec_rint, num_part_rint, 
		sbVj, coag_on, siz_str, wall_on)
	
	print('end of coagulation kernel testing part, now moving onto Smoluchowski analytical solution test')
		
	# end of coagulation kernel test -------------------------------------------

	# start of Smoluchowski analytical solution test ----------------------------

	RH = 0.5 # relative humidity
	T = 298. # temperature (K)
	PInit = 1.e5 # pressure in chamber (Pa)

	
	# using an initially monodisperse population attempt to replicate Smoluchowski's
	# analytical solution
	# prepare required inputs
	nsb = 20 # number of size bins for semiimplicit approach
	Dpb = np.logspace(np.log10(5.1e-3), -1, num=(nsb+1)) # particle diameters bounds (um)
	Dpbw = Dpb[1::]-Dpb[0:-1] # size bin distances (um)
	Dpbc = Dpb[0:-1]+Dpbw/2. # diameters at size bin centres
	
	# radius of single particles (um)
	sbr = Dpbc/2.
	# volumes of single particles (um3)
	sbVi = ((4./3.)*np.pi*sbr**3.).reshape(-1)
	sbbound = ((4./3.)*np.pi*(Dpb/2.)**3.).reshape(1, -1) # volume boundaries per size bin (um3)
	# remember initial radius bounds (um)
	rbou00 = [((3./(4.*np.pi))*sbbound[0, 0])**(1./3.), ((3./(4.*np.pi))*sbbound[0, -1])**(1./3.)]
	# ensure upper bound on final bin very high
	sbbound[0, -1] = sbbound[0, -1]*1.e6
	# ensure lower bound on lower bin very low
	sbbound[0, 0] = 0.
	rbou = ((3./(4.*np.pi))*sbbound)**(1./3.) # radius boundaries per size bin
	
	num_comp = 1 # number of components
	M = np.ones((1, 1))*100. # molecular weight of components (g/mol)
	rho = np.array((1.)) # component density (g/cm3)
	# particle number concentration per size bin (# particles/cm3)
	num_part = np.zeros((1, sbVi.shape[0]))
	num_part[0, 0] = 1.e6 # (# particles/cm3)
	# molecular concentration (# molecules/cm3) per size bin, assuming 
	# one component present
	y = ((rho*1.e-12)/M)*(sbVi*num_part)*si.N_A
	tint = 600. # time to coagulate over (s)
	testf = 3 # testing flag
	coag_on = 1 # coagulation turned on
	siz_str = 1 # size structure
	wall_on = 0 # flag for whether to consider wall
	vdWon = 0 # flag for Van Der Waals kernel
	ytot0 = y.sum() # starting total concentration of components (molecules/cc)
	ttot = 0. # cumulative time through simulation (s)
	tot_time = 12.*3600. # total simulation time
	#while (ttot < tot_time):
			
	#	[num_part, y, sbr, Gi, eta_ai, sbVi, sbbound, rbou] = coag.coag(RH, T, 
	#		sbr*1.e-6, sbVi.reshape(1, -1)*1.e-18, M, sbr*1.e-6, y, num_part, tint, sbbound*1.e-18, rbou, 
	#		num_comp, vdWon, rho, sbVi, sbr, PInit, testf, y, num_part,
	#		sbVi.reshape(1, -1)*1.e-18, coag_on, siz_str, wall_on)
	#	ttot += tint # total time through simulation (s)
	#	num_part = num_part.reshape(1, -1)
	#	y = y.reshape(1, -1)
	#	sbbound = sbbound.reshape(1, -1)
	#	rbou = rbou.reshape(1, -1)
		
	# check for mass conservation from semiimplicit
	#print(str('change (%) in total mass of particles between start and end of simulation using semi-implicit number conserving approach, where a positive means an increase in mass: ' + str(((y.sum()-ytot0)/ytot0)*100.) + ' %'))
	
	# integrated solution based on a mass-conserving approach ---------------
	from coag_integ2 import integ_coag
	
	nsb = 8 # number of size bins for integrated solution
	num_comp = 1 # number of components
	# particle number concentration per size bin (# particles/cm3)
	N_perbin = np.zeros((1, nsb))
	N_perbin[0, 0] = 1.e6 # (# particles/cm3)
	
	# molar volume of single component (um3/mol)
	y_MV = M*(1./(rho*1.e-12))
	
	x = np.logspace(np.log10(3.e-3), np.log10(2.e-2), num = nsb) # radius of particles per size bin (um)
	# size bin radius bounds (um)
	xb = np.zeros((nsb+1))
	if (nsb > 1):
		
		xb[1:-1] = x[0:-1]+(x[1::]-x[0:-1])/2.
		xb[-1] = x[-1]+(x[-1]-xb[-2])
		xb[-1] = x[-1]*1.e6
	else:
		xb[0] = 0.
		xb[1] = x[-1]*1.e6
	
	# single particle volume (um3) per size bin
	V0 = (4./3.)*np.pi*(x**3.)
	
	# component concentration (# molecules/cm3)
	y = ((np.tile((V0*N_perbin).reshape(-1, 1), [1, num_comp]))/np.tile(y_MV.reshape(1, -1), [nsb, 1]))*si.constants.N_A
	y = (y.flatten().reshape(-1, 1))
	y = y[:, 0]
	
	tint = tot_time/1.e1 # integration time interval (s)
	ttot = 0. # cumulative time through simulation (s) 
	ybfe = sum(y)
	
	msdf = 1 # monomer size distribution flag
	
	while (ttot < tot_time):
		
		[N_perbin, y, x] = integ_coag(N_perbin, T, x, xb, nsb, tint, y, y_MV, num_comp, msdf)
		#N_perbin = integ_coag(N_perbin, T, x, xb, nsb, tint)
		ttot += tint # total time (s)
	yaft = sum(y)
	
	y_diff = ((yaft-ybfe)/ybfe)*100.
	print(str('% change in total particle-phase molecular concentration from value prior to integration solution, where positive means an increase during the solution: ' + str(y_diff) + ' %'))

	xb[-1] = x[-1]*2. # reverse amplification of final size bin bound
	xb[0] = x[0]/2. # attain a lowermost size bin bound greater than 0
	dNdlogDp = N_perbin/((np.diff(np.log10(xb*2))).reshape(1, -1))
	
	# -----------------------------------------------

	# Smoluchowski analytical solution
	V0 = 1.e-3 # initial volume in first size bin (assuming monodisperse population)
	Beta = 3.e-4 # the coagulation kernel used in Fig. 4 Smoluchowski (1916)
	t = np.linspace(0., 4.*3600.) # time steps (s)
	nsb = 210 # number of size bins
	# array for storing results of volume concentration per size bin
	Vres = np.zeros((len(t), nsb+1))
	Vres[0, 0] = Vres[0, -1] = V0
	exp_lo = np.linspace(0, Vres.shape[1]-1, nsb) # exponent for numerator
	exp_hi = np.linspace(2, Vres.shape[1]-1+2, nsb) # exponent for denominator
	
	for it in range(1, len(t)): # loop through times
		Vres[it, 0:-1] = (V0*(Beta*t[it])**exp_lo)/((1.+Beta*t[it])**exp_hi)
		Vres[it, -1] = V0/(1.+Beta*t[it])
	#plt.plot(t/3600., Vres[:, 0]/V0, '-g')
	#plt.plot(t/3600., Vres[:, 1]/V0, '--g')
	#plt.plot(t/3600., Vres[:, 2]/V0, ':g')
	#plt.plot(t/3600., Vres[:, 3]/V0, '.g')
	#plt.plot(t/3600., Vres[:, 0:-1].sum(axis=1)/V0, '-g')
	#plt.plot(t/3600., Vres[:, -1]/V0, '--r')
	#plt.show()
	
	h = np.linspace(0, (tot_time)) # time steps (s)
	npb = np.zeros((len(h), nsb)) # particle number concentration per size bin (# particles/cm3)
	npb[0, 0] = 1.e6 # initial monodisperse particle population (# particles/cm3)
	r0 = 3.e-3 # particle radius in first size bin (um)
	V0 = (4./3.)*np.pi*(r0)**3. # particle volume in first size bin (um3)	
	km1 = np.linspace(0, nsb-1, nsb) # exponents
	
	# suggested kernel calculation from Eq. 15.16 Jacobson 2005 ------
	# dynamic viscosity of air (g/m.s) (Eq. 4.54 Jacobson 2005)
	# note the conversion of the universal gas constant from 
	# kg.m2/s2.mol.K to g.m2/s2.mol.K
	na = 5./(16.*si.N_A*3.673e-10**2.)*((28.966*(si.R*1.e3)*T/np.pi))**0.5
	# convert to kg/m.s
	na = na*1.e-3 
	Beta = (8.*si.k*T/(3.*na)) # coagulation kernel (m3/particle.s)
	# convert to cm3/particle.s
	Beta = Beta*1.e6
	
	Beta = Beta*npb[0, 0] # from Eq. 15.15 Jacobson (2005)
		
	# ----------------------------------------------------------------
	
	# loop through times
	for ti in range(1, len(h)):
		num = npb[0, 0]*((0.5*h[ti]*Beta)**km1)
		den = (1.+(0.5*h[ti]*Beta))**(km1+2.)
		npb[ti, :] = num/den

	# create plot, note when using the Beta value of the Smoluchowski analytical solution
	# the first four hours of this plot should agree with Fig. 4 of the 1916 Smoluchowski
	# paper (http://eudml.org/doc/215805)
	#plt.plot(h/3600., npb[:, 0]/npb[0, 0], '-k', label = 'sb1')
	#plt.plot(h/3600., npb[:, 1]/npb[0, 0], '--k', label = 'sb2')
	#plt.plot(h/3600., npb[:, 2]/npb[0, 0], ':k', label = 'sb3')
	#plt.plot(h/3600., npb[:, 3]/npb[0, 0], '.k', label = 'sb4')
	#plt.plot(h/3600., (npb[:, :]/npb[0, 0]).sum(axis=1), '-b', label = '$\Sigma N$')
	#plt.xlabel('Time (hours)')
	#plt.ylabel('$N/N_{0}$')
	#plt.legend()
	#plt.show()
	
	# prepare for normalising Smoluchowski output number concentrations to width of size bins ---------------
	# volume of single particle per size bin (um3) - see Figure 15.1 of Jacobson (2005) 
	# for monomer size distribution
	sbVi = np.arange(1, nsb+1)*V0
	# radius of single particles per size bin (um3)
	sbri = ((3./(4.*np.pi))*sbVi)**(1./3.)
	# volume bounds (um3)
	sbVb = sbVi[0:-1] + (sbVi[1::]-sbVi[0:-1])/2.
	# radius bounds
	rb = ((3./(4.*np.pi))*sbVb)**(1./3.)
	# concatenate lower and upper bounds
	rb = np.concatenate([(sbri[0]-(rb[0]-sbri[0])).reshape(1), rb, (sbri[-1]+(sbri[-1]-rb[-1])).reshape(1)])
	
	# check for mass conservation from analytical solution
	#print((sbVi*npb[0, :]).sum())
	#print((sbVi*npb[-1, :]).sum())
	
	# normalised number size distribution of particles (# particle/cm3)
	# analytical results
	npb = npb/(np.log10(rb[1::]*2.)-np.log10(rb[0:-1]*2))
	# semiimplicit results
	rbou[0, 0] = rbou00[0]
	rbou[0, -1] = rbou00[1]
	npbs = num_part/(np.log10(rbou[0, 1::]*2.)-np.log10(rbou[0, 0:-1]*2))
	# two-point moving average
	npbs = (npbs[0, 1::]+npbs[0, 0:-1])/2.
	sbDp = (sbr[1::]*2.+sbr[0:-1]*2)/2.
		
	# plot results
	plt.loglog(sbri*2, npb[0, :], '--.', label = 'initial') # initial	
	plt.loglog(sbri*2, npb[-1, :], '-+', label = 'Smoluchowski') # final
	plt.loglog(sbDp, npbs, ':x', label = 'semi-implicit number conserving') # final
	plt.loglog(x*2., dNdlogDp[0, :], '-+', label = 'integrated volume conserving')
	
	ax0 = plt.gca()	
	ax0.set_ylim([1, 1.e8])
	ax0.set_xlim([5.e-3, 1.e-1])
	ax0.legend()
	
	locmaj = matplotlib.ticker.LogLocator(base=10., subs=(1., ))
	ax0.yaxis.set_major_locator(locmaj)
	locmin = matplotlib.ticker.LogLocator(base=10., subs=np.arange(2, 10)*.1) 
	ax0.yaxis.set_minor_locator(locmin)
	ax0.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
	
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	plt.xlabel('$D_{p}\, (\mathrm{\mu m})$')
	plt.ylabel('$dN$ (# particles$\, \mathrm{cm^{-3}}$)/$d$log$_{10}D_p (\mathrm{\mu m})$')
	plt.show()
	# ----------------------------------------------------------------------

test_coag() # call on test	
