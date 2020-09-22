'''module to model particle loss to walls, called on by ode_gen'''
# the process of particle deposition to walls is reproduced here using either
# the McMurry and Rader (1985) model, or the user-defined parameters

import numpy as np
import scipy.constants as si
from scipy import integrate
import matplotlib.pyplot as plt 

def wallloss(Pn, Cn, Gi, eta_ai, Dp, MW, Varr, sbn, nc, TEMP, t, 
			inflectDp, pwl_xpre, pwl_xpro, inflectk, ChamR, Rader, testf, p_char, 
			e_field, num_asb):

	# inputs:----------------------------------------------------------
	
	# Pn - particle number concentration per size bin before time step
	# (# particle/cc (air))
	# Cn - particle phase concentration per component per size bin
	# after a time step over which gas phase reaction and 
	# gas-particle partitioning have occurred (molecules/cc (air))
	# Gi - Cunningham slip-correction factor (dimensionless)
	# eta_ai - dynamic viscosity of air (g/m.s)
	# Dp - particle diameters (m)
	# MW - component molecular weight (g/mol)
	# Varr - volume of single particles per size bin (m3)
	# sbn - number of size bins
	# nc - number of species
	# TEMP - system temperature (K)
	# t - time that wall loss occurs over (s)
	# inflectDp - particle diameter at which wall loss inflection occurs (m)
	# pwl_xpre - x value preceding inflection point
	# pwl_xpro - x value proceeding inflection point
	# inflectk - deposition rate at inflection (/s)
	# ChamR - spherical equivalent radius of chamber (below eq. 2 Charan (2018)) (m)
	# Rader - flag of whether or not to use Rader and McMurry approach
	# testf - flag of whether in testing mode or not (1 for one 0 for normal mode)
	# p_char - average number of charges per particle (/particle)
	# e_field - average electric field inside chamber (g.m/A.s3)
	# num_asb - number of actual particle size bins
	# ----------------------------------------------------------------
	if Rader == 0: # manual input of wall loss rate
		
		Beta = np.zeros((Pn.shape))
		Beta[Dp<inflectDp, 0] = 10**((np.log10(inflectDp)-
								np.log10(Dp[Dp<inflectDp]))*pwl_xpre+np.log10(inflectk))
		Beta[Dp>=inflectDp, 0] = 10**((np.log10(Dp[Dp>=inflectDp])-
								np.log10(inflectDp))*pwl_xpro+np.log10(inflectk))
		
	if Rader == 1:
		# -------------------------------------------------------------------------
		# McMurry & Rader option McMurry 1985, DOI: 10.1080/02786828508959054
		# average number of charges per particle (/particle)
		n = p_char
		# elementary charge (C==A.s) (Charan et al. 2018)
		e = 1.602e-19
		# average magnitude of electric field (g.m/A.s3) (note: V/m == kg.m/A.s3) 
		E = e_field
		
		# electrostatic migration velocity/deposition velocity 
		# (eq. 11 McMurry and Rader (1985), eq. 5 Charan (2018)) (m/s)
		ve = np.abs((n*e*Gi*E)/(3.0*np.pi*eta_ai*Dp))
		
		
		# mass of components in particles (g)
		ish = np.squeeze(Pn<1.0e-10) # only use size bins where particles available
		ve[ish] = 0.0
		mass = np.zeros((sbn-1, nc))
		ish = (Pn>1.0e-20)[:, 0] # size bins where particles present
		
		mass[ish, :] = ((Cn.reshape(sbn-1, nc)[ish, :]/Pn[ish, 0].reshape(-1, 1))/
						si.N_A)*MW.reshape(1, nc)
		# density of particles (g/m3)
		rho = np.ones((sbn-1))
		rho[ish] = np.sum(mass[ish, :], 1)/Varr[ish]
		
		
		# terminal particle settling velocity (m/s) (eq. 4 Charan (2018) also 
		# eq. 20.4 Jacobson 2005)
		# gravitational constant has units m/s2
		vs = (Dp**2.0*rho*si.g*Gi)/(18.0*eta_ai)
	
		# particle diffusion coefficient (15.29 of Jacobson (2005) and 
		# 8.73 of Seinfeld and pandis (1998)) (m2/s)
		# multiply eta_a by 1.0e-3 to convert from g/m.s to kg/m.s
		# this makes it consistent with the units of Boltzmann constant
		Dpi = (((si.k*TEMP)/(3.0*np.pi*Dp*(eta_ai*1.0e-3)))*Gi)
		
		# eddy diffusion coefficient (/s) (scalar)
		if testf == 1:
			Ke = 6.4e-3 
		else:
			Ke = 1.0
			
		# x and y terms in eq. 2 of Charan (2018)
		# x (eq. 13 McMurry (1985)) (dimensionless)
		x = (np.pi*vs)/(2.0*((Ke*Dpi)**(1.0/2.0)))
		# y (eq. 14 McMurry (1985)) (dimensionless)
		y = (np.pi*ve)/(2.0*((Ke*Dpi)**(1.0/2.0)))
		
		z1 = (x+y)
		z2 = (x-y)
		
		# empty array for Debye function results
		D1 = np.zeros(sbn-1)
		D11 = np.zeros(sbn-1)
		
		for sbi in range(sbn-1): # size bin loop
			
			if Pn[sbi]<1.0e-10: # only consider size bins with particles inside
				continue 
			
			# state integration function
			Dint = lambda t: t/(np.exp(t)-1.0)
			# integrate for first order Debye function (eq. 3 Charan (2018))
			# there is an upper limit to the variable, above which an overflow occurs in 
			# the integration, the integration reaches a threshold result for high 
			# positive numbers (>~700), so just set this threshold as the result if input 
			# too great
			if z1[sbi] > 700.0:
				a = np.array(([1.644934066848228]))
			else:
				a = integrate.quad(Dint, 0.0, z1[sbi])
			if z2[sbi] > 700.0:
				a1 = np.array(([1.644934066848228]))
			else:
				a1 = integrate.quad(Dint, 0.0, z2[sbi])
	
			D1[sbi] = (1.0/z1[sbi])*a[0]
			D11[sbi] = (1.0/z2[sbi])*a1[0]
		
		# first bit of Beta (loss rate to walls (/s)) calculation (eq. 2 Charan (2018))
		Beta1 = np.zeros((sbn-1))
		Beta1[ish] = (3.0*((Ke*Dpi[ish])**(0.5)))/(np.pi*ChamR*x[ish])
		Beta2 = ((x+y)**2.0)/2.0+(x+y)*D1+(x-y)*D11 # second bit
	
		# Beta (loss rate to walls (/s)) calculation (eq. 2 Charan (2018) and eq. 15 
		# McMurry and Rader 1985) - value
		# represents fraction of particles lost to walls every second (Beta meaning is 
		# just above eq. 1 of McMurry (1985) DOI: 10.1080/02786828508959054)
		Beta = (Beta1*Beta2).reshape(-1, 1)
	
	if (Beta<0).sum()>0:
		Beta[Beta<0] = 0.0
		print('Warning Beta in walloss.py estimated below 0, which is not possible, so value forced to 0.0')
	
	if testf == 1: # if in test mode
		return(Beta)
	
	# integrate this fraction over the time step interval to give total 
	# fraction lost over interval	
	Beta = Beta*t
	
	# find where loss of particles exceeds available number of particles and change
	# Beta to the realistic maximum
	ish = Pn<(Beta*Pn)
	Beta[ish] = 1.0
		
	# new particle number concentration and particle-phase concentrations of components
	Pn -= (Beta*Pn)
	Cn -= (Beta*Cn.reshape(sbn-1, nc)).flatten(order='C')
			
	# remove particles and their components if particle number negative
	ish = np.array((np.where(Pn<1.0e-8)))
	for i in ish[0, :]:
		Pn[i] = 0.
		Cn[nc*(i):nc*(i+1)] = 0.
	
	# prepare output
	Pn.reshape(num_asb, 1)
	
	return(Pn, Cn)
