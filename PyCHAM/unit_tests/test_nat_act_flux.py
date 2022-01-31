'''module for testing the two-stream radiative transfer model contained in nat_act_flux'''
# unit test for the two-stream radiative transfer model contained in nat_act_flux module,
# which is based on Shettle and Weinman (1970)
# (doi.org/10.1175/1520-0469(1970)027<1048:TTOSIT>2.0.CO;2)
#. Reproduction of the Shettle and Weinman (1970) figures is the basis of this unit
# test
# assumes calling from the PyCHAM home folder
print('testing the two-stream radiative transfer model contained in nat_act_flux by reproducing the Shettle and Weinman (1970) (doi.org/10.1175/1520-0469(1970)027<1048:TTOSIT>2.0.CO;2) figures')

import os
import sys
dir_path = os.getcwd() # current working directory
# temporarily add the PyCHAM folder to path
sys.path.append(str(dir_path+'/PyCHAM'))
import nat_act_flux
import numpy as np
import matplotlib.pyplot as plt # for plots


# define function
def test_nat_act_flux():

	# inputs ------------------------------------------
	# ---------------------------------------------------
	
	F0 = 1370. # incoming solar radiation flux (top of atmosphere) (W/m2)
	theta = np.pi/4. # solar zenith angle (radians)
	tau = np.log(2.)
	
	# calling function flag (<=0 for unit testing: <=-1 for reproducing 
	# Shettle and Weinman results, 0 for reproducing Jacobson (2005) results)
	# is variable callf


	# ----------------------------------------------------------
	# testing against Fig. 1 Shettle and Weinman (1970)
	# albedo as a function of solar zenith angle
	# note that we only consider reproducing the
	# Eddington approximation of Fig.1, not the Irvine (1968)
	# solution
	
	callf = -1 # identity of caller
	# total optical depth of atmosphere text below 
	# Eq. 19 of Shettle and Weinman (1970)
	tau = np.array((16.)).reshape(1)
	# scattering asymmetry factor text below 
	# Eq. 19 of Shettle and Weinman (1970)
	g = np.array((0.786)).reshape(1)
	
	# Figure 1 of Shettle and Weinman (1970) reproduction,
	# note that Eq. 20 used to estimate 
	# surface albedo (A) for Figure 1 and single 
	# scattering albedo less than 1, whilst Eq. 22 used for a
	# single scattering albedo equal to 1.
	a = np.array((0.95, 1.))
	# solar zenith angle (radians) to evaluate at
	theta0 = np.arange(0., 1.57*1.01, 1.57/20.)
	# filler for phase function calculation
	Pfunc_text = []
	
	# prepare results matrix
	Ares = np.zeros((len(a), len(theta0))) 
	ain = 0 # counter through single scattering albedo (a)
	# surface albedo, note that when testing against Fig.1,
	# my understanding is that the albedo found in Eqs. 
	# 20 and 22 of Shettle and Weinman (1970) represents
	# the upward diffuse radiation of the atmosphere,
	# with zero contribution by the surface.
	A = 0.
	for ai in a: # single scattering albedo loop
		ainow = np.array((ai)).reshape(1)
		
		tin = 0 # reset counter through incident sun angle (theta)
		for ti in theta0:
			
			mu0 = np.cos(ti) # cosine of incident sun angle
			
			# call function for results
			Ares[ain, tin]  = nat_act_flux.nat_act_flux(A, ainow, F0, ti, tau, callf, mu0, 1, g, Pfunc_text)
		
			tin += 1 # keep count
		ain += 1 # keep count
	
	# plot figure 1
	fig, ax0 = plt.subplots() # prepare plot
	ax0.plot(theta0, Ares[0, :], label = str('a = ' + str(a[0])))
	ax0.plot(theta0, Ares[1, :], label = str('a = ' + str(a[1])))
	ax0.legend(loc = 'lower right')
	ax0.set_title('Reproduction of Eddington approximation in Figure 1 of \nShettle and Weinman (1970) to test single layer atmosphere case')
	ax0.set_ylabel(r'Albedo')
	ax0.set_xlabel(r'$\mu_{0}$ (RAD)')
	ax0.tick_params(direction = 'in', which = 'both')
	plt.show()
	
	# ----------------------------------------------------------
	# Figure 2 of Shettle and Weinman (1970), 
	# note we only address the Eddington approximation here
	A = 0. # surface albedo
	ai = np.array((1.)).reshape(1) # single scattering albedo
	# optical depths
	tau_all = np.array((1., 2., 4., 8., 16., 32., 64.))
	# scattering asymmetry factor from Fig. 2 of Shettle and Weinman (1970)
	g = 0.75
	Pfunc_text = [] # filler
	# prepare results matrix
	Ares = np.zeros((len(tau_all), len(theta0))) 
	# caller identity
	callfi = -2
	
	tauin = 0 # resent tau count
	
	for taui in tau_all: # loop through optical depth
		tin = 0 # reset theta count
		for ti in theta0: # loop through incident sun angle
		
			mu0 = np.cos(ti) # cosine of incident sun angle
			
			# call function for results
			Ares[tauin, tin]  = nat_act_flux.nat_act_flux(A, ai, F0, ti, taui, callfi, mu0, 1, g, Pfunc_text)
			
			tin += 1 # keep count
		tauin += 1 # keep count
	
	# plot figure 2
	from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
	fig, ax0 = plt.subplots() # prepare plot
	ax0.plot(np.cos(theta0), Ares[0, :], label = str(r'$\tau$ = ' + str(tau_all[0])))
	ax0.plot(np.cos(theta0), Ares[1, :], label = str(r'$\tau$ = ' + str(tau_all[1])))
	ax0.plot(np.cos(theta0), Ares[2, :], label = str(r'$\tau$ = ' + str(tau_all[2])))
	ax0.plot(np.cos(theta0), Ares[3, :], label = str(r'$\tau$ = ' + str(tau_all[3])))
	ax0.plot(np.cos(theta0), Ares[4, :], label = str(r'$\tau$ = ' + str(tau_all[4])))
	ax0.plot(np.cos(theta0), Ares[5, :], label = str(r'$\tau$ = ' + str(tau_all[5])))
	ax0.plot(np.cos(theta0), Ares[6, :], label = str(r'$\tau$ = ' + str(tau_all[6])))
	ax0.set_ylabel(r'Albedo')
	ax0.set_xlabel(r'COSINE OF SOLAR ZENITH ANGLE ($\mu_{0}$)')
	ax0.legend(loc = 'lower right')
	ax0.set_title('Reproduction of Eddington approximation from Figure 2 of \nShettle and Weinman (1970) to test single layer atmosphere case')
	ax0.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f')) # x axis tick labels at two decimal places
	ax0.invert_xaxis() # reverse x axis
	ax0.tick_params(direction = 'in', which = 'both')
	plt.show()
	
	# -------------------------------------------------------------------------------
	# Fig. 3 of Shettle and Weinman (1970)
	# Atmospheric albedo as a function of asymmetry factor
	# and the effective optical thickness for a single layer
	# atmosphere
	
	A = 0. # surface albedo
	ai = np.array((1.)).reshape(1) # single scattering albedo
	theta0 = 0. # solar zenith angle
	mu0 = np.cos(theta0) # cosine of solar zenith angle
	tau_all =  np.array((2., 16.)) # the two total optical thicknesses
	g = np.arange(-1., 1.*1.001, 0.1) # range of asymmetry factor used in function
	Ares = np.zeros((len(tau_all)+5, len(g)))
	
	gcnt = 0 # keep count on asymmetry factor
	callfi = -3 # caller identity
	
	# count through tau (atmosphere optical depth)
	for gi in g: # loop through asymmetry factors
		tau_cnt = 0 # keep count on changing optical depth
		for taui in tau_all:

			Ares[tau_cnt, gcnt] = nat_act_flux.nat_act_flux(A, ai, F0, theta0, taui, callfi, mu0, 1, gi, Pfunc_text)
			tau_cnt += 1 # keep count on changing optical depth
		gcnt += 1 # keep count on asymmetry factor
		
	# now try keeping T constant by varying atmosphere optical depth
	T_all = np.array((1., 2., 4., 8., 16.))
	callfi = -3 # caller identity
	
	callfi_all = np.array((-4, -5, -6, -7, -8))
	gcnt = 0 # keep count on asymmetry factor

	for gi in g: # loop through asymmetry factors
		tau_cnt = 2 # keep count on changing optical depth
		for Ti in T_all: # loop through effective optical thickness
		
			# optical thickness
			taui = Ti/(1.-gi)
			# call function for results
			Ares[tau_cnt, gcnt]  = nat_act_flux.nat_act_flux(A, ai, F0, ti, taui, callfi, mu0, 1, gi, Pfunc_text)
		
			tau_cnt += 1 # keep count on T
		gcnt += 1 # keep count on asymmetry factor
		
	ig, ax0 = plt.subplots() # prepare plot
	ax0.plot(g, Ares[0, :], label = str(r'$\tau$ = ' + str(tau_all[0])))
	ax0.plot(g, Ares[1, :], label = str(r'$\tau$ = ' + str(tau_all[1])))
	ax0.plot(g, Ares[2, :], label = str(r'T = ' + str(T_all[0])))
	ax0.plot(g, Ares[3, :], label = str(r'T = ' + str(T_all[1])))
	ax0.plot(g, Ares[4, :], label = str(r'T = ' + str(T_all[2])))
	ax0.plot(g, Ares[5, :], label = str(r'T = ' + str(T_all[3])))
	ax0.plot(g, Ares[6, :], label = str(r'T = ' + str(T_all[4])))
	ax0.legend(loc = 'lower right')
	ax0.set_title('Reproduction of Eddington approximation in Figure 3 of \nShettle and Weinman (1970) to test single layer atmosphere case')
	ax0.set_ylabel(r'ALBEDO')
	ax0.set_xlabel(r'ASYMMETRY FACTOR, g')
	ax0.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f')) # x axis tick labels at two decimal places	
	ax0.tick_params(direction = 'in', which = 'both')
	plt.show()
	
	
	# ---------------------------------------
	# comparison with Table 1 of Shettle and Weinman (1970)
	theta0 = 0. # solar zenith angle
	mu0 = np.cos(theta0)
	# atmosphere optical depth
	tau_all = np.array((1., 5.))
	A_all = np.array((0., 0.25, 0.8, 1.)) # surface albedo
	callfi = -9
	# conservative atmosphere assumption (no absorption)
	ai = np.array((1.)).reshape(1)
	tau_cnt = 0 # keep count on atmosphere optical depths
	# for this case, which is comparing against Kahle (1968), doi.org/10.1086/149463
	F0 = 1.
	gi = 0. # scattering asymmetry factor for Rayleigh atmosphere
	# Rayleigh (gas) scattering, Eq. 9.85 of Jacobson (2005)
	Pfunc_text = ['(3./4.)*(1.+(np.cos(ThetaP))**2.)']
	#g = (1./2.)*(P[0]*((1.**2.)/2.)-P[-1]*((1.**2.)/2.))
	
	
	res = np.zeros((len(tau_all), len(A_all)*2)) # results matrix
	
	for taui in tau_all: # loop through optical depths
	
		A_cnt = 0 # keep count on surface albedos
		for A in A_all: # surface albedo loop
			[Fdown, Fup] = nat_act_flux.nat_act_flux(A, ai, F0, theta0, taui, callfi, mu0, 1, gi, Pfunc_text)
			res[tau_cnt, A_cnt] = Fdown
			res[tau_cnt, 4+A_cnt] = Fup
			A_cnt += 1 # keep count on surface albedos
		
		tau_cnt += 1 # keep count on atmosphere optical depths
	
	res = np.around(res, decimals=2) # round to two decimal points, ready for preparation
	
	ig, ax0 = plt.subplots() # prepare plot
	ax0.plot()
	ax0.set_title('Reproduction of Table 1 of Shettle and Weinman (1970)\nto test single layer atmosphere case')
	# results in a text box
	str0 = str(r'F down using Eddington approximation, $\tau$=1 and')
	str0 += '\n' 
	str0 += str(r'surface albedos = 0, 0.25, 0.8, 1.: ' +str(res[0, 0:4]))
	ax0.text(0.1, 0.8, str0, color='m')
	
	str0 = str(r'F up using Eddington approximation, $\tau$=1 and')
	str0 += '\n' 
	str0 += str(r'surface albedos = 0, 0.25, 0.8, 1.: ' +str(res[0, 4::]))
	ax0.text(0.1, 0.7, str0, color='m')
	
	str0 = str(r'F down using Eddington approximation, $\tau$=5 and')
	str0 += '\n' 
	str0 += str(r'surface albedos = 0, 0.25, 0.8, 1.: ' +str(res[1, 0:4]))
	ax0.text(0.1, 0.5, str0, color='m')
	
	str0 = str(r'F up using Eddington approximation, $\tau$=5 and')
	str0 += '\n' 
	str0 += str(r'surface albedos = 0, 0.25, 0.8, 1.: ' +str(res[1, 4::]))
	ax0.text(0.1, 0.4, str0, color='m')
	
	ax0.set_xlim(0, 3.5)
	ax0.set_ylim(0, 1)
	plt.show()
	
	# Figure 4 of Shettle and Weinman (1970) ------------------------------
	callfi = -10 # caller identity
	A = 0.02 # first wavelength (0.76mu) ocean albedo Table 2 Shettle and Weinman (1970)
	gi = np.array((0.553, 0.848, 0.708))
	ai = np.array((0.124, 0.997, 0.434))
	F0 = 1370 # incoming solar radiance (Wm-2)
	NL = 3 # number of vertical atmosphere layer boundaries
	tau = np.array((0.492, 20., 0.370)) # contribution of each layer to optical depth
	
	# solar zenith angle from Figure 4 of Shettle and Weinman (1970)
	# note, this should be changed to the corresponding value in Fig. 4
	theta0 = np.arccos(0.2)
	mu0 = np.cos(theta0)
	
	# call function to get downward and upward diffuse irradiances at each vertical layer
	# and the total global downward directed irradiance
	[Fdown02, Fup02, tau02] = nat_act_flux.nat_act_flux(A, ai, F0, theta0, tau, callfi, mu0, NL, gi, Pfunc_text)
	
	# Eq. 26 of Shettle and Weinman (1970)
	Gdown02 = Fdown02+mu0*np.pi*F0*np.exp(-tau02/mu0)
	
	# normalisation factor
	norm = mu0*np.pi*F0

	# convert optical depth to altitude
	alt = np.zeros((len(tau02)))
	indx = tau02<0.492 # top layer index
	alt[indx] = np.interp(tau02[indx], [0., 0.492], [5., 4.])
	indx = tau02>=0.492 # top layer index
	indx = indx*(tau02<0.492+20.)
	alt[indx] = np.interp(tau02[indx], [0.492, 0.492+20.], [4., 3.])
	indx = (tau02>=0.492+20.)
	alt[indx] = np.interp(tau02[indx], [0.492+20., 0.492+20.+0.370], [3., 0.])
	
	ig, ax0 = plt.subplots() # prepare plot
	ax0.plot(Fdown02/norm, alt, label = str(r'$F\downarrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
	ax0.plot(Gdown02/norm, alt, '--', label = str(r'$G\downarrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
	ax0.plot(Fup02/norm, alt, '--', label = str(r'$F\uparrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
	
	# next solar zenith angle for Figure 4
	# solar zenith angle from Figure 4 of Shettle and Weinman (1970)
	# note, this should be changed to the corresponding value in Fig. 4
	theta0 = np.arccos(0.6)
	mu0 = np.cos(theta0)
	
	# call function to get downward and upward diffuse irradiances at each vertical layer
	# and the total global downward directed irradiance
	[Fdown06, Fup06, tau06] = nat_act_flux.nat_act_flux(A, ai, F0, theta0, tau, callfi, mu0, NL, gi, Pfunc_text)
	
	# Eq. 26 of Shettle and Weinman (1970)
	Gdown06 = Fdown06+mu0*np.pi*F0*np.exp(-tau02/mu0)
	
	norm = mu0*np.pi*F0 # normalisation factor
	
	ax0.plot(Fdown06/norm, alt, label = str(r'$F\downarrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
	ax0.plot(Gdown06/norm, alt, '--', label = str(r'$G\downarrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
	ax0.plot(Fup06/norm, alt, '--', label = str(r'$F\uparrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))

	# next solar zenith angle for Figure 4
	# solar zenith angle from Figure 4 of Shettle and Weinman (1970)
	# note, this should be changed to the corresponding value in Fig. 4
	theta0 = np.arccos(1.)
	mu0 = np.cos(theta0)
	
	# call function to get downward and upward diffuse irradiances at each vertical layer
	# and the total global downward directed irradiance
	[Fdown1, Fup1, tau1] = nat_act_flux.nat_act_flux(A, ai, F0, theta0, tau, callfi, mu0, NL, gi, Pfunc_text)
	
	# Eq. 26 of Shettle and Weinman (1970)
	Gdown1 = Fdown1+mu0*np.pi*F0*np.exp(-tau02/mu0)
	
	norm = mu0*np.pi*F0 # normalisation factor
	
	ax0.plot(Fdown1/norm, alt, label = str(r'$F\downarrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
	ax0.plot(Gdown1/norm, alt, '--', label = str(r'$G\downarrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
	ax0.plot(Fup1/norm, alt, '--', label = str(r'$F\uparrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))

	ax0.legend(loc = 'lower right')
	ax0.set_title('Reproduction of Figure 4 of Shettle and Weinman (1970)\nto test non-conservative atmosphere case')
	ax0.set_ylabel(r'Altitude (km)')
	ax0.set_xlabel(r'$F\downarrow/\mu_{0} \pi F_{0}$ and $G\downarrow/\mu_{0} \pi F_{0}$ and $F\uparrow/\mu_{0} \pi F_{0}$')
	ax0.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f')) # x axis tick labels at two decimal places	
	ax0.tick_params(direction = 'in', which = 'both')
	plt.show() # Fig.4 of Shettle and Weinman (1970) reproduction
	
	
	
	# Figure 5 of Shettle and Weinman (1970) ---------------------------------------------
	
	callfi = -10 # caller identity
	A = 0.75 # first wavelength (0.76mu) snow albedo Table 2 Shettle and Weinman (1970)
	gi = np.array((0.553, 0.848, 0.708))
	ai = np.array((0.124, 0.997, 0.434))
	F0 = 1370 # incoming solar radiance (Wm-2)
	NL = 3 # number of vertical atmosphere layer boundaries
	tau = np.array((0.492, 20., 0.370)) # contribution of each layer to optical depth
	
	# solar zenith angle from Figure 4 of Shettle and Weinman (1970)
	# note, this should be changed to the corresponding value in Fig. 4
	theta0 = np.arccos(0.2)
	mu0 = np.cos(theta0)
	
	# call function to get downward and upward diffuse irradiances at each vertical layer
	# and the total global downward directed irradiance
	[Fdown02, Fup02, tau02] = nat_act_flux.nat_act_flux(A, ai, F0, theta0, tau, callfi, mu0, NL, gi, Pfunc_text)
	
	# Eq. 26 of Shettle and Weinman (1970)
	Gdown02 = Fdown02+mu0*np.pi*F0*np.exp(-tau02/mu0)
	
	# normalisation factor
	norm = mu0*np.pi*F0

	# convert optical depth to altitude
	alt = np.zeros((len(tau02)))
	indx = tau02<0.492 # top layer index
	alt[indx] = np.interp(tau02[indx], [0., 0.492], [5., 4.])
	indx = tau02>=0.492 # top layer index
	indx = indx*(tau02<0.492+20.)
	alt[indx] = np.interp(tau02[indx], [0.492, 0.492+20.], [4., 3.])
	indx = (tau02>=0.492+20.)
	alt[indx] = np.interp(tau02[indx], [0.492+20., 0.492+20.+0.370], [3., 0.])
	
	ig, ax0 = plt.subplots() # prepare plot
	ax0.plot(Fdown02/norm, alt, label = str(r'$F\downarrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
	ax0.plot(Gdown02/norm, alt, '--', label = str(r'$G\downarrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
	ax0.plot(Fup02/norm, alt, '--', label = str(r'$F\uparrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
	
	# next solar zenith angle for Figure 5
	# solar zenith angle from Figure 5 of Shettle and Weinman (1970)
	# note, this should be changed to the corresponding value in Fig. 4
	theta0 = np.arccos(0.6)
	mu0 = np.cos(theta0)
	
	# call function to get downward and upward diffuse irradiances at each vertical layer
	# and the total global downward directed irradiance
	[Fdown06, Fup06, tau06] = nat_act_flux.nat_act_flux(A, ai, F0, theta0, tau, callfi, mu0, NL, gi, Pfunc_text)
	
	# Eq. 26 of Shettle and Weinman (1970)
	Gdown06 = Fdown06+mu0*np.pi*F0*np.exp(-tau06/mu0)
	
	norm = mu0*np.pi*F0 # normalisation factor
	
	ax0.plot(Fdown06/norm, alt, label = str(r'$F\downarrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
	ax0.plot(Gdown06/norm, alt, '--', label = str(r'$G\downarrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
	ax0.plot(Fup06/norm, alt, '--', label = str(r'$F\uparrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))

	# next solar zenith angle for Figure 5
	# solar zenith angle from Figure 5 of Shettle and Weinman (1970)
	# note, this should be changed to the corresponding value in Fig. 4
	theta0 = np.arccos(1.)
	mu0 = np.cos(theta0)
	
	# call function to get downward and upward diffuse irradiances at each vertical layer
	# and the total global downward directed irradiance
	[Fdown1, Fup1, tau1] = nat_act_flux.nat_act_flux(A, ai, F0, theta0, tau, callfi, mu0, NL, gi, Pfunc_text)
	
	# Eq. 26 of Shettle and Weinman (1970)
	Gdown1 = Fdown1+mu0*np.pi*F0*np.exp(-tau1/mu0)
	
	norm = mu0*np.pi*F0 # normalisation factor
	
	ax0.plot(Fdown1/norm, alt, label = str(r'$F\downarrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
	ax0.plot(Gdown1/norm, alt, '--', label = str(r'$G\downarrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
	ax0.plot(Fup1/norm, alt, '--', label = str(r'$F\uparrow(\lambda$ = 0.76$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))

	ax0.legend(loc = 'lower right')
	ax0.set_title('Reproduction of Figure 5 of Shettle and Weinman (1970)\nto test non-conservative atmosphere case')
	ax0.set_ylabel(r'Altitude (km)')
	ax0.set_xlabel(r'$F\downarrow/\mu_{0} \pi F_{0}$ and $G\downarrow/\mu_{0} \pi F_{0}$ and $F\uparrow/\mu_{0} \pi F_{0}$')
	ax0.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f')) # x axis tick labels at two decimal places	
	ax0.tick_params(direction = 'in', which = 'both')
	plt.show() # Fig. 5 of Shettle and Weinman (1970) reproduction


		# Figure 6 of Shettle and Weinman (1970) blue (0.4mu) wavelength -------------------------------------------
	
	theta0s= [np.arccos(1.), np.arccos(0.6), np.arccos(0.2)]
	ig, ax0 = plt.subplots() # prepare plot (Fig. 6 Shettle and Weinman (1970) reproduction) for blue wavelength (red wavelength after)
	
	for theta0 in theta0s: # loop through solar zenith angles given in Figure 6
	
		callfi = -12 # caller identity
		A = 0.09 # 0.4mu wavelength ocean albedo Table 2 Shettle and Weinman (1970)
		ai = np.array((1., 1., 1.))
		#0.4mu column of Table 2 Shettle and Weinman (1970)
		gi = np.array((0.180, 0.848, 0.507))
		F0 = 1370 # incoming solar radiance (Wm-2)
		tau = np.array((0.291, 20., 0.346)) # contribution of each layer to optical depth
		mu0 = np.cos(theta0)
	
		NL = 3 # number of vertical atmosphere layer boundaries
		
		# increase vertical resolution
		# new number of vertical atmosphere layers
		NLn = 300
		gin = np.zeros((NLn))
		ain = np.zeros((NLn))
		# new optical depth per vertical layer
		taun = np.diff(np.arange(0., np.sum(tau)*1.00001, np.sum(tau)/NLn))
	
		indx = np.cumsum(taun)<tau[0] # original top layer
		gin[indx] = gi[0]
		ain[indx] = ai[0]
		indx = np.cumsum(taun)>=tau[0] # original middle layer
		indx = indx*(np.cumsum(taun)<(tau[0]+tau[1])) # original middle layer
		gin[indx] = gi[1]
		ain[indx] = ai[1]
		indx = np.cumsum(taun)>=(tau[0]+tau[1]) # original bottom layer
		gin[indx] = gi[2]
		ain[indx] = ai[2]
		
		# call function to get downward and upward diffuse irradiances at each vertical layer
		# and the total global downward directed irradiance
		[Fdown, Fup, taun] = nat_act_flux.nat_act_flux(A, ain, F0, theta0, taun, callfi, mu0, NLn, gin, Pfunc_text)
		
		# Eq. 26 of Shettle and Weinman (1970)
		Gdown = Fdown+mu0*np.pi*F0*np.exp(-taun/mu0)
		
		# normalisation factor
		norm = mu0*np.pi*F0
	
		# convert optical depth to altitude
		alt = np.zeros((len(taun)))
		indx = taun<0.291 # top layer index
		alt[indx] = np.interp(taun[indx], [0., 0.291], [5., 4.])
		indx = taun>=0.291 # top layer index
		indx = indx*(taun<20.291)
		alt[indx] = np.interp(taun[indx], [0.291, 20.291], [4., 3.])
		indx = (taun>=20.291)
		alt[indx] = np.interp(taun[indx], [20.291, 20.291+0.346], [3., 0.])
	
		ax0.plot(Fdown/norm, alt, label = str(r'$F\downarrow(\lambda$ = 0.4$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
		ax0.plot(Gdown/norm, alt, '--', label = str(r'$G\downarrow(\lambda$ = 0.4$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
		
	ax0.legend(loc = 'lower right')
	ax0.set_title('Reproduction of Figure 6 of Shettle and Weinman (1970)\nto test three layer atmosphere case for blue wavelength,\nred wavelength in next figure')
	ax0.set_ylabel(r'Altitude (km)')
	ax0.set_xlabel(r'$F\downarrow/\mu_{0} \pi F_{0}$ and $G\downarrow/\mu_{0} \pi F_{0}$')
	ax0.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f')) # x axis tick labels at two decimal places	
	ax0.tick_params(direction = 'in', which = 'both')
		
	# show reproduction of Fig. 6 of Shettle and Weinman (1970) for blue wavelength
	plt.show()


	# Figure 6 of Shettle and Weinman (1970) red (0.7mu) wavelength -------------------------------------------
	
	theta0s= [np.arccos(1.), np.arccos(0.6), np.arccos(0.2)]
	ig, ax0 = plt.subplots() # prepare plot (Fig. 6 Shettle and Weinman (1970) reproduction) for red wavelength (blue wavelength above)
	
	for theta0 in theta0s: # loop through solar zenith angles given in Figure 6
	
		callfi = -12 # caller identity
		A = 0.02 # 0.7mu wavelength ocean albedo Table 2 Shettle and Weinman (1970)
		ai = np.array((1., 1., 1.))
		#0.4mu column of Table 2 Shettle and Weinman (1970)
		gi = np.array((0.511, 0.848, 0.701))
		F0 = 1370 # incoming solar radiance (Wm-2)
		tau = np.array((0.069, 20., 0.169)) # contribution of each layer to optical depth
		mu0 = np.cos(theta0)
	
		NL = 3 # number of vertical atmosphere layer boundaries
		
		# increase vertical resolution
		# new number of vertical atmosphere layers
		NLn = 300
		gin = np.zeros((NLn))
		ain = np.zeros((NLn))
		# new optical depth per vertical layer
		taun = np.diff(np.arange(0., np.sum(tau)*1.00001, np.sum(tau)/NLn))
		
		indx = np.cumsum(taun)<tau[0] # original top layer
		gin[indx] = gi[0]
		ain[indx] = ai[0]
		indx = np.cumsum(taun)>=tau[0] # original middle layer
		indx = indx*(np.cumsum(taun)<(tau[0]+tau[1])) # original middle layer
		gin[indx] = gi[1]
		ain[indx] = ai[1]
		indx = np.cumsum(taun)>=(tau[0]+tau[1]) # original bottom layer
		gin[indx] = gi[2]
		ain[indx] = ai[2]
		
		# call function to get downward and upward diffuse irradiances at each vertical layer
		# and the total global downward directed irradiance
		[Fdown, Fup, taun] = nat_act_flux.nat_act_flux(A, ain, F0, theta0, taun, callfi, mu0, NLn, gin, Pfunc_text)
		
		# Eq. 26 of Shettle and Weinman (1970)
		Gdown = Fdown+mu0*np.pi*F0*np.exp(-taun/mu0)
		
		# normalisation factor
		norm = mu0*np.pi*F0
	
		# convert optical depth to altitude
		alt = np.zeros((len(taun)))
		indx = taun<0.069 # top layer index
		alt[indx] = np.interp(taun[indx], [0., 0.069], [5., 4.])
		indx = taun>=0.069 # middle layer index
		indx = indx*(taun<20.069)
		alt[indx] = np.interp(taun[indx], [0.069, 20.069], [4., 3.])
		indx = (taun>=20.069) # bottom layer index
		alt[indx] = np.interp(taun[indx], [20.069, 20.069+0.169], [3., 0.])
	
		ax0.plot(Fdown/norm, alt, label = str(r'$F\downarrow(\lambda$ = 0.7$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
		ax0.plot(Gdown/norm, alt, '--', label = str(r'$G\downarrow(\lambda$ = 0.7$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
		
	ax0.legend(loc = 'lower right')
	ax0.set_title('Reproduction of Figure 6 of Shettle and Weinman (1970)\nto test three layer atmosphere case for red wavelength,\nblue wavelength in previous figure')
	ax0.set_ylabel(r'Altitude (km)')
	ax0.set_xlabel(r'$F\downarrow/\mu_{0} \pi F_{0}$ and $G\downarrow/\mu_{0} \pi F_{0}$')
	ax0.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f')) # x axis tick labels at two decimal places	
	ax0.tick_params(direction = 'in', which = 'both')
		
	# show reproduction of Fig. 6 of Shettle and Weinman (1970) for red wavelength
	plt.show()

	
	# Figure 7 of Shettle and Weinman (1970) blue (0.4mu) wavelength -------------------------------------------
	
	theta0s= [np.arccos(1.), np.arccos(0.6), np.arccos(0.2)]
	ig, ax0 = plt.subplots() # prepare plot (Fig. 7 Shettle and Weinman (1970) reproduction) for blue wavelength (red wavelength below)
	
	for theta0 in theta0s: # loop through solar zenith angles given in Figure 7
	
		callfi = -12 # caller identity
		A = 0.78 # 0.4mu wavelength snow albedo Table 2 Shettle and Weinman (1970)
		ai = np.array((1., 1., 1.))
		#0.4mu column of Table 2 Shettle and Weinman (1970)
		gi = np.array((0.180, 0.848, 0.507))
		F0 = 1370 # incoming solar radiance (Wm-2)
		tau = np.array((0.291, 20., 0.346)) # contribution of each layer to optical depth
		mu0 = np.cos(theta0)
	
		NL = 3 # number of vertical atmosphere layer boundaries
		
		# increase vertical resolution
		# new number of vertical atmosphere layers
		NLn = 300
		gin = np.zeros((NLn))
		ain = np.zeros((NLn))
		# new optical depth per vertical layer
		taun = np.diff(np.arange(0., np.sum(tau)*1.00001, np.sum(tau)/NLn))
		
		indx = np.cumsum(taun)<tau[0] # original top layer
		gin[indx] = gi[0]
		ain[indx] = ai[0]
		indx = np.cumsum(taun)>=tau[0] # original middle layer
		indx = indx*(np.cumsum(taun)<(tau[0]+tau[1])) # original middle layer
		gin[indx] = gi[1]
		ain[indx] = ai[1]
		indx = np.cumsum(taun)>=(tau[0]+tau[1]) # original bottom layer
		gin[indx] = gi[2]
		ain[indx] = ai[2]
		
		# call function to get downward and upward diffuse irradiances at each vertical layer
		# and the total global downward directed irradiance
		[Fdown, Fup, taun] = nat_act_flux.nat_act_flux(A, ain, F0, theta0, taun, callfi, mu0, NLn, gin, Pfunc_text)
		
		# Eq. 26 of Shettle and Weinman (1970)
		Gdown = Fdown+mu0*np.pi*F0*np.exp(-taun/mu0)
		
		# normalisation factor
		norm = mu0*np.pi*F0
	
		# convert optical depth to altitude
		alt = np.zeros((len(taun)))
		indx = taun<0.291 # top layer index
		alt[indx] = np.interp(taun[indx], [0., 0.291], [5., 4.])
		indx = taun>=0.291 # middle layer index
		indx = indx*(taun<20.291)
		alt[indx] = np.interp(taun[indx], [0.291, 20.291], [4., 3.])
		indx = (taun>=20.291) # bottom layer index
		alt[indx] = np.interp(taun[indx], [20.291, 20.291+0.346], [3., 0.])
	
		ax0.plot(Fdown/norm, alt, label = str(r'$F\downarrow(\lambda$ = 0.4$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
		ax0.plot(Gdown/norm, alt, '--', label = str(r'$G\downarrow(\lambda$ = 0.4$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
		
	ax0.legend(loc = 'lower right')
	ax0.set_title('Reproduction of Figure 7 of Shettle and Weinman (1970)\nto test three layer atmosphere case for blue wavelength,\nred wavelength in next figure')
	ax0.set_ylabel(r'Altitude (km)')
	ax0.set_xlabel(r'$F\downarrow/\mu_{0} \pi F_{0}$ and $G\downarrow/\mu_{0} \pi F_{0}$')
	ax0.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f')) # x axis tick labels at two decimal places	
	ax0.tick_params(direction = 'in', which = 'both')
		
	# show reproduction of Fig. 7 of Shettle and Weinman (1970) for blue wavelength
	plt.show()
	
	
	# Figure 7 of Shettle and Weinman (1970) red (0.7mu) wavelength -------------------------------------------
	
	theta0s= [np.arccos(1.), np.arccos(0.6), np.arccos(0.2)]
	ig, ax0 = plt.subplots() # prepare plot (Fig. 7 Shettle and Weinman (1970) reproduction) for red wavelength (blue wavelength above)
	
	for theta0 in theta0s: # loop through solar zenith angles given in Figure 7
	
		callfi = -12 # caller identity
		A = 0.78 # 0.4mu wavelength snow albedo Table 2 Shettle and Weinman (1970)
		ai = np.array((1., 1., 1.))
		#0.4mu column of Table 2 Shettle and Weinman (1970)
		gi = np.array((0.511, 0.848, 0.701))
		F0 = 1370 # incoming solar radiance (Wm-2)
		tau = np.array((0.069, 20., 0.169)) # contribution of each layer to optical depth
		mu0 = np.cos(theta0)
	
		NL = 3 # number of vertical atmosphere layer boundaries
		
		# increase vertical resolution
		# new number of vertical atmosphere layers
		NLn = 300
		gin = np.zeros((NLn))
		ain = np.zeros((NLn))
		# new optical depth per vertical layer
		taun = np.diff(np.arange(0., np.sum(tau)*1.00001, np.sum(tau)/NLn))
		
		indx = np.cumsum(taun)<tau[0] # original top layer
		gin[indx] = gi[0]
		ain[indx] = ai[0]
		indx = np.cumsum(taun)>=tau[0] # original middle layer
		indx = indx*(np.cumsum(taun)<(tau[0]+tau[1])) # original middle layer
		gin[indx] = gi[1]
		ain[indx] = ai[1]
		indx = np.cumsum(taun)>=(tau[0]+tau[1]) # original bottom layer
		gin[indx] = gi[2]
		ain[indx] = ai[2]
		
		# call function to get downward and upward diffuse irradiances at each vertical layer
		# and the total global downward directed irradiance
		[Fdown, Fup, taun] = nat_act_flux.nat_act_flux(A, ain, F0, theta0, taun, callfi, mu0, NLn, gin, Pfunc_text)
		
		# Eq. 26 of Shettle and Weinman (1970)
		Gdown = Fdown+mu0*np.pi*F0*np.exp(-taun/mu0)
		
		# normalisation factor
		norm = mu0*np.pi*F0
	
		# convert optical depth to altitude
		alt = np.zeros((len(taun)))
		indx = taun<0.069 # top layer index
		alt[indx] = np.interp(taun[indx], [0., 0.069], [5., 4.])
		indx = taun>=0.069 # middle layer index
		indx = indx*(taun<20.069)
		alt[indx] = np.interp(taun[indx], [0.069, 20.069], [4., 3.])
		indx = (taun>=20.069) # bottom layer index
		alt[indx] = np.interp(taun[indx], [20.069, 20.069+0.169], [3., 0.])
	
		ax0.plot(Fdown/norm, alt, label = str(r'$F\downarrow(\lambda$ = 0.7$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
		ax0.plot(Gdown/norm, alt, '--', label = str(r'$G\downarrow(\lambda$ = 0.7$\mu,\, \mu_0=$' +str(np.around(mu0, decimals=1)) + ')' ))
		
	ax0.legend(loc = 'lower right')
	ax0.set_title('Reproduction of Figure 7 of Shettle and Weinman (1970)\nto test three layer atmosphere case for red wavelength,\nblue wavelength in previous figure')
	ax0.set_ylabel(r'Altitude (km)')
	ax0.set_xlabel(r'$F\downarrow/\mu_{0} \pi F_{0}$ and $G\downarrow/\mu_{0} \pi F_{0}$')
	ax0.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f')) # x axis tick labels at two decimal places	
	ax0.tick_params(direction = 'in', which = 'both')
		
	# show reproduction of Fig. 7 of Shettle and Weinman (1970) for red wavelength
	plt.show()
	
	

test_nat_act_flux() # call function