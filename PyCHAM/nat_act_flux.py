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
'''module for actinic flux by natural light'''
# uses Eq. 1 of Madronich (1987) (doi.org/10.1029/JD092iD08p09740),
# with radiance given by the two-stream method of Shettle and 
# Weinman (1970), 
# (doi.org/10.1175/1520-0469(1970)027<1048:TTOSIT>2.0.CO;2)
# which is recommended by 
# Joseph et al.
# (1976) (doi.org/10.1175/1520-0469(1976)033<2452:TDEAFR>2.0.CO;2)
# Note the delta-Eddington phase function for radiance scattering 
# given by Joseph et al.
# (1976) (doi.org/10.1175/1520-0469(1976)033<2452:TDEAFR>2.0.CO;2)
# is recommended by Madronich (1987) (doi.org/10.1029/JD092iD08p09740)

# For the two-stream method Shettle and 
# Weinman (1970), the simultaneous equations given by their 
# Non-conservative case equations (12a-15b) are
# solved using a numpy iterative solver 

import numpy as np
import datetime
try:
	import scatt_Pfunc # the function for scattering phase function
except:
	import os
	if os.path.exists('scatt_Pfunc'): # remove any bad functions
		os.remove(scatt_Pfunc)
import importlib # for reloading scattering phase function

def nat_act_flux(A, a, F0, theta, tau, callf, mu0, NL, g, Pfunc_text):

	# inputs: ---------------------------------------------
	# A - surface (ground) albedo (fraction of radiation (reflected 0-1))
	# a - single scattering albedo (probability that radiation scattered rather than absorbed) (fraction 0-1)
	# F0 - solar irradiance at top of atmosphere (W/m2)
	# theta - solar zenith angle (angle between the perpendicular to planet surface 
	#	and incident solar radiance at top of atmosphere) (radians)
	# tau - optical depth of each vertical layer of atmosphere (not cumulative) (natural logarithm 
	# 	of ratio of incident radiation/transmitted radiation)
	# callf - flag for the calling function (0 for unit testing)
	# mu0 - cosine of the angle between the perpendicular of the 
	# surface of interest and the beam being investigated (the incident radiation angle)
	# equals the solar zenith angle when the surface of interest is Earth surface
	# NL - number of vertical layers in atmosphere
	# g - the asymmetry factor for scattering (below Eq. 4 
	# 	Shettle and Weinman (1970)) g>0 means forward scattering favoured,
	# 	g=0 means forward and backward scattering equally favoured and
	# 	g<0 means backward scattering favoured
	# Pfunc_text - text for writing a function for calculating the scattering phase function
	# -------------------------------------------------------
	
	# function for calculating C1 and C2 of Shettle and Weinman (1970) Eqs. 12-15
	def C_calc(ai, gi, tau, A, mu0, NL, F0, P_all): # function definition
		
		# inputs: -------------
		# ai - single scattering albedo per vertical layer
		# gi - scattering asymmetry factor per vertical layer
		# tau - optical depth
		# A - surface albedo
		# mu0 - cosine of the solar zenith angle
		# NL - number of vertical atmospheric layers
		# F0 - top-of-the-atmosphere flux (W/m2)
		# P_all - phase function for downward scattering per vertical layer
		# -----------------------
		
		I0 = I1 = ki = Fdown = Fup = B2 = 0. # fillers for return values
		x = np.array((0,0)).reshape(1,2) # fillers for return values
		# if single scattering albedo below 1 
		# (i.e. some absorption occurs in addition to scattering)
		if (np.any(a < 1.)):
		
			# empty A array for Ax=B problem
			Amat = np.zeros((NL*2, NL*2))
			# empty B array for Ax=B problem
			Barr = np.zeros((NL*2))
			k_all = np.zeros((NL)) # record for k values per vertical layer
			alp_all = np.zeros((NL)) # record for alpha values per vertical layer
			bet_all = np.zeros((NL)) # record for alpha values per vertical layer
		
			# loop through vertical layers
			for NLi in range(NL):
				
				# constants below Eq. 12 of Shettle and Weinman (1970)
				# this vertical atmospheric layer
				ki = (3.*(1.-ai[NLi])*(1.-ai[NLi]*gi[NLi]))**(1./2.)
				pi = (3.*(1.-ai[NLi])/(1.-ai[NLi]*gi[NLi]))**(1./2.)
				alp = 3.*ai[NLi]*F0*mu0**2.*(1.+gi[NLi]*(1.-ai[NLi]))/(4.*(1.-ki**2.*mu0**2.))
				bet = 3.*ai[NLi]*F0*mu0*(1.+3.*gi[NLi]*(1.-ai[NLi])*mu0**2)/(4.*(1.-ki**2.*mu0**2.))
				
				k_all[NLi] = ki
				alp_all[NLi] = alp
				bet_all[NLi] = bet
				
				if (NLi > 0): # if below top layer
				
					# neighbouring vertical atmospheric layer above this one
					k2i = (3.*(1.-ai[NLi-1])*(1.-ai[NLi-1]*gi[NLi-1]))**(1./2.)
					p2i = (3.*(1.-ai[NLi-1])/(1.-ai[NLi-1]*gi[NLi-1]))**(1./2.)
					alp2 = 3.*ai[NLi-1]*F0*mu0**2.*(1.+gi[NLi-1]*(1.-ai[NLi-1]))/(4.*(1.-k2i**2.*mu0**2.))
					bet2 = 3.*ai[NLi-1]*F0*mu0*(1.+3.*gi[NLi-1]*(1.-ai[NLi-1])*mu0**2)/(4.*(1.-k2i**2.*mu0**2.))
				
				if (NLi < NL-1): # if above bottom layer
				
					# neighbouring vertical layer below this one
					k3i = (3.*(1.-ai[NLi+1])*(1.-ai[NLi+1]*gi[NLi+1]))**(1./2.)
					p3i = (3.*(1.-ai[NLi+1])/(1.-ai[NLi+1]*gi[NLi+1]))**(1./2.)
					alp3 = 3.*ai[NLi+1]*F0*mu0**2.*(1.+gi[NLi+1]*(1.-ai[NLi+1]))/(4.*(1.-k3i**2.*mu0**2.))
					bet3 = 3.*ai[NLi+1]*F0*mu0*(1.+3.*gi[NLi+1]*(1.-ai[NLi+1])*mu0**2)/(4.*(1.-k3i**2.*mu0**2.))

				
				if (NLi == 0): # top vertical layer of atmosphere
					# Eq. 13 Shettle and Weinman (1970)
					Amat[0, 0:2] = [1.+2.*pi/3., 1.-2.*pi/3.]
					Barr[0] = (alp+2.*bet/3.)
					
					if (NL > 1): # if more than one layer
						# Eq. 12b and 15a Shettle and Weinman (1970)
						Amat[1, 0:4] = [P_all[NLi]*np.exp(-ki*tau[NLi]), -P_all[NLi]*np.exp(ki*tau[NLi]), -P_all[NLi+1]*np.exp(-k3i*tau[NLi]), P_all[NLi+1]*np.exp(k3i*tau[NLi])]
						Barr[1] = bet*np.exp(-tau[NLi]/mu0)-bet3*np.exp(-tau[NLi]/mu0)
				
				if (NLi != 0 and NLi != NL-1): # middle vertical layers of atmosphere
					# Eq. 12a and 15a Shettle and Weinman (1970)
					Amat[NLi*2, NLi*2-2:NLi*2+2] = [np.exp(-k2i*tau[NLi-1]), np.exp(k2i*tau[NLi-1]), -np.exp(-ki*tau[NLi-1]), -np.exp(ki*tau[NLi-1])]
					Barr[NLi*2] = alp2*np.exp(-tau[NLi-1]/mu0)-alp*np.exp(-tau[NLi-1]/mu0)
					# Eq. 12b and 15a Shettle and Weinman (1970)
					Amat[NLi*2+1, NLi*2:NLi*2+4] = [P_all[NLi]*np.exp(-ki*tau[NLi]), -P_all[NLi]*np.exp(ki*tau[NLi]), -P_all[NLi+1]*np.exp(-k3i*tau[NLi]), P_all[NLi+1]*np.exp(k3i*tau[NLi])]
					Barr[NLi*2+1] = bet*np.exp(-tau[NLi]/mu0)-bet3*np.exp(-tau[NLi]/mu0)
				
				if (NLi == NL-1): # bottom vertical layer of atmosphere
					# Eq. 14 Shettle and Weinman (1970)
					Amat[-1, -2::] = [(1.-A-2.*(1.+A)*pi/3.)*np.exp(-ki*tau[NLi]), (1.-A+2.*(1.+A)*pi/3.)*np.exp(ki*tau[NLi])]
					Barr[-1] = ((1.-A)*alp-2.*(1.+A)*bet/3.+A*mu0*F0)*np.exp(-tau[-1]/mu0)
					
					if (NL > 1): # if more than one layer
						# Eq. 12a and 15a Shettle and Weinman (1970) (assume tau = 1)
						Amat[-2, -4::] = [np.exp(-k2i*tau[NLi-1]), np.exp(k2i*tau[NLi-1]), -np.exp(-ki*tau[NLi-1]), -np.exp(ki*tau[NLi-1])]					
						Barr[-2] = alp2*np.exp(-tau[NLi-1]/mu0)-alp*np.exp(-tau[NLi-1]/mu0)
				
			
			# get C1 and C2 values for each vertical layer of atmosphere
			x = np.linalg.solve(Amat, Barr) # solve simultaneous equations of the form Ax = b
			# reshape so that C1 and C2 for each vertical layer are on the same row
			x = x.reshape(NL, 2)
			
			# get radiance I1 and I0 from Eqs. 12a and 12b Shettle and Weinman (1970)
			I0 = x[:, 0]*np.exp(-k_all*tau)+x[:, 1]*np.exp(k_all*tau)-alp_all*np.exp(-tau/mu0)
			I1 = P_all*(x[:, 0]*np.exp(-k_all*tau)-x[:, 1]*np.exp(k_all*tau))-bet_all*np.exp(-tau/mu0)
			
			# Eq. 8 of Shettle and Weinman (1970) for downward irradiance
			# at each vertical layer
			Fdown = np.pi*(I0+2./3.*I1)
			Fup = np.pi*(I0-2./3.*I1)
		
		if (np.any(a == 1.)): # conservative atmosphere
		
			tau_sum = np.cumsum(tau) # cumulative optical depth through atmosphere
			T = np.sum((1.-g)*tau) # Eq. 18 of Shettle and Weinman (1970)
			B2 = (3.*mu0*F0*(1.-A)*(2.+3.*mu0+(2.-3.*mu0)*np.exp(-tau_sum[-1]/mu0)))/(4.*(4.+3.*(1.-A)*T))
			# Eq. 17b of Shettle and Weinmann (1970)
			B1 = (3.*mu0**2./4.+mu0/2.)*F0-2.*B2/3.
			# Eq. 16a of Shettle and Weinmann (1970) for radiance
			I0 = B1-(3./4.)*mu0**2.*F0*np.exp(-tau/mu0)-B2*T
			
			# Eq. 16b of Shettle and Weinmann (1970) for radiance
			I1 = B2-(3./4.)*mu0*F0*np.exp(-tau/mu0)

			# downward diffuse irradiance at Earth surface from Eq. 8 
			# of Shettle and Weinman (1970)
			Fdown = np.pi*(I0+2./3.*I1)
			# downward diffuse irradiance at Earth surface from Eq. 8 
			# of Shettle and Weinman (1970)
			Fup = np.pi*(I0-2./3.*I1)
		
		return(I0, I1, x[:, 0], x[:, 1], ki, Fdown, Fup, B2)
	
	# write the function for scattering phase function
	f = open('PyCHAM/scatt_Pfunc.py', mode='w')
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
	f.write('\'\'\'module for calculating scattering phase function\'\'\'\n')
	f.write('# generated by nat_act_flux, with content provided by user \n')
	f.write('\n')
	f.write('# File Created at %s\n' %(datetime.datetime.now()))
	f.write('\n')
	f.write('import numpy as np \n')
	f.write('\n')
	# following part is the function (there should be an indent at the start of each line)
	# suggest 1 Tab
	f.write('def Pfunc(theta):\n')
	f.write('	\n')
	f.write('	# inputs -------------------\n')
	f.write('		theta - the downward scattering angle - equal to the solar zenith angle\n')
	f.write('	-------------------------------\n')
	f.write('	\n')
	if (Pfunc_text != []): # if not an empty list
		# looping through list that contains the lines for the scattering phase function module
		for line in Pfunc_text:
			f.write(	str(line + '\n'))
		f.write('	return(Pall)')
	# if an empty list, then default to delta-Eddington scattering phase 
	# function of Joseph et al. (1976) 
	# (doi.org/10.1175/1520-0469(1976)033<2452:TDEAFR>2.0.CO;2) 
	# which is recommended by Madronich (1987) (doi.org/10.1029/JD092iD08p09740)
	else:
		f.write('	# phase function for scattering as a function of scattering angle \n')
		f.write('	# Eq. 3 of Shettle and Weinman (1970) \n')
		f.write('	# note that \'omega\' is in place of \'a\' in Eq. 3 but \'a\' is used \n')
		f.write('	# elsewhere in Shettle and Weinman (1970) \n')
		f.write('	# backward scattering \n')
		f.write('	P0 = 1.+a[0]*np.cos(np.pi)\n')
		f.write('	# forward scattering \n')
		f.write('	P1 = 1.+a[0]*np.cos(0.)\n')
		f.write('	\n')
		f.write('	# g, the asymmetry factor for scattering (below Eq. 4\n')
		f.write('	# Shettle and Weinman (1970)): \n')
		f.write('	# g (asymmetry factor of the original phase function\n')
		f.write('	# Eq. 2a of Joseph et al. (1976))\n')
		f.write('	g = (1./2.)*(P1*((1.**2.)/2.)-P0*((1.**2.)/2.))\n')
		f.write('	\n')
		f.write('	# fractional scattering in the forward peak\n')
		f.write('	# Eq. 5a of Joseph et al. (1976)\n')
		f.write('	f = g**2.\n')
		f.write('	\n')
		f.write('	# g\' for the delta-Eddington phase function (phase function \n')
		f.write('	# approximation by a Dirac delta function forward scatter \n')
		f.write('	# peak and a two-term expansion of the phase function)\n')
		f.write('	# Eq. 2b of Joseph et al. (1976)\n')
		f.write('	gprim = (g-f)/(1.-f)\n')
		f.write('	\n')
		f.write('	# note further discussion on the role of the \n')
		f.write('	# Dirac delta function is given in Crosbie \n')
		f.write('	# and Davidson (1985)\n')
		f.write('	# doi.org/10.1016/0022-4073(85)90200-6\n')
		f.write('	\n')
		f.write('	if (theta < np.pi/5.):\n')
		f.write('		Ddel = 1.\n')
		f.write('	\n')
		f.write('	# the phase function for scattering of light (probability \n')
		f.write('	# of light being scattered in a particular direction)	\n')
		f.write('	\n')
		f.write('	P = 2.*f*Ddel*(1.-np.cos(theta))+(1.-f)*(1.+3.*gprim*np.cos(theta))\n')
		f.write('	\n')
		f.write('	return(P, g)\n')
		f.close()
	
	try:
		
		importlib.reload(scatt_Pfunc) # ensure latest version uploaded

	except:
		import os
		if os.path.exists('scatt_Pfunc'):
			os.remove(scatt_Pfunc) # remove bad function
		erf = 1
		err_mess = 'Error: bad phase function calculation, please check relevant model variable inputs and see README for further guidance'
	



	# if reproducing Fig. 1 of Shettle and Weinman (1970) ---------------
	if (callf == -1):
		
		# note that we only consider reproducing the
		# Eddington approximation of Fig.1, not the Irvine (1968)
		# solution
		
		NL = 1 # number of vertical layers in atmosphere
		
		# to reproduce Fig. 1 of Shettle and Weinmann (1970)
		# note that g is an input to this module
		if (a[0] < 1.):
		
			# get the values of C1 and C2
			[I0, I1, C1, C2, ki, Fdown, Fup, B2] = C_calc(a, g, tau, A, mu0, NL, F0, 0)

			# calculate atmospheric albedo using their Eq. 20
			atmosA = 2.*(C1+C2)/(mu0*F0)-3.*a*mu0*(1.+g-a*g)/(2.*(1.-ki**2.*mu0**2))
			
		if (a[0] == 1.):
			# get the B2 value
			[I0, I1, C1, C2, ki, Fdown, Fup, B2] = C_calc(a, g, tau, A, mu0, NL, F0, 0)
			# calculate atmospheric albedo using their Eq. 22
			atmosA = 1.-(4.*B2/(3.*mu0*F0))
		
		return(atmosA)
	
	# if reproducing Fig. 2 of Shettle and Weinman (1970), note we 
	# only address the Eddington approximation here
	if (callf == -2 or callf == -3):		
		
		# get the B2 value for a conservative atmosphere
		[I0, I1, C1, C2, ki, Fdown, Fup, B2] = C_calc(a, g, tau, A, mu0, NL, F0, 0)
		# calculate atmospheric albedo using their Eq. 22
		atmosA = 1.-(4.*B2/(3.*mu0*F0))
		
		return(atmosA)
	
	# if reproducing Table 1 of Shettle and Weinman (1970), which 
	# compares against the results of Kahle (1968):
	# doi.org/10.1086/149463
	if (callf == -9):
	
		[I0, I1, C1, C2, ki, Fdown, Fup, B2] = C_calc(a, g, tau, A, mu0, NL, F0, 0)
			
		# upward  irradiance at top of atmosphere - the amount
		# of incoming irradiance that does not reach Earth
		# surface directly or by forward scattering (diffusely)
		Fup = np.pi*(F0-(I0+2./3.*I1+mu0*F0*np.exp(-tau/mu0)))
			
		# upward irradiance at top of atmosphere - the amount
		# of incoming irradiance that does reach Earth
		# surface directly or by forward scattering (diffusely)
		# and is then reflected back up by Earth surface
		#(Eq. 10 Shettle and Weinman (1970))
		Fup += A*np.pi*(I0+2./3.*I1+mu0*F0*np.exp(-tau/mu0))
			
		return(Fdown, Fup)

		
	# if reproducing Jacobson (2005) phase function figure
	if (callf == 0):
		import matplotlib.pyplot as plt
		
		fig, ax0 = plt.subplots() # prepare plot
		
		Pt = phase_func(ThetaP)
		x = np.cos(ThetaP)*Pt
		y = np.sin(ThetaP)*Pt
		ax0.plot(x, y, '-c', label = r'$P(\Theta)$ being used by model') # plot the phase function above the x=0 line
		ax0.plot(x, -y, '-c') # plot the phase function below the x=0 line
		
		Pt = 1.+a*np.cos(ThetaP)
		x = np.cos(ThetaP)*Pt
		y = np.sin(ThetaP)*Pt
		ax0.plot(x, y, '--y', label = r'$P(\Theta) = 1+a(cos(\Theta))$') # plot the phase function above the x=0 line
		ax0.plot(x, -y, '--y') # plot the phase function below the x=0 line
		
		Pt = (3./4.)*(1.+(np.cos(ThetaP))**2.) # Eq. 9.85 Jacobson (2005) - Rayleigh (gas) scattering
		x = np.cos(ThetaP)*Pt
		y = np.sin(ThetaP)*Pt
		ax0.plot(x, y, '--k', label = r'$P(\Theta) = \frac{3}{4}(1+(cos\Theta)^{2})$') # plot the phase function above the x=0 line
		ax0.plot(x, -y, '--k') # plot the phase function below the x=0 line
		
		Pt = np.ones((len(ThetaP))) # Eq. 9.84 Jacobson (2005)
		x = np.cos(ThetaP)*Pt
		y = np.sin(ThetaP)*Pt
		ax0.plot(x, y, '--b', label = r'$P(\Theta) = 1$') # plot the phase function above the x=0 line
		ax0.plot(x, -y, '--b') # plot the phase function below the x=0 line
		
		# scattering angle text box
		ax0.text(0.1, 0.05, '$\Theta$', color='m')
		# add lines for scattering angle
		ax0.plot([0., 0.1], [0., 0.1], 'm')
		ax0.plot([0., 0.1], [0., 0.], 'm')
		ax0.plot((np.cos(ThetaP[ThetaP<np.pi/4.]))*0.1, (np.sin(ThetaP[ThetaP<np.pi/4.]))*0.1, 'm')
		ax0.legend(loc = 'lower left')
		ax0.set_title('Reproduction of Figure 9.29 of Jacobson (2005)\nto test scattering phase function')
		ax0.tick_params(direction = 'in', which = 'both')
		plt.show()
	 # --------------------------------------------------------

	mu = np.cos(theta) # page 1 of Shettle and Weinman (1970)

		
	if (callf == -10): # Figure 4 and 5 Shettle and Weinman (1970)
		
		NL = 300 # increase resolution of vertical layers
		# remember original scattering asymmetry factor
		g0 = np.zeros((len(g)))
		g0[:] = g[:]
		# remember original contribution of tau from each layer
		tau0 = np.zeros((len(tau)))
		tau0[:] = tau[:]
		# remember original single scattering albedo for each layer
		ai0 = np.zeros((len(tau)))
		ai0[:] = a[:]
		# generate higher resolution optical depth, note this is cumulative
		# (increasing from top of atmosphere down)
		tau = np.arange(0., sum(tau)*1.000001, sum(tau)/(NL-1))
		gi = np.zeros((NL))
		ai = np.zeros((NL))
		for zindx in range(len(tau0)): # loop through original values
			if (zindx == 0): # top original layer
				indx = tau<np.sum(tau0[0:zindx+1])
			if (zindx != 0 and zindx != (len(tau0)-1)): # middle original layers
				indx = tau<np.sum(tau0[0:zindx+1])
				indx = indx*tau>=np.sum(tau0[0:zindx])
			if (zindx == (len(tau0)-1)): # bottom original layer
				indx = tau>=np.sum(tau0[0:zindx])
						
			gi[indx] = g0[zindx]
			ai[indx] = ai0[zindx] # single scattering albedo for this layer
		
		# phase function for scattering as a function of scattering angle
		# Eq. 3 of Shettle and Weinman (1970),
		# note we use the solar zenith angle as this equals the 
		# downward angle
		# note also that Shettle and Weinman (1970) don't state these
		# phase functions directly, so the values here for each layer
		# have been derived by fitting to Fig. 4 of Shettle and Weinman (1970)
		P_all = np.zeros((len(tau)))
		indx =  (gi == 0.553)
		P_all[indx] = 1.7
		indx =  (gi == 0.848)
		P_all[indx] = 0.3
		indx =  (gi == 0.708)
		P_all[indx] = 1.+ai[indx]*np.cos(theta)
		
		# get fluxes
		[I0, I1, C1, C2, ki, Fdown, Fup, B2] = C_calc(ai, gi, tau, A, mu0, NL, F0, P_all)
			
		return(Fdown, Fup, tau)
	
	# if single scattering albedo equal to 1 
	# (i.e. no absorption occurs, only refraction)
	# note that using single scattering albebo = 1
	# for the above equations for I0 and I1 gives
	# a singular matrix that cannot be inverted
	if (np.any(a == 1)):
	
		tau_sum = np.cumsum(tau) # cumulative optical depth through atmosphere
		
		T = ((1.-g)*tau) # Eq. 18 of Shettle and Weinman (1970)
		T = np.cumsum(T) # integrate over atmosphere depth
		
		# Eq. 17a of Shettle and Weinman (1970)
		B2 = (3.*mu0*F0*(1.-A)*(2.+3.*mu0+(2.-3.*mu0)*np.exp(-tau_sum[-1]/mu0)))/(4.*(4.+3.*(1.-A)*T[-1]))
		
		# Eq. 17b of Shettle and Weinman (1970)
		B1 = ((3.*mu0**2.)/4.+mu0/2.)*F0-2.*B2/3.
		
		# Eq. 16a of Shettle and Weinmann (1970) for radiance
		I0 = B1-(3./4.)*mu0**2.*F0*np.exp(-tau_sum/mu0)-B2*T
			
		# Eq. 16b of Shettle and Weinmann (1970) for radiance
		I1 = B2-(3./4.)*mu0*F0*np.exp(-tau_sum/mu0)
		
		# Eq. 8 of Shettle and Weinmann (1970)
		# total downward diffuse irradiance
		Fdown = np.pi*(I0+2.*I1/3.)
		
		# downward diffuse irradiance at Earth surface from Eq. 8 
		# of Shettle and Weinman (1970)
		Fup = np.pi*(I0-2./3.*I1)
		
		return(Fdown, Fup, tau_sum)
		
		
		
			
	# diffuse radiance at Earth surface (W/m2) using the Eddington 
	# approximation (Eq. 2 of Shettle and Weinman (1970))
	L = I0 + I1*mu

	
	if (callf == -5): # Figure 3 of Shettle and Weinman (1970)
		T = 2. # effective optical thickness
		g_all = np.arange(-1., 1.*1.00001, 0.1)
		gi  = 0 # keep count on asymmetry factors
		atmosA = np.zeros((len(g_all)))
		for g in g_all: # loop through asymmetry factors
			tau = T/(1.-g) # Eq. 18 of Shettle and Weinmann (1970)
			# Eq. 17a of Shettle and Weinmann (1970)
			B2 = (3.*mu0*F0*(1.-A)*(2.+3.*mu0+(2.-3.*mu0)*np.exp(-tau/mu0)))/(4.*(4.+3.*(1.-A)*T))
			# calculate atmospheric albedo using their Eq. 22
			atmosA[gi] = 1.-(4.*B2/(3.*mu0*F0))
			gi  += 1 # keep count on asymmetry factors
	
	if (callf == -6): # Figure 3 of Shettle and Weinman (1970)
		T = 4. # effective optical thickness
		g_all = np.arange(-1., 1.*1.00001, 0.1)
		gi  = 0 # keep count on asymmetry factors
		atmosA = np.zeros((len(g_all)))
		for g in g_all: # loop through asymmetry factors
			tau = T/(1.-g) # Eq. 18 of Shettle and Weinmann (1970)
			# Eq. 17a of Shettle and Weinmann (1970)
			B2 = (3.*mu0*F0*(1.-A)*(2.+3.*mu0+(2.-3.*mu0)*np.exp(-tau/mu0)))/(4.*(4.+3.*(1.-A)*T))
			# calculate atmospheric albedo using their Eq. 22
			atmosA[gi] = 1.-(4.*B2/(3.*mu0*F0))
			gi  += 1 # keep count on asymmetry factors
	
	if (callf == -7): # Figure 3 of Shettle and Weinman (1970)
		T = 8. # effective optical thickness
		g_all = np.arange(-1., 1.*1.00001, 0.1)
		gi  = 0 # keep count on asymmetry factors
		atmosA = np.zeros((len(g_all)))
		for g in g_all: # loop through asymmetry factors
			tau = T/(1.-g) # Eq. 18 of Shettle and Weinmann (1970)
			# Eq. 17a of Shettle and Weinmann (1970)
			B2 = (3.*mu0*F0*(1.-A)*(2.+3.*mu0+(2.-3.*mu0)*np.exp(-tau/mu0)))/(4.*(4.+3.*(1.-A)*T))
			# calculate atmospheric albedo using their Eq. 22
			atmosA[gi] = 1.-(4.*B2/(3.*mu0*F0))
			gi  += 1 # keep count on asymmetry factors
	
	if (callf == -8): # Figure 3 of Shettle and Weinman (1970)
		T = 16. # effective optical thickness
		g_all = np.arange(-1., 1.*1.00001, 0.1)
		gi  = 0 # keep count on asymmetry factors
		atmosA = np.zeros((len(g_all)))
		for g in g_all: # loop through asymmetry factors
			tau = T/(1.-g) # Eq. 18 of Shettle and Weinmann (1970)
			# Eq. 17a of Shettle and Weinmann (1970)
			B2 = (3.*mu0*F0*(1.-A)*(2.+3.*mu0+(2.-3.*mu0)*np.exp(-tau/mu0)))/(4.*(4.+3.*(1.-A)*T))
			# calculate atmospheric albedo using their Eq. 22
			atmosA[gi] = 1.-(4.*B2/(3.*mu0*F0))
			gi  += 1 # keep count on asymmetry factors
	
	if (callf == -10): # Figure 4 of Shettle and Weinman (1970)
		# total global downward directed irradiance (Eq. 26 Shettle and Weinman (1970))
		Gdown = Fdown+mu0*np.pi*F0*np.exp(-sum(tau)/mu0)
		
		
		
		return(Fdown, Fup, Gdown)
	
	return(L, atmosA)