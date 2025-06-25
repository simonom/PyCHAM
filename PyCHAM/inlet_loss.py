########################################################################
#								       #
# Copyright (C) 2018-2025					       #
# Simon O'Meara : simon.omeara@manchester.ac.uk			       #
#								       #
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
'''module for inlet loss correction'''
# when particles pass through an inlet system this module accounts for size-dependent loss

import numpy as np
import datetime
import importlib
import matplotlib.pyplot as plt

def inlet_loss(call, Nwet, xn, yp, losst, num_comp, self):

	# inputs: -----------------------------------------------
	# call - the calling module flag
	# Nwet - number concentration of particles entering the inlet (# particles/cm3)
	# xn - radii of particles, corrected for relative humidity in inlet (um)
	# yp - concentrations of components in the particle phase, corrected for 
	#	relative humidity in inlet (molecules/cm3)
	# losst - time of passage through inlet (s)
	# num_comp - number of components
	# self - PyCHAM series of variables
	# ---------------------------------------------------------

	# calculate fractional loss rate of particles (fraction /s)
	# diameter array to solve over (um), note this calculation 
	# follows that in wallloss.py
	Dp = xn*2.
	Beta = np.zeros((Dp.shape[0], Dp.shape[1]))
	Beta[Dp<self.inflectDp] = (10**((np.log10(self.inflectDp)-np.log10(
		Dp[Dp<self.inflectDp]))*self.pwl_xpre+np.log10(self.inflectk)))
	Beta[Dp>=self.inflectDp] = (10**((np.log10(Dp[Dp>=self.inflectDp])-np.log10(
		self.inflectDp))*self.pwl_xpro+np.log10(self.inflectk)))
	
	if (call >= 3): # if button pressed to plot loss rate as a function of diameter

		if isinstance(Beta, str):
			if (Beta[0:5] == 'Error'):
				return(Beta)		

		plt.ion()
		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
		# plot fractional loss rate of particles (/s)
		if (any(Beta[0, :]==0)):
			ax0.semilogx(Dp[0, :], Beta[0, :], '-k')
		else:
			ax0.loglog(Dp[0, :], Beta[0, :], '-k')

		# plot details
		if (call == 3):
			ax0.set_title(str('Loss rate of particles during passage ' +
				'through instrument inlet'))
		if (call == 4):
			ax0.set_title(str('Loss rate of particles to walls'))
		ax0.set_ylabel('Loss rate (fraction of particles (0-1) s$^{-1}$)', size = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
		ax0.set_xlabel('Particle diameter ($\mu\,$m)', fontsize=14)
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
		plt.show()
		return([])

	# integrate over entire time in inlet (fraction)
	sd_lrate = Beta*losst
	
	Nwet -= Nwet*sd_lrate # correct for inlet loss (# particles/cm3)
	# repeat over components and correct for loss of components (molecules/cm3)
	yp -= yp*np.repeat(sd_lrate, num_comp, axis=1)

	return(Nwet, yp)

# function for correcting for wall loss of particles and converting particle
# volume to mass for estimation of SOA mass concentration with time 
def wl_correct(self):

	# inputs: ----
	# self - contains PyCHAM variables
	# ------------

	# get required variables from self

	Ndry = np.zeros((self.ro_obj.Nrec_dry.shape[0], self.ro_obj.Nrec_dry.shape[1]))
	Ndry[:, :] = self.ro_obj.Nrec_dry[:, :]
	timehr = self.ro_obj.thr
	x_um = self.ro_obj.cen_size

	# first, calculate fractional loss rate of particles (fraction /s)
	# diameter array to solve over (um), note this calculation 
	# follows that in wallloss.py
	Dp = x_um*2.
	Beta = np.zeros((Dp.shape[0], Dp.shape[1]))
	Beta[Dp<self.inflectDp] = (10**((np.log10(self.inflectDp)-np.log10(
		Dp[Dp<self.inflectDp]))*self.pwl_xpre+np.log10(self.inflectk)))
	Beta[Dp>=self.inflectDp] = (10**((np.log10(Dp[Dp>=self.inflectDp])-np.log10(
		self.inflectDp))*self.pwl_xpro+np.log10(self.inflectk)))
	
	# prepare to hold particle wall loss corrected number concentration (particles/cm^3)
	Ndry_corr = np.zeros((Ndry.shape))
	
	# flag for seeing injection of seed particles
	seed_view_flag = -1

	# prepare to hold cumulative (over time) number concentration of particles 
	# lost per sze bin (# particles/cm^3)
	n_loss = np.zeros((1, Ndry.shape[1]))

	# now loop through times of results to correct particle loss to wall at
	# each time
	for it in range(1, len(timehr)):

		# skip this time if no seed particles yet injected
		if (sum(Ndry[it, :]) == 0.):
			continue
		
		if (seed_view_flag == -1):	
			# get particle number concentration at seed injection 
			# (prior to any particle loss to
			# wall) (# particles/cm^3)
			Ndry_corr[it, :] = Ndry[it, :]
			# remember that we have registered seed particle injection
			seed_view_flag = it

		# time step (s)
		tn = (timehr[it]-timehr[it-1])*3600.

		# get the arithmetics mean particle diameter (um) per size bin at this 
		# time and the preceding time step
		Dpn = ((x_um[it, :]+x_um[it-1, :])/2.)*2.
		
		# get fractional loss rates (fraction of particles /s)
		Betan = Beta[it, :]
		# multiply by time step to get fraction of particles lost over time step (0-1)
		# per size bin
		frac_loss = Betan*tn
		# multiply by arithmetic mean of particle number concentration per size bin
		# to get the number concentration of particles lost (particles/cm^3)
		n_loss += ((Ndry[it, :]+Ndry[it-1, :])/2.)*frac_loss
		# add this number to the number concentration to correct for wall loss of
		# particles (# particles/cm^3)
		Ndry_corr[it, :] = Ndry[it, :]+n_loss

	# estimate volume concentration of seed particles (um^3/cm^3), summing over size bins
	seedv = np.sum((
		Ndry_corr[seed_view_flag, :]*((4./3.)*np.pi*x_um[seed_view_flag, :]**3.)))

	# estimate mass concentration of seed particles (ug/m^3), note conversion of seedv to
	# cm^3/cm^3 and conversion of result from g/cm^3 to ug/m^3. Note that we use the 
	# density for SOA rather than the density of the seed material, this allows for
	# subtraction of the seed material later on when we don't know the fraction of
	# particle number concentration that is seed vs. that is SOA 
	seedm = ((seedv*1.e-12)*self.soa_rho)*1.e12

	# estimate volume concentration of particles at all times (um^3/cm^3), summing over size
	# bins
	pv = np.sum((Ndry_corr[seed_view_flag::, :]*((4./3.)*
		np.pi*x_um[seed_view_flag::, :]**3.)), axis= 1)

	# estimate mass concentration of all particles (ug/m^3), note conversion of seedv to
	# cm^3/cm^3 and conversion of result from g/cm^3 to ug/m^3.
	pm = ((pv*1.e-12)*self.soa_rho)*1.e12

	# subtract mass of seed (ug/m^3)
	soam = pm-seedm
	# plot SOA mass concentration against time (ug/m^3)
	plt.ion()
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
	ax0.plot(timehr[seed_view_flag::], soam)

	ax0.set_ylabel(str('Particle wall loss corrected SOA mass concentration ' +
		'($\mathrm{\mu g\,m^{-3}}$)'), size = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	ax0.set_xlabel('Time through simulation (hours)', fontsize=14)
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')

	plt.show()

	return()