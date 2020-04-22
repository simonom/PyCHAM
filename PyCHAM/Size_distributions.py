'''module to create number size distributions'''
##########################################################################################
#                                                                                        #
#    Module to create log-normal size distribution. Code copied/modified from            #
#    http://all-geo.org/volcan01010/2013/09/how-to-use-lognormal-distributions-in-python/#
#                                                                                        #
#                                                                                        #
#    Copyright (C) 2017  David Topping : david.topping@manchester.ac.uk                  #
#                                      : davetopp80@gmail.com                            #
#    Personal website: davetoppingsci.com                                                #
#                                                                                        #
#    This program is free software: you can redistribute it and/or modify                #
#    it under the terms of the GNU Affero General Public License as published            #
#    by the Free Software Foundation, either version 3 of the License, or                #
#    (at your option) any later version.                                                 # 
#                                                                                        #   
#    This program is distributed in the hope that it will be useful,                     #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                      # 
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                       #
#    GNU Affero General Public License for more details.                                 #
#                                                                                        #
#    You should have received a copy of the GNU Affero General Public License            # 
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.               #
#                                                                                        #
#                                                                                        #
##########################################################################################

from scipy import stats # import the scipy.stats module
import numpy as np
import matplotlib.pyplot as plt
import ipdb

def lognormal(num_bins, pconc, std, lowersize, uppersize, loc, scale, space_mode):

	# ---------------------------------
	# Inputs:
	
	# num_bins - number of size bins (not including wall)
	# pconc - starting number concentration of particles (# particle/cc (air)), if this
	# is scalar, the number concentration will be split between size bins, and if it has
	# the same length as number of size bins, each element will be allocated to its
	# corresponding size bin
	# loc - shift of lognormal probability distribution function for seed particles 
	# number-size distribution (um)
	# scale - scaling factor of lognormal probability distribution function for seed 
	# particles (dimensionless) 
	# std - geometric standard deviation (dimensionless)
	# lowersize - lower radius of particles (um)
	# uppersize - upper radius of particles (um)
	# space_mode - string saying whether to space size bins logarithmically or linearly
	# ---------------------------------
    
	# volume in largest and smallest size bin (um3)
	vNb = (4.0/3.0)*np.pi*uppersize**3.0
	v1 = (4.0/3.0)*np.pi*lowersize**3.0
	
	
		
		
	# constant volume ratio method --------------
	# constant volume ratio (13.3 Jacobson 2005)
# 		Vrat = (vNb/v1)**(1/(num_bins-1))
	# volume at centre of each size bin (um3) (13.2 Jacobson 2005)
# 		Varr = v1*(Vrat**np.arange(0, num_bins))
# 		Varr = np.append(Varr, ((4.0*np.pi)/3.0)*(1.0**3.0)) # wall volume (um3)
# 		# radius at centre of each size bin (um)
# 		x_output = ((3.0/(4.0*np.pi))*Varr[0:-1])**(1.0/3.0)
	# ---------------------------------------
	# logarithmic method
	if space_mode == 'log':
		rad_bounds = 10.0**(np.linspace(np.log10(lowersize), 
						np.log10(uppersize), num=(num_bins+1)))
		rwid = (rad_bounds[1::]-rad_bounds[0:-1]) # width of size bins (um)
		x_output = rad_bounds[0:-1]+rwid/2.0 # particle radius (um)
	
	# ---------------------------------------
	# linear method
	if space_mode == 'lin' or space_mode == 'none':
		rad_bounds = np.linspace(lowersize, uppersize, (num_bins+1))
		rwid = (rad_bounds[1]-rad_bounds[0]) # width of size bins (um)
		x_output = rad_bounds[0:-1]+rwid/2.0 # particle radius (um)



	# ---------------------------------------
	# enhance upper radius bound (um) to reduce possibility of particles growing beyond 
	# this (reversed in saving.py)
	upper_bin_rad_amp = 1.0e6
	rad_bounds[-1] = rad_bounds[-1]*upper_bin_rad_amp
	
	if len(pconc)==1 and sum(pconc)>0.0: # for calculating the number size distribution
		# number fraction-size distribution - enforce high resolution to ensure size
		# distribution of seed particles fully captured
		hires = 10**(np.linspace(np.log10(lowersize), np.log10(uppersize), num_bins*1.0e2))
		pdf_output = stats.lognorm.pdf(hires, std, loc, scale)
		pdf_out = np.interp(x_output, hires, pdf_output)	
		# number concentration of all size bins (# particle/cc (air))
		Nperbin = (pdf_out/sum(pdf_out))*pconc
	
	# if number concentration (#/cc (air)) explicitly stated in inputs
	if len(pconc)>1 and sum(pconc)>0.0:
		Nperbin = np.array((pconc))
	if sum(pconc)==0.0:
		Nperbin = np.zeros((num_bins))
	
		# a plot of initial number size distribution (dN/dlog10Dp) 
# 	print(x_output*2.0)
# 	print(Nperbin/(np.log10(rad_bounds[1::]*2.0)-np.log10(rad_bounds[0:-1]*2.0)))
# 	plt.loglog(x_output*2.0, Nperbin/(np.log10(rad_bounds[1::]*2.0)-np.log10(rad_bounds[0:-1]*2.0)), '-x')
# 	plt.show()
# 	ipdb.set_trace()
	
	# volume of single particles per size bin (um3) - use with lognormal method
	Varr = ((4.0*np.pi)/3.0)*(x_output**3.0)
	
	# volume bounds (um3)
	V_bounds = ((4.0*np.pi)/3.0)*(rad_bounds**3.0)


	# return number of particles per size bin (the number size distribution)
	return(Nperbin, x_output, rad_bounds, V_bounds, Varr, upper_bin_rad_amp)