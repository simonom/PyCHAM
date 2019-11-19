'''module to create log-normal size distribution'''
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

# Last modification 28/5/17

from scipy import stats # Import the scipy.stats module
import numpy as np
import matplotlib.pyplot as plt

def lognormal(num_bins, total_conc, std, lowersize, uppersize, loc, scale):

	# ---------------------------------
	# Inputs:
	
	# num_bins - number of size bins (not including wall)
	# loc - shift of lognormal probability distribution function for seed particles 
	# number-size distribution (um)
	# scale - scaling factor of lognormal probability distribution function for seed 
	# particles (dimensionless) 
	# std - geometric standard deviation (dimensionless)
	# lowersize - lower radius of particles (um)
	# uppersize - upper radius of particles (um)
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
	# lognormal method
# 		x_output = np.exp(np.linspace(np.log(lowersize), 
# 						np.log(uppersize), num=(num_bins)))
	# radius bounds on size bins (um)
# 		rad_bounds = np.exp((np.log(x_output[1:-1])+np.log(x_output[0:-2]))/2.0)
	
	# ---------------------------------------
	# linear method
	rad_bounds = np.linspace(lowersize, uppersize, (num_bins+1))
	rwid = (rad_bounds[1]-rad_bounds[0]) # width of size bins (um)
	x_output = rad_bounds[0:-1]+rwid/2.0 # particle radius (um)

	# extreme radius bound (um) to reduce possibility of particles growing beyond this
	# (reversed in res_plot.py)
	rad_bounds[-1] = rad_bounds[-1]*1.0e2
	# ---------------------------------------
	
	if total_conc>0.0:
		# number fraction-size distribution - enforce high resolution to ensure size
		# distribution of seed particles fully captured
		hires = np.linspace(lowersize, uppersize, 1000)
		pdf_output = stats.lognorm.pdf(hires, std, loc, scale)
		pdf_out = np.zeros((num_bins))
		# now loop through actual bins and sum probabilities bounded by them
		for i in range(num_bins):
			hires_ind = (hires>=rad_bounds[i])*(hires<rad_bounds[i+1])
			pdf_out[i] = pdf_output[hires_ind].sum()
		# number concentration of all size bins (# particle/cc (air))
		Nperbin = (pdf_out/sum(pdf_out))*total_conc
	else:
		Nperbin = np.zeros((num_bins))
	
	# volume of particles per size bin (um3) - use with lognormal method
	Varr = ((4.0*np.pi)/3.0)*(x_output**3.0)
	
	# volume bounds (um3)
	V_bounds = ((4.0*np.pi)/3.0)*(rad_bounds**3.0)
	V_bounds[-1] = V_bounds[-1]*1.0e1 # set upper bound very high


	# return number of particles per size bin (the number size distribution)
	return(Nperbin, x_output, rad_bounds, rwid, V_bounds, Varr)