'''module to estimate component volatilities and liquid densities'''
import numpy as np
import sys
import os
from git import Repo
import shutil
import scipy.constants as si
import errno
import stat

# download latest version of umansysprop
cwd = os.getcwd() # address of current working directory
if os.path.isdir(cwd + '/umansysprop'): # check if there is an existing umansysprop folder
	def handleRemoveReadonly(func, path, exc):
		excvalue = exc[1]
		if not os.access(path, os.W_OK):
			# Is the error an access error ?
			os.chmod(path, stat.S_IWUSR)
			func(path)
		else:
			raise
	shutil.rmtree(cwd + '/umansysprop', ignore_errors=False, onerror=handleRemoveReadonly) # remove existing folder, onerror will change permission of directory if needed. 
sys.path.insert(1, (cwd + '/umansysprop')) # address for updated version
git_url = 'https://github.com/loftytopping/UManSysProp_public.git'
Repo.clone_from(git_url, (cwd + '/umansysprop'))

from umansysprop import boiling_points
from umansysprop import vapour_pressures
from umansysprop import liquid_densities

def volat_calc(spec_list, Pybel_objects, TEMP, H2Oi, num_speci, Psat_water, voli, volP, 
				testf, corei):

	# ------------------------------------------------------------
	# inputs:
	# spec_list - array of SMILE strings for components 
	# (omitting water and core, if present)
	# Pybel_objects - list of Pybel objects representing the species in spec_list
	# (omitting water and core, if present)
	# testf - flag for whether in normal mode (0) or testing mode (1)
	# corei - index of seed particle component
	# ------------------------------------------------------------
	
	if testf==1:
		return(0,0,0) # return dummies

	# if voli is in relative index (-n), change to absolute
	for i in range(len(voli)):
		if voli[i]<0:
			voli[i] = num_speci+voli[i]

	NA = si.Avogadro # Avogadro's number (molecules/mol)
	y_dens = np.zeros((num_speci, 1)) # components' liquid density (kg/m3)
	Psat = np.zeros((num_speci, 1)) # species' vapour pressure
	
	for i in range (num_speci):
		
		# omit estimation for water as it's value is already given in Psat_water 
		# (log10(atm))
		if i == H2Oi:
			Psat[i] = Psat_water
			y_dens[i] = 1.0*1.0E3 # (kg/m3 (particle))
			continue 
		if i == corei: # core properties (only used if seed particles present)
			y_dens[i] = 1.77*1.0E3 # core density (kg/m3 (particle))
			continue 
		if spec_list[i] == '[HH]': # omit H2 as unliked by liquid density code
			# liquid density code does not like H2, so manually input kg/m3
			y_dens[i] = 1.0e3
		else:
			# density (convert from g/cc to kg/m3)
			y_dens[i] = liquid_densities.girolami(Pybel_objects[i])*1.0E3
		
		# vapour pressure (log10 atm) (# eqn. 6 of Nannoolal et al. (2008), with dB of 
		# that equation given by eq. 7)
		Psat[i] = ((vapour_pressures.nannoolal(Pybel_objects[i], TEMP, 
						boiling_points.nannoolal(Pybel_objects[i]))))

	Psat = (np.power(10.0, Psat)*101325.0) # convert to Pa
	# manually assigned vapour pressures (Pa) (including seed component (if applicable))
	if len(voli)>0:
		Psat[voli, 0] = volP
	Psat_Pa = np.zeros((len(Psat), 1)) # for storing vapour pressures in Pa (Pa)
	Psat_Pa[:, 0] = Psat[:, 0]
    	
    # convert saturation vapour pressures from Pa to molecules/cc (air) using ideal
    # gas law, R has units cc.Pa/K.mol
	Psat = Psat*(NA/(8.3144598e6*TEMP))

	return Psat, y_dens, Psat_Pa