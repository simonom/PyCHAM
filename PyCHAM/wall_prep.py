'''module to provide PyCHAM with required wall properties'''

import numpy as np
import scipy.constants as si

def wall_prep(TEMP, Psat_Pa, y_mw, num_speci, wall_accom, testf):

	# ------------------------------------------------------------
	# inputs: 
	
	# TEMP - temperature (K) (1)
	# Psat_Pa - component vapour pressures in Pa (num_speci, 1)
	# y_mw - molecular weight (g/mol) (num_speci, 1)
	# num_speci - number of components
	# wall_accom - wall accommodation coefficient for all components (1)
	# testf - flag for whether in normal mode (0) or testing mode (1)
	# ------------------------------------------------------------
	
	if testf==1:
		return(0,0) # return dummies
	# accommodation coefficient of components with chamber wall (dimensionless)
	# (eq. 14 of Zhang et al. 2015) (value given in table 2 of Zhang)
	wall_accom = np.ones((num_speci, 1))*wall_accom
	
	# activity coefficient of components on wall
	act_coeff_wall = 1.0
	
	# vapour-wall partitioning coefficient (m3/g) (Zhang et al. 2015 eq. 12)
	Kw = np.zeros((len(Psat_Pa), 1))
	
	Kw[:, 0] = (8.314e0*TEMP)/(Psat_Pa[:, 0]*act_coeff_wall*y_mw[:, 0])
	
	
	return(wall_accom, Kw)