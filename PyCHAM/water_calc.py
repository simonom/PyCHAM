'''module to calculate saturation vapour pressure of water'''

import numpy as np

def water_calc(TEMP, RH, NA):

	# saturation vapour pressure (Pa) (TEMP in K.) (need a reference - have emailed Dave)
	Psat_water = (np.exp((-0.58002206E4/TEMP)+0.13914993E1-(0.48640239E-1*TEMP) 
		+ (0.41764768E-4*(TEMP**2.0E0))-(0.14452093E-7*(TEMP**3.0E0))+ 
		(0.65459673E1*np.log(TEMP))))
	
	# convert saturation vapour pressures from Pa to molecules/cc (air) using ideal
    # gas law, R has units cc.Pa/K.mol
	H2Oc = (RH*Psat_water)*(NA/(8.3144598e6*TEMP))
	
	# convert Psat_water to log10(atm) from Pa to be consistent with vapour pressures 
	# in volat_calc later on
	Psat_water = np.log10(Psat_water/101325.0)
	H2O_mw = 18.0 # state molecular weight of water (g/mol)
	
	return H2Oc, Psat_water, H2O_mw