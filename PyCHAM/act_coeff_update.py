'''module to update particle-phase activity coefficients of components'''
# here the particle-phase activity coefficients of components can be updated, 
# for example due to deliquescence and efflorescence with respect 
# to water

import numpy as np
import scipy.constants as si
import importlib
import hyst_eq

# define function
def ac_up(y, H2Oi, y_H2O0, TEMP, wat_hist, act_coeff):

	# inputs: -------------------------------
	# y - concentrations of components (molecules/cm3)
	# H2Oi - index of water
	# y_H2O0 - gas-phase concentration of water on previous 
	#	integration step (molecules/cm3)
	# TEMP - temperature inside chamber now (K)
	# wat_hist - flag for history of water (0 if dried and therefore 
	#	depends on deliquescence relative humidity, 1if wet 
	#	and therefore dependent on efflorescence relative humidity)
	# act_coeff - original activity coefficients (fraction)
	# -----------------------------------------
	
	
	if (y[H2Oi] != y_H2O0): # check whether relative humidity has changed
		
		# saturation vapour pressure now (Pa) (TEMP in K.)
		# (from water_calc module)
		# R has units cc.Pa/K.mol
		Psat_watern = (np.exp((-0.58002206E4/TEMP)+0.13914993e1-(0.48640239e-1*TEMP) 
		+ (0.41764768e-4*(TEMP**2.e0))-(0.14452093E-7*(TEMP**3.e0))+ 
		(0.65459673e1*np.log(TEMP))))
		
		# relative humidity now
		# (rearranged equation from water_calc module)
		# R has units cc.Pa/K.mol
		RHn = y[H2Oi]*(8.3144598e6*TEMP)/(Psat_watern*si.N_A)
	
		importlib.reload(hyst_eq) # import most recent version
	
		# if history of relative humidity is below the deliquescence RH
		if (wat_hist == 0):
			# deliquescence relative humidity at this temperature
			DRH = hyst_eq.drh(TEMP)
			if (RHn >= DRH):
				act_coeff[:, H2Oi] = 1.
				wat_hist = 1
		# if history of relative humidity is above the  the deliquescence RH
		if (wat_hist == 1):
			# efflorescence relative humidity at this temperature
			ERH = hyst_eq.erh(TEMP)
			if (RHn < ERH):
				act_coeff[:, H2Oi] = 1.e30
				wat_hist = 0

	return(act_coeff, wat_hist)