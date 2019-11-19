'''function to estimate the initial concentration of water on particles and walls'''

import numpy as np
from kimt_calc import kimt_calc
import ipdb
from mov_cen_main import mov_cen_main as movcen # moving centre method for rebinning

def init_water_partit(x, y, H2Oi, Psat, mfp, num_sb, num_speci, 
						accom_coeff, y_mw, surfT, R_gas, TEMP, NA, y_dens, 
						N_perbin, DStar_org, RH, core_diss, Varr, Vbou, Vol0, tmax, MV,
						therm_sp, Cw, total_pconc, kgwt):
						
	# --------------------------------------------------------------
	# inputs:
	
	# Psat - saturation vapour pressure of components (molecules/cc (air))
	# therm_sp - thermal speed of components (m/s) (num_speci)
	# Cw - concentration of wall (g/m3 (air))
	# total_pconc - total initial particle concentration (#/cc (air))
	# kgwt - mass transfer coefficient for vapour-wall partitioning (cm3/molecule.s)
	# --------------------------------------------------------------
	if total_pconc>0.0: # if seed particles present
		sbstep = 0 # count on size bins
		
		for radius in x:
			
			# get core properties
			# concentration (molecules/cc (air))
			ycore = y[num_speci*(sbstep+2)-1]*core_diss
			
			# first guess of particle-phase water concentration based on RH = mole 
			# fraction (molecules/cc (air))
			y[num_speci*(sbstep+1)+H2Oi] = ycore*RH
			
			# gas phase concentration of water (molecules/cc (air)), will stay constant
			# because RH is constant
			Cgit = y[H2Oi]
			
			# concentration of water at particle surface in gas phase (molecules/cc (air))
			Csit = y[num_speci*(sbstep+1)+H2Oi]
			
			y0 = y # initial concentration (molecules/cc (air))
			
			while np.abs(Cgit-Csit)>1.0e8:
				
				# call on the moving centre method for redistributing particles of 
				# different sizes
				(N_perbin, Varr, y[num_speci::], radius, redt, blank, tnew) = movcen(N_perbin, 
					Vbou, 
					np.transpose(y[num_speci::].reshape(num_sb, num_speci)), 
					(np.squeeze(y_dens*1.0e-3)), num_sb, num_speci, y_mw, x, Vol0, 0.0,
					tmax, 0, y0[num_speci::], MV)
					
				y0 = y # update initial concentrations (molecules/cc (air))
				
				# update partitioning coefficients
				[kimt, kelv_fac] = kimt_calc(y, mfp, num_sb, num_speci, accom_coeff, y_mw,   
								surfT, R_gas, TEMP, NA, y_dens, N_perbin, DStar_org, 
								x.reshape(1, -1)*1.0e-6, Psat, therm_sp, 
								H2Oi)
				
				# concentration of water at particle surface in gas phase (molecules/cc (air))
				Csit = y[num_speci*(sbstep+1)+H2Oi]
				
				# total particle surface gas-phase concentration of all species 
				# (molecules/cc (air)):
				conc_sum = (np.sum(y[num_speci*(sbstep+1):num_speci*(sbstep+2)-1])+
									y[num_speci*(sbstep+2)-1]*core_diss)
				
				Csit = (Csit/conc_sum)*Psat[H2Oi]*kelv_fac[sbstep]
				
				
				dydt = kimt[H2Oi, sbstep]*(Cgit-Csit) # partitioning rate (molecules/cc.s)
				# new estimate of concentration of condensed water (molecules/cc (air))
				y[num_speci*(sbstep+1)+H2Oi] += dydt*(1.0e-4*radius[sbstep])
				
			sbstep += 1
	else:
		if num_sb>1: 
			radius = np.zeros((len(x)))
			radius[:] = x[:]
		else:
			radius = 0.0
	
	# first guess of wall-phase water concentration based on RH = mole 
	# fraction (molecules/cc (air))
	y[num_speci*num_sb+H2Oi] = Cw*RH
	
	# gas phase concentration of water (molecules/cc (air)) responds to loss to walls
	y[H2Oi] = y[H2Oi]-y[num_speci*(num_sb)+H2Oi]
	
	scaling = 10**(19.0-np.log10(Cw)) # scaling for new estimate
	
	# allow gas-phase water to equilibrate with walls
	while np.abs(y[H2Oi]-Psat[H2Oi,0]*((y[num_speci*(num_sb)+H2Oi]/Cw)))>1.0e8:
		
												
		# concentration of water at wall surface in gas phase (molecules/cc (air))
		Csit = Psat[H2Oi,0]*((y[num_speci*(num_sb)+H2Oi]/Cw))
		
		dydt = (kgwt*Cw)*(y[H2Oi]-Csit) # partitioning rate (molecules/cc.s)
		
		# new estimate of concentration of condensed water (molecules/cc (air))
		y[num_speci*(num_sb)+H2Oi] += dydt*(scaling)
		y[H2Oi] -= dydt*(scaling)
		
	print('finished initiating water condensation')
	return(y, Varr, radius)