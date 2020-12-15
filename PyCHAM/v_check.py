'''module to check whether ode solver estimates partitioning that leads to unacceptable change in volume of particles'''
# code called by the mov_cen_main.py module to ensure that the volume change of particles
# estimated by the ode solver is within acceptable bounds

import numpy as np
from compl_evap import compl_evap as compl_evap

def Vchange_check(res, MV, sbb, sbn, NA, n0, nc, solv_time, ts0, bc_red, Vol0, Psat):

	# inputs: ---------------------------------
	
	# res - concentrations (molecules/cc (air)) of all components in all phases (columns) 
	#		at the adaptive time steps set by the ode solver
	# MV - molar volume of components (cc/mol)
	# sbb - fixed volume bounds (um3) of size bins
	# sbn - number of size bins
	# NA - Avogadro's number (molecules/mol)
	# n0 - number concentration of particles per size bin (# particle/cc (air)) before 
	#		integration
	# nc - number of components
	# solv_time - adaptive time step (s) used by ode solver
	# ts0 - original time step passed to integrator (s)
	# bc_red - flag for whether boundary conditions have caused a change to time step
	# Vol0 - default volume per size bin (um3), volume at centre of size bin bounds
	# Psat - saturation vapour pressure of components (molecules/cm3 (air))
	# -----------------------------------------
	
	# flag for breaking condition of volume change, updated once volumes checked below
	Vchang_flag = 2
	# count on adaptive time steps - work backwards through the adaptive time steps of the
	# ode solver to check whether volume condition met
	tsi = 1
	# acceptable number of size bins to change by
	acc_sb_chng = 1000
	Vnew = np.zeros((sbn))
	
	while Vchang_flag>=1:
		
		# estimated concentrations (moleculecs/cc (air)) at this time step
		ytest = res
		
		# loop through size bins to check whether volume condition met
		for sbi in range(sbn):
			
			# calculating new volumes of particles in size bins --------------------------
			if n0[sbi]<1.0e-10: # if no particles, assign the default volume (um3)
				Vnew[sbi] = Vol0[sbi]
				continue # onto next size bin
			# concentration of components in this size bin (molecules/cc (air))
			Cnow = ytest[(sbi+1)*nc:(sbi+2)*nc]
			# molar concentration of components in one particle (mol/cc (air))
			Cnow = Cnow/(NA*n0[sbi])
			# new volume (um3) of single particles, note MV has units cc/mol, so needs 
			# conversion to um3/mol
			Vnew[sbi] = sum(Cnow*(MV[:, 0]*1.e12))
			
			# comparing new volumes against condition for volume change ------------------
			# check for volume conditions 
			if (Vnew[sbi]<0.0): # shrunk to unrealistic negative volume
				# allow complete evaporation of smallest size bin if only volatiles 
				# present and time step already relatively small
				if sbi == 0.0 and sum(Cnow[Psat<1.0e-20]<1.0e-20) and ts0<1.0e-3:
					(ytest, n0, Vnew) = compl_evap(ytest, n0, Vnew, Vol0, nc, sbn)
				else: # change integration time step
					Vchang_flag = 1
# 					print('negative volume seen in Vchange_check.py, size bin affected: ' + str(sbi) + ', concentration of components: ' + str(Cnow))
			if (Vnew[sbi]>sbb[-1]): # grown to unpractical positive volume
				Vchang_flag = 1
# 				print('volume above uppermost size bin bound seen in Vchange_check.py, size bin affected: ' + str(sbi) + ', concentration of components: ' + str(Cnow))
			if ((sbi-acc_sb_chng) >= 1):
				if (Vnew[sbi]<sbb[sbi-acc_sb_chng]): # excessive shrink
					Vchang_flag = 1
					print('volume below acceptable change in Vchange_check.py, size bin affected: ' + str(sbi) + ', concentration of components: ' + str(Cnow) + ', percentage away from allowed change bound: ' + str(((sbb[sbi-acc_sb_chng]-Vnew[sbi])/sbb[sbi-acc_sb_chng])*100.0))
			if (((sbi+1)+acc_sb_chng) < sbn+1):
				if (Vnew[sbi]>sbb[(sbi+1)+acc_sb_chng]): # excessive growth
					Vchang_flag = 1
					print('volume above acceptable change in Vchange_check.py, size bin affected: ' + str(sbi) + ', concentration of components: ' + str(Cnow) + ', percentage away from allowed change bound: ' + str(((Vnew[sbi]-sbb[(sbi+1)+acc_sb_chng])/sbb[(sbi+1)+acc_sb_chng])*100.0))
					
			# decision based on volume change condition ----------------------------------
			if (Vchang_flag==1): # if volume condition not met
				tsi += 1
				# break out of size bin loop and change adaptive time step index
				if (tsi < res.shape[0]):
					Vchang_flag = 3 # update flag and try next integrator time step
# 					print('trying earlier ode solver result, number of results from the final: ', tsi)
					break # stop size bin loop
				else: # need to reduce integration time step (s) below the minimum here
					redt = 1 # flag for time step reduction due to volume change
					t = solv_time[0]/2.0
					bc_red = 0 # flag for time step reduction due to boundary conditions
					Vchang_flag = 0 # will exit while loop for volume condition not met
					break # stop size bin loop
			else: # if volume condition met
				Vchang_flag = 2
		
		# applies if size bin loop finished and volume condition met ---------------------		
		# if volume condition met after all size bins checked adopt the results at this 
		# time
		if (Vchang_flag==2):
			
			if (tsi==1): # integration time step unchanged
				t = ts0 # time for integration on next step (s)
				redt = 0 # flag for no time step reduction due to volume change
			if (tsi>1): # if a sub-step used, then get associated time
				t = solv_time[-tsi] # new time for integration on next step (s)
				redt = 2 # flag for time step reduction due to volume change
			Vchang_flag = 0 # will exit while loop for volume condition not met
			
	return(redt, t, bc_red, Vnew, tsi)
