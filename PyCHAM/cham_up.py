'''update the chamber variables depending on time'''
# based on the user inputs and time through simulation,
# chamber variables are updated before calling the 
# ODE solver

import numpy as np
import os
import pp_dursim
import volat_calc

# define function
def cham_up(sumt, temp, tempt, Pnow, light_stat, light_time, 
	light_time_cnt, light_ad, tnew, nuc_ad, nucv1, nucv2, nucv3, 
	new_part_sum1, update_stp, update_count, lat, lon, dayOfYear,
	photo_par_file, act_flux_path, injectt, gasinj_cnt, inj_indx, 
	Ct, pmode, pconc, pconct, seedt_cnt, num_comp, y, N_perbin, 
	mean_rad, corei, seedVr, seed_name, lowsize, uppsize, num_sb, MV, rad0, radn, std, 
	y_dens, H2Oi, rbou, const_infl_t, infx_cnt, Cinfl, wall_on, Cfactor, seedi, diff_vol, DStar_org):

	# inputs: ------------------------------------------------
	# sumt - cumulative time through simulation (s)
	# temp - temperature in chamber (K)
	# tempt - times that temperatures reached (s)
	# Pnow - pressure in chamber (Pa)
	# light_stat - status of lights
	# light_time - times that light attain status (s)
	# light_time_cnt - light status counter
	# light_ad - marker for whether to change time interval 
	#	in response to changing natural light intensity
	# tnew - time interval between chamber updates (s)
	# nuc_ad - flag for whether user wants time step adapted 
	# to nucleation
	# nucv1 - nucleation parameter one
	# nucv2 - nucleation parameter two
	# nucv3 - nucleation parameter three
	# new_part_sum1 - total number concentration of new 
	#	particles so far (#/cc (air))
	# update_stp - time interval between operator-split 
	#	updates (s)
	# update_count - count since operator-split last 
	#	updated (s)
	# lat - latitude (degrees)
	# lon - longitude (degrees)
	# dayOfYear - number of days since 31st December
	# photo_par_file - photochemistry parameter file
	# act_flux_path - actinic flux file
	# injectt - time of instantaneous injections of 
	#	components (s)
	# gasinj_cnt - count on injection times of component(s)
	# inj_indx - index of component(s) being injected after 
	#	experiment start
	# Ct - concentration(s) (ppb) of component(s) injected 
	#	instantaneously after experiment start
	# pmode - whether particle number size distributions stated explicitly or by mode
	# pconc - concentration of injected 
	#	particles (#/cc (air))
	# pconct - times of particle injection (s)
	# seedt_cnt - count on injections of particles
	# num_comp - number of components
	# y - concentration of components (molecules/cc (air))
	# N_perbin - concentration of particles (#/cc (air))
	# mean_rad - mean radius for particle number size 
	#	distribution (um)
	# corei - index of core component
	# seedVr - volume ratio of component(s) comprising seed particles
	# seed_name - name(s) of component(s) comprising seed 
	#	particles
	# lowsize - lower size bin boundary (um)
	# uppsize - upper size bin boundary (um)
	# num_sb - number of size bins
	# MV - molar volume of components (cc/mol)
	# rad0 - initial radius at size bin centres (um)
	# radn - current radius at size bin centres (um)
	# std - standard deviation for injected particle number size 
	#	distributions
	# y_dens - component densities (g/cm3)
	# H2Oi - index of water
	# rbou - size bin radius bounds (um)
	# const_infl_t - times for constant influxes (s)
	# infx_cnt - count on constant influx occurrences
	# Cinfl - influx rate for components with constant 
	# 	influx (ppb/s)
	# wall_on - marker for whether wall is on
	# Cfactor - conversion factor from ppb to molecules/cc (air)
	# seedi - index of seed component(s)
	# diff_vol - diffusion volumes of components according to 
	#	Fuller et al. (1969)
	# DStar_org - gas-phase diffusion coefficients of components (cm2/s)
	# -----------------------------------------------------------------------

	# check on change of light setting --------------------------------------

	# begin by assuming no change to time interval required due to chamber 
	# condition/nucleation
	bc_red = 0
	
	if (len(light_time))>0:
	
		# whether lights on (1) or off (0) during this step
		lightm = light_stat[int(sum(light_time<=sumt)-1)]
		
		# check whether changes occur at start of this time step
		if (sumt == light_time[light_time_cnt] and light_time_cnt>-1):
			
			if (light_time_cnt<(len(light_stat)-1)):
				light_time_cnt += 1 # keep count of light setting index
			else:
				light_time_cnt = -1 # reached end
			# reset flag for time step reduction due to chamber condition
			bc_red = 0
				
		
		# check whether light on/off changes during proposed integration time step
		if (sumt+tnew > light_time[light_time_cnt] and light_time_cnt != -1):
			# if yes, then reset integration time step so that next step coincides 
			# with change
			tnew = light_time[light_time_cnt]-sumt
			bc_red = 1 # flag for time step reduction due to boundary conditions
			
		# if reached final status of lights, then keep this status
		if light_time_cnt == -1:
			lightm = light_stat[light_time_cnt]
 
	# if lights are on during this step and lighting is natural, then check whether
	# proposed time step needs reducing to limit change to light intensity, if this
	# time interval adaption is requested
	# if using natural light
	cwd = os.getcwd() # address of current working directory

	if (lightm == 1 and photo_par_file == str(cwd + '/PyCHAM/photofiles/MCMv3.2') and act_flux_path == 'no' and light_ad == 1):	
		# check time step required to limit change to rate of 
		# MCM photochemical equation number 6 (python index 5), 
		# which the unit test for
		# the zenith module shows to be most photosensitive (compared to
		# other photochemical equations)
		import zenith
		# photochemical rate now
		(secxn, cosxn) = zenith.zenith(sumt, lat, lon, dayOfYear)
		Jn =1.747E-01*cosxn**(0.155)*np.exp(-1.0*0.125*secxn)
		# photochemical rate after proposed time step
		(secxt, cosxt) = zenith.zenith(sumt+tnew, lat, lon, dayOfYear)
		Jt =1.747E-01*cosxn**(0.155)*np.exp(-1.0*0.125*secxn)
		# iteratively reduce proposed time interval until photochemical
		# rate changes by acceptable amount
		while (abs(Jt-Jn) > 5.e-3):
			tnew = tnew*0.9
			# photochemical rate after proposed time step
			(secxt, cosxt) = zenith.zenith(sumt+tnew, lat, lon, dayOfYear)
			Jt =1.747E-01*cosxt**(0.155)*np.exp(-1.0*0.125*secxt)
			bc_red = 1

	# check on updates to temperature (K) --------------------------------------	
	if len(temp)>1: # because a temperature must be given for experiment start
	
		# check whether changes occur at start of this time step
		if (sumt == tempt[temp_count]):

			# new temperature (K)
			temp_now = temp[temp_count]
			
			# update vapour pressure of water (log10(atm)), but don't change 
			# gas-phase concentration because we assume RH allowed to change with
			# varying temperature
			[_, Psat_water, _] = water_calc(temp_now, RH, 6.02214129e+23)
			
			# update vapour pressures of all components (molecules/cc and Pa), 
			# ignore density output
			[Psat, _, Psat_Pa] = volat_calc.volat_calc(0, Pybel_objects, temp_now, H2Oi,   
							num_comp, Psat_water, [], [], 0, corei, seed_name, 
							pconc, 0, 0.0, [], 1)

			# note, assume that air pressure inside chamber stays constant despite
			# varying temperature, therefore total molecular concentration must vary
			# (molecules/cc (air))
			ntot = Pnow*(NA/(si.R*temp_now))
			# remember number of molecules in one billionth of total number prior
			# to update
			Cfactor0 = Cfactor
			# update number of molecules in one billionth of this
			Cfactor = ntot*1.e-9 # ppb-to-molecules/cc
			
			# dynamic viscosity of air (kg/m.s), eq. 4.54 of Jacobson 2005
			dyn_visc = 1.8325e-5*((416.16/(temp_now+120.))*(temp_now/296.16)**1.5)
	
			ma = 28.966e-3 # molecular weight of air (kg/mol) (Eq. 16.17 Jacobson 2005)
			
			# air density (kg/m3 (air)), ideal gas law
			rho_a =  (Pnow*ma)/((si.R)*temp_now)
						
			# update mean free path and thermal speed
			# mean thermal speed of each molecule (m/s) (11.151 Jacobson 2005)
			# note that we need the weight of one molecule, which is why y_mw is divided by
			# Avogadro's constant, and we need it in kg, which is why we multiply by 1e-3
			therm_sp = ((8.*si.k*temp_now)/(np.pi*(y_mw/si.N_A)*1.e-3))**0.5
			
			# mean free path (m) for each component (15.24 of Jacobson 2005)
			# molecular weight of air (28.966 g/mol taken from table 16.1 Jacobson 2005)
			mfp = (2.*dyn_visc/(rho_a*therm_sp)).reshape(-1, 1)
			
			# diffusion coefficient (m2/s) of components in gas phase (air), eq 4.1.4 of
			# the Taylor (1993) textbook 
			# Multicomponent Mass Transfer, ISBN: 0-471-57417-1, note diffusion 
			# volume for air (19.7) taken from Table 4.1 of Taylor (1993) and mw of 
			# air converted to g/mol from kg/mol.  This is a replication of the original method 			# from Fuller et al. (1969): doi.org/10.1021/j100845a020
			DStar_org = 1.013e-2*temp_now**1.75*(((y_mw+ma*1.e3)/(y_mw*ma*1.e3))**0.5)/(Pnow*(diff_vol**(1./3.)+19.7**(1./3.))**2.)
			# convert to cm2/s
			DStar_org = DStar_org*1.e4
			
			
			# alter constant concentration (molecules/cc) of any components
			# with constant gas-phase concentration (ppb)
			if num_const_compi>0:
				
				y[const_compi[:]] = y[const_compi[:]]*(Cfactor/Cfactor0)

			if (temp_count<(len(tempt)-1)):
				temp_count += 1 # keep count of temperature setting index
			else:
				temp_count = -1 # reached end
			bc_red = 0 # reset flag for time step reduction due to boundary conditions
			
		# check whether temperature changes during proposed integration time step
		if (sumt+tnew > tempt[temp_count] and temp_count!=-1):
			# if yes, then reset integration time step so that next step coincides 
			# with change
			tnew = tempt[temp_count]-sumt
			bc_red = 1 # flag for time step reduction due to boundary conditions

	if (len(temp) == 1):
		temp_now = temp[0] # temperature constant if only one value given

	# check on instantaneous injection of components ---------------------------------------
	if len(injectt)>0 and gasinj_cnt>-1: # if any injections occur
	
		# check whether changes occur at start of this time step
		if (sumt == injectt[gasinj_cnt]):
			# account for change in gas-phase concentration,
			# convert from ppb/s to molecules/cc.s (air)
			y[inj_indx] += Ct[:, gasinj_cnt]*Cfactor-y[inj_indx]
			if (gasinj_cnt<(Ct.shape[1]-1)):
				gasinj_cnt += 1 # update count on injections
			else:
				gasinj_cnt = -1 # reached end
			bc_red = 0 # reset flag for time step reduction due to boundary conditions
				
		# check whether changes occur during proposed integration time step
		if (sumt+tnew > injectt[gasinj_cnt] and gasinj_cnt != -1):
			# if yes, then reset integration time step so that next step coincides 
			# with change
			tnew = injectt[gasinj_cnt]-sumt
			bc_red = 1 # flag for time step reduction due to boundary conditions
	
	# check on instantaneous injection of particles --------------------------------------	
	if ((sum(pconct[0, :]) > 0) and (seedt_cnt > -1) and (num_sb-wall_on > 0)): # if constant influx occurs
	
		# check whether changes occur at start of this time step
		if (sumt == pconct[0, seedt_cnt]):
			
			# account for change in seed particles
			[y[num_comp:num_comp*(num_sb-wall_on+1)], N_perbin, x, 
					Varr] = pp_dursim.pp_dursim(y[num_comp:num_comp*(num_sb-wall_on+1)], 
					N_perbin, 
					mean_rad[:, seedt_cnt], pmode, pconc[:, seedt_cnt], seedi, seedVr, lowsize, 
					uppsize, num_comp, (num_sb-wall_on), MV, rad0, radn, 
					std[:, seedt_cnt], y_dens, H2Oi, rbou)
			if (seedt_cnt<(pconct.shape[1]-1)):
				seedt_cnt += 1
			else:
				seedt_cnt = -1 # reached end
			bc_red = 0 # reset flag for time step reduction due to boundary conditions
			
		# check whether changes occur during proposed integration time step
		if (sumt+tnew > pconct[0, seedt_cnt] and seedt_cnt!=-1): 
			# if yes, then reset integration time step so that next step coincides 
			# with change
			tnew = pconct[0, seedt_cnt]-sumt
			bc_red = 1 # flag for time step reduction due to boundary conditions

	# check on continuous influxes of components ----------------------------------------------
	if (len(const_infl_t) > 0): # if influx occurs
	
		# in case influxes begin after simulation start create a zero array of correct shape
		if (sumt == 0. and const_infl_t[infx_cnt] != 0.):
			Cinfl_now = np.zeros((con_infl_C.shape[0], 1))
		
		# if the final input for influxes reached
		if (infx_cnt == -1):
			# influx of components now, convert from ppb/s to molecules/cc.s (air)
			Cinfl_now = (Cinfl[:, infx_cnt]*Cfactor).reshape(-1, 1)
			
		# check whether changes occur at start of this time step
		if (sumt == const_infl_t[infx_cnt] and (infx_cnt != -1)):
			
			# influx of components now, convert from ppb/s to molecules/cc.s (air)
			Cinfl_now = (Cinfl[:, infx_cnt]*Cfactor).reshape(-1, 1)
			
			# update index counter for constant influxes - used in integrator below
			if (infx_cnt<(Cinfl.shape[1]-1)):
				infx_cnt += 1
			else:
				infx_cnt = -1 # reached end
			bc_red = 0 # reset flag for time step reduction due to boundary conditions
			
		# check whether changes occur during proposed integration time step
		if (sumt+tnew > const_infl_t[infx_cnt] and (infx_cnt != -1)):
			# if yes, then reset integration time step so that next step coincides 
			# with change
			tnew = const_infl_t[infx_cnt]-sumt
			bc_red = 1 # flag for time step reduction due to boundary conditions
			
	else: # if no influxes, provide filler
		Cinfl_now = np.zeros((1, 1))
	
	# check on nucleation ---------------------------------------------------------
	# if automatic time step adaption to nucleation requested, check whether number of new particles
	# exceeds 10 % of total number formed during nucleation event.  Second part of condition is that
	# the specifiec nucleation event has not yet reached its defined finishing particle number
	# concentration (#/cc (air))
	if ((nuc_ad == 1) and (new_part_sum1 < nucv1*0.9) and ((num_sb-wall_on) > 0)):
		# the time step (s) needed to increase number concentration of nucleated particles by 10 %
		t_need = (0.1*nucv1+new_part_sum1)
		t_need = np.log(t_need/nucv1)
		t_need = np.log(t_need/nucv2)
		t_need = t_need*nucv3*-1.-sumt
	
		if tnew>t_need: # if suggested time step exceeds this, then reduce to required time step 
			tnew = t_need
			update_stp = t_need
			update_count = 0.
			bc_red = 1
			
	# nucleation check end -------------------------------------------------------------------------	


	return(temp_now, Pnow, lightm, light_time_cnt, tnew, bc_red, update_stp, update_count, 
		Cinfl_now, seedt_cnt, Cfactor, infx_cnt, gasinj_cnt, DStar_org)
