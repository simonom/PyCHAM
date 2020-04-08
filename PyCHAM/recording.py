'''function to record ode results in model real-time, but not to save a hard copy of results'''
import numpy as np
import ipdb

def recording(y, N_perbin, x, step, sumt, y_mat, Nresult_dry, Nresult_wet, x2, t_out, 
				tot_stps, num_speci, num_sb, MW, rho, Cn, Vbou, rindx, rstoi, 
				pindx, nprod, dydt_vst, RO2_indices, H2Oi, TEMP, lightm, nreac, 
				pconc, core_diss, Psat, kelv_fac, kimt, kwgt, Cw, timeoday, lat, lon, 
				act_flux_path, DayOfYear, act_coeff, PInit, photo_par_file, Jlen, 
				reac_coef):

	# -------------------------------------------------		
	# inputs:
	# step - number of time step results recorded so far
	# Nresult_dry - stored number concentrations assuming dry distribution
	#				(#particles/cc (air))
	# Nresult_wet - stored number concentrations assuming wet distribution
	#				(#particles/cc (air))
	# MW - molecular weight of components (g/mol)
	# rho - density of components (g/cc)
	# Cn - number of molecules per component (rows) in size bins (columns) (/cc (air))
	# Vbou - volume bounds on the size bins (um3)
	# dydt_vst - dictionary for tracking the tendency to change of user-specified 
	#			components
	# H2Oi - index of water
	# act_coeff - activity coefficients of components
	# PInit - chamber pressure (Pa)
	# photo_par_file - name of file containing photochemical reaction parameters
	# Jlen - number of photochemical reactions
	# reac_coef - reaction rate coefficients during this time step (/s)
	# -------------------------------------------------
	
    
    # set up recording matrices on the first step
	if step == 0:
		y_mat = np.zeros((tot_stps+1, num_speci+num_sb*num_speci))
		if num_sb>0:
			x2 = np.zeros((tot_stps+1, num_sb-1))
			Nresult_dry = np.zeros((tot_stps+1, num_sb-1))
			Nresult_wet = np.zeros((tot_stps+1, num_sb-1))
			CSres = np.zeros((tot_stps+1, 2))
		else:
			x2 = 0.0
			Nresult_dry = 0.0
			Nresult_wet = 0.0
		t_out = np.zeros((tot_stps+1))
		
    
	# note that radius is saved as part of mov_cen_main.py
	y_mat[step, :] = y
	if num_sb>1:
	
		MV = (MW/rho).reshape(num_speci, 1) # molar volume (cc/mol)
	
		Vnew = np.zeros((num_sb-1))
		ish = N_perbin>1.0e-10

		if ish.sum()==0:
			N_perbin_no_wat = np.zeros((num_sb-1))
		else:
			# new volume of single particle per size bin (um3) excluding volume of water
			Vnew[ish] = (np.sum((Cn[:, ish]/(6.0221409e+23*N_perbin[ish]))*MV*1.0e12, 0)-
					((Cn[H2Oi, ish]/(6.0221409e+23*N_perbin[ish]))*MV[H2Oi]*1.0e12))
		
			# prepare for new number of particles per size bin (# particle/cc (air))
			N_perbin_no_wat = np.zeros((N_perbin.shape[0]))
			# loop through size bins to find number of particles in each 
			# (# particle/cc (air))
			for Ni in range(0, num_sb-1):
				ish = (Vnew>=Vbou[Ni])*(Vnew<Vbou[Ni+1])
				N_perbin_no_wat[Ni] = N_perbin[ish].sum()
	
	
		Nresult_wet[step, :] = N_perbin # record with water
		Nresult_dry[step, :] = N_perbin_no_wat # record with water removed
		
		if (N_perbin_no_wat<0).sum()>0:
			print('negative particle concentration')
			print(N_perbin)
			print(N_perbin_no_wat)
			ipdb.set_trace()
		# note, this saves the single particle radius (um) at size bin centre 
		# including contribution of water
		x2[step, :] = x
	
	else:
		x2 = 0.0
		Nresult = 0.0
	
	t_out[step] = sumt
	
	
	if len(dydt_vst)>0:
		# note, using __import__ rather than import allows opening in run time, thereby using
		# updated module
		dydt_rec = __import__('dydt_rec')
		
		# update fraction loss estimate of user-specified species by calling on function
		# generated automatically in eqn_parser
		dydt_vst = dydt_rec.dydt_rec(y, rindx, rstoi, reac_coef, pindx, nprod, step, 
					dydt_vst, nreac, num_sb, num_speci, pconc, core_diss, Psat, kelv_fac, 
					kimt, kwgt, Cw, act_coeff)

	return(t_out, y_mat, Nresult_dry, Nresult_wet, x2, dydt_vst)