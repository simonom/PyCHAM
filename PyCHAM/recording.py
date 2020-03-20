'''function to record ode results in model real-time, but not to save a hard copy of results'''
import numpy as np
import ipdb

def recording(y, N_perbin, x, step, sumt, y_mat, Nresult, x2, t_out, tot_stps, 
				num_speci, num_sb, MW, rho, Cn, Vbou):

	# -------------------------------------------------		
	# inputs:
	# MW - molecular weight of components (g/mol)
	# rho - density of components (g/cc)
	# Cn - number of molecules per component (rows) in size bins (columns) (/cc (air))
	# Vbou - volume bounds on the size bins (um3)
	# -------------------------------------------------
	
    
    # set up recording matrices on the first step
	if step == 0:
		y_mat = np.zeros((tot_stps+1, num_speci+num_sb*num_speci))
		if num_sb>0:
			x2 = np.zeros((tot_stps+1, num_sb-1))
			Nresult = np.zeros((tot_stps+1, num_sb-1))
			CSres = np.zeros((tot_stps+1, 2))
		else:
			x2 = 0.0
			Nresult = 0.0
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
					((Cn[-2, ish]/(6.0221409e+23*N_perbin[ish]))*MV[-2]*1.0e12))
		
			# prepare for new number of particles per size bin (# particle/cc (air))
			N_perbin_no_wat = np.zeros((N_perbin.shape[0]))
			# loop through size bins to find number of particles in each 
			# (# particle/cc (air))
			for Ni in range(0, num_sb-1):
				ish = (Vnew>=Vbou[Ni])*(Vnew<Vbou[Ni+1])
				N_perbin_no_wat[Ni] = N_perbin[ish].sum()
	
	
		Nresult[step, :] = N_perbin_no_wat
		if (N_perbin_no_wat<0).sum()>0:
			print('negative particle concentration')
			print(N_perbin)
			print(N_perbin_no_wat)
			ipdb.set_trace()
		x2[step, :] = x # note, this saves the size including contribution of water (um)
	
	else:
		x2 = 0.0
		Nresult = 0.0
	
	t_out[step] = sumt

	return(t_out, y_mat, Nresult, x2)