'''module for calculating and recording change tendency (molecules/cm3/s) of components'''
# changes due to gas-phase photochemistry and partitioning are included; 
# generated in eqn_pars and treats loss from gas-phase as negative

# File Created at 2021-08-03 11:40:25.917648

import numpy as np 

def dydt_rec(y, rindx, rstoi, reac_coef, pindx, pstoi, nprod, step, dydt_vst, nreac, num_sb, num_comp, pconc, core_diss, Psat, kelv_fac, kimt, kw, Cw, act_coeff, seedi, dydt_erh_flag, H2Oi, wat_hist):
	
	# loop through components to record the tendency of change 
	for compi in dydt_vst.get('comp_index'): 
		
		# open relevant dictionary value 
		dydt_rec = dydt_vst.get(compi) 
		# keep count on relevant reactions 
		reac_count = 0 
		# loop through relevant reactions 
		for i in dydt_rec[0, 0:-2]: # final two rows for particle- and wall-partitioning 
			i = int(i) # ensure reaction index is integer - this necessary because the dydt_rec array is float (the tendency to change records beneath its first row are float) 
			# estimate gas-phase change tendency for each reaction involving this component 
			gprate = ((y[rindx[i, 0:nreac[i]]]**rstoi[i, 0:nreac[i]]).prod())*reac_coef[i] 
			# identify whether this component reacted and/or produced
			# and get its stoichiometric index
			stoi_indx = (np.where(rindx[i, 0:nreac[i]]==compi))[0]
			if len(stoi_indx)>0: 
				dydt_rec[step+1, reac_count] -= rstoi[i, stoi_indx]*((gprate))
			stoi_indx = (np.where(pindx[i, 0:nprod[i]]==compi))[0]
			if len(stoi_indx)>0: 
				dydt_rec[step+1, reac_count] += pstoi[i, stoi_indx]*((gprate))
			reac_count += 1 
			
		# now estimate and record tendency to change due to particle- and wall-partitioning  
		# particle-partitioning 
		# if efflorescence has occurred (modelled inside act_coeff_update), then need
		# to account for the immediate transfer of water from particles to gas
		if (wat_hist == 0):
			if (dydt_erh_flag != 0): 
				dydt_rec[step+1, reac_count] += dydt_erh_flag
			else: 
				dydt_rec[step+1, reac_count] = 0
		else:
			for ibin in range(num_sb-1): # size bin loop
				Csit = y[num_comp*(ibin+1):num_comp*(ibin+2)]
				conc_sum = np.zeros((1)) 
				if any(pconc > 0.): # if seed particles present 
					conc_sum[0] = ((Csit.sum()-sum(Csit[seedi]))+sum(Csit[seedi]*core_diss))
				else: 
					conc_sum[0] = Csit.sum() 
				# prevent numerical error due to division by zero 
				ish = (conc_sum == 0.) 
				conc_sum[ish] = 1.e-40 
				# particle surface gas-phase concentration (molecules/cc (air)) 
				Csit = (Csit[compi]/conc_sum)*Psat[0, compi]*kelv_fac[ibin, 0]*act_coeff[0, compi] 
				# partitioning rate (molecules/cc.s) 
				dydt_all = kimt[ibin, compi]*(y[compi]-Csit) 
				# gas-phase change (molecules/cc/s) 
				dydt_rec[step+1, reac_count] -= dydt_all 
				
		# wall-partitioning 
		if (kw > 1.e-10): 
			# concentration at wall (molecules/cc (air)) 
			Csit = y[num_comp*num_sb:num_comp*(num_sb+1)] 
			Csit = (Psat[0, :]*(Csit/Cw)*act_coeff[0, compi])
			dydt_all = (kw)*(y[compi]-Csit[compi]) 
			# gas-phase change (molecules/cc/s) 
			dydt_rec[step+1, reac_count+1] -= dydt_all 
		
	return(dydt_vst) 
