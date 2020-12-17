'''module for calculating and recording change tendency of components'''
# changes due to gas-phase photochemistry and partitioning are included

# File Created at 2020-12-15 18:13:17.839200

import numpy as np 

def dydt_rec(y, rindx, rstoi, reac_coef, pindx, pstoi, nprod, step, dydt_vst, nreac, num_sb, num_speci, pconc, core_diss, Psat, kelv_fac, kimt, kwgt, Cw, act_coeff, seedi):
	
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
		for ibin in range(num_sb-1): # size bin loop
			Csit = y[num_speci*(ibin+1):num_speci*(ibin+2)]
			conc_sum = np.zeros((1)) 
			if any(pconc > 0.): # if seed particles present 
				conc_sum[0] = ((Csit.sum()-Csit[seedi])+Csit[seedi]*core_diss)
			else: 
				conc_sum[0] = Csit.sum() 
			# prevent numerical error due to division by zero 
			ish = conc_sum==0.0 
			conc_sum[ish] = 1.0e-40 
			# particle surface gas-phase concentration (molecules/cc (air)) 
			Csit = (Csit[compi]/conc_sum)*Psat[0, compi]*kelv_fac[ibin, 0]*act_coeff[0, compi] 
			# partitioning rate (molecules/cc.s) 
			dydt_all = kimt[ibin, compi]*(y[compi]-Csit) 
			# gas-phase change (molecules/cc/s) 
			dydt_rec[step+1, reac_count] -= dydt_all 
		# wall-partitioning 
		if (kwgt)>1.0e-10: 
			# concentration at wall (molecules/cc (air)) 
			Csit = y[num_speci*num_sb:num_speci*(num_sb+1)] 
			Csit = (Psat[0, :]*(Csit/Cw)*act_coeff[0, compi])
			dydt_all = (kwgt)*(y[compi]-Csit[compi]) 
			# gas-phase change (molecules/cc/s) 
			dydt_rec[step+1, reac_count+1] -= dydt_all 
		
	return(dydt_vst) 
