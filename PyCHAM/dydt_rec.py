##########################################################################################
#                                                                                        #
#    Copyright (C) 2018-2025 Simon O'Meara : simon.omeara@manchester.ac.uk               #
#                                                                                        #
#    All Rights Reserved.                                                                #
#    This file is part of PyCHAM                                                         #
#                                                                                        #
#    PyCHAM is free software: you can redistribute it and/or modify it under             #
#    the terms of the GNU General Public License as published by the Free Software       #
#    Foundation, either version 3 of the License, or (at your option) any later          #
#    version.                                                                            #
#                                                                                        #
#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT               #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              #
#    details.                                                                            #
#                                                                                        #
#    You should have received a copy of the GNU General Public License along with        #
#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                #
#                                                                                        #
##########################################################################################
'''module for calculating and recording change tendency (# molecules/cm3/s)
 of components'''
# changes due to gas-phase photochemistry and partitioning are included; 
# generated in init_conc and treats loss from gas-phase as negative

# File Created at 2025-08-13 10:44:11.661383

import numpy as np 

def dydt_rec(y, reac_coef, step, num_sb, num_comp,
kelv_fac, kimt, act_coeff, dydt_erh_flag, H2Oi, wat_hist, self):
	
	# number of particle size bins excluding wall
	num_asb = num_sb-self.wall_on
	
	# loop through components to record the tendency of
	#change, note that components can be grouped, e.g. RO2 for non-HOM-RO2 
	dydtnames = self.dydt_vst['comp_names'] 
	for comp_name in dydtnames: # get name of this component
		
		key_name = str(str(comp_name) + '_comp_indx') # get index of this component
		compi = self.dydt_vst[key_name] 
		# open relevant dictionary value containing reaction numbers and results 
		key_name = str(str(comp_name) + '_res')
		dydt_rec = self.dydt_vst[key_name]
		if hasattr(self, 'sim_ci_file'):
			# open relevant dictionary value containing estimated continuous influx
			key_name = str(str(comp_name) + '_ci')
			ci_array = self.dydt_vst[key_name]
		# open relevant dictionary value containing flag
		# for whether component is reactant (1) or product (0) in each reaction 
		key_name = str(str(comp_name) + '_reac_sign')
		reac_sign = self.dydt_vst[key_name] 
		# keep count on relevant reactions 
		reac_count = 0 
		# loop through relevant reactions 
		# note that final three rows are for
		# particle- and wall-partitioning and dilution 
		reac_indx = np.zeros((len(dydt_rec[0, 0:-3]))).astype('int') 
		reac_indx[:] = (dydt_rec[0, 0:-3]).astype('int') 
		for indx in reac_indx: 
			# estimate gas-phase change tendency for each reaction involving this component 
			gprate = ((y[self.rindx_g[indx, 0:self.nreac_g[indx]]]**self.rstoi_g[indx, 0:self.nreac_g[indx]]).prod())*reac_coef[indx] 
			dydt_rec[step+1, reac_count] += reac_sign[reac_count]*((gprate))
			reac_count += 1 # keep count on relevant reactions 
			
		# now estimate and record tendency to change due
		# to particle- and wall-partitioning  
		# particle-partitioning 
		# if efflorescence has occurred (modelled inside act_coeff_update), then need
		# to account for the immediate transfer of water from particles to gas
		if (wat_hist == 0):
			if (dydt_erh_flag != 0): 
				dydt_rec[step+1, reac_count] += dydt_erh_flag
			else: 
				dydt_rec[step+1, reac_count] = 0
		else:
			for ibin in range(num_asb-1): # particle bin loop
				Csit = y[num_comp*(ibin+1):num_comp*(ibin+2)]
				conc_sum = np.zeros((1)) 
				if (sum(sum(self.pconc > 0.)) > 0): # if seed particles present 
					conc_sum[0] = ((Csit.sum()-sum(Csit[self.seedi]))+sum(Csit[self.seedi]*self.core_diss))
				else: 
					conc_sum[0] = Csit.sum() 
				# prevent numerical error due to division by zero 
				ish = (conc_sum == 0.) 
				conc_sum[ish] = 1.e-40 
				# particle surface gas-phase concentration (# molecules/cm3 (air)) 
				Csit = (Csit[compi]/conc_sum)*self.Psat[0, compi]*kelv_fac[ibin, 0]*act_coeff[0, compi] 
				# partitioning rate (# molecules/cc.s) 
				dydt_all = kimt[ibin, compi]*(y[compi]-Csit) 
				# gas-phase change (# molecules/cm3/s) 
				dydt_rec[step+1, reac_count] -= np.sum(dydt_all) 
				
		# wall-partitioning 
		if (self.wall_on > 0): 
			# concentration at walls (# molecules/cm3 (air)) 
			Csit = y[num_comp*(num_asb+1)::].reshape(self.wall_on, num_comp) 
			Csit = (self.Psat[-self.wall_on::, :]*(Csit/self.Cw)*act_coeff[-self.wall_on::, :])
			dydt_all = (self.kw)*(y[0:num_comp]-Csit) 
			# gas-phase change due to gas-wall partitioning(# molecules/cm3/s), 
			# note sum over all wall bins 
			dydt_rec[step+1, reac_count+1] -= np.sum(dydt_all[:, compi]) 
		
		# dilution
		if (self.dil_fac_now > 0):
			dydt_rec[step+1, reac_count+2] -= y[compi]*self.dil_fac_now
	
		# estimated continuous influx
		if (hasattr(self, 'sim_ci_file') and step < ci_array.shape[0]):
			ci_array[step, 0] = -1*np.sum(dydt_rec[step+1, :])
	
	return(self) 
