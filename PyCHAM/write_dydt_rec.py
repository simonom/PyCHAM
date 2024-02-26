########################################################################
#								       #
# Copyright (C) 2018-2024					       #
# Simon O'Meara : simon.omeara@manchester.ac.uk			       #
#								       #
# All Rights Reserved.                                                 #
# This file is part of PyCHAM                                          #
#                                                                      #
# PyCHAM is free software: you can redistribute it and/or modify it    #
# under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or    #
# (at  your option) any later version.                                 #
#                                                                      #
# PyCHAM is distributed in the hope that it will be useful, but        #
# WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# General Public License for more details.                             #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with PyCHAM.  If not, see <http://www.gnu.org/licenses/>.      #
#                                                                      #
########################################################################
'''generates module that records change tendencies'''
# function to automatically generate a module that is used to record the 
# tendency
# of components (components specified by the user, their index given by 
# the tracked_comp user input) 
# to change in response to box model 
# mechanisms - with the resulting function called on each time step

import datetime

def write_dydt_rec(self): # define function

	# inputs -------
	# self - the PyCHAM class
	# --------------

	f = open(self.PyCHAM_path + '/PyCHAM/dydt_rec.py', mode='w')
	f.write('##########################################################################################\n')
	f.write('#                                                                                        #\n')
	f.write('#    Copyright (C) 2018-2024 Simon O\'Meara : simon.omeara@manchester.ac.uk               #\n')
	f.write('#                                                                                        #\n')
	f.write('#    All Rights Reserved.                                                                #\n')
	f.write('#    This file is part of PyCHAM                                                         #\n')
	f.write('#                                                                                        #\n')
	f.write('#    PyCHAM is free software: you can redistribute it and/or modify it under             #\n')
	f.write('#    the terms of the GNU General Public License as published by the Free Software       #\n')
	f.write('#    Foundation, either version 3 of the License, or (at your option) any later          #\n')
	f.write('#    version.                                                                            #\n')
	f.write('#                                                                                        #\n')
	f.write('#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT               #\n')
	f.write('#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       #\n')
	f.write('#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              #\n')
	f.write('#    details.                                                                            #\n')
	f.write('#                                                                                        #\n')
	f.write('#    You should have received a copy of the GNU General Public License along with        #\n')
	f.write('#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                #\n')
	f.write('#                                                                                        #\n')
	f.write('##########################################################################################\n')
	f.write('\'\'\'module for calculating and recording change tendency (# molecules/cm3/s)\n')
	f.write(' of components\'\'\'\n')
	f.write('# changes due to gas-phase photochemistry and partitioning are included; \n')
	f.write('# generated in init_conc and treats loss from gas-phase as negative\n')
	f.write('\n')
	f.write('# File Created at %s\n' %(datetime.datetime.now()))
	f.write('\n')
	f.write('import numpy as np \n')
	f.write('\n')
	# following part is the function (there should be an indent at the start of each line)
	# suggest 1 Tab
	f.write('def dydt_rec(y, reac_coef, step, num_sb, num_comp, core_diss, kelv_fac, kimt, act_coeff, dydt_erh_flag, H2Oi, wat_hist, pconc, self):\n')
	f.write('	\n')
	f.write('	# number of particle size bins excluding wall\n')
	f.write('	num_asb = num_sb-self.wall_on\n')
	f.write('	\n')
	f.write('	# loop through components to record the tendency of change, note that components can be grouped, e.g. RO2 for non-HOM-RO2 \n')
	f.write('	dydtnames = self.dydt_vst[\'comp_names\'] \n')
	f.write('	for comp_name in dydtnames: # get name of this component\n')
	f.write('		\n')
	f.write('		key_name = str(str(comp_name) + \'_comp_indx\') # get index of this component\n')
	f.write('		compi = self.dydt_vst[key_name] \n')
	f.write('		# open relevant dictionary value containing reaction numbers and results \n')
	f.write('		key_name = str(str(comp_name) + \'_res\')\n')
	f.write('		dydt_rec = self.dydt_vst[key_name]\n')
	f.write('		if hasattr(self, \'sim_ci_file\'):\n')
	f.write('			# open relevant dictionary value containing estimated continuous influx\n')
	f.write('			key_name = str(str(comp_name) + \'_ci\')\n')
	f.write('			ci_array = self.dydt_vst[key_name]\n')

	f.write('		# open relevant dictionary value containing flag for whether component is reactant (1) or product (0) in each reaction \n')
	f.write('		key_name = str(str(comp_name) + \'_reac_sign\')\n')
	f.write('		reac_sign = self.dydt_vst[key_name] \n')
	f.write('		# keep count on relevant reactions \n')
	f.write('		reac_count = 0 \n')
	f.write('		# loop through relevant reactions \n')
	f.write('		for i in dydt_rec[0, 0:-3]: # final three rows for particle- and wall-partitioning and dilution \n')
	f.write('			i = int(i) # ensure reaction index is integer - this necessary because the dydt_rec array is float (the tendency to change records beneath its first row are float) \n')
	f.write('			# estimate gas-phase change tendency for each reaction involving this component \n')
	f.write('			gprate = ((y[self.rindx_g[i, 0:self.nreac_g[i]]]**self.rstoi_g[i, 0:self.nreac_g[i]]).prod())*reac_coef[i] \n')
	f.write('			dydt_rec[step+1, reac_count] += reac_sign[reac_count]*((gprate))\n')
	f.write('			reac_count += 1 # keep count on relevant reactions \n')
	f.write('			\n')
	f.write('		# now estimate and record tendency to change due to particle- and wall-partitioning  \n')
	f.write('		# particle-partitioning \n')
	f.write('		# if efflorescence has occurred (modelled inside act_coeff_update), then need\n')
	f.write('		# to account for the immediate transfer of water from particles to gas\n')
	f.write('		if (wat_hist == 0):\n')
	f.write('			if (dydt_erh_flag != 0): \n')
	f.write('				dydt_rec[step+1, reac_count] += dydt_erh_flag\n')
	f.write('			else: \n')
	f.write('				dydt_rec[step+1, reac_count] = 0\n')
	f.write('		else:\n')
	f.write('			for ibin in range(num_asb-1): # particle bin loop\n')
	f.write('				Csit = y[num_comp*(ibin+1):num_comp*(ibin+2)]\n')
	f.write('				conc_sum = np.zeros((1)) \n')
	f.write('				if (sum(sum(pconc > 0.)) > 0): # if seed particles present \n')
	f.write('					conc_sum[0] = ((Csit.sum()-sum(Csit[self.seedi]))+sum(Csit[self.seedi]*core_diss))\n')
	f.write('				else: \n')
	f.write('					conc_sum[0] = Csit.sum() \n')
	f.write('				# prevent numerical error due to division by zero \n')
	f.write('				ish = (conc_sum == 0.) \n')
	f.write('				conc_sum[ish] = 1.e-40 \n')
	f.write('				# particle surface gas-phase concentration (# molecules/cm3 (air)) \n')
	f.write('				Csit = (Csit[compi]/conc_sum)*self.Psat[0, compi]*kelv_fac[ibin, 0]*act_coeff[0, compi] \n')
	f.write('				# partitioning rate (# molecules/cc.s) \n')
	f.write('				dydt_all = kimt[ibin, compi]*(y[compi]-Csit) \n')
	f.write('				# gas-phase change (# molecules/cm3/s) \n')
	f.write('				dydt_rec[step+1, reac_count] -= np.sum(dydt_all) \n')
	f.write('				\n')
	f.write('		# wall-partitioning \n')
	f.write('		if (self.wall_on > 0): \n')
	f.write('			# concentration at walls (# molecules/cm3 (air)) \n')
	f.write('			Csit = y[num_comp*(num_asb+1)::].reshape(self.wall_on, num_comp) \n')
	f.write('			Csit = (self.Psat[-self.wall_on::, :]*(Csit/self.Cw)*act_coeff[-self.wall_on::, :])\n')
	f.write('			dydt_all = (self.kw)*(y[0:num_comp]-Csit) \n')
	f.write('			# gas-phase change due to gas-wall partitioning(# molecules/cm3/s), \n')
	f.write('			# note sum over all wall bins \n')
	f.write('			dydt_rec[step+1, reac_count+1] -= np.sum(dydt_all[:, compi]) \n')
	f.write('		\n')
	f.write('		# dilution\n')
	f.write('		if (self.dil_fac_now > 0):\n')
	f.write('			dydt_rec[step+1, reac_count+2] -= y[compi]*self.dil_fac_now\n')
	f.write('	\n')
	f.write('		# estimated continuous influx\n')
	f.write('		if hasattr(self, \'sim_ci_file\'):\n')
	f.write('			ci_array[step, 0] = -1*np.sum(dydt_rec[step+1, :])\n')
	f.write('	\n')
	f.write('	return(self) \n')
	f.close()
