##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2023 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
#                                                                                       											 #
#    All Rights Reserved.                                                                									 #
#    This file is part of PyCHAM                                                         									 #
#                                                                                        											 #
#    PyCHAM is free software: you can redistribute it and/or modify it under              						 #
#    the terms of the GNU General Public License as published by the Free Software       					 #
#    Foundation, either version 3 of the License, or (at your option) any later          						 #
#    version.                                                                            										 #
#                                                                                        											 #
#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT                						 #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       			 #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              				 #
#    details.                                                                            										 #
#                                                                                        											 #
#    You should have received a copy of the GNU General Public License along with        					 #
#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                 							 #
#                                                                                        											 #
##########################################################################################
'''generates an array of component indices for the components that constitute a particular component type'''
# for peroxy radicals makes a two column array, with the first column giving the index of components
# included in the peroxy radical list and the
# second column containing the index based on where the RO2 occurs in comp_namelist (i.e., relative
# to all other components from the chemical scheme)

import numpy as np

def group_indices(Hcount, SMILES, i, self):

	# inputs: -----------------------------------------------------------------------------------------
	# Hcount - number of hydrogens in component
	# SMILES - SMILES being considered
	# i - index of component being considered
	# self - reference to PyCHAM
	# ---------------------------------------------------------------------------------------------------

	# carbon and oxygen numbers of this component
	Cn = SMILES.count('C') + SMILES.count('c')
	On = SMILES.count('O') + SMILES.count('o')

	# check for peroxy radical
	if ('O[O]' in SMILES or '[O]O' in SMILES):
		if (Cn > 0):
			if (self.RO2_indices.shape[0] == 0): # if first installment
				self.RO2_indices = np.array((1, i)).reshape(1, 2)
			else:
				self.RO2_indices = np.concatenate((self.RO2_indices, np.array((self.RO2_indices[-1, 0]+1, i))), axis = 0)
			if (SMILES.count('O') + SMILES.count('o') >= 6):
				# HOMs peroxy radicals
				self.HOM_RO2_index.append(i) 

	# check if alkoxy radical present in this component and that component is organic
	if ('[O]' in SMILES):
		if (Cn > 0):
			# if it is an alkoxy radical (rather than alkyl peroxy radical) add its index to list
			if ('O[O]' not in SMILES and '[O]O' not in SMILES): # ensure it's not alkyl peroxy radical
				self.RO_indx.append(comp_num)					

	# check for HOMs
	if (SMILES.count('O') + SMILES.count('o') >= 6):
		# HOMs radicals
		self.HOMs_indx.append(i)


	# check on hydroperoxides and HOM-OOH. Note that HOM-OOH may 
	# not be named identifiably in HOM extensions, so need to 
	# rely on reactants and SMILES	
	# get equation number where this component formed	
	if (Cn >= 1 and On >=2 and Hn >= 1):
		print('Cn check', self.comp_namelist[i])
		# get the equations where this component formed
		produci = np.sum(self.pindx_g == i, axis = 1) == 1
	
		# reactions with both HO2 as reactant and this component as 
		# product
		eqn_of_inter = produci*(np.sum(self.rindx_g == self.HO2i, axis=1) == 1)  

		# check if reactant in any of these equations is HO2
		if (sum(eqn_of_inter) > 0): 
			# loop through the reactions of interest
			for eqni in np.where(eqn_of_inter == 1)[0]:

				print('HO2 check', self.comp_namelist[i])
				# get the index of the precursor peroxy radical
				pre_RO2i = self.rindx_g[eqni, :] != self.HO2i
				
				# get SMILE string of precursor RO2
				import ipdb; ipdb.set_trace()
				pre_RO2_SMILE = self.rel_SMILES[pre_RO2i]
				# carbon count
				Cn_pre = pre_RO2_SMILE.count('C')+pre_RO2_SMILE.count('c')			
				# oxygen count
				On_pre = pre_RO2_SMILE.count('O')+pre_RO2_SMILE.count('o')

				# hydroperoxides are formed when a peroxy radical reacts with HO2,
				# and the hydroperoxide has one H atom more than the orginal 
				# peroxy radical
				
				if (Hcount-1 == self.Hn_list[pre_RO2i] and Cn == Cn_pre and On == On_pre):
					print('count check', self.comp_namelist[i])
					import ipdb; ipdb.set_trace()
					self.OOH.append(int(i))
					if (On >= 6): # HOMs hydroperoxides
						self.HOM_OOH.append(int(i))

	# alcohols and carbonyls
	# alcohols form from peroxy radical reaction with the RO2 pool,
	# leading to an alcohol with one less oxygen and one more
	# hydrogen than the peroxy radical
	# carbonyls form from peroxy radical reaction with the RO2 pool,
	# leading to a carbonyl with one less oxygen and one less
	# hydrogen that the precursor peroxy radical
	if (Cn >= 1):

		# get the equations where this component formed
		produci = np.sum(self.pindx_g == i, axis = 1)
		
		# get equations where RO2 in reaction rate coefficient
		RO2_rrci = 'RO2' in self.reac_coef_g

		# get equations where both this component produced and RO2 in reaction
		# rate coefficient
		eq_of_inter = (produci*RO2_rrci)

		# check if RO2 present in the reaction rate coefficient in any of these equations
		if (sum(eq_of_inter) > 0): 
		
			# cumulative sum of equation indices where this component formed
			produci_cs = (np.cumsum(eq_of_inter))[eq_of_inter == 0] = 0.
			# get the reaction number where RO2 present
			for pi in range(max(produci_cs)):
				# reaction index
				reac_now = np.where(produci_cs == pi)[0][0]
				
				# get the precursor index
				prei = self.rindx_g[reac_now, 0]
				# get the precursor carbon and oxygen count
				pCn = self.rel_SMILES.count('C')+self.rel_SMILES.count('c')
				pOn = self.rel_SMILES.count('O')+self.rel_SMILES.count('o')
				pHn = self.Hn_list[prei]

				if (Hcount-1 == pHn and Cn == pCn and On+1 == pOn):
					print('OH count check', self.comp_namelist[i])
					import ipdb; ipdb.set_trace()
					self.OH.append(int(i))
					if (On >= 6): # HOMs alcohols
						self.HOM_OH.append(int(i))
				
				if (Hcount+1 == pHn and Cn == pCn and On+1 == pOn):
					print('=O count check', self.comp_namelist[i])
					import ipdb; ipdb.set_trace()
					self.carbonyl.append(int(i))
					if (On >= 6): # HOMs carbonyls
						self.HOM_carbonyl.append(int(i))

	# nitrates
	if (name_only[-3::] == 'NO3'): # MCM and PRAM carbonyl
		if (Cn >= 1):
			self.NO3.append(int(comp_num))
			if (On >= 6):
				self.HOM_NO3.append(int(comp_num)) 

	return(self)
