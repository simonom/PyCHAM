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
'''generates an array of component indices for the components that constitute a particular component type'''
# for peroxy radicals makes a two column array, with the first column giving the 
# index of components
# included in the peroxy radical list and the
# second column containing the index based on where the RO2 occurs in comp_namelist 
# (i.e., relative
# to all other components from the chemical scheme)

import numpy as np

def group_indices(Hcount, SMILES, i, self):

	# inputs: ----------------------------------------------------------------
	# Hcount - number of hydrogens in component
	# SMILES - SMILES being considered
	# i - index of component being considered
	# self - reference to PyCHAM
	# ------------------------------------------------------------------------

	# carbon and oxygen numbers of this component
	Cn = SMILES.count('C') + SMILES.count('c')
	On = SMILES.count('O') + SMILES.count('o')

	# if name string of this component matches a string section in the list of
	# reactive RO2 names
	if (self.comp_namelist[i] in self.RO2_names):
		
		# check that this name exactly matches (so not just a string subsection)
		# one of those in the list of reactive RO2 names. Note that 
		# self.reac_RO2_indx is initiated in sch_interr.py
		for RO2_namei in self.RO2_names:
			if (RO2_namei == self.comp_namelist[i]):
				self.reac_RO2_indx = np.concatenate((self.reac_RO2_indx, 
					np.array((i)).reshape(1)), axis=0)
				self.reac_RO2_indx = self.reac_RO2_indx.astype('int')
	
	# check if alkoxy radical present in this component 
	# and that component is organic
	if ('[O]' in SMILES):
		if ('C' in SMILES or 'c' in SMILES):
			# if it is an alkoxy radical (rather than alkyl peroxy radical) 
			# add its index to list
			# ensure it's not alkyl peroxy radical
			if ('O[O]' not in SMILES and '[O]O' not in SMILES):
				self.RO_indx.append(i)


	# check for peroxy radical, note that this list of peroxy radicals
	# is based on the SMILES string, and so may be differenet to the list of
	# reactive peroxy radicals in self.reac_RO2_indx
	if ('O[O]' in SMILES or '[O]O' in SMILES):
		if (Cn > 0):
			if (self.RO2_indices.shape[0] == 0): # if first installment
				self.RO2_indices = np.array((1, i)).reshape(1, 2)
			else:
				self.RO2_indices = np.concatenate((self.RO2_indices, 
					np.array((self.RO2_indices[-1, 0]+1, i)).reshape(
					1, 2)), axis = 0)
			if (On >= 6):
				# HOMs peroxy radicals
				self.HOM_RO2_indx.append(i) 

	# check if alkoxy radical present in this component and that component is organic
	if ('[O]' in SMILES):
		if (Cn > 0):
			# if it is an alkoxy radical (rather than alkyl peroxy radical) add 
			# its index to list
			# ensure it's not alkyl peroxy radical
			if ('O[O]' not in SMILES and '[O]O' not in SMILES):
				self.RO_indx.append(i)					

	# check for HOMs and accretion products
	if (On >= 6):
		# highly oxidised molecules
		self.HOMs_indx.append(i)
		if (Cn > 10):
			# HOMs accretion products
			self.ROOR_indx.append(i)

	# check for HOMs present in PRAM, first find the index of the
	# first PRAM species

	# fillers
	C10H15O2O2_indx = 1.e10
	BZo_RO2_O7_indx = 1.e10
		
	if 'C10H15O2O2' in self.comp_namelist:
		C10H15O2O2_indx = self.comp_namelist.index('C10H15O2O2')
		mon_Cnum = 10
	if 'BZo_RO2_O7' in self.comp_namelist:
		BZo_RO2_O7_indx = self.comp_namelist.index('BZo_RO2_O7')
		mon_Cnum = 6

	PRAM_start_indx = min([C10H15O2O2_indx, BZo_RO2_O7_indx])
	
	if (i >= PRAM_start_indx): # if a PRAM component
		self.PRAM_indx.append(i)

		# if a closed-shell HOM PRAM species
		if ('[O]' not in SMILES and On >= 6):
			if (Cn == mon_Cnum):
				# index of closed-shell PRAM monomer
				self.PRAMcsmon_indx.append(i)
			if (Cn < mon_Cnum):
				# index of closed-shell PRAM fragmentation, but
				# note that in 'Methods/The PRAM section' of
				# Roldin et al. 2019 (doi.org/10.1038/s41467-019-12338-8)
				# they explain, that no matter how oxidised an alkoxy radical
				# is, they assume that it fragments (decomposes) to 
				# the MCM species C717O2 and CH3COCH3
				self.PRAMcsfrag_indx.append(i)
				
			if (Cn > mon_Cnum):
				# index of closed-shell PRAM accretion product
				self.PRAMcsacc_indx.append(i)

		# if a peroxy radical HOM PRAM species
		if (On >= 6):
			if ('O[O]' in SMILES or '[O]O' in SMILES ):
				self.PRAMpr_indx.append(i)

	# Ademipo start------------------------------------------------------
	# check for product classes of  HOMs as defined by Baker et al. 2024
	# doi.org/10.5194/acp-24-4789-2024
	if (Cn > 10):
		if ('[O]' not in SMILES and '[o]' not in SMILES):
			if (On >= 4):
				# HOMs accretion products as defined by Baker et al. 2024
				# i is component index
				# append is a python command
				# ROORBaker__index is the list containing the indices
				# of components that fit the HOM accretion product 
				# definition
				# self is a place to store the ROORBaker_indx variable
				# and self.ROORBaker_indx is defined in 
				self.ROORBaker_indx.append(i)
				
        # check for closed-shell fragment product HOMs as defined by Baker et al. 2024
	if (5 <= Cn < 10):
		if ('[O]' not in SMILES and '[o]' not in SMILES):
			if (On >= 4):
                       	 # HOMs fragment products as defined by Baker at al. 2024
				# i is component index
				# append is a python command
				# HOMFragBaker__index is the list containing the indices
				# of components that fit the HOM fragment product
				# defintion
				# self is a place to store the HOMFragBaker_indx variable
				# and self.HOMFragBaker_indx is defined in 
				self.HOMFragBaker_indx.append(i)

        # check for peroxy radical product HOMs as defined by Baker et al. 2024
	if (Cn == 10 and On >=4):
		if ('O[O]' in SMILES or '[O]O' in SMILES or 'o[o]' in SMILES or 
		'[o]o' in SMILES):
			# HOMs peroxy radical products  as defined  by Baker et al. 2024
			# i is component index
			# append is a python command
			# HOMRO2Baker__index is the list containing the indices 
			# of components that fit the HOM fragment product
			# definition
			# self is a place to store the HOMRO2Baker_indx variable
			# and self.HOMRO2Baker_indx is defined in 
			self.HOMRO2Baker_indx.append(i)

	# check for closed-shell monomer product HOMs as defined by Baker et al. 2024
	if (Cn == 10 and On >= 4):
		if ('[O]' not in SMILES and '[o]' not in SMILES):
			# HOMs monomer products as defined by Baker et al. 2024
			# i is component index
			# append is a python command
			# HOMRO2Baker__index is the list containing the indices
			# of components that fit the HOM fragment product
			# definition
			# self is a place to store the HOMMonBaker_indx variable
			# and self.HOMMonBaker_indx is defined in 
			self.HOMMonBaker_indx.append(i)

	# Ademipo finish-----------------------------------------------------


	# check on hydroperoxides and HOM-OOH. Note that HOM-OOH may 
	# not be named identifiably in HOM extensions, so need to 
	# rely on reactants and SMILES	
	# get equation number where this component formed	
	if (Cn >= 1 and On >=2 and Hcount >= 1):
		
		# get the equations where this component formed
		produci = np.sum(self.pindx_g == i, axis = 1) == 1
	
		# reactions with both HO2 as reactant and this component as 
		# product
		if self.HO2i != []:
		
			if (self.rindx_g.ndim == 2):
				eqn_of_inter = produci*(np.sum(self.rindx_g == 
					self.HO2i, axis=1) == 1)  
			if (self.rindx_g.ndim == 1):
				eqn_of_inter = produci*(np.sum(self.rindx_g == self.HO2i) == 1)  

		else:
			eqn_of_inter = np.zeros((self.rindx_g.shape[0]))

		# check if reactant in any of these equations is HO2
		if (sum(eqn_of_inter) > 0): 
			# loop through the reactions of interest
			for eqni in np.where(eqn_of_inter == 1)[0]:

				# get the index of the precursor peroxy radical
				pre_RO2i = self.rindx_g[eqni, :][
					self.rindx_g[eqni, :] != self.HO2i]
				# get SMILE string of precursor RO2
				
				pre_RO2_SMILE = self.rel_SMILES[pre_RO2i[0]]
				# carbon count
				Cn_pre = pre_RO2_SMILE.count('C')+pre_RO2_SMILE.count('c')			
				# oxygen count
				On_pre = pre_RO2_SMILE.count('O')+pre_RO2_SMILE.count('o')

				# check on whether hydrogen number given in product 
				# component name
				if 'H' in self.comp_namelist[i]:
					try: # check if number comes after H letter
						# get index of H character
						Hindx = self.comp_namelist[i].index('H')	
						Hcount = float(self.comp_namelist[i][Hindx+1])
						# in case two digit number after H character
						try:
							Hcount = float(self.comp_namelist[
								i][Hindx+1:Hindx+3])	
						except: # just one number
							Hcount = Hcount
					except:
						Hcount = Hcount

				# do the same hydrogen number in name check for precursor	
				try: # check if number comes after H letter
					# get index of H character
					Hindx = self.comp_namelist[pre_RO2i[0]].index('H')	
					Hn_pre = float(self.comp_namelist[pre_RO2i[0]][Hindx+1])
					try: # in case two digit number after H character
						Hn_pre = float(self.comp_namelist[pre_RO2i[0]][Hindx+1:Hindx+3])	
					except: # just one number
						Hn_pre = Hn_pre
				except:

					# get hydrogen number of precursor component	
					try: # if H count already found for this precursor component
						Hn_pre = self.Hn_list[pre_RO2i[0]]
					except: # if H count not already found for this precursor component
						# if hydrogen is present in this molecule
						if ('H' in self.Pybel_objects[pre_RO2i[0]].formula):

							Hindx_start = self.Pybel_objects[pre_RO2i[0]].formula.index('H')+1
							Hindx_end = Hindx_start
							for Hnum_test in self.Pybel_objects[pre_RO2i[0]].formula[Hindx_start::]:
								try:
									float(Hnum_test) # only continue if this character is a number
									Hindx_end += 1
									if (Hindx_end == len(self.Pybel_objects[pre_RO2i[0]].formula)):
										Hn_pre = float(self.Pybel_objects[pre_RO2i[0]].formula[Hindx_start:Hindx_end])
								except:
									if (Hindx_end != Hindx_start):
										Hn_pre = float(self.Pybel_objects[pre_RO2i[0]].formula[Hindx_start:Hindx_end])
									else:
										Hn_pre = 1. # if no number then Hydrogen must be alone
									break

						else: # if no hydrocarbons
							Hn_pre = 0.
								
				# hydrogen peroxide molecule
				if (Hcount-1 == Hn_pre and Cn == Cn_pre and On == On_pre):	
					self.OOH.append(int(i))
					if (On >= 6): # HOMs hydroperoxides
						self.HOM_OOH.append(int(i))

	# alcohols and carbonyls
	# alcohols form from peroxy radical reaction with the RO2 pool,
	# leading to an alcohol with one less oxygen and one more
	# hydrogen than the peroxy radical
	# carbonyls form from peroxy radical reaction with the RO2 pool,
	# leading to a carbonyl with one less oxygen and one less
	# hydrogen thaN the precursor peroxy radical
	if (Cn >= 1 and On >= 1):

		# get the equations where this component formed
		produci = np.sum(self.pindx_g == i, axis = 1)	

		# get the equations where just one reactant present
		eqi_single_reac = (self.nreac_g == 1)

		# get equations where both this component produced and RO2 in reaction
		# rate coefficient and no other reactants beside the precursor
		eq_of_inter = (produci*eqi_single_reac*self.RO2_in_rrc)

		# check if RO2 present in the reaction rate coefficient in 
		# any of these equations
		if (sum(eq_of_inter) > 0): 	
			# get index of zeros
			eoi_zindx = (eq_of_inter == 0)
			# cumulative sum of equation indices where this component formed
			produci_cs = (np.cumsum(eq_of_inter))
			# ensure zeros remain as zeros
			produci_cs[eoi_zindx] = 0.
			
			# get the reaction number where RO2 present
			for pi in range(int(max(produci_cs))):
				# reaction index
				reac_now = np.where(produci_cs == pi+1)[0][0]
				
				# get the precursor index
				pre_RO2i = self.rindx_g[reac_now, 0]
				# get the precursor carbon and oxygen count
				pCn = (self.rel_SMILES[pre_RO2i].count('C')+
					self.rel_SMILES[pre_RO2i].count('c'))
				pOn = (self.rel_SMILES[pre_RO2i].count('O')+
					self.rel_SMILES[pre_RO2i].count('o'))
				
				# check on whether hydrogen number given in product 
				# component name
				if 'H' in self.comp_namelist[i]:
					try: # check if number comes after H letter
						# get index of H character
						Hindx = self.comp_namelist[i].index('H')	
						Hcount = float(self.comp_namelist[i][Hindx+1])
						# in case two digit number after H character
						try:
							Hcount = float(self.comp_namelist[
								i][Hindx+1:Hindx+3])	
						except: # just one number
							Hcount = Hcount
					except:
						Hcount = Hcount

				# do the same hydrogen number in name check for precursor	
				try: # check if number comes after H letter
					# get index of H character
					Hindx = self.comp_namelist[pre_RO2i].index('H')	
					Hn_pre = float(self.comp_namelist[pre_RO2i][Hindx+1])
					try: # in case two digit number after H character
						Hn_pre = float(self.comp_namelist[
							pre_RO2i][Hindx+1:Hindx+3])	
					except: # just one number
						Hn_pre = Hn_pre
				except:

					# get hydrogen number of precursor component	
					try: # if H count already found for this precursor component
						Hn_pre = self.Hn_list[pre_RO2i]
					except: # if H count not already found for this precursor component
						# if hydrogen is present in this molecule
						if ('H' in self.Pybel_objects[pre_RO2i].formula):

							Hindx_start = self.Pybel_objects[pre_RO2i].formula.index('H')+1
							Hindx_end = Hindx_start
							for Hnum_test in self.Pybel_objects[pre_RO2i].formula[Hindx_start::]:
								try:
									float(Hnum_test) # only continue if this character is a number
									Hindx_end += 1
									if (Hindx_end == len(self.Pybel_objects[pre_RO2i].formula)):
										Hn_pre = float(self.Pybel_objects[pre_RO2i].formula[Hindx_start:Hindx_end])
								except:
									if (Hindx_end != Hindx_start):
										Hn_pre = float(self.Pybel_objects[pre_RO2i].formula[Hindx_start:Hindx_end])
									else:
										Hn_pre = 1. # if no number then Hydrogen must be alone
									break

						else: # if no hydrocarbons
							Hn_pre = 0.
				
				if (Hcount-1 == Hn_pre and Cn == pCn and On+1 == pOn):
					self.OH.append(int(i))
					if (On >= 6): # HOMs alcohols
						self.HOM_OH.append(int(i))
					break # don't loop through anymore reactions	
				if (Hcount+1 == Hn_pre and Cn == pCn and On+1 == pOn):	
					self.carbonyl.append(int(i))
					if (On >= 6): # HOMs carbonyls
						self.HOM_carbonyl.append(int(i))
					break # don't loop through anymore reactions

	# nitrates
	if (self.comp_namelist[i][-3::] == 'NO3'): # MCM and PRAM carbonyl
		if (Cn >= 1):
			self.NO3.append(int(i))
			if (On >= 6):
				self.HOM_NO3.append(int(i)) 

	return(self)
