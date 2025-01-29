########################################################################
#								       #
# Copyright (C) 2018-2025					       #
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
'''isolate sections of chemical scheme'''
# using the chemical scheme markers, sections of the 
# chemical scheme are separated

import re
import numpy as np
import formatting
import openbabel.pybel as pybel

# define function
def sch_interr(total_list_eqn, self):

	# inputs: ------------------------------------------------------------------
	# total_list_eqn - all lines from the chemical scheme file
	# self - reference to PyCHAM
	# self.chem_sch_mrk - markers to denote different section of the 
	#	chemical scheme
	# --------------------------------------------------------------------------
	
	# preparatory part ---------------------------------------------------------
	self.eqn_list = [] # empty list for gas-phase reaction equation
	self.aqeqn_list = [] # empty list for particle-phase reaction equation
	self.sueqn_list = [] # empty list for surface (e.g. wall) reaction equation
	self.RO2_names = [] # empty list for peroxy radicals
	# initiate array of indices of reactive RO2 components
	self.reac_RO2_indx = np.zeros((0)) # initiate with empty array
	self.rrc = [] # empty list for reaction rate coefficients
	self.rrc_name = [] # empty list for reaction rate coefficient labels
	eqn_flag = 0 # don't collate reaction equations until seen
	pr_flag = 0 # don't collate peroxy radicals until seen
	RO2_count = 0 # count on number of lines considered in peroxy radical list
	# if thirteenth marker missing, assume this is due to old inputs
	# when only 12 markers were needed
	if (len(self.chem_sch_mrk) == 12):
		self.chem_sch_mrk.append('')
	# begin by not looking out for component atomic composition
	ac_flag = 0
	# the atom name register
	atom_reg = ['H', 'C', 'N', 'O', 'S', 'Cl']
	# get the molar masses of these components according to pybel (g/mol)
	self.atom_reg_mm = np.zeros((len(atom_reg)))

	# get the molar mass of a single hydrogen atom according to pybel
	Hnum = (pybel.readstring('smi', 'F'))

	# difference between molwt, exactmass functions in pybel given 
	# here (24/01/2025): 
	# https://openbabel.org/api/3.0/
	# classOpenBabel_1_1OBMol.shtml#a7cac960f30506aa53d083983845032df,
	# and I chose molwt (for molecules), which is equivalent to atomic mass
	# for atoms, because this represents the average isotopic abundance,
	# whereas exactmass should only be used when the isotopic abundance
	# has been purposefully set in an experiment

	# molar mass (g/mol) of whole FH molecule
	MM = Hnum.molwt
	for atom in Hnum:
		# subtract the mass of fluorine get the mass
		# of a single hydrogen atom
		MM -= atom.atomicmass

	for atomi in range(len(atom_reg)):
		
		if atom_reg[atomi] == 'H':
			self.atom_reg_mm[atomi] = MM
		else:
			for atom in (pybel.readstring('smi', atom_reg[atomi])):
				self.atom_reg_mm[atomi] = atom.atomicmass

	# if a generic rate constant file provided, 
	# then read in the lines from this
	if (self.rate_cons_name != []):

		# prepare to hold codes and numbers of photolysis rates
		jcode_list = []
		jnum_list = []

		# open the rate constants file
		try:
			f_open_eqn = open(self.rate_cons_name, mode='r')
		# in case in same folder as model variables file
		except:
			pd_indx = self.inname[::-1].index('/')
			pd = self.inname[0:-pd_indx]
			self.rate_cons_path = str(pd + self.rate_cons_name)
			f_open_grc = open(self.rate_cons_path, mode='r')
		
		# read the file and store everything into a list
		total_list_grc = f_open_grc.readlines()
		f_open_grc.close()

		for line in total_list_grc:

			# based on the constants_mcm.f90 file from
			# the mcm website when kpp formatting for chemical scheme
			# file selected
			if (line.count('=') == 0): # if no rate constant
				continue

			# remove white space
			line = line.replace(' ', '')

			# if a photolysis rate constant
			if (line.count('=') == 2 and line.count('::') == 1):
				# get the code for this rate constant
				jcode_now = line[line.index('::')+2:line.index('=')]
			
				# get the photolysis rate number
				try: # in case two digits
					jnum_now = int(line[line.index('!MCM')+6:
						line.index('!MCM')+8])
				except: # in case one digit
					jnum_now = int(line[line.index('!MCM')+6:
						line.index('!MCM')+7])

				jcode_list.append(jcode_now)
				jnum_list.append(jnum_now)

			# if a generic rate constant
			if (line.count('=') == 1 and line.count('!') == 0):

				# convert fortran-type 
				# scientific notation to 
				# python type
				line = formatting.SN_conversion(line)
				# ensure rate coefficient 
				# is python readable
				line = formatting.convert_rate_mcm(line)
				self.rrc.append(line)

	# -------------------------------------------------------------------------
	
	# obtain lists for reaction rate coefficients, peroxy radicals 
	# and equation reactions using markers for separating chemical 
	# scheme elements
	for line in total_list_eqn:
		
		line1 = line.strip() # remove bounding white space
		
	
		# ----------------------------------------------------------------
		# atomic composition per component, only use if not using
		# SMILES to define components
		if ('#DEFVAR' in line1 and self.ac_by_cs == 1):
			# invoke looking out for component atomic composition
			ac_flag = 1
			# prepare to hold (in the following order per component):
			# hydrogen number, carbon number, nitrogen
			# number, oxygen number, sulphur number, chlorine number
			self.ac_dic = {}

		# if this is a line containing atomic composition
		# of component
		if (ac_flag == 1 and '=' in line):

			# get component name
			comp_name = str(line.split('=')[0].replace(' ', ''))
			
			# create new dictionary key-value pair
			self.ac_dic[comp_name] = [0]*len(atom_reg)

			# get atom numbers
			atomn = line.split('=')[1].replace(' ', '')

			# these components currently (23/01/2025) not
			# assigned atomic numbers in MCM kpp v3.3.1 file
			if ('IGNORE' in atomn):
				if (comp_name == 'NA'):
					self.ac_dic[comp_name][atom_reg.index('H')] = 1 
					self.ac_dic[comp_name][atom_reg.index('N')] = 1
					self.ac_dic[comp_name][atom_reg.index('O')] = 3
					continue
				if (comp_name == 'SA'):
					self.ac_dic[comp_name][atom_reg.index('H')] = 2 
					self.ac_dic[comp_name][atom_reg.index('S')] = 1
					self.ac_dic[comp_name][atom_reg.index('O')] = 4
					continue
			
			# attempt getting individual atom numbers
			# starting index for considering
			istart = 0
			# work through string elements
			for atomni in range(1, len(atomn)):

				# if not ready to look for a number, continue to
				# next element of string
				if (atomni <= istart):
					continue
				try:
					anum = int(atomn[istart:atomni])
				# if can't make into an integer, then register atom
				# number
				except:
					try:
						# atom number
						anum = int(atomn[istart:atomni-1])
					except:
						anum = 1
					# get index where atom name ends
					try:
						atomnend = atomni+atomn[
							atomni::].index('+')
					except:
						atomnend = atomni+atomn[
							atomni::].index(';')
					at_name = str(atomn[atomni-1:atomnend])
					# reference index of this atom
					at_indx = atom_reg.index(at_name)
						
					self.ac_dic[comp_name][at_indx] = anum

					# move index for starting to look for atom
					# numbers up
					istart = atomnend+1

		if (ac_flag == 1 and '#INLINE F90_RCONST' in line):
			# stop looking for components and their atomic composition
			ac_flag = 0
			

		# ----------------------------------------------------------
		# generic reaction rate coefficients part
		# marker at end of generic reaction rate coefficient lines
		# the first \ allows python to interpret the second \ as a dash
		# to use in regex which means an escape in case the marker is a
		# regex special character
		# the $ means occurs at end of string
		end_mark = str('\\' + self.chem_sch_mrk[7]+'$')
		# look out for start of generic reaction rate coefficients
		# could be generic reaction coefficient if just one = in line
		
		if (len(line1.split('=')) == 2 and ac_flag == 0):
			rrc_flag = 1
			# don't record if nothing preceding '=' 
			# (can occur in KPP files, e.g.
			# =IGNORE)
			if (len((line1.split('=')[0]).strip()) == 0):
				rrc_flag = 0
			# don't record if this just an IGNORE command
			if len((line1.split('=')[1]).strip()) >= 6:
				if (line1.split('=')[1]).strip()[0:6] == 'IGNORE':
					rrc_flag = 0
			# don't record if marker (if one present) 
			# for end of generic reaction rate
			# coefficient lines not present
			if (len(self.chem_sch_mrk[7]) > 0):
				if re.search(end_mark, line1.strip()) == None:
					rrc_flag = 0

			# don't record if a separate file path is provided for 
			# rate constants
			if (self.rate_cons_name != []):
				rrc_flag = 0			

			if (rrc_flag == 1): 
				
				# dont consider if start of peroxy radical list
				if (line1.split('=')[0]).strip() != self.chem_sch_mrk[1]:
					# don't consider if a gas-phase chemical scheme reaction
					if (((line1.split('=')[0]).strip())[0] != 
						self.chem_sch_mrk[0]):
						# don't consider if an aqueous-phase 
						# chemical scheme reaction
						if (((line1.split('=')[0]).strip())[0] != 
							self.chem_sch_mrk[8]):
							# don't consider if a surface 
							# (e.g. wall) chemical scheme 
							# reaction
							if (((line1.split('=')[0]).strip())[0] 
								!= self.chem_sch_mrk[12]):
							
								# remove end characters
								line2 = line1.replace(str(
								self.chem_sch_mrk[7]), '')
								# remove all white space
								line2 = line2.replace(' ', '')
								# convert fortran-type 
								# scientific notation to 
								# python type
								line2 = \
								formatting.SN_conversion(line2)
								# ensure rate coefficient 
								# is python readable
								line2 = \
								formatting.convert_rate_mcm(
								line2)
								self.rrc.append(line2.strip())
								# get just name of generic
								# reaction rate coefficient
								self.rrc_name.append((
								line2.split('=')[0]).strip())		
			
		# ---------------------------------------------------------------------------
		# peroxy radical part
		# start logging peroxy radicals
		RO2_start_mark = str('^' + self.chem_sch_mrk[1])
		
		# if starting marker for peroxy radical list seen, flag that recording needed
		if (re.match(RO2_start_mark, line1) != None):
			
			# to double check that recording needed for peroxy radicals (in case 
			# self.chem_sch_mrk[1] is not unique)
			
			# first check whether the RO2 list comprises just one line, as this will
			# mean its end marker is present
			if (len(self.chem_sch_mrk[5].strip()) > 0):
				# .* allows search across all elements of 
				# line, \\ ensures marker is
				# recognised as string
				mark = str('.*\\' + self.chem_sch_mrk[5])
				if (re.match(mark, line1) != None):
					pr_flag = 1
			
			# look for presence of marker for RO2 list 
			# continuing onto next line, which
			# confirms this is the RO2 list when it covers more than one line
			# .* allows search across all elements of line, \\ ensures marker is
			# recognised as string
			mark = str('.*\\' + self.chem_sch_mrk[6])
			if (re.match(mark, line1) != None):
				pr_flag = 1
			
			# if line end or continuation marker not supplied then 
			# assume the RO2 start
			# marker is unique
			if (len(self.chem_sch_mrk[5].strip()) == 0 and 
				len(self.chem_sch_mrk[6].strip()) == 0):
				pr_flag = 1
					
		if (pr_flag == 1):
			# get the elements in line separated by peroxy radical separator
			line2 = line1.split(self.chem_sch_mrk[2])
			
			RO2_count += 1 # count on number of lines considered in peroxy radical list
			
			for line3 in line2: # loop through elements in line
				if len(line3.split('='))>1: # in case of RO2 = ...
					line3 = (line3.split('='))[1]
				if len(line3.split(';'))>1: # in case of RO2 list finishing with ...;
					line3 = (line3.split(';'))[0]
				if len(line3.split('&'))>1: # in case of RO2 list finishing with &
					line3 = (line3.split('&'))[0]
				
				# remove any white space
				line3 = line3.strip()
				# don't include white space or ampersands
				if (line3 == '' or line3 == '&'):
					continue
				
				else: # if not these, then first strip surrounding marks
					if line3[0:len(self.chem_sch_mrk[3])] == self.chem_sch_mrk[3]:
						line3 = line3[len(self.chem_sch_mrk[3])::]
					if line3[-len(self.chem_sch_mrk[4])::] == self.chem_sch_mrk[4]:
						line3 = line3[0:-len(self.chem_sch_mrk[4])]
					
					self.RO2_names.append(line3)
			# check for end of RO2 list - given either by marker for end or absence of
			# marker for continuation onto next line of RO2
			# check for marker for end of RO2 list
			if (len(self.chem_sch_mrk[5].strip()) > 0):
				# .* allows search across all elements of line, \\ ensures marker is
				# recognised as string
				mark = str('.*\\' + self.chem_sch_mrk[5])
				if (re.match(mark, line1) != None):
					pr_flag = 0
			
			else: # look for absence of marker for RO2 list continuing onto next line
				# .* allows search across all elements of line, \\ ensures marker is
				# recognised as string
				mark = str('.*\\' + self.chem_sch_mrk[6])
				if (re.match(mark, line1) == None):
					pr_flag = 0
		
		# -------------------------------------------------
		# gas-phase reaction equation part
		# ^ means occurs at start of line and, first \ means 
		# second \ can be interpreted 
		# and second \ ensures recognition of marker
		marker = str('^\\' +  self.chem_sch_mrk[0])
		
		# first check is whether equation start marker is present
		if (re.match(marker, line1) != None):
			
			# second check is whether markers for 
			# starting reaction rate coefficients
			# part, and markers for end of equation lines, are present
			eqn_markers = [str('.*\\' +  self.chem_sch_mrk[9]), 
				str('.*\\' +  self.chem_sch_mrk[11])]
			if (re.match(eqn_markers[0], line1) != None and 
				re.match(eqn_markers[1], line1) != None):

				# if generic rate constant file provided
				if (self.rate_cons_name != []):
				
					# if photolysis rate constant present
					if 'J(' in line1:

						# prepare to hold updated line
						line2 = str(line1[:])

						# loop through line characters
						for si in range(len(line1)):

							if (line2[si:si+2] == 'J('):
								
								# get end index
								ei = si+line2[si::].index(')')	

								# the photolysis rate code
								jcode_now = line2[si+2:ei]
								# get the photolysis rate number 
								# associated with
								# this code
								prn_now = jnum_list[
									jcode_list.index(
									jcode_now)]
								line2 = line2.replace(
									jcode_now, 
									str(prn_now))

						# updated photolysis rate constants
						line1 = str(line2[:])
						
							
				self.eqn_list.append(line1) # store reaction equations

		# aqueous-phase reaction equation part
		# ^ means occurs at start of line and, first \ means second \ can be interpreted 
		# and second \ ensures recognition of marker
		# first, check if a marker given, if not bypass
		if self.chem_sch_mrk[8] == '':
			continue
		else:
			
			marker = str('^\\' +  self.chem_sch_mrk[8])
		
			if (re.match(marker, line1) != None):
				# second check is whether markers for 
				# starting reaction rate coefficients
				# part, and markers for end of equation lines, 
				# are present
				eqn_markers = [str('.*\\' +  self.chem_sch_mrk[9]), 
					str('.*\\' +  self.chem_sch_mrk[11])]
				if (re.match(eqn_markers[0], line1) != None and 
					re.match(eqn_markers[1], line1) != None):
					self.aqeqn_list.append(line1) # store reaction equations

		# surface (e.g. wall) reaction equation part
		# ^ means occurs at start of line and, first \ means 
		# second \ can be interpreted 
		# and second \ ensures recognition of marker
		# first, check if a marker given, if not bypass
		if (self.chem_sch_mrk[12] == ''):
			continue
		else:
			
			marker = str('^\\' +  self.chem_sch_mrk[12])
		
			if (re.match(marker, line1) != None):
				# second check is whether markers for starting reaction 
				# rate coefficients
				# part, and markers for end of equation lines, are present
				eqn_markers = [str('.*\\' +  self.chem_sch_mrk[9]), 
					str('.*\\' +  self.chem_sch_mrk[11])]
				if (re.match(eqn_markers[0], line1) != None and 
					re.match(eqn_markers[1], line1) != None):
					# store reaction equations
					self.sueqn_list.append(line1)

	# number of equations
	self.eqn_num = np.array((len(self.eqn_list), len(self.aqeqn_list), 
		len(self.sueqn_list)))


	return(self)
