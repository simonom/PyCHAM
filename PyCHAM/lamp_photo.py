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
'''module for calculating photolysis rates''' 
# photolysis rates calculated for specified absorption cross-section 
# and quantum yield files using actinic flux file'''

import numpy as np
import os

def lamp_photo(J, TEMP, self):

	# --------------------------------------------------------------
	# inputs
	# self - reference to program
	# --------------------------------------------------------------
	
	if (self.af_path != 'nat_act_flux'): # if using actinic fluxes from file
		# open wavelengths (nm) we have total actinic flux for (photon/cm2/nm/s)
		# from chamber
		f = open(self.af_path, 'r')
		wl_chm = np.empty(0) # wavelengths (nm)
		act_chm = np.empty(0) # actinic flux
	
		for line in f: # loop through line
		
			try: # omit any headers
				wl = float((line.strip()).split(',')[0])
				act = float((line.strip()).split(',')[1])
				wl_chm = np.append(wl_chm, (np.array(wl)).reshape(1), axis=0)
				act_chm = np.append(act_chm, (np.array(act)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		
		# ensure that wavelengths and actinic flux have a resolution of 1 nm wavelength
		# to ensure correct integration of photolysis rate over full spectrum
		wl_chm_ref = np.arange(np.min(wl_chm), np.max(wl_chm)+1)
		# unit-resolution wavelength actinic flux (photon/cm2/nm/s)
		act_chm = np.interp(wl_chm_ref, wl_chm, act_chm)
		# unit-resolution wavelength wavelength (nm)
		wl_chm = wl_chm_ref
		
		# apply any wavelength-dependent transmission factor
		if (self.tf_range == 1): # if wavelength-dependent transmission factor supplied
			for wl_cond in self.tf:
				# index where underscore separates wavelength from transmission factor
				us_indx = wl_cond.index('_')
				# dissociate the wavelength from the transmission factor
				wln = float(wl_cond[0:us_indx])
				tfn = float(wl_cond[us_indx+1::])
				
				if 'wlo' in locals(): # if there's an old range
					gt_indx = wl_chm > wln # the greater than index
					lt_indx = wl_chm <= wlo # the lesser than index
					act_chm[gt_indx*lt_indx] = act_chm[gt_indx*lt_indx]*tfn
				else: # if first range
					act_chm[wl_chm > wln] = act_chm[wl_chm>wln]*tfn
					
				# remember old range
				wlo = wln
		
		# if using a sinusoidal diurnal variation in actinic flux
		if (self.light_stat_now == 2):
			# fraction through day
			ftd = self.sumt/(24.*3600.)
			# distance through 2*pi radians
			dist = ftd*2.*np.pi
			# shift, so that peak is at midday
			dist -= np.pi/2.
			# find the sinusoidal function
			act_fact = np.sin(dist)
			# adjust so that peak has value 1, minimum has value 0
			act_fact = (act_fact+1.)/2.
			# apply factor to actinic flux	
			act_chm = act_chm*act_fact

	if (self.af_path == 'nat_act_flux'): # if using modelled solar actinic flux
		
		import nat_act_flux
		import zenith
		
		# get required solar properties from zenith
		(secx, cosx) = zenith.zenith(self)
		
		naff = 4 # flag for calling module
		
		# cosine of the solar zenith angle (radians), note that solar declination angle 
		# (self.dec), latitude (in radians) (self.lat_rad) and local hour angle (lha) are 
		# from zenith.py
		self.mu0 = np.sin(self.dec)*np.sin(self.lat_rad)+np.cos(self.dec)*np.cos(self.lat_rad)*np.cos(self.lha)
		
		# get wavelengths (nm) and their associated actinic fluxes (photon/cm2/nm/s)
		
		# loop through top-of-the-atmosphere downward actinic flux 
		# (photon/cm2/nm/s) per wavelength (nm)
		
		for F0 in F_dep_wl:
			[act_chm] = nat_act_flux.nat_act_flux(A, a, F0, theta, tau, naff, mu0, NL)
		
	
	# get UV-C transmission factor now
	tf_UVCn = self.tf_UVC[(np.sum(self.tf_UVCt<=self.sumt)-1)]

	if (tf_UVCn != 1.):
		# get indices of UVC region (100-280nm)
		UVCindx = np.zeros((wl_chm.shape[0]))
		UVCindx[wl_chm>100.] += 1
		UVCindx[wl_chm<280.] += 1
		UVCindx = UVCindx == 2
		act_chm[UVCindx] = act_chm[UVCindx]*tf_UVCn

	# --------------------------------------------------------------
	# in the below sections photolysis rates are calculated using either the
	# user-supplied file of absorption cross-sections and quantum yields
	# or using those files recommended by MCM 
	# (http://mcm.leeds.ac.uk/MCMv3.3.1/parameters/photolysis.htt)
	
	cwd = os.getcwd() # address of current working directory
	
	
	# determine whether to use MCM-recommended wavelength-dependent absorption 
	# cross-sections (cm2/molecule) and quantum yields (fraction) or other
	if self.photo_path != str(cwd+'/PyCHAM/photofiles/MCMv3.2'): # using user-supplied estimates
		 
		# open file to read
		f = open(self.photo_path, 'r')
		
		# keep count on photolysis reactions, note Fortran indexing to be consistent with 
		# MCM photochemical reaction numbers
		Ji = 1
		
		wlxs = np.empty(0) # will contain wavelengths for cross-sections (nm)
		all_xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		wlqy = np.empty(0) # will contain wavelengths for quantum yields (nm)
		all_qy = np.empty(0) # will contain quantum yields (fraction)
		
		# flags for when to record absorption cross section and quantum yields
		xs_rec = 0
		qy_rec = 0
		
		# loop through lines of file containing absorption cross-sections 
		# and quantum yields
		for line in f:

			# know when end reached - see below this section for the arranging
			# of absorption cross-sections and quantum yields
			if line.strip() == str('J_'+str(Ji+1) + '_axs') or line.strip() == str('J_'+str(Ji+1) + '_qy') or line.strip() =='J_end':

				# absorption cross section (cm2/molecule) interpolation to wavelengths
				# given in actinic flux file
				all_xs = np.interp(wl_chm, wlxs, all_xs)
				# quantum yield (fraction) interpolation to wavelengths
				# given in actinic flux file
				all_qy = np.interp(wl_chm, wlqy, all_qy)
				# note, length of J set in eqn_parser.py
				
				J[Ji] = sum(all_xs*all_qy*act_chm)
				
				# reset to empty
				wlxs = np.empty(0) # wavelengths for cross-sections (nm)
				all_xs = np.empty(0) # absorption cross-sections (cm2/molecule)
				wlqy = np.empty(0) # wavelengths for quantum yields (nm)
				all_qy = np.empty(0) # quantum yields (fraction)
				
				if line.strip() =='J_end':
					continue
				
				Ji += 1 # prepare for next photochemical reaction
				xs_rec = 0
				qy_rec = 0
		
			# absorption cross sections for this reaction starts
			if line.strip() == str('J_'+str(Ji) + '_axs'):
				xs_rec = 1
				qy_rec = 0
				continue
				
			# quantum yields for this reaction starts
			if line.strip() == str('J_'+str(Ji) + '_qy'):
				xs_rec = 0
				qy_rec = 1
				continue
			
			if xs_rec == 1: # when looking for absorption cross-sections
				wl = float(line.split(',')[0])
				xs = float(line.split(',')[1])
				wlxs = np.append(wlxs, (np.array(wl)).reshape(1), axis=0)
				all_xs = np.append(all_xs, (np.array(xs)).reshape(1), axis=0)
				
			
			if qy_rec == 1: # when looking for quantum yields
				wl = float(line.split(',')[0]) # get the wavelength part
				qy = float(line.split(',')[1]) # get the quantum yield at this wavelength
				wlqy = np.append(wlqy, (np.array(wl)).reshape(1), axis=0)
				all_qy = np.append(all_qy, (np.array(xs)).reshape(1), axis=0)
			
			
				
				
		f.close() # close file
	
	
	# Using MCM recommended component absorption cross-sections 
	# and quantum yields (http://mcm.leeds.ac.uk/MCMv3.3.1/parameters/photolysis.htt)
	# in combination with the actinic flux found above
	if (self.photo_path == str(cwd+'/PyCHAM/photofiles/MCMv3.2')):
	
		# --------------------------------------------------------------
		# J<1> and J<2> for O3 (ozone) photolysis
		# cross-section file
		f = open(str(self.photo_path+'/O3/o3_molina86_cs.txt'), 'r')
		wlO3xs = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xsO3 = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
	
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wl = float(line.split('	')[0])
				xs = float(line.split('	')[1])*np.exp(float(line.split('	')[2])/TEMP)
				wlO3xs = np.append(wlO3xs, (np.array(wl)).reshape(1), axis=0)
				xsO3 = np.append(xsO3, (np.array(xs)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xsO3 = np.interp(wl_chm, wlO3xs, xsO3)
		
		# same for J<1> quantum yield O3=O(1D)
		f = open(str(self.photo_path+'/O3/o3_o1d_matsumi02_qy_298.txt'), 'r')
		wlO3qy = np.empty(0) # will contain wavelengths for qy (nm)
		qyO3 = np.empty(0) # will contain qy (dimensionless fraction (0-1))
		
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wl = float(line.split('	')[0])
				qy = float(line.split('	')[1])
				wlO3qy = np.append(wlO3qy, (np.array(wl)).reshape(1), axis=0)
				qyO3 = np.append(qyO3, (np.array(qy)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# quantum yield (J<1>) interpolation
		qyO3 = np.interp(wl_chm, wlO3qy, qyO3)
		
		
		# same for J<2> quantum yield O3=O
		f = open(str(self.photo_path+'/O3/o3_o3p_matsumi02_qy_298.txt'), 'r')
		wlO3Pqy = np.empty(0) # will contain wavelengths for qy (nm)
		qyO3P = np.empty(0) # will contain qy (dimensionless fraction (0-1))
		
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wl = float(line.split('	')[0])
				qy = float(line.split('	')[1])
				wlO3Pqy = np.append(wlO3Pqy, (np.array(wl)).reshape(1), axis=0)
				qyO3P = np.append(qyO3P, (np.array(qy)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# quantum yield (J<2>) interpolation
		qyO3P = np.interp(wl_chm, wlO3Pqy, qyO3P)
		
		
		# photolysis rate for J<1> and J<2> (/s)
		J[1] = sum(xsO3*qyO3*act_chm)
		J[2] = sum(xsO3*qyO3P*act_chm)
		indx = wl_chm<264.

		# --------------------------------------------------------------
		# J<3> for H2O2 (hydrogen peroxide) photolysis: H2O2 = OH + OH
		# cross-section file
		f = open(str(self.photo_path+'/H2O2/h2o2_iupac2003_cs_298.txt'), 'r')
		wlH2O2xs = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xsH2O2 = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
	
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wl = float(line.split('	')[0])
				xs = float(line.split('	')[1])
				wlH2O2xs = np.append(wlH2O2xs, (np.array(wl)).reshape(1), axis=0)
				xsH2O2 = np.append(xsH2O2, (np.array(xs)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xsH2O2 = np.interp(wl_chm, wlH2O2xs, xsH2O2)
		
		# for quantum yield of H2O2 (J<3>), MCM site links to IUPAC recommendation
		# which is 1.0 above a wavelength of 230 nm and which states uncertainty below this.
		# This is true on 04/12/2019.
		# Therefore, here assume quantum yield of one for all wavelengths.
		# photolysis rate for J<3> (/s)
		J[3] = sum(xsH2O2*1.0*act_chm)
		
		# --------------------------------------------------------------
		# J<4> for NO2 (nitrogen dioxide) photolysis: NO2 = NO + O
		# cross-section file
		f = open(str(self.photo_path+'/NO2/no2_iupac03_cs_298.txt'), 'r')
		wlNO2xs = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xsNO2 = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split(' ')[0])
				wl = float(line.split(' ')[0])
				xs = float(line.split(' ')[1])*1.0e-20
				wlNO2xs = np.append(wlNO2xs, (np.array(wl)).reshape(1), axis=0)
				xsNO2 = np.append(xsNO2, (np.array(xs)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xsNO2 = np.interp(wl_chm, wlNO2xs, xsNO2)
		
		# same for J<4> quantum yield: NO2 = NO + O
		f = open(str(self.photo_path+'/NO2/no2_iupac03_qy_298.txt'), 'r')
		wlNO2qy = np.empty(0) # will contain wavelengths for cross-sections (nm)
		qyNO2 = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wl = float(line.split('	')[0])
				qy = float(line.split('	')[1])
				wlNO2qy = np.append(wlNO2qy, (np.array(wl)).reshape(1), axis=0)
				qyNO2 = np.append(qyNO2, (np.array(qy)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		qyNO2 = np.interp(wl_chm, wlNO2qy, qyNO2)
		
		# photolysis rate for J<4>
		J[4] = sum(xsNO2*qyNO2*act_chm)

		# --------------------------------------------------------------
		# J<5> for NO3 (nitrate radical) photolysis: NO3 = NO ;
		# cross-section file
		f = open(str(self.photo_path+'/NO3/no3_iupac03_cs_298.txt'), 'r')
		wlNO3xs = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xsNO3 = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wl = float(line.split('	')[0])
				xs = float(line.split('	')[1])*1.0e-19
				wlNO3xs = np.append(wlNO3xs, (np.array(wl)).reshape(1), axis=0)
				xsNO3 = np.append(xsNO3, (np.array(xs)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xsNO3 = np.interp(wl_chm, wlNO3xs, xsNO3)
		
		# J<5> quantum yield: NO3 = NO
		f = open(str(self.photo_path+'/NO3/no3_no_o2_johnson96_qy_298.txt'), 'r')
		wlNO3qy = np.empty(0) # will contain wavelengths for cross-sections (nm)
		qyNO3 = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wl = float(line.split('	')[0])
				qy = float(line.split('	')[1])
				wlNO3qy = np.append(wlNO3qy, (np.array(wl)).reshape(1), axis=0)
				qyNO3 = np.append(qyNO3, (np.array(qy)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		qyNO3 = np.interp(wl_chm, wlNO3qy, qyNO3)
		
		# photolysis rate for J<5>
		J[5] = sum(xsNO3*qyNO3*act_chm)
		
		# J<6> quantum yield: NO3 = NO2 + O
		f = open(str(self.photo_path+'/NO3/no3_no2_o_johnson96_qy_298.txt'), 'r')
		wlNO3qy = np.empty(0) # will contain wavelengths for cross-sections (nm)
		qyNO3 = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wl = float(line.split('	')[0])
				qy = float(line.split('	')[1])
				wlNO3qy = np.append(wlNO3qy, (np.array(wl)).reshape(1), axis=0)
				qyNO3 = np.append(qyNO3, (np.array(qy)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		qyNO3 = np.interp(wl_chm, wlNO3qy, qyNO3)
		
		# photolysis rate for J<6>
		J[6] = sum(xsNO3*qyNO3*act_chm)
		
		
		# --------------------------------------------------------------
		# J<7> for HONO (nitrous acid) photolysis: HONO = OH + NO ;
		# cross-section file
		f = open(str(self.photo_path+'/HONO/hono_bongartz91_94_cs_298.txt'), 'r')
		wlHONOxs = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xsHONO = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wl = float(line.split('	')[0])
				xs = float(line.split('	')[1])
				wlHONOxs = np.append(wlHONOxs, (np.array(wl)).reshape(1), axis=0)
				xsHONO = np.append(xsHONO, (np.array(xs)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xsHONO = np.interp(wl_chm, wlHONOxs, xsHONO)
		
		# J<7> quantum yield: HONO = OH + NO
		# assume quantum yield of 1.0 for J<7> following recommendation of MCM website on
		# 4/12/2019
		
		# photolysis rate for J<7>
		J[7] = sum(xsHONO*1.0*act_chm)
		
		# --------------------------------------------------------------
		# J<8> for HNO3 (nitric acid) photolysis: HNO3 = OH + NO2;
		# cross-section file
		f = open(str(self.photo_path+'/HNO3/hno3_burkholder93_cs.txt'), 'r')
		wl_xs = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wl = float(line.split('	')[0])
				xsA = float(line.split('	')[1])*1.0e-20
				xsB = float(line.split('	')[2])
				xss = xsA*np.exp(xsB*(TEMP-298.0))
				wl_xs = np.append(wl_xs, (np.array(wl)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl_xs, xs)
		
		# J<8> quantum yield: HNO3 = OH + NO2
		# following MCM website recommendation (true on 4/12/2019), assume 1.0
		
		# photolysis rate
		J[8] = sum(xs*1.0*act_chm)
		
		# --------------------------------------------------------------
		# J<11> for HCHO (formaldehyde) photolysis: HCHO = CO + HO2 + HO2;
		# cross-section file
		f = open(str(self.photo_path+'/HCHO/hcho_meller00_cs.txt'), 'r')
		wl_xs = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wl = float(line.split('	')[0])
				xsA = float(line.split('	')[1])*1.0e-21
				xsB = float(line.split('	')[2])*1.0e-24
				xss = xsA+xsB*(TEMP-298.0)
				wl_xs = np.append(wl_xs, (np.array(wl)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl_xs, xs)
		
		# J<11> and J<12> quantum yield: HCHO = CO + HO2 + HO2 and HCHO = H2 + CO
		f = open(str(self.photo_path+'/HCHO/hcho_iupac03_qy_298.txt'), 'r')
		wl_qy = np.empty(0) # will contain wavelengths for cross-sections (nm)
		qy = np.empty(0) # will contain quantum yields (cm2/molecule)
		qy2 = np.empty(0) # will contain quantum yields (cm2/molecule)
		
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wl = float(line.split('	')[0])
				qys = float(line.split('	')[1])
				qys2 = float(line.split('	')[2])
				wl_qy = np.append(wl_qy, (np.array(wl)).reshape(1), axis=0)
				qy = np.append(qy, (np.array(qys)).reshape(1), axis=0)
				qy2 = np.append(qy2, (np.array(qys2)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# quantum yield interpolation for J<11> and J<12>
		qy = np.interp(wl_chm, wl_qy, qy)
		qy2 = np.interp(wl_chm, wl_qy, qy2)
		
		# photolysis rate
		J[11] = sum(xs*qy*act_chm)
		J[12] = sum(xs*qy2*act_chm)
		
		# -------------------------------------------------------------------
		# J<13> for → CH3 + HCO
		# cross-section file
		f = open(str(self.photo_path+'/CH3CHO/ch3cho_iupac03_cs_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		
		for line in f: # loop through line
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xss = float(line.split('	')[1])
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)	
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# quantum yield file
		f = open(str(self.photo_path+'/CH3CHO/ch3cho_iupac03_qy_298.txt'), 'r')
		wl_qy = np.empty(0) # will contain wavelengths for cross-sections (nm)
		qy = np.empty(0) # will contain quantum yields (cm2/molecule)
		
		for line in f: # loop through line
			try: # omit headers
				float(line.split('	')[0])
				wl = float(line.split('	')[0])
				qys = float(line.split('	')[3])
				wl_qy = np.append(wl_qy, (np.array(wl)).reshape(1), axis=0)
				qy = np.append(qy, (np.array(qys)).reshape(1), axis=0)
			except:
				continue
		f.close() # close file
		# quantum yield interpolation
		qy = np.interp(wl_chm, wl_qy, qy)
		# photolysis rate
		J[13] = sum(xs*qy*act_chm)
		
		# -------------------------------------------------------------------
		# J<14> for → C2H5 + HCO
		# cross-section file
		f = open(str(self.photo_path+'/C2H5CHO/c2h5cho_iupac03_cs_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		
		for line in f: # loop through line
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xss = float(line.split('	')[1])
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)	
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# quantum yield file
		f = open(str(self.photo_path+'/C2H5CHO/c2h5cho_chen&zhu01_qy_298.txt'), 'r')
		wl_qy = np.empty(0) # will contain wavelengths for cross-sections (nm)
		qy = np.empty(0) # will contain quantum yields (cm2/molecule)
		
		for line in f: # loop through line
			try: # omit headers
				float(line.split('	')[0])
				wl = float(line.split('	')[0])
				qys = float(line.split('	')[1])
				wl_qy = np.append(wl_qy, (np.array(wl)).reshape(1), axis=0)
				qy = np.append(qy, (np.array(qys)).reshape(1), axis=0)
			except:
				continue
		f.close() # close file
		# quantum yield interpolation
		qy = np.interp(wl_chm, wl_qy, qy)
		# photolysis rate
		J[14] = sum(xs*qy*act_chm)
		
		# --------------------------------------------------------------
		# J<15> for multiple photolysis
		# cross-section file
		f = open(str(self.photo_path+'/n_C3H7CHO/n_c3h7cho_iupac05_cs_qy_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xss = float(line.split('	')[1])*1.0e-21
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# using recommendation inside n_c3h7cho_iupac05_cs_qy_298 (true on 4/12/2019)
		qy = 0.21
		
		# photolysis rate
		J[15] = sum(xs*qy*act_chm)
		
		# -------------------------------------------------------------------
		# J<16> for → C2H4 + CH3CHO
		# cross-section file
		# use the same cross-section as used in J[15] above
		
		# using recommendation inside n_c3h7cho_iupac05_cs_qy_298 (true on 4/12/2019)
		qy = 0.10
		
		# photolysis rate
		J[16] = sum(xs*qy*act_chm)
		
		# --------------------------------------------------------------
		# J<17> for multiple photolysis
		# cross-section file
		f = open(str(self.photo_path+'/i_C3H7CHO/i-c3h7cho_martinez92_cs_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xss = float(line.split('	')[1])*1.0e-21
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# quantum yield file
		f = open(str(self.photo_path+'/i_C3H7CHO/i_c3h7cho_chen02_hco _qy_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths (nm)
		qy = np.empty(0) # will contain quantum yields (fractions)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				qys = float(line.split('	')[1])
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				qy = np.append(qy, (np.array(qys)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# quantum yield interpolation
		qy = np.interp(wl_chm, wl, qy)
		
		# photolysis rate
		J[17] = sum(xs*qy*act_chm)
		
		# --------------------------------------------------------------
		# J<18> for → CH2=CCH3 + HCO and J<19> for → CH2=C(CH3)CO + H
		# cross-section file
		f = open(str(self.photo_path+'/MACR/macr_iupac05_cs_qy_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		for line in f: # loop through line
	
			try: # omit headers
				float(line.split('	')[0])
				float(line.split('	')[1])
				wls = float(line.split('	')[0])
				xss = float(line.split('	')[1])
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# quantum yield
		# for J<18> and J<19> we take the recommendation of 'macr_iupac05_cs_qy_298.txt' and 
		# use 1.95e-3 (true on 20/12/2019)
		qy = 1.95e-3
		
		# photolysis rate
		J[18] = sum(xs*qy*act_chm)
		J[19] = sum(xs*qy*act_chm)
		
		# --------------------------------------------------------------
		# J<20> for → CH3C(CHO)=CHCH2O + OH
		# cross-section file
		f = open(str(self.photo_path+'/C5HPALD1/C5HPALD_cs_qy_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xss = float(line.split('	')[1])
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# quantum yield
		# for J<20> take the recommendation of 'C5HPALD_cs_qy_298.txt' and 
		# use 1.0 (true on 20/12/2019)
		qy = 1.0
		
		# photolysis rate
		J[20] = sum(xs*qy*act_chm)
		
		# --------------------------------------------------------------
		# J<21> for CH3COCH3 (acetone) photolysis: CH3COCH3 = CH3CO3 + CH3O2;
		# cross-section file
		f = open(str(self.photo_path+'/CH3COCH3/ch3coch3_iupac05_cs.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xsA = float(line.split('	')[1])
				xsB = float(line.split('	')[2])
				xsC = float(line.split('	')[3])
				xsD = float(line.split('	')[4])
				xss = xsA*(1.0+xsB*TEMP+xsC*TEMP**2.0+xsD*TEMP**3.0)
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		f = open(str(self.photo_path+'/CH3COCH3/ch3coch3_iupac05_qy_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths (nm)
		qy = np.empty(0) # will contain quantum yields
		for line in f: # loop through lines
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				qys = float(line.split('	')[3])
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				qy = np.append(qy, (np.array(qys)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# quantum yield interpolation
		qy = np.interp(wl_chm, wl, qy)
		
		# photolysis rate
		J[21] = sum(xs*qy*act_chm)
		
		# --------------------------------------------------------------
		# J<22> for → CH3CO + C2H5
		# cross-section file
		f = open(str(self.photo_path+'/MEK/ch3coc2h5_iupac05_cs_qy_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xss = float(line.split('	')[1])
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# quantum yield, following recommendation from ch3coc2h5_iupac05_cs_qy_298
		# (true on 4/12/2019) to use 0.16
		qy = 0.16
		
		# photolysis rate
		J[22] = sum(xs*qy*act_chm)
		
		# --------------------------------------------------------------
		# J<23> for → CH3CH=CH2 + CO and J<24> for → CH3CO + CH2=CH
		# cross-section file
		f = open(str(self.photo_path+'/MVK/mvk_iupac05_cs_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xss = float(line.split('	')[1])
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# quantum yield file
		f = open(str(self.photo_path+'/MVK/mvk_iupac05_qy_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		qy = np.empty(0) # will contain quantum yields for J<23> (fraction)
		qy24 = np.empty(0) # will contain quantum yields for J<24> (fraction)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				qys = float(line.split('	')[1])
				qys24 = float(line.split('	')[2])
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				qy = np.append(qy, (np.array(qys)).reshape(1), axis=0)
				qy24 = np.append(qy24, (np.array(qys24)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
	
		# quantum yield interpolation for J<23>
		qy = np.interp(wl_chm, wl, qy)
		# quantum yield interpolation for J<24>
		qy24 = np.interp(wl_chm, wl, qy24)
		
		# photolysis rate
		J[23] = sum(xs*qy*act_chm)
		J[24] = sum(xs*qy24*act_chm)
		
		# --------------------------------------------------------------
		# J<31>, J<32>, J<33> for glyoxal photolysis
		# cross-section file
		f = open(str(self.photo_path+'/CHOCHO/chocho_volkamer05_cs_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xss = float(line.split('	')[1])*1.0e-20
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		f = open(str(self.photo_path+'/CHOCHO/chocho_tadic06_qy_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths (nm)
		qy = np.empty(0) # will contain quantum yields
		qy2 = np.empty(0) # will contain quantum yields
		qy3 = np.empty(0) # will contain quantum yields
		for line in f: # loop through lines
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				qys = float(line.split('	')[1])
				qys2 = float(line.split('	')[3])
				qys3 = float(line.split('	')[2])
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				qy = np.append(qy, (np.array(qys)).reshape(1), axis=0)
				qy2 = np.append(qy2, (np.array(qys2)).reshape(1), axis=0)
				qy3 = np.append(qy3, (np.array(qys3)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# quantum yield interpolation
		qy = np.interp(wl_chm, wl, qy)
		qy2 = np.interp(wl_chm, wl, qy2)
		qy3 = np.interp(wl_chm, wl, qy3)
		
		# photolysis rate
		J[31] = sum(xs*qy*act_chm)
		J[32] = sum(xs*qy2*act_chm)
		J[33] = sum(xs*qy3*act_chm)
		
		# --------------------------------------------------------------
		# J<34> for → CH3CO + HCO
		# cross-section file
		f = open(str(self.photo_path+'/CH3COCHO/ch3cocho_jpl_iupac05_cs_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xss = float(line.split('	')[1])*1.0e-20
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# quantum yield
		f = open(str(self.photo_path+'/CH3COCHO/ch3cocho_iupac05_qy_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths (nm)
		qy = np.empty(0) # will contain quantum yields
	
		for line in f: # loop through lines
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				qys = float(line.split('	')[2])
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				qy = np.append(qy, (np.array(qys)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# quantum yield interpolation
		qy = np.interp(wl_chm, wl, qy)
		
		# photolysis rate
		J[34] = sum(xs*qy*act_chm)
		
		# --------------------------------------------------------------
		# J<35> for → CH3CO + CH3CO
		# cross-section file
		f = open(str(self.photo_path+'/BIACET/biacet_horowitz01_cs_plum83_qy_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xss = float(line.split('	')[1])*1.0e-20
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# following quantum yield value given in biacet_horowitz01_cs_plum83_qy_298.txt
		# (true on 4/12/2019) 
		qy = 0.158
		
		# photolysis rate
		J[35] = sum(xs*qy*act_chm)
		
		# --------------------------------------------------------------
		# J<41> for → CH3O + OH
		# cross-section file
		f = open(str(self.photo_path+'/CH3OOH/ch3ooh_iupac05_cs_qy_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xss = float(line.split('	')[1])*1.0e-20
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# follow quantum yield value given in ch3ooh_iupac05_cs_qy_298.txt
		# (true on 4/12/2019)
		qy = 1.0
		
		# photolysis rate
		J[41] = sum(xs*qy*act_chm)
		
		# --------------------------------------------------------------
		# J<51> for → CH3O + NO2
		# cross-section file
		f = open(str(self.photo_path+'/CH3ONO2/ch3ono2_iupac05_cs_qy.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xss = float(line.split('	')[1])*1.0e-20
				xssB = float(line.split('	')[2])*1.0e-3
				xss = xss*np.exp(xssB*(TEMP-298.0))
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# follow quantum yield value given in ch3ono2_iupac05_cs_qy.txt
		# (true on 4/12/2019)
		qy = 1.0
		
		# photolysis rate
		J[51] = sum(xs*qy*act_chm)
	
		# --------------------------------------------------------------
		# J<52> for → C2H5O + NO2
		# cross-section file
		f = open(str(self.photo_path+'/CH3CH2ONO2/ch3ch2ono2_iupac05_cs_qy.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xss = float(line.split('	')[1])
				xssB = float(line.split('	')[2])
				xss = xss*np.exp(xssB*(TEMP-298.0))
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# follow quantum yield value given in ch3ch2ono2_iupac05_cs_qy.txt
		# (true on 4/12/2019)
		qy = 1.0
		
		# photolysis rate
		J[52] = sum(xs*qy*act_chm)
	
		# --------------------------------------------------------------
		# J<53> for → n-C3H7O + NO2
		# cross-section file
		f = open(str(self.photo_path+'/NOA/noa_barnes93_cs_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xss = float(line.split('	')[1])
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# quantum yield
		f = open(str(self.photo_path+'/NOA/noa_estimated_qy_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths (nm)
		qy = np.empty(0) # will contain quantum yields
	
		for line in f: # loop through lines
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				qys = float(line.split('	')[1])
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				qy = np.append(qy, (np.array(qys)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# quantum yield interpolation
		qy = np.interp(wl_chm, wl, qy)
		
		# photolysis rate
		J[53] = sum(xs*qy*act_chm)
		
		# --------------------------------------------------------------
		# J<54> for → CH3C(O.)CH3 + NO2
		# cross-section file
		f = open(str(self.photo_path+'/CH3CH2ONO2/ch3ch2ono2_iupac05_cs_qy.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xss1 = float(line.split('	')[1])
				xss2 = float(line.split('	')[2])
				xss = xss1 #+xss2*(TEMP-298.0)
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# quantum yield uses value given in ch3ch2ono2_iupac05_cs_qy.txt
		# (true on 4/12/2019)
		qy=1.0
		
		# photolysis rate
		J[54] = sum(xs*qy*act_chm)
	
		# --------------------------------------------------------------
		# J<55> for → t-C4H9O + NO2
		# comparing the link given in the MCM site: http://mcm.leeds.ac.uk/MCMv3.3.1/parameters/photolysis/t_C4H9ONO2/t_c4h9ono2_Roberts&Fajar89_cs_298.txt
		# for tert-butyl nitrate to the download folder of photolysis rates,
		# this process does not seem to be included in the download, so manual input provided
		wl = np.array(([270.0, 275.0, 280.0, 285.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0]))
		xs = np.array(([4.3e-20, 4.0e-20, 3.7e-20, 3.1e-20, 2.6e-20, 2.0e-20, 1.5e-20, 1.0e-20, 7.0e-21, 4.5e-21, 2.7e-21, 1.50e-21, 8.60e-22]))
		
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# quantum yield uses value given in http://mcm.leeds.ac.uk/MCMv3.3.1/parameters/photolysis/t_C4H9ONO2/t_c4h9ono2_Roberts&Fajar89_cs_298.txt
		# (true on 4/12/2019)
		qy=1.0
		
		# photolysis rate
		J[55] = sum(xs*qy*act_chm)
		
		# --------------------------------------------------------------
		# J<56> for → CH3C(O)CH2(O.) + NO2 and → CH3CO + HCHO + NO2
		# cross-section file
		f = open(str(self.photo_path+'/NOA/noa_barnes93_cs_298.txt'), 'r')
		wl = np.empty(0) # will contain wavelengths for cross-sections (nm)
		xs = np.empty(0) # will contain absorption cross-sections (cm2/molecule)
		for line in f: # loop through line
			
			try: # omit headers
				float(line.split('	')[0])
				wls = float(line.split('	')[0])
				xss = float(line.split('	')[1])
				wl = np.append(wl, (np.array(wls)).reshape(1), axis=0)
				xs = np.append(xs, (np.array(xss)).reshape(1), axis=0)
				
			except:
				continue
		f.close() # close file
		# absorption cross section (cm2/molecule) interpolation
		xs = np.interp(wl_chm, wl, xs)
		
		# quantum yield uses value given in MCM website table (http://mcm.leeds.ac.uk/MCMv3.3.1/parameters/photolysis.htt)
		# (true on 4/12/2019)
		qy=0.9
		
		# photolysis rate
		J[56] = sum(xs*qy*act_chm)
	
	return J
