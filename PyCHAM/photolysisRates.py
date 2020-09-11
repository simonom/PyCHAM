'''module for calculating photolysis rates'''
# can use ambient sunlight expression or artificial chamber lights for calculating a
# photolysis rate that is fixed when the lights are on

import scipy
import os
import numpy as np
import requests, zipfile, io # for downloading
import shutil
from lamp_photo import lamp_photo
import zenith

def PhotolysisCalculation(time, lat, lon, TEMP, act_flux_path, DayOfYear, photo_par_file,
							Jlen, tf):

	# inputs:-----------------------------------------------------------------------------	
	# act_flux_path - name of path to file containing known actinic flux (only used if 
	#				lights on inside chamber)
	# DayOfYear - number of days through the calendar year
	# photo_par_file - name of file containing estimates for wavelength-dependent
	# 					absorption cross-sections and quantum yields
	# Jlen - number of photolysis reactions
	# tf - the transmission factor (for natural light intensity)
	# ------------------------------------------------------------------------------------
	
	(secx, cosx) = zenith.zenith(time, lat, lon, DayOfYear)
	
	J = np.zeros((Jlen)) # prepare output
    
	# if using MCM and natural light
	cwd = os.getcwd() # address of current working directory

	if (photo_par_file == str(cwd + '/PyCHAM/photofiles/MCMv3.2') and act_flux_path == 'no'):
		#J          L           M          N
		J[1]=6.073E-05*cosx**(1.743)*np.exp(-1.0*0.474*secx)
		J[2]=4.775E-04*cosx**(0.298)*np.exp(-1.0*0.080*secx)
		J[3]=1.041E-05*cosx**(0.723)*np.exp(-1.0*0.279*secx)
		J[4]=1.165E-02*cosx**(0.244)*np.exp(-1.0*0.267*secx)
		J[5]=2.485E-02*cosx**(0.168)*np.exp(-1.0*0.108*secx)
		J[6]=1.747E-01*cosx**(0.155)*np.exp(-1.0*0.125*secx)
		J[7]=2.644E-03*cosx**(0.261)*np.exp(-1.0*0.288*secx)
		J[8]=9.312E-07*cosx**(1.230)*np.exp(-1.0*0.307*secx)
		J[11]=4.642E-05*cosx**(0.762)*np.exp(-1.0*0.353*secx)
		J[12]=6.853E-05*cosx**(0.477)*np.exp(-1.0*0.323*secx)
		J[13]=7.344E-06*cosx**(1.202)*np.exp(-1.0*0.417*secx)
		J[14]=2.879E-05*cosx**(1.067)*np.exp(-1.0*0.358*secx)
		J[15]=2.792E-05*cosx**(0.805)*np.exp(-1.0*0.338*secx)
		J[16]=1.675E-05*cosx**(0.805)*np.exp(-1.0*0.338*secx)
		J[17]=7.914E-05*cosx**(0.764)*np.exp(-1.0*0.364*secx)
		J[18]=1.140E-05*cosx**(0.396)*np.exp(-1.0*0.298*secx)
		J[19]=1.140E-05*cosx**(0.396)*np.exp(-1.0*0.298*secx)
		J[21]=7.992E-07*cosx**(1.578)*np.exp(-1.0*0.271*secx)
		J[22]=5.804E-06*cosx**(1.092)*np.exp(-1.0*0.377*secx)
		J[23]=1.836E-05*cosx**(0.395)*np.exp(-1.0*0.296*secx)
		J[24]=1.836E-05*cosx**(0.395)*np.exp(-1.0*0.296*secx)
		J[31]=6.845E-05*cosx**(0.130)*np.exp(-1.0*0.201*secx)
		J[32]=1.032E-05*cosx**(0.130)*np.exp(-1.0*0.201*secx)
		J[33]=3.802E-05*cosx**(0.644)*np.exp(-1.0*0.312*secx)
		J[34]=1.537E-04*cosx**(0.170)*np.exp(-1.0*0.208*secx)
		J[35]=3.326E-04*cosx**(0.148)*np.exp(-1.0*0.215*secx)
		J[41]=7.649E-06*cosx**(0.682)*np.exp(-1.0*0.279*secx)
		J[51]=1.588E-06*cosx**(1.154)*np.exp(-1.0*0.318*secx)
		J[52]=1.907E-06*cosx**(1.244)*np.exp(-1.0*0.335*secx)
		J[53]=2.485E-06*cosx**(1.196)*np.exp(-1.0*0.328*secx)
		J[54]=4.095E-06*cosx**(1.111)*np.exp(-1.0*0.316*secx)
		J[55]=1.135E-05*cosx**(0.974)*np.exp(-1.0*0.309*secx)
		J[56]=7.549E-06*cosx**(1.015)*np.exp(-1.0*0.324*secx)
		J[57]=3.363E-06*cosx**(1.296)*np.exp(-1.0*0.322*secx)
		J[61]=7.537E-04*cosx**(0.499)*np.exp(-1.0*0.266*secx)

		J = J*tf

	# from MAC spectral analysis and Mainz database (xsproc.py)
# 	J[1] = 2.3706768705670786e-05
# 	J[2] = 7.128976494635883e-05
# 	J[3] = 1.26339378672e-05
# 	J[4] = 0.0011216654096221574
# 	J[5] = 0.0024080663555999995
# 	J[6] = 0.02122656504799999
# 	J[7] = 0.0002830012299999995
# 	J[8] = 2.8579387129999973e-07
# 	J[11] = 1.3984589959570004e-05
# 	J[12] = 6.155370761783531e-06
# 	J[15] = 3.657979548045003e-05
# 	J[21] = 2.4121216268999992e-06
# 	J[22] = 1.0067723552799998e-05
# 	J[31] = 8.998612475770001e-07
# 	J[32] = 3.710891288701262e-06
# 	J[33] = 8.468352055021745e-06
# 	J[34] = 1.039182783128049e-05
# 	J[35] = 3.1229494794723276e-05
# 	J[41] = 0.014555553853016991
	
	# if a file path for the chamber's actinic flux is supplied, then use the MCM 
	# equations for photolysis rates as a function of wavelength-dependent actinic flux
	# this is done if artificial lights are turned on
	if (act_flux_path != 'no'):
		# using MCM recommended values
# 		folder_url = "http://mcm.leeds.ac.uk/MCMv3.3.1/parameters/photolysis/MCMv3.2photolysis.zip"
# 		r = requests.get(folder_url) # create HTTP response object
# 		if r.ok==0: # check object is ready
# 			print('Photolysis information not found at the folowing URL, please check and update in PhotolysisRates.py: http://mcm.leeds.ac.uk/MCMv3.3.1/parameters/photolysis/MCMv3.2photolysis.zip')
# 		z = zipfile.ZipFile(io.BytesIO(r.content)) #  download object
# 		z.extractall('PyCHAM/MCMphotofiles') # creates folder and stores files there
		
		# call on MCM_photo module to process photolysis files and estimate J values
		J = lamp_photo(photo_par_file, J, TEMP, act_flux_path)
	
		# remove photolysis information folder
# 		cwd = os.getcwd() # address of current working directory
# 		if os.path.isdir(cwd + '/PyCHAM/MCMphotofiles/' + Photo_name[0]): # check if there is an existing Photo_name folder
# 			def handleRemoveReadonly(func, path, exc):
# 				excvalue = exc[1]
# 				if not os.access(path, os.W_OK):
# 					# Is the error an access error ?
# 					os.chmod(path, stat.S_IWUSR)
# 					func(path)
# 				else:
# 					raise
# 			shutil.rmtree(cwd + '/PyCHAM/MCMphotofiles/' + Photo_name[0], ignore_errors=False, onerror=handleRemoveReadonly) # remove existing folder, onerror will change permission of directory if needed

	
	return(J)
