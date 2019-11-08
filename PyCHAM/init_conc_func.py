'''function to initiate concentrations of components, obtain MCM reaction rate constants and produce the reaction coefficient file'''

import numpy as np
from water_calc import water_calc
import MCM_const
import eqn_parser
import scipy.constants as si

def init_conc_func(num_speci, init_SMIL, smiles_array, init_conc, TEMP, RH, M, N2, 
					O2, reac_coef, filename, PInit, time, lat, lon, Pybel_objects,
					testf, pconc):
		
	# -----------------------------------------------------------
	# inputs:
	
	# M - initial concentration of M (molecules/cc (air))
	# N2 - initial concentration of N2 (molecules/cc (air))
	# O2 - initial concentration of O2 (molecules/cc (air))
	# PInit - initial pressure (Pa)
	# init_SMIL - SMILES of components present at start of experiment (whose 
	# concentrations are given in init_conc)
	# testf - flag for whether in normal mode (0) or testing mode (1/2)
	# pconc - initial concentration of particles (# particles/cc (air))
	# -----------------------------------------------------------

	if testf==1: # testing mode
		# return dummies
		return(0,0,0,0,0,0,0,0)

	NA = si.Avogadro # Avogadro's number (molecules/mol)
	# empty array for storing species' concentrations, must be an array
	y = np.zeros((num_speci))
	y_mw = np.zeros((num_speci, 1)) # species' molecular weight (g/mol)
	# empty array for storing index of interesting gas-phase components
	y_indx_plot = []
	
	# convert concentrations
	# total number of molecules in 1 cc air using ideal gas law.  R has units cc.Pa/K.mol
	ntot = PInit*(NA/(8.3144598e6*TEMP))
	# number of molecules in one billionth of this
	Cfactor = ntot*1.0e-9 # ppb-to-molecules/cc

	# insert initial concentrations where appropriate (init_SMIL contains only the
	# species stated in the .def files)
	for i in range (len(init_SMIL)):

    	# index of where initial species occurs in SMILE string
		y_indx = smiles_array.index(init_SMIL[i])
		y[y_indx] = init_conc[i]*Cfactor
		# remember index for plotting gas-phase concentrations later
		y_indx_plot.append(y_indx) 
	
	for i in range(num_speci): # loop through all species
		y_mw[i] = Pybel_objects[i].molwt # molecular weight (g/mol)
	
	# account for water's properties
	H2Oi = num_speci # index for water
	num_speci += 1 # update number of species to account for water
	# append empty element to y and y_mw to hold water values
	y = np.append(y, 1.0e-40)
	y_mw = (np.append(y_mw, np.zeros((1)))).reshape(-1, 1)
	if pconc>0.0: # if seed particles present
		# append core gas-phase concentration (molecules/cc (air)) and molecular weight 
		# (g/mol) (needs to have a 1 length in second dimension for the kimt calculations)
		y = np.append(y, 1.0e-40) 
		y_mw = (np.append(y_mw, 132.14)).reshape(-1, 1)
		corei = num_speci # index of core component
		num_speci += 1 # update number of species to account for core material
	else:
		corei = 1e40 # dummy
	
	
	# concentration (molecules/cc (air))
	[y[H2Oi], Psat_water, y_mw[H2Oi]] = water_calc(TEMP, RH, NA)
	
	# obtain a list of MCM constants and their equations
	# get photolysis rates at initial time
	(mcm_constants, MCMConstNameList) = MCM_const.mcm_constants(TEMP, 
											y[H2Oi], M, N2, O2, time, lat, lon)
	
	# now create reaction rate file (reaction rates are set up to have units /s)
	eqn_parser.write_rate_file(filename, reac_coef, mcm_constants, MCMConstNameList, 
								testf, M, N2, O2)
	
	return (y, H2Oi, Psat_water, y_mw, num_speci, Cfactor, y_indx_plot, corei)