'''module to test kimt_calc, for comparison with Fig. 16.1 on p. 537 of Jacobson (2005)'''
print('function to test kimt_calc.py and ode_gen.py, please call once in the Unit_Testing folder')

def test_kimt_calc(): # define testing function

	# get Kelvin factor (curvature effect)
	from kimt_calc import kimt_calc # import function
	# call on function to calculate Kelvin factor
	[kimt, kelv_fac] = kimt_calc(y, mfp, num_sb, num_speci, accom_coeff, y_mw,   
							surfT, R_gas, TEMP, NA, y_dens, N_perbin, DStar_org, 
							x.reshape(1, -1)*1.0e-6, Psat, therm_sp, H2Oi)
	
	from ode_gen import ode_gen # import function						
	# call on ode_gen to get solute effect
	sol_eff = ode_gen(1, y, num_speci, 1, 1, 1, 1, 1, H2Oi, 
			TEMP, 1, num_sb, 
			1, mfp, accom_coeff, surfT, y_dens, 
			N_perbin, DStar_org, y_mw, x, core_diss, 1, 1, 1, 1, 
			1, 1, pconc, 
			1, 1, therm_sp,
			1, 1, 1, 1, 
			1, 1, 1, 1, 1, 1, 1, 1, 1, 
			1, 1, 1, 1, 1, 1, 1, 2, 1,
			1, 1, 1, 1, 1, 1)
	
	fig, (ax0) = plt.subplots(1, 1, figsize=(6,5))
	ax0.semilogx(x, kelv_fac, label='Curvature effect')
	ax0.semilogx(x, sol_eff, label='Solute effect')
	# note that the Raoult (soulte) effect is calculated differently to eq. 16.39 of
	# Jacobson, the PyCHAM calculation is the mole fraction.  Following eq. 3 of 
	# Riipinen et al. (2010): doi:10.1016/j.atmosenv.2009.11.022, we get the equilibrium
	# saturation ratio by the product of the two effects
	# Whereas the exact inputs used for Fig. 16.1 of Jacobson (2005) are unknown, the 
	# plot here should have similar features despite probably not giving an exact
	# replication
	ax0.semilogx(x, (kelv_fac*sol_eff), '--', label='Equilibrium saturation ratio')
	ax0.set_ylim([0.98, 1.02e0])
	ax0.set_xlabel(r'Particle radius ($\rm{\mu}$m)', fontsize=10)
	ax0.set_ylabel(r'Saturation ratio', fontsize=10)
	ax0.legend()
	plt.show()			
	
# ----------------------------------------------------------------------------------------
# inputs

import os
import sys
import scipy.constants as si
import numpy as np
import matplotlib.pyplot as plt

dirpath = os.getcwd() # get current path
sys.path.append(os.path.split(dirpath)[0]) # add path to system path

pconc = 1 # whether core material present as seed particles (1 for yes)
mfp = np.ones((2, 1)) # components' mean free path (dummy values)
num_sb = 101 # number of size bins (including wall)
num_speci = 2 # number of components
accom_coeff = np.ones((num_speci, 1)) # accommodation coefficient (dummy value)
y_mw = np.array((18.015, 100.0)).reshape(-1, 1) # molecular weight (g/mol)
surfT = 72.0 # surface tension (assume surface tension of water (g/s2==mN/m==dyn/cm))
R_gas = si.R # ideal gas constant (kg.m2.s-2.K-1.mol-1)
TEMP = 298.15 # temperature (K)
NA = si.Avogadro # Avogadro's number (molecules/mol)
y_dens = np.array((1.0e3, 1.0e3)).reshape(-1, 1) # component density (kg/m3)
N_perbin = np.ones((num_sb-1)) # number of particles per size bin (#particles/cc (air))
DStar_org = np.ones((2, 1)) # dummy gas-phase diffusion coefficient of components (m2/s)
x = np.logspace(-2.0, 1.0, num=(num_sb-1), endpoint=True, base=10.0) # particle radius (um)
Vsol = ((4.0/3.0)*np.pi*x[0]**3.0)
Vall = (4.0/3.0)*np.pi*x**3.0
# concentrations of solute and solvent in gas-, particle- and wall-phase (ten size bins)
# assume same volume of solute in each
y = np.zeros((num_speci*(num_sb-1)))
for i in range(num_sb-1):
	y[i*num_speci+1] = Vsol/Vall[i]
	y[i*num_speci] = (Vall[i]-Vsol)/Vall[i]
y = np.append(np.array((0.0, 0.0)), y) # append dummy gas-phase concentration
y = np.append(y, np.array((0.0, 0.0))) # append dummy wall-phase concentration
core_diss = 4.0 # dissociation constant

Psat = np.ones((2)) # dummy saturation vapour pressures (molecules/cc (air))
therm_sp = np.ones((2)) # dummy thermal speed
H2Oi = 0 # index of water

# call testing function with inputs
test_kimt_calc()