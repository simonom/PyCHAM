'''test for the Kelvin and Raoult effects, for comparison with Fig. 16.1 on p. 537 of Jacobson (2005)'''
print('function to test the Kelvin and Raoult effects, please call from the PyCHAM home folder')

import sys
import os
# ensure modules can be seen 
# (assumes calling function is in the home folder)
sys.path.append(str(os.getcwd() + '/PyCHAM'))
import write_ode_solv

def test_partit_var(): # define testing function

	# import functions
	import partit_var
	import ode_solv

	# get Kelvin factor (curvature effect)
	# call on function to calculate Kelvin factor
	[kimt, kelv_fac] = partit_var.kimt_calc(y, mfp, num_sb, num_comp, accom_coeff, y_mw,   
		surfT, R_gas, TEMP, NA, y_dens, N_perbin, DStar_org, 
		x.reshape(1, -1)*1.0e-6, Psat, therm_sp, H2Oi, act_coeff, wall_on, 1)
	
	# call on ode_gen to get solute (Raoult) effect (mole fraction of water in particles)
	mf = ode_solv.ode_solv(y, 1., [], [], [], [],
		[], [], [], [], [], [], [],
		[], [], [], [], [], [], 
		[], [], [], 
		[], [], [], [], [], num_comp, 
		num_sb, wall_on, Psat, [], act_coeff, [], [],
		(num_comp-1), core_diss, kelv_fac, kimt, (num_sb-wall_on), 
		[],
		[], [], [], [],
		[], [], [], [], [], [], 
		[], [], [], [], [], 
		[], [], [], 
		[], [], [], [])

	
	fig, (ax0) = plt.subplots(1, 1, figsize=(6,5))
	ax0.semilogx(x, kelv_fac, label='Curvature (Kelvin) effect')
	ax0.semilogx(x, mf[:, 0].reshape(-1, 1), label='Solute (Raoult) effect')
	# note that the Raoult (solute) effect is calculated differently to eq. 16.39 of
	# Jacobson, the PyCHAM calculation is the mole fraction.  Following eq. 3 of 
	# Riipinen et al. (2010): doi:10.1016/j.atmosenv.2009.11.022, we get the equilibrium
	# saturation ratio by the product of the two effects
	# Whereas the exact inputs used for Fig. 16.1 of Jacobson (2005) are unknown, the 
	# plot here should have similar features despite probably not giving an exact
	# replication
	ax0.semilogx(x, (kelv_fac*(mf[:, 0].reshape(-1, 1))), '--', label='Equilibrium saturation ratio')
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
num_sb = 100 # number of size bins (excluding wall)
num_comp = 2 # number of components (including water and core)
accom_coeff = np.ones((num_comp, 1)) # accommodation coefficient
act_coeff = np.ones((num_comp, 1)) # activity coefficient
y_mw = np.array((18.015, 100.0)).reshape(-1, 1) # molecular weight (g/mol)
surfT = 72.0 # surface tension (assume surface tension of water (g/s2==mN/m==dyn/cm))
R_gas = si.R # ideal gas constant (kg.m2.s-2.K-1.mol-1)
TEMP = 298.15 # temperature (K)
NA = si.Avogadro # Avogadro's number (molecules/mol)
y_dens = np.array((1.0e3, 1.0e3)).reshape(-1, 1) # component density (kg/m3)
N_perbin = np.ones((num_sb, 1)) # number of particles per size bin (#particles/cc (air))
DStar_org = np.ones((2, 1)) # dummy gas-phase diffusion coefficient of components (m2/s)
x = np.logspace(-2.0, 1.0, num=(num_sb), endpoint=True, base=10.0) # particle radius (um)
Vsol = ((4.0/3.0)*np.pi*x[0]**3.0)
Vall = (4.0/3.0)*np.pi*x**3.0
# concentrations of solute and solvent in particle-phase
# assume same volume of solute in each
y = np.zeros((num_comp*(num_sb)))
for i in range(num_sb):
	y[i*num_comp+1] = Vsol/Vall[i] # solute concentration
	y[i*num_comp] = (Vall[i]-Vsol)/Vall[i] # solvent concentration
y = np.append(np.array((0.0, 0.0)), y) # append dummy gas-phase concentration

core_diss = 4.0 # dissociation constant

Psat = np.ones((2)) # dummy saturation vapour pressures (molecules/cc (air))
therm_sp = np.ones((2)) # dummy thermal speed
H2Oi = 0 # index of water
wall_on = 0 # forget about wall

# call function to generate ordinary differential equation (ODE)
# solver module, add two to comp_num to account for water and seed material
write_ode_solv.ode_gen([], [], [], wall_on, num_comp+2, 
		(num_sb-wall_on), [], 2)

# call testing function with inputs
test_partit_var()
