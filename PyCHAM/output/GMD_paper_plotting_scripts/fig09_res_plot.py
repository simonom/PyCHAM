'''module to test wallloss, for comparison with Fig. 2 on p. 254 of McMurry & Rader (1985), DOI: 10.1080/02786828508959054'''
print('function to test wallloss.py, can be called from the PyCHAM home directory or GMD paper Results folder')

def test_wallloss(Pn, Cn, Gi, eta_a, Dp, MW, Varr, sbn, nc, TEMP, t, 
			inflectDp, pwl_xpre, pwl_xpro, inflectk, ChamR, Rader, testf, p_char, 
			e_field, fig, ax0, num_asb): # define testing function

	import sys
	import os
	
	cwd = os.getcwd() # get current working directory
	import sys
	# ensure modules can be seen 
	# (assumes calling from the home folder, but will also work when in the GMD paper Results folder)
	sys.path.append(str(os.getcwd() + '/PyCHAM'))

	
	from wallloss import wallloss # call on wall loss module
	
	Beta = wallloss(Pn, Cn, Gi, eta_a, Dp, MW, Varr, sbn, nc, TEMP, t, 
			inflectDp, pwl_xpre, pwl_xpro, inflectk, ChamR, Rader, testf, p_char, 
			e_field, num_asb)
	
	if p_char == 0:
		
		ax0.loglog(Dp*1.0e6, Beta, 'k')
# 		ax0.set_ylim([8.0e-6, 1.0e0])
# 		ax0.set_title('Wall Loss Rate of Neutral and Charged Particles')
		ax0.set_xlim([1.0e-2, 1.0e1])
		ax0.set_xlabel(r'$D_p$ ($\rm{\mu}$m)', fontsize=12)
		ax0.set_ylabel(r'$\beta$ $(s^{-1})$', fontsize=12)
		ax0.text(0.5, 0.05, r'R = 35 cm')
		ax0.text(0.5, 0.03, r'E = 45 volts/cm')
		ax0.text(0.5, 0.018, r'$\rm{k_{e}=6.4x10^{-3}\, s^{-1}}$')
		ax0.text(0.025, 5e-5, 'n=0')
	if p_char == 1:
		ax0.loglog(Dp*1.0e6, Beta, 'k')
		ax0.text(0.03, 3e-3, 'n=1')
	if p_char == 2:
		ax0.loglog(Dp*1.0e6, Beta, 'k')
		ax0.text(0.05, 4e-3, '2')
	if p_char == 3:
		ax0.loglog(Dp*1.0e6, Beta, 'k')
		ax0.text(0.085, 5e-3, '3')
	if p_char == -1e6:
		ax0.loglog(Dp*1.0e6, Beta, '--k', label = r'$D_{p,flec}=1\, \mathrm{\mu m},\, \nabla_{pre}=1,\, \nabla_{pro}=1,\, \beta_{flec}=\mathrm{2x10^{-5} s^{-1}}$')
	if p_char == -2e6:
		ax0.loglog(Dp*1.0e6, Beta, '-.k', label = r'$D_{p,flec}=0.8\, \mathrm{\mu m},\, \nabla_{pre}=0,\, \nabla_{pro}=-1,\, \beta_{flec}=\mathrm{1x10^{-3} s^{-1}}$')
		ax0.legend()
		

# ----------------------------------------------------------------------------------------
# inputs

import os
import sys
import numpy as np
import scipy.constants as si
import matplotlib.pyplot as plt

import sys
# ensure modules can be seen 
# (assumes calling from the home folder, but will also work when in the GMD paper Results folder)
sys.path.append(str(os.getcwd() + '/PyCHAM'))

from fl_reg_determ import reg_determ

# number of particle size bins, +1 to replicate wall
sbn = 100+1
# particle number concentration per size bin before time step (# particle/cc (air))
Pn  = np.ones((sbn-1, 1))
# relative humidity
RH = 0.5
TEMP = 298.15 # temperature (K)
sbr = np.logspace(-8.4, -5.7, num=(sbn-1), endpoint=True, base=10.0) # particle radius (m)
PInit = 1.0e5 # pressure inside chamber (Pa)
# Knudsen number (dimensionless) and dynamic viscosity ((eta_a) eq. 4.54 Jacobson 2005) 
# (g/m.s)
[Kn, eta_a, rho_ai, kin_visc] = reg_determ(RH, TEMP, sbr, PInit)
# particle diameters (m)
Dp = sbr*2.0
# Cunningham slip-flow correction (eq. 15.30 of Jacobson (2005)) with constant taken
# from text below textbook equation (dimensionless)
Gi = 1.0+Kn*(1.249+0.42*(np.exp(-0.87/Kn)))
# volume of single particles per size bin (m3)
Varr = (4.0/3.0)*np.pi*sbr**3.0
nc = 1 # number of components
MW = np.ones((1, nc))*100.0 # component molecular weight (g/mol)
# particle phase concentration per component per size bin (molecules/cc (air))
Cn = Varr*(1.0e6)*(1.0/MW)*si.Avogadro
# time that wall loss occurs over (s)
t = 60.0
# dummy values for manual setting of wall loss equation
inflectDp = 1.0
pwl_xpre = 1.0
pwl_xpro = 1.0 
inflectk = 1.0
# spherical equivalent radius of chamber (Fig. 2 legend) (m)
ChamR = 0.35
Rader = 1 # flag to use the McMurry and Rader (1985) model
testf = 1
p_char = 0 # average number of charges per particle
e_field = 4500000 # average electric field in chamber, (g.m/A.s3) (note: V/m == kg.m/A.s3) 

fig, (ax0) = plt.subplots(1, 1, figsize=(6,7.0))

test_wallloss(Pn, Cn, Gi, eta_a, Dp, MW, Varr, sbn, nc, TEMP, t, 
			inflectDp, pwl_xpre, pwl_xpro, inflectk, ChamR, Rader, testf, p_char, 
			e_field, fig, ax0, sbn) # call testing function

p_char = 1 # average number of charges per particle

test_wallloss(Pn, Cn, Gi, eta_a, Dp, MW, Varr, sbn, nc, TEMP, t, 
			inflectDp, pwl_xpre, pwl_xpro, inflectk, ChamR, Rader, testf, p_char, 
			e_field, fig, ax0, sbn) # call testing function
	
p_char = 2 # average number of charges per particle

test_wallloss(Pn, Cn, Gi, eta_a, Dp, MW, Varr, sbn, nc, TEMP, t, 
			inflectDp, pwl_xpre, pwl_xpro, inflectk, ChamR, Rader, testf, p_char, 
			e_field, fig, ax0, sbn) # call testing function

p_char = 3 # average number of charges per particle

test_wallloss(Pn, Cn, Gi, eta_a, Dp, MW, Varr, sbn, nc, TEMP, t, 
			inflectDp, pwl_xpre, pwl_xpro, inflectk, ChamR, Rader, testf, p_char, 
			e_field, fig, ax0, sbn) # call testing function
			
# ----------------------------------------------------------------------------------------
# now for the customised example cases
inflectDp = 1.0e-6 # diameter at which inflection in deposition rate occurs (m)
pwl_xpre = 1.0
pwl_xpro = 1.0 
inflectk = 2.0e-5
p_char = -1e6
Rader = 0 # flag to use the McMurry and Rader (1985) model
test_wallloss(Pn, Cn, Gi, eta_a, Dp, MW, Varr, sbn, nc, TEMP, t, 
			inflectDp, pwl_xpre, pwl_xpro, inflectk, ChamR, Rader, testf, p_char, 
			e_field, fig, ax0, sbn) # call testing function

inflectDp = 8.0e-7 # diameter at which inflection in deposition rate occurs
pwl_xpre = 0.0
pwl_xpro = -1.0 
inflectk = 1.0e-3
p_char = -2e6
Rader = 0 # flag to use the McMurry and Rader (1985) model
test_wallloss(Pn, Cn, Gi, eta_a, Dp, MW, Varr, sbn, nc, TEMP, t, 
			inflectDp, pwl_xpre, pwl_xpro, inflectk, ChamR, Rader, testf, p_char, 
			e_field, fig, ax0, sbn) # call testing function

plt.show() # show figure
fig.savefig('fig09.png')
