'''module to unit test the coag.py module, by producing plot of coagulation kernel as function of particle diameter, for comparison with Fig. 15.7 on p. 512 of Jacobson (2005)'''
print('function to test coag.py, please call once in the Unit_Testing folder')
import os
import sys
import scipy.constants as si
import numpy as np

dirpath = os.getcwd() # get current path
sys.path.append(os.path.split(dirpath)[0]) # add path to system path


# define function
def test_coag(RH, T, sbr, sbVi, M, rint, num_molec, num_part, tint, sbbound,
			num_comp, vdWon, rho, rad0, PInit, testf, num_molec_rint, num_part_rint,
			sbVj):

	from coag import coag # import coagulation module
	
	# call on coagulation module
	coag(RH, T, sbr, sbVi, M, rint, num_molec, num_part, tint, sbbound,
			num_comp, vdWon, rho, rad0, PInit, testf, num_molec_rint, num_part_rint, 
			sbVj)


# ----------------------------------------------------------------------------------------
# state inputs

RH = 0.50 # relative humidity (fraction (0-1))
T = 298.0 # temperature (K)
sbr = np.logspace(-8.0, -4.6, num=50, endpoint=True, base=10.0) # size bin radius (m)
sbVi = (4.0/3.0)*np.pi*sbr**3.0 # single particle volume for i sizes (m3)
M = 200.0 # molecular weight of components (g/mol) 
rho = 1.0 # components densities (g/cm3) in a 1D array
rint = np.array((1.0e-8, 1.0e-5)) # size of interest (m)
# concentration of particles per size bin (sbr) (particle/cc (air))
num_part = np.ones((1, len(sbr)))
# molecular concentration (molecules/cc (air)), arranged by component in rows and size 
# bins (sbr) in columns
num_molec = ((rho*1.0e6)*((1.0/M)*si.N_A))*sbVi*num_part
# concentration of particles per size bin (sbr) (particle/cc (air))
num_part_rint = np.ones((1, len(rint)))
# molecular concentration (molecules/cc (air)), arranged by component in rows and size 
# bins (sbr) in columns
num_molec_rint = ((rho*1.0e6)*((1.0/M)*si.N_A))*((4.0/3.0)*np.pi*rint**3.0)*num_part_rint
tint = 1.0 # time interval coagulation occurs over (s)
sbbound = sbVi[0:-1]+(sbVi[1::]-sbVi[0:-1])/2.0 # size bin volume boundaries (m3)
sbbound = np.append(0.0, sbbound)
sbbound = np.append(sbbound, sbVi[-1]*2.0)
num_comp = 1.0 # number of components
vdWon = 0 # saying whether the van der Waals kernel should be calculated or ignored
rad0 = sbr # original radius at size bin centre (um)
PInit = 1.0e5 # pressure inside chamber (Pa)
testf = 1 # unit testing flag on
sbVj = (4.0/3.0)*np.pi*rint**3.0 # single particle volume for j sizes (m3)

# call on test function
test_coag(RH, T, sbr, sbVi, M, rint, num_molec, num_part, tint, sbbound,
			num_comp, vdWon, rho, rad0, PInit, testf, num_molec_rint, num_part_rint, 
			sbVj)