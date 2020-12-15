'''module to estimate coagulation kernel using p. 508 of Jacobson (2005), called on by ode_gen'''

# number- but not mass-conserving coagulation code, called from ode_gen and operator-split
# from ode solver

import numpy as np 
import scipy.constants as si
from fl_reg_determ import reg_determ
from reyn_num import Reyn_num
from wk_int import W_k_int
import scipy.integrate as integ
import matplotlib.pyplot as plt
import scipy.constants as si
from mov_cen_water_eq import mov_cen_main as movcen # moving centre method for rebinning
import fullmov

def coag(RH, T, sbr, sbVi, M, rint, num_molec, num_part, tint, sbbound, rbou,
			num_comp, vdWon, rho, V0, rad0, PInit, testf, num_molec_rint, num_part_rint, 
			sbVj, coag_on, siz_str, wall_on):

	# inputs:---------------------------------------------------------
	
	# RH - relative humidity (fraction)
	# T - temperature (K)
	# sbr - size bin radius (m)
	# sbVi - single particle volume for i sizes (relating to sbr) (m3)
	# M - molecular weight of components (g/mol)	
	# rint - size(s) of interest (m)
	# num_molec - molecular concentration (molecules/cc (air)),
	# arranged by component in rows and size bins in columns, for particles in sbr
	# num_part - concentration of particles per size bin 
	# (particle/cc (air)) (columns) (excluding walls), for particles in sbr
	# tint - time interval coagulation occurs over (s)
	# sbbound - size bin volume boundaries (m3)
	# rbou - size bin radius boundaries (um)
	# num_comp - number of components
	# vdWon - flagging whether the van der Waals kernel should be calculated or ignored (0
	# for ignore, 1 for calculate)
	# rho - component densities (g/cm3) in a 1D array
	# V0 - original volume at size bin centre (um3)
	# rad0 - original radius at size bin centre (um)
	# PInit - pressure inside chamber (Pa)
	# testf - unit testing flag (0 for off, 1 for kernel testing, 2 for coagulation 
	#	testing with fixed kernel)
	# num_molec_rint - molecular concentration (molecules/cc (air)),
	# arranged by component in rows and size bins in columns, for particles in rint
	# num_part_rint - concentration of particles per size bin 
	# (particle/cc (air)) (columns) (excluding walls), for particles in rint
	# sbVj -  - single particle volume for j sizes (relating to rint) (m3)
	# coag_on - whether to allow coagulation to occur (1) or not (0)
	# siz_str - the size structure to use
	# wall_on - marker for whether wall being considered
	# --------------------------------------------------------------
	
	num_part = num_part.reshape(1, -1)
	# volume concentration of particles (m3/cc (air))
	vol_part = (sbVi*num_part).reshape(1, -1)

	# ensure sbn is integer
	sbrn = np.int(np.max(sbr.shape))
	sbn = np.int(np.max(rint.shape)) # number of size bins
	
	# call on function to determine the Knudsen no. and therefore flow 
	# regime of each size bin
	[Kni, eta_ai, rho_ai, kin_visc] = reg_determ(RH, T, sbr, PInit)
	[Knj, eta_aj, rho_aj, kin_visc] = reg_determ(RH, T, rint, PInit)

	# Reynold number and terminal fall velocity for each size bin
	[Rei, Vfi] = Reyn_num(sbr, eta_ai, rho_ai, kin_visc, 1.0e6, Kni)
	[Rej, Vfj] = Reyn_num(rint, eta_aj, rho_aj, kin_visc, 1.0e6, Knj)
	
	# repeat Knudsen number over number of size bins
	Kni_m = np.tile(Kni, (sbn, 1))
	
	# Cunningham slip-flow correction (15.30) with constant taken
	# from text below textbook equation (dimensionless)
	Gi = 1.0+Kni*(1.249+0.42*(np.exp(-0.87/Kni)))
	Gj = 1.0+Knj*(1.249+0.42*(np.exp(-0.87/Knj)))
	# particle diffusion coefficient (15.29) (m2/s)
	# multiply eta_a by 1.0e-3 to convert from g/m.s to kg/m.s
	# this makes it consistent with the units of Boltzmann constant
	Dpi = (((si.k*T)/(6.0*np.pi*sbr*(eta_ai*1.0e-3)))*Gi).reshape(sbrn, 1)
	Dpj = (((si.k*T)/(6.0*np.pi*rint*(eta_aj*1.0e-3)))*Gj).reshape(1, sbn)
	
	# repeat Dpi over size bins of j (m2/s)
	Dp_mi = np.repeat(Dpi, sbn, 1)
	# repeat Dpj over size bins of i (m2/s)
	Dp_mj = np.repeat(Dpj, sbrn, 0)
	# matrix for sums of Dp in size bins (m2/s)
	Dp_sum = Dp_mi+Dp_mj

	# zero the upper triangle part of the matrix as only the lower is needed
# 	Dp_sum[np.triu_indices(sbrn, 1, m=sbn)] = 0.0
	
	# Brownian collision kernel (m3/particle.s) (p. 508) 
	K_B = np.zeros((sbrn, sbn))

	# ensure sbr is a column array and rint is a row array
	sbr2 = np.zeros((np.max(sbr.shape), 1))
	rint2 = np.zeros((1, np.max(rint.shape)))
	sbr2[:,0] = sbr[:]
	rint2[0,:] = rint[:]
	
	# repeat size bin radii over size bins (m)
	rint2 = np.repeat(rint2, sbrn, 0)
	sbr_m = np.repeat(sbr2, sbn, 1)
	# zero the upper triangle part of the matrix as only the lower is needed
# 	rint2[np.triu_indices(sbrn, 1, m=sbn)] = 0.0
# 	sbr_m[np.triu_indices(sbrn, 1, m=sbn)] = 0.0
	
	# matrix for size bins sums (m)
	sbr_sum = sbr_m+rint2
	
	
	# size bins in continuum regime
	i = ((Kni<1.0).reshape(sbrn,1))
	# spread across rint (j)
	i = np.repeat(i, sbn, 1)
	
	# Brownian collision kernel (15.28) (m3/particle.s) in the continuum regime
	K_B[i] = 4.0*np.pi*(sbr_sum[i])*(Dp_sum[i])
	
	# size bins in free-molecular regime	
	i = ((Kni>10.0).reshape(sbrn,1))
	# spread across rint (j)
	i = np.repeat(i, sbn, 1)
	
	# single particle mass (g):
	# first, number of moles per component in a single particle (relating to sbr)
	num_mol_single_sbr = np.zeros((num_molec.shape[0], (num_part.shape[1])))
	ish = num_part[0, :]>0.
	num_mol_single_sbr[:, ish] = (num_molec[:, ish]/num_part[0, ish])/si.N_A
	# second product of number of moles and molecular weight
	weight_compon = num_mol_single_sbr*M
	# final sum mass across components for single particle mass (g)
	Mpi = np.sum(weight_compon, 0)
	# number of moles per component in a single particle (relating to rint)
	num_mol_single_rint = np.zeros((num_molec_rint.shape[0], (num_part_rint.shape[1])))
	ish = num_part_rint[0, :]>0.	
	num_mol_single_rint[:, ish] = (num_molec_rint[:, ish]/num_part_rint[0, ish])/si.N_A
	# second product of number of moles and molecular weight
	weight_compon = num_mol_single_rint*M
	Mpj = np.sum(weight_compon, 0)

	# thermal speed of particle (15.32) (m/s) (multiply mass by 1.0e-3 to 
	# convert from g to kg and therefore be consistent with Boltzmann's 
	# constant (1.380658e-23kgm2/s2.K.molec))
	Mpi = np.pi*(Mpi*1.0e-3).reshape(sbrn, 1).repeat(sbn, 1)
	nu_pi = np.zeros((sbrn, sbn))
	nu_pi[Mpi>0.] = ((8.0*si.k*T)/Mpi[Mpi>0.])**0.5
	
	Mpj = np.pi*(Mpj*1.0e-3).reshape(1, sbn).repeat(sbrn, 0)
	nu_pj = np.zeros((sbrn, sbn))
	nu_pj[Mpj>0.] = ((8.0*si.k*T)/Mpj[Mpj>0.])**0.5	

	# sum of squares of speeds (m2) (15.31)
	nu_p_sum = nu_pi**2.0+nu_pj**2.0
			
	# Brownian collision kernel (15.31) (m3/particle.s)
	K_B[i] = np.pi*(sbr_sum[i]**2.0)*(nu_p_sum[i]**0.5)
	
	# size bins in transition regime, where the transition regime corresponds to particles
	# of size i, not size j
	i = (1.0<=Kni) 
	j = (Kni<=10.0) # size bins in transition regime		
	i = (i*j).reshape(sbrn, 1)	# size bins in transition regime
	# spread across rint (j)
	i = np.repeat(i, sbn, 1)
	
	# particle mean free path (15.34) (m)
	lam_pi = np.zeros((len(nu_pi[:,0])))
	ish = nu_pi[:, 0]>0.
	lam_pi[ish] = ((8.0*Dpi[ish, 0])/(np.pi*nu_pi[ish, 0]))
	lam_pj = np.zeros((len(nu_pj[0,:])))
	ish = nu_pj[0, :]>0.
	lam_pj[ish] = ((8.0*Dp_mj[0, ish])/(np.pi*nu_pj[0, ish]))

	# mean distance from centre of a sphere travelled by particles
	# leaving sphere's surface and travelling lam_p (m) (15.34)
	num = (2.0*sbr+lam_pi)**3.0-(4.0*sbr**2.0+lam_pi**2.0)**1.5
	den = (6.0*sbr*lam_pi)-2.0*sbr
	sig_pi = np.zeros((len(sbr), 1))
	
	ish = den > 0.
	
	sig_pi[ish, 0] = num[ish]/den[ish]	
	sig_pi = sig_pi.repeat(sbn, 1)
	
	num = (2.0*rint+lam_pj)**3.0-(4.0*rint**2.0+lam_pj**2.0)**1.5
	den = (6.0*rint*lam_pj)-2.0*rint
	sig_pj = np.zeros((1, len(rint)))
	ish = den>0.
	sig_pj[0, ish] = num[ish]/den[ish]
	sig_pj = sig_pj.repeat(sbrn, 0)

	# sum mean distances (m)
	sig_p_sum = sig_pi**2.0+sig_pj**2.0
		
	# kernel numerator
	K_Bnum = 4.0*np.pi*sbr_sum[i]*Dp_sum[i]

	# left term kernel denominator
	K_Blden = (sbr_sum[i]/(sbr_sum[i]+(sig_p_sum[i])**0.5))
	
	# right term kernel denominator
	num = 4.0*Dp_sum[i]
	den = ((nu_p_sum[i])**0.5)*sbr_sum[i]
	K_Brden = np.zeros((len(num)))
	ish = den>0.
	K_Brden[ish] = num[ish]/den[ish] 
	
	# collision kernel (15.33) (m3/particle.s)
	K_B[i] = (K_Bnum/(K_Blden+K_Brden))
	

	
	if (testf == 1): # testing flag on
		fig, (ax0,ax1) = plt.subplots(1, 2, figsize=(12,6))
		ax0.loglog(sbr*10**6, K_B[:,0]*10**6, label='Brownian')
		ax0.set_xlabel(r'Radius of second particle ($\rm{\mu}$m)', fontsize=10)
		ax0.set_ylabel(r'Coagulation kernel ($\rm{cm^{3}particle^{-1}s^{-1}}$)', fontsize=10)
		ax1.loglog(sbr*10**6, K_B[:,1]*10**6, label='Brownian')
		ax1.set_xlabel(r'Radius of second particle ($\rm{\mu}$m)', fontsize=10)
		ax1.set_ylabel(r'Coagulation kernel ($\rm{cm^{3}particle^{-1}s^{-1}}$)', fontsize=10)
		
	
	# Convective Brownian Diffusion Enhancement kernel:

	
	# particle Schmidt number for i and jsize bins (dimensionless) (15.36)
	Scpi = kin_visc/Dpi
	Scpj = kin_visc/Dpj

	# repeat Rej over sbrn and Rei over sbn
	Re_jm = np.repeat(Rej.reshape(1,sbn), sbrn, 0)
	Re_im = np.repeat(Rei.reshape(sbrn,1), sbn, 1)
	Scpj_m = np.repeat(Scpj.reshape(1,sbn), sbrn, 0)
	Scpi_m = np.repeat(Scpi.reshape(sbrn,1), sbn, 1)

	# repeat Reynold number and Schmidt number arrays over size bins
	i = (Re_jm<=1.0)
	j = rint2>=sbr2
	i2 = i*j # Rej less than or equal to 1, and rj greater than or equal to ri
	i3 = (Re_im<=1.0)
	j2 = rint2<sbr2
	j3 = i3*j2 # Rei less than or equal to 1 and ri greater than rj
	
	i = (Re_jm>1.0)
	i3 = i*j # Rej greater than 1 and rj greater than or equal to ri
	i4 = (Re_im>1.0)
	j4 = i4*j2 # Rei greater than 1 and ri greater than rj
	

	# convective Brownian diffusion enhancement kernel (15.35)
	K_DE = np.zeros((sbrn, sbn))
	
	# convective Brownian diffusion enhancement kernel (15.35)
	# condition for both K_DE equation is r_j>=r_i
	K_DE[i2] = (K_B[i2]*0.45*Re_jm[i2]**(1.0/3.0)*Scpi_m[i2]**(1.0/3.0))
	K_DE[j3] = (K_B[j3]*0.45*Re_im[j3]**(1.0/3.0)*Scpj_m[j3]**(1.0/3.0))
	K_DE[i3] = (K_B[i3]*0.45*Re_jm[i3]**(1.0/2.0)*Scpi_m[i3]**(1.0/3.0))
	K_DE[j4] = (K_B[j4]*0.45*Re_im[j4]**(1.0/2.0)*Scpj_m[j4]**(1.0/3.0))
	
	if testf==1:
		ax0.loglog(sbr*10**6, K_DE[:,0]*10**6, label='Diff. Enhancement')
		ax1.loglog(sbr*10**6, K_DE[:,1]*10**6, label='Diff. Enhancement')
	

	
	# Gravitational Collection Kernel:
	Ecoll = np.zeros((sbrn, sbn))
	j = rint2>=sbr2
	Ecoll[j] = ((sbr2.repeat(sbn, 1))[j])**2.0/((sbr_sum[j])**2.0)
	j = rint2<sbr2
	Ecoll[j] = (rint2[j])**2.0/((sbr_sum[j])**2.0)
	
	# Gravitational collection kernel (15.37)
	# difference in terminal fall velocities
	del_Vf = np.abs(np.repeat(Vfj.reshape(1, sbn), sbrn, 0)-Vfi.reshape(sbrn,1))
	K_GC = Ecoll*np.pi*((sbr_sum)**2.0)*del_Vf
	
	if testf==1:
		ax0.loglog(sbr*10**6, K_GC[:,0]*10**6, label='Settling')
		ax1.loglog(sbr*10**6, K_GC[:,1]*10**6, label='Settling')
	
	
	# Kernel for Turbulent Inertial Motion:

	# rate of dissipation of turbulent kinetic energy per gram of medium 
	# (m2/s3) (8.4 for a typical value (which is apparently taken 
	# from Pruppacher and Klett 1997 (p. 511 Jacobson (2005)))
	epsilon	= 5.0e-4
	# kernel for turbulent inertial motion (15.40)
	K_TI = (((np.pi*epsilon**(3.0/4.0))/(si.g*kin_visc**(1.0/4.0)))*(
		sbr_sum**2.0)*del_Vf)
		
	if testf==1:
		ax0.loglog(sbr*10**6, K_TI[:,0]*10**6, label='Turb. inertia')
		ax1.loglog(sbr*10**6, K_TI[:,1]*10**6, label='Turb. inertia')
	
	# kernel for Turbulent Shear (15.41)
	K_TS = ((8.0*np.pi*epsilon)/(15.0*kin_visc))**0.5*(sbr_sum**3.0)
	
	if testf==1:
		ax0.loglog(sbr*10**6, K_TS[:,0]*10**6, label='Turb. shear')
		ax1.loglog(sbr*10**6, K_TS[:,1]*10**6, label='Turb. shear')
	
	# -----------------------------------------------------------------
	# Van der Waals/viscous collision kernel:	
	
	# particle radius (m)
	radi = sbr2 # radii array (m)
	
	# van der Waals factor results array and particle pair 
	# Knudsen number array: radii of particle i
	# in 1st dim., radii of particle j in 2nd
	# ratios of second particle radius to first particle - this 
	# approach commented out but useful for comparing with Fig. 15.8 
	# Jacobson (2005)
	res_all = np.zeros((radi.shape[0], rint.shape[0]))
	res_Knp = np.zeros((radi.shape[0], rint.shape[0]))
	
	A_H = 200.0*(si.k*1.0e3*T) # Hamaker constant
	
	for j in range(0, rint.shape[0]):
		if vdWon == 0: # when omitting van der Waals calculation for expediency
			break 
		for i in range(0,radi.shape[0]):
	
			ri = radi[i]
			rj = rint[j]
			a = (-A_H/(6.0*si.k*1.0e3*T))*2.0*ri*rj
			a1 = (-A_H/6.0)*2.0*ri*rj
			b = (-A_H/(6.0*si.k*1.0e3*T))
			b1 = (-A_H/6.0)
			# square of sum of size bin radii (m2)
			c = (ri+rj)**2.0
			# square of difference in size bin radii (m2)
			d = (ri-rj)**2.0
			# product of radii (m2)
			e = ri*rj
			# sum of radii (m)
			f = ri+rj
			# difference of radii (m)
			g = ri-rj
	
			# define the integration in 15.44
			def integrand(x,a,b,c,d,e,f,g):
				Dterm = (1.0+((2.6*e)/(c))*((e/(f*(x-g)))**0.5)+e/(f*(x-g)))
				Ep0_1 = a/(x**2.0-c)
				Ep0_2 = a/(x**2.0-d)
				Ep0_3 = b*np.log((x**2.0-c)/(x**2.0-d))
				rterm = 1.0/(x**2.0)
				return Dterm*np.exp(Ep0_1+Ep0_2+Ep0_3)*rterm
			# define the integration in 15.43
			def integrand2(x,a1,b1,c,d,e,f,g,T):
				# terms of Ep0
				Ep0_1 = a1/(x**2.0-c)
				Ep0_2 = a1/(x**2.0-d)
				Ep0_3 = b1*np.log((x**2.0-c)/(x**2.0-d))
				Ep0 = Ep0_1+Ep0_2+Ep0_3
				# terms of first differential Ep0
				Ep1_1 = (-2.0*a1*x)/((x**2.0-c)**2.0)
				Ep1_2 = (-2.0*a1*x)/((x**2.0-d)**2.0)
				Ep1_3 = (2.0*b1*x)/(x**2.0-c)
				Ep1_4 = (-2.0*b1*x)/(x**2.0-d)
				Ep1 = Ep1_1+Ep1_2+Ep1_3+Ep1_4
				# terms of second differential Ep0
				Ep2_1 = (6.0*a1*x**4.0-4.0*a1*c*x**2.0-2.0*a1*c**2.0)/((x**2.0-c)**4.0)	
				Ep2_2 = (6.0*a1*x**4.0-4.0*a1*d*x**2.0-2.0*a1*d**2.0)/((x**2.0-d)**4.0)
				Ep2_3 = (-2.0*b1*x**2.0-2.0*b1*c)/((x**2.0-c)**2.0)
				Ep2_4 = (2.0*b1*x**2.0+2.0*b1*d)/((x**2.0-d)**2.0)
				Ep2 = Ep2_1+Ep2_2+Ep2_3+Ep2_4
		
				return (Ep1+x*Ep2)*np.exp((-1.0/(si.k*1.0e3*T))*((x/2.0)*Ep1+Ep0))*(x**2.0)
	
			# integration bounds - note both integral functions
			# fall to negligible values after (ri+rj)*1.0e2 and if
			# infinity used as the upper bound numerical issues 
			# arise, therefore use (ri+rj)*1.0e2 for upper bound
			ilu = (ri+rj)*1.0e2 # upper
			ill = (ri+rj) # lower
			# integration in 15.44
			res = integ.quad(integrand, ill, ilu,args=(a,b,c,d,e,f,g),
				points=([ill*2.0]), limit=1000)
			# integration in 15.43
			res2 = integ.quad(integrand2,ill,ilu,
				args=(a1,b1,c,d,e,f,g,T),points=([ill*2.0]),limit=1000)
			# 15.44 and 15.43
			W_c = 1.0/(f*res[0])
			W_k = (-1.0/(2.0*c*(si.k*1.0e3)*T))*res2[0]
	
			# -----------------------------------------
			# particle Knudsen number calculated above
			
			# Cunningham slip-correction factor (dimensionless)
			Gi = 1.0+Kni[i]*(1.249+0.42*(np.exp(-0.87/Kni[i])))
			Gj = 1.0+Knj[j]*(1.249+0.42*(np.exp(-0.87/Knj[j])))
			# particle diffusion coefficient (15.29) (m2/s) (note the 
			# Boltzmann constant has units (kg m2)/(s2 K), so *1e3 to
			# convert to g from kg)
			Dpi = (((si.k*1.0e3)*T)/(6.0*np.pi*ri*eta_ai))*Gi
			Dpj = (((si.k*1.0e3)*T)/(6.0*np.pi*rj*eta_aj))*Gj	
			Mi = ((4.0/3.0)*np.pi*ri**3.0)*1.0e6
			Mj = ((4.0/3.0)*np.pi*rj**3.0)*1.0e6
			vbari = ((8.0*si.k*1.0e3*T)/(np.pi*Mi))**0.5 #15.32
			vbarj = ((8.0*si.k*1.0e3*T)/(np.pi*Mj))**0.5 #15.32
			# mean free path (15.34)
			lami = (8.0*Dpi)/(np.pi*vbari)
			lamj = (8.0*Dpj)/(np.pi*vbarj)
	
			Knp = ((lami**2.0+lamj**2.0)**0.5)/(ri+rj)
	
			# -----------------------------------------
			# get the van der Waals/viscous collision correction 
			# factor (15.42):		
			fac = 4.0*(Dpi+Dpj)/(((vbari**2.0+vbarj**2.0)**0.5)*(ri+rj))
			V_E = (W_c*(1.0+fac))/(1.0+(W_c/W_k)*fac)
			
			# fill in results
			res_Knp[i, j] = Knp
			res_all[i, j] = V_E
		
	
	# Van der Waals/collision coagulation kernel
	if vdWon == 0:
		K_V = K_B*(1.0-1.0) # when omitting van der Waals correction for expediency
	else:
		K_V = K_B*(res_all-1.0)

	# -------------------------------------------------------------
	
	# total coagulation kernel (m3/particle.s), with sbr in rows and rint in columns
	# eq. 15.27 of Jacobson (2005) says that the sum of kernels should be multiplied 
	# by a dimensionless coalescence efficiency.  For particles under 2um this should be
	# close to unity it says in the coalescence efficiency section.
	Beta = (K_B+K_DE+K_GC+K_TI+K_TS+K_V)
	
	if (testf == 1): # plot to compare to Fig. 15.7 of Jacobson (2005)
		ax0.loglog(sbr*10**6, Beta[:,0]*10**6, label='Total')
		ax1.loglog(sbr*10**6, Beta[:,1]*10**6, label='Total')
		plt.legend()
		plt.show()
		return()
	
	# zero beta for any size bins that have low number of particles
	ish = np.squeeze(num_part<1.0e-10)
	Beta[ish, :] = 0.
	ish = np.squeeze(num_part_rint<1.e-10)
	Beta[:, ish] = 0.
	
	# Perform coagulation, using the implicit approach of Jacobson (2005), given in 
	# eq. 15.5 

	# scale up by 1e6 to convert to cm3/particle.s from m3/particle.s
	# and therefore be consistent with particle concentrations which are 
	# (# particles/cm3 (air)).  For the production term due to coagulation below, we use
	# only the k,j coordinates in beta, where j goes as high as k-1, therefore, production
	# only uses the lower left triangle in beta.  However, the loss term below uses
	# k,j coordinates where j goes from 1 to the number of size bins
	if (coag_on == 1):
		Beta = Beta*1.e6
	if (coag_on == 0):
		Beta = Beta*0.
	if (testf == 2): # if testing the coagulation equations below, fix the kernels
		Beta = np.array(([[0.2, 0.4], [0.4, 0.5]]))
	if (testf == 3): # if comparing against Smoluchowski (1916 Fig. 4), fix kernel
		Beta = np.ones((sbn, sbn))*5.969331112806462e-10

	# matrix with volume of coagulated particles resulting from pairing of k and j 
	# particles single particle volumes of k and j (m3) 
	sbVmatk = (sbVi.reshape(-1, 1)).repeat(sbn, 1)
	sbVmatj = (sbVj.reshape(1, -1)).repeat(sbrn, 0)
	
	# combined volume of single coagulated particles (m3)
	coagV = sbVmatk+sbVmatj
	
	# matrix for number concentration at t-h in size bins repeated across rows
	num_partj = np.tile(num_part.reshape(1, -1), (sbn, 1))
	# matrix for volume concentration at t-h in size bins repeated across rows
	vol_partj = np.tile(vol_part.reshape(1, -1), (sbn, 1))
	
	# molecular concentration for k and j particles (# particles /cc (air)), components in
	# rows and size bins in columns.  j will represent the original concentration,
	# whereas molec_k will be updated during size bin loop below
	molec_k = np.zeros((num_molec.shape))
	molec_j = np.zeros((num_molec.shape))
	molec_k[:,:] = num_molec[:,:]
	molec_j[:,:] = num_molec[:,:]
	
	# use eq. 15.8 of Jacobson (2005) to estimate n_{k,t}
	# size bin loop, starting (importantly) with smallest size bin
	for sbi in range(sbn):
	
		# matrix for updated number concentrations (all those with index<sbi) are
		# concentrations at t, not t-h (concentrations spread across columns)
		# (# particles/cc (air))
		num_partk = np.tile(num_part.reshape(-1, 1), (1, sbrn))
		
		# full-moving structure - where coagulated particles are distributed between existing
		# size bins so that particle volumes at size bin centres are maintained
		#if (siz_str == 1):
		#	volind = np.zeros((coagV.shape[0], coagV.shape[1]))
			# only the coagV between this size bin and the one below, note that coagV equal
			# to this size bin accounted for below (coagV between this size bin and 
			# the one above)
		#	if (sbi > 0):
		#		coagV_sub_ind = (coagV>sbVi[0, sbi-1])*(coagV<sbVi[0, sbi])
				# fraction going to this size bin rather than one smaller
		#		frac_sub_in = (coagV-sbVi[0, sbi-1])/(sbVi[0, sbi]-sbVi[0, sbi-1])
		#		volind[coagV_sub_ind] = frac_sub_in[coagV_sub_ind] 

			# only the coagV between this size bin and the one above
		#	if (sbi < sbn-1):
		#		coagV_exc_ind = (coagV >= sbVi[0, sbi])*(coagV < sbVi[0, sbi+1])
				# fraction going to this size bin rather than one larger
		#		frac_exc_in = (sbVi[0, sbi+1]-coagV)/(sbVi[0, sbi+1]-sbVi[0, sbi])				
		#		volind[coagV_exc_ind] = frac_exc_in[coagV_exc_ind]
		#	if (sbi == sbn-1): # uppermost size bin
		#		coagV_exc_ind = (coagV >= sbVi[0, sbi])
		#		volind[coagV_exc_ind] = 1.

		if (siz_str == 0 or siz_str == 1): # moving-centre structure
			# index of k,j size bin pairs that can coagulate to give a volume that fits 
			# into k
			volind = (coagV>=sbbound[0, sbi])*(coagV<sbbound[0, sbi+1])

		# disallow number concentration gain when coagulating with itself 
		volind[sbi, :] = 0.
		volind[:, sbi] = 0.

		# using only the relevant k-j pairs in beta, number of particles from k and j
		# size bins coagulating to give new particle in sbi (# particles/cc.s).
		# This accounts for both the upper and lower diagonal of Beta, so considers
		# k-j pairs as well as j-k
		numsum_ind0 = (Beta*volind)*(num_partk*num_partj)
		# integrate over time
		numsum_ind0 = tint*numsum_ind0
		# multiply by 0.5 since 2 particles make 1
		numsum_ind0 = 0.5*numsum_ind0
		# sum to get the total rate of number of particles coagulating (from k and j) 
		# (# particles/cc.s)
		numsum0 = (numsum_ind0.sum()).sum()
		# numerator, note we use num_partj as want n_{k,t-h}
		num = num_partj[0, sbi]+numsum0
		
		#if (siz_str == 1): # full-moving structure
		#	if (sbi < sbn-1):
				# index of sbi pairs that give a volume beyond next size bin up
		#		volind = coagV[sbi, :] >= sbVi[0, sbi+1]
		#		volind = volind.astype(float)
				# lose a proportion of the number of particles from sbi coagulating 
				# with smaller particles to give a particle between sbi and sbi+1
		#		indx = coagV[sbi, :] < sbVi[0, sbi+1]
		#		volind[indx] = (coagV[sbi, indx]-sbVi[0, sbi])/(sbVi[0, sbi+1]-sbVi[0, sbi])
		#	if (sbi == sbn-1): # uppermost size bin
		#		volind = np.zeros((coagV.shape[1]))
				# lose half the number of particles of sbi coagulating with itself to 
				# give a volume in sbi
		#		volind[sbi] = 0.5

		if (siz_str == 0 or siz_str == 1): # moving-centre structure
			# particle number only from sbi when newly coagulated particles give a volume
			# outside the volume bounds of sbi
			volind = (coagV[sbi, :] >= sbbound[0, sbi+1])
			# lose half the number of particles of sbi coagulating with itself to give a 
			# volume within sbi.  Note, if uncertain about whether this should be 0.5 or
			# whether a factor 2 should be applied when self-coagulation gives a 
			# particle that fits a larger size bin, then run the test_coag to ensure 
			# that the loss from here
			# is number conservative with the gain in larger size bins
			if (coagV[sbi, sbi] < sbbound[0, sbi+1]):
				volind = volind.astype(float)
				volind[sbi] = 0.5

		# denominator, representing particle number loss from this size bin
		den_parts = Beta[sbi, :]*volind*num_partj[0, :]	
		den = 1.0+tint*(den_parts.sum())
		# updated number concentration (# particles/cc (air)) eq. 15.8 Jacobson (2005)
		num_part[0, sbi] = num/den
		
		# --------------------------------------------------------------------------------
		# change to molecular concentration part (molecules/cc (air))

		#if (siz_str == 1): # full-moving structure

		#	volind = np.zeros((coagV.shape[0], coagV.shape[1]))
			# only the coagV between this size bin and the one below, note that coagV equal
			# to this size bin accounted for below (coagV between this size bin and 
			# the one above)
		#	if (sbi > 0):
		#		coagV_sub_ind = (coagV>sbVi[0, sbi-1])*(coagV<sbVi[0, sbi])
				# fraction going to this size bin rather than one smaller
		#		frac_sub_in = (coagV-sbVi[0, sbi-1])/(sbVi[0, sbi]-sbVi[0, sbi-1])
		#		volind[coagV_sub_ind] = frac_sub_in[coagV_sub_ind] 

			# the coagV between this size bin and the one above
		#	if (sbi < sbn-1):
		#		coagV_exc_ind = (coagV >= sbVi[0, sbi])*(coagV < sbVi[0, sbi+1])
				# fraction going to this size bin rather than one larger
		#		frac_exc_in = (sbVi[0, sbi+1]-coagV)/(sbVi[0, sbi+1]-sbVi[0, sbi])				
		#		volind[coagV_exc_ind] = frac_exc_in[coagV_exc_ind]
		#	if (sbi == sbn-1): # uppermost size bin
		#		coagV_exc_ind = (coagV >= sbVi[0, sbi])
		#		volind[coagV_exc_ind] = 1.
		

		if (siz_str == 0 or siz_str == 1): # moving-centre structure
			# index of k,j size bin pairs that can coagulate to give a volume that fits 
			# into k
			volind = (coagV>=sbbound[0, sbi])*(coagV<sbbound[0, sbi+1])
			volind = volind.astype(float)		

		# consistency with particle number gain calculation above - note that we cannot
		# use N_{k,t} as this not yet known
		volind[sbi, :] = 0.
		
		# to find the number of molecules transferring from other size bins to this 
		# one do not multiply by 0.5 as needed to for number concentration
		numsum_ind = (Beta*volind)*(num_partk*num_partj)
		numsumk = ((numsum_ind.sum(axis=1))).reshape(-1, 1)
		
		# particle number concentration (# particles/cc (air)) gained by sbi from each bin
		num_contr = ((numsumk)*tint)[:, 0]
		# molecular concentration gained by sbi from each bin, note that if 
		# molec_k (new concentration) used instead of molec_j (old concentration), 
		# mass conservation issues can arise when two neighbouring size bins with very
		# different original number concentration coagulate to 
		# give a particle in a larger size bin
		ish = (num_partj[0, :] > 0.) # index of bins containing particles
		molec_contr = ((num_contr[ish]/num_partj[0, ish]).reshape(1, -1)*molec_j[:, ish]).sum(axis=1)	

		# molecular concentration loss part (molecules/cc (air)) ------------------------------
		#if (siz_str == 1): # full-moving structure
		#	volind = np.zeros((coagV.shape[1]))
		#	if (sbi < sbn-1): # can't lose molecules from largest size bin due to coagulation
				# index of sbi pairs that give a volume beyond next size bin up
		#		volind = coagV[sbi, :]>=sbVi[0, sbi+1]
				# lose a fraction of the molecular concentration of sbi coagulating 
				# with itself to give a volume in sbi
		#		if (coagV[sbi, sbi] < sbVi[0, sbi+1]):
		#			volind[sbi] = 1.-(sbVi[0, sbi+1]-coagV[sbi, sbi])/(sbVi[0, sbi+1]-sbVi[0, sbi])
		
		if (siz_str == 0 or siz_str == 1): # moving-centre structure
			# particle number only from sbi when newly coagulated particles give a volume
			# outside the current bounds, note size bin bound on uppermost size bin set very
			# high, so should be excluded from molecular concentration loss
			volind = (coagV[sbi, :] >= sbbound[0, sbi+1])
			# half loss of particle number when self-coagulation produces a particle within
			# sbi - note we are not yet at accounting for molecules lost, as still considering
			# particle number
			if (coagV[sbi, sbi] < sbbound[0, sbi+1]):
				volind = volind.astype(float)
				volind[sbi] = 0.

		num_lost_corr = num_partj[0, sbi]-(num_partj[0, sbi]/(1.+tint*sum(Beta[sbi, :]*volind*num_partj[0, :])))

		# number concentration of particle number lost from sbi, 
		# Eq. 15.8 with production term removed
		#num_lost =  num_partj[0, sbi] - num_partj[0, sbi]/(1.0+tint*((Beta[sbi, :]*volind*num_partj[0, :]).sum()))
		
		# fraction of particles lost via coagulation with each size bin
		# identify size bins where coagulation occurs
		#ish = ((Beta[sbi, :]*volind[:]*num_partj[0, :]) > 0.)
		#frac_loss = np.zeros((sbn))
		
		#frac_loss[ish] = ((Beta[sbi, ish]*volind[ish]*num_partj[0, ish])/
		#	(Beta[sbi, ish]*volind[ish]*num_partj[0, ish]).sum())
		
		# ignore any particles lost due to self-coagulation that gives a particle in sbi
		#if (coagV[sbi, sbi] < sbbound[0, sbi+1]):
		#	num_lost_corr = num_lost*(1.-frac_loss[sbi])
		#else:
		#	num_lost_corr = num_lost
		# molecular concentration represented by this loss
		# note that loss equation in GMD paper is equal to the num_lost_corr equation 
		# above and this
		# equation combined
		molec_loss = 0.
		if (num_partj[0, sbi]>0.):
			molec_loss = molec_j[:, sbi]*(num_lost_corr/num_partj[0, sbi])
		
		# new molecular concentration in sbi
		molec_k[:, sbi] += (molec_contr-molec_loss)

		if (testf == 2):
			# number of particles (not molecules) being lost
			num_part_lost = num_partj[0, sbi]-num_partj[0, sbi]/den
			# those lost from each bin
			num_part_lost_sb = num_part_lost*(den_parts/den_parts.sum())	
			print(str('total particle number concentration lost from size bin ' +str(sbi)+': '+str(num_part_lost)))
			print(str('particle number concentration lost from size bin ' +str(sbi)+' by coagulation with each size bin (in ascending order of size bin): '+str(num_part_lost_sb)))
			print(str('total particle number concentration gained by size bin '+str(sbi)+': '+str(numsum0)))
			print(str('particle number concentration gained by size bin ' +str(sbi)+' from each coagulating pair (in ascending order of size bin): '+str(numsum_ind0)))
			print(str('new particle concentrations per size bin : '+str(num_part) ))
			print(str('new molecular concentration per size bin : ' + str(molec_k)))
			print('type c return at the command line to continue')
			if (sbi == 1):
				return()
	# using new molecular concentration, calculate new dimensions per size bin 
	MV = (M[:, 0]/(rho)).reshape(num_comp, 1) # molar volume (cc/mol)
	# new volume of single particle per size bin (um3)
	ish = (num_part[0, :] > 1.0e-20) # only use size bins where particles reside
	Vnew = np.zeros((sbrn))
	Vnew[ish] = np.sum(((molec_k[:, ish]/(si.N_A*num_part[0, ish]))*MV*1.e12), 0)

	# new combined volume of all particles per size bin (um3)
	Vtot = Vnew*num_part[0, :]
	# remove particles and their corresponding component concentrations if their volume
	# is negligibly small
	negl_indx = (Vtot == 0.) # index of size bins with zero particle concentration
	num_part[0, negl_indx] = 0.
	Vnew[negl_indx] = 0.
	molec_k[:, negl_indx] = 0.
	
	# new radius per size bin (um)
	rad = ((3.0*Vnew)/(4.0*np.pi))**(1.0/3.0)
	# size bins with no particle assigned central radius
	ish = num_part[0, :] <= 1.e-20
	rad[ish] = rad0[ish]
	Vnew[ish] = V0[ish]
	
	# just want particle number concentration as an array with one dimension
	if num_part.ndim>1:
		num_part = num_part[0, :]
	
	# molecular concentrations (molecules/cc (air))
	y = molec_k.flatten(order='F')

	if (siz_str == 0): # moving centre
		(num_part, Vnew, y, rad, redt, blank, tnew) = movcen(num_part.reshape(-1, 1), 
		sbbound[0, :]*1.e18, 
		np.transpose(y.reshape(sbn, num_comp)), 
		sbn, num_comp, sbVi[0, :]*1.e18, 0.0, 0, MV*1.e12)
		# revert volume bounds to um3 from m3 before returning
		sbbound = sbbound[0, :]*1.e18

	if (siz_str == 1): # full-moving
		(Vnew, rad, y, 
		num_part, sbbound, rbou) = fullmov.fullmov(sbn, num_part.reshape(-1, 1),
 		num_comp, y, MV*1.e12, V0, sbbound[0, :]*1.e18, rbou)

	sbbound = sbbound.reshape(-1) # return to 1D array

	return(num_part, y, rad, Gi, eta_ai, Vnew, sbbound, rbou)
