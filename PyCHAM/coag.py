'''module to estimate coagulation kernel using p. 508 of Jacobson (2005), called on by ode_gen'''

import numpy as np 
import scipy.constants as si
from fl_reg_determ import reg_determ
from Reyn_num import Reyn_num
from W_k_int import W_k_int
import scipy.integrate as integ
import matplotlib.pyplot as plt
import ipdb

def coag(RH, T, sbr, sbVi, M, rint, num_molec, num_part, tint, sbbound,
			num_comp, vdWon, rho, rad0, PInit):

	# --------------------------------------------------------------
	# inputs:
	
	# RH - relative humidity (fraction)
	# T - temperature (K)
	# sbr - size bin radius (m)
	# sbVi - single particle volume for i sizes (m3)
	# M - molecular weight of components (g/mol)	
	# rint - size of interest (m)
	# num_molec - molecular concentration (molecules/cc (air)),
	# arranged by component in rows and size bins in columns
	# num_part - concentration of particles per size bin 
	# (particle/cc (air)) (columns) (excluding walls)
	# tint - time interval coagulation occurs over (s)
	# sbbound - size bin volume boundaries (m3)
	# num_comp - number of components
	# vdWon - saying whether the van der Waals kernel should be calculated or ignored
	# rho - components densities (g/cm3) in a 1D array
	# rad0 - original radius at size bin centre (um)
	# PInit - pressure inside chamber (Pa)
	
	# --------------------------------------------------------------
	# outputs:
	
	# Beta - sum of coagulation kernels
	# --------------------------------------------------------------

	# wall parts
	num_molec = num_molec[:, 0:-1]
	num_part = num_part.reshape(1, -1)
	
	# ensure sbn is integer
	sbn = np.int(np.max(rint.shape))
	sbrn = np.int(np.max(sbr.shape))

	# call on function to determine the Knudsen no. and therefore flow 
	# regime of each size bin
	[Kni, eta_ai, rho_ai, kin_visc] = reg_determ(RH, T, sbr, PInit)
	[Knj, eta_aj, rho_aj, kin_visc] = reg_determ(RH, T, rint, PInit)
	# Reynold number and terminal fall velocity for each size bin
	[Rei, Vfi] = Reyn_num(sbr, eta_ai, rho_ai, kin_visc, 1.0e6, Kni)
	[Rej, Vfj] = Reyn_num(sbr, eta_aj, rho_aj, kin_visc, 1.0e6, Knj)
	
	# repeat Knudsen number over number of size bins
	Kni_m = np.tile(Kni, (sbn, 1))
	
	# Cunningham slip-flow correction (15.30) with constant taken
	# from text below (dimensionless)
	Gi = 1.0+Kni*(1.249+0.42*(np.exp(-0.87/Kni)))
	Gj = 1.0+Knj*(1.249+0.42*(np.exp(-0.87/Knj)))
	# particle diffusion coefficient (15.29) (m2/s)
	# multiply eta_a by 1.0e-3 to convert from g/m.s to kg/m.s
	# this makes it consistent with the units of Boltzmann constant
	Dpi = (((si.k*T)/(6.0*np.pi*sbr*(eta_ai*1.0e-3)))*Gi).reshape(sbrn, 1)
	Dpj = (((si.k*T)/(6.0*np.pi*rint*(eta_aj*1.0e-3)))*Gj).reshape(1, sbn)

	# repeat Dpj over size bins of i (m)
	Dp_mj = np.repeat(Dpj, sbrn, 0)
	# matrix for sums of Dp in size bins (m2/s)
	Dp_sum = Dpi+Dp_mj
	# zero the upper triangle part of the matrix as only the lower is needed
	Dp_sum[np.triu_indices(sbrn, 1)] = 0.0
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
	rint2[np.triu_indices(sbrn, 1)] = 0.0
	sbr_m[np.triu_indices(sbrn, 1)] = 0.0
	
	# matrix for size bins sums (m)
	sbr_sum = sbr2+rint2
	
	
	# size bins in continuum regime
	i = ((Kni<1.0).reshape(sbrn,1))
	# spread across rint (j)
	i = np.repeat(i, sbn, 1)
	
	# collision kernel (15.28) (m3/particle.s)
	K_B[i] = 4.0*np.pi*(sbr_sum[i])*(Dp_sum[i])
	
	
	
	# size bins in free-molecular regime	
	i = ((Kni>10.0).reshape(sbrn,1))
	# spread across rint (j)
	i = np.repeat(i, sbn, 1)
	
	# single particle mass (g):
	# first, number of moles per component in a single particle
	num_mol_single = (num_molec/num_part)/si.N_A
	# second product of number of moles and molecular weight
	weight_compon = num_mol_single*M
	# final sum mass across components for single particle mass (g)
	Mpi = np.sum(weight_compon, 0)
	Mpj = np.sum(weight_compon, 0)

	# thermal speed of particle (15.32) (m/s) (multiply mass by 1.0e-3 to 
	# convert from g to kg and therefore be consistent with Boltzmann's 
	# constant (1.380658e-23kgm2/s2.K.molec))
	nu_pi = (((8.0*si.k*T)/(np.pi*(Mpi*1.0e-3)))**0.5).reshape(sbrn, 1)
	nu_pj = np.zeros((sbrn, sbn))
	ind = Mpj>0 # prevent division by zero
	nu_pj[ind] = ((8.0*si.k*T)/(np.pi*(Mpj[ind]*1.0e-3)))**0.5
	
	# sum of squares of speeds (m2) (15.31)
	nu_p_sum = nu_pi**2.0+nu_pj**2.0
			
	# collision kernel (15.31)
	K_B[i] = np.pi*(sbr_sum[i]**2.0)*((nu_p_sum[i])**(0.5))
	
	# size bins in transition regime
	i = (1.0<=Kni) 
	j = (Kni<=10.0) # size bins in transition regime		
	i = (i*j).reshape(sbrn,1)	# size bins in transition regime
	# spread across rint (j)
	i = np.repeat(i, sbn, 1)
	
	# particle mean free path (15.34) (m)
	lam_pi = ((8.0*Dpi)/(np.pi*nu_pi))
	lam_pj = np.zeros((sbrn, sbn))
	ind = nu_pj>0 # prevent division by zero
	lam_pj[ind] = ((8.0*Dp_mj[ind])/(np.pi*nu_pj[ind]))
	# mean distance from centre of a sphere travelled by particles
	# leaving sphere's surface and travelling lam_p (m) (15.34)
	sig_pi = (((2.0*sbr+lam_pi)**3.0-(4.0*sbr**2.0+lam_pi**2.0)**1.5)/(6.0*sbr*lam_pi)-2.0*sbr)
	sig_pj = np.zeros((sbrn, sbn))
	ind = rint2>0 # prevent division by zero
	sig_pj[ind] = (((2.0*rint2[ind]+lam_pj[ind])**3.0-(4.0*rint2[ind]**2.0+lam_pj[ind]**2.0)**1.5)/(6.0*rint2[ind]*lam_pj[ind])-2.0*np.repeat(sbr2, sbn, 1)[ind])

	# sum mean distances (m)
	sig_p_sum = sig_pi**2.0+sig_pj**2.0
		
	# kernel numerator
	K_Bnum = 4.0*np.pi*sbr_sum[i]*Dp_sum[i]

	# left term kernel denominator
	K_Blden = (sbr_sum[i]/(sbr_sum[i]+(sig_p_sum[i])**0.5))
	
	# right term kernel denominator
	K_Brden = ((4.0*Dp_sum[i])/(((nu_p_sum[i])**0.5)*sbr_sum[i]))
	
	# collision kernel (15.33)	
	K_B[i] = (K_Bnum/(K_Blden+K_Brden))
	
	
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
	
	# Kernel for Turbulent Inertial Motion:

	# rate of dissipation of turbulent kinetic energy per gram of medium 
	# (m2/s3) (8.4 for a typical value (which is apparently taken 
	# from Pruppacher and Klett 1997 (p. 511 Jacobson (2005)))
	epsilon	= 5.0e-4
	# kernel for turbulent inertial motion (15.40)
	K_TI = (((np.pi*epsilon**(3.0/4.0))/(si.g*kin_visc**(1.0/4.0)))*(
		sbr_sum**2.0)*del_Vf)
	
	
	# kernel for Turbulent Shear (15.41)
	K_TS = ((8.0*np.pi*epsilon)/(15.0*kin_visc))**0.5*(sbr_sum**3.0)
	
	# -----------------------------------------------------------------
	# van der Waals/viscous collision kernel:	
	
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
		
	
	# van der Waals/collision coagulation kernel
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
	# zero beta for any size bins that don't have particles
	ish = np.squeeze(num_part<1.0e-10)
	Beta[ish, :] = 0.0
	Beta[:, ish] = 0.0
	
	# Perform coagulation, using the implicit approach of Jacobson (2005), given in 
	# eq. 15.5 

	# product of Beta and time interval (m3/particle), scale up by 1e6 to convert to 
	# cm3/particle and therefore consistent with particle concentrations which are 
	# (# particles/cm3 (air))
	Beta = Beta*tint*1.0e6
	
	# now produce matrices for number concentration in size bins repeated across rows
	# and columns
	num_partk = np.tile(num_part.reshape(-1, 1), (1, sbrn))
	num_partj = np.tile(num_part.reshape(1, -1), (sbn, 1))
	
	# matrix with volume of coagulated particles resulting from pairing of k and j 
	# particles
	# single particle volumes of i and j (m3) 
	sbVmati = (sbVi.reshape(sbrn, 1)).repeat(sbrn, 1)
	sbVmatj = sbVi.repeat(sbrn, 0)
	
	# combined volume of single coagulated particles (m3)
	coagV = sbVmati+sbVmatj
	# set as lower triangular matrices, as the upper triangle is just a repetition of the 
	# lower and therefore would cause duplication
	coagV[np.triu_indices(sbrn, k=1, m=None)] = -1.0
	
	MV = ((M[:, 0]/(rho)).reshape(num_comp, 1)) # molar volume (cc/mol)
	# total volume of molecules in each size bin (m3)
	V0 = np.sum(((num_molec[:, :]/(6.0221409e+23*num_part[0, :]))*MV*1.0e-6), 0)
	
	# for recording new number of molecules following coagulation
	num_molec2 = np.zeros((num_comp, sbn))
	num_molec2[:, :] = num_molec[:, :]
	
	
	for sbi in range(sbn): # loop through size bins
	
		# index of where new particles for this size bin come from
		new_ind = (coagV>=sbbound[0, sbi])*(coagV<sbbound[0, sbi+1])
		Beta2 = np.zeros((sbn, sbrn))
		Beta2[:, :] = Beta*new_ind # Beta with only the relevant kernels left
		
		# number of particles coagulating for every k-j pair, eq. 15.5 of Jacobson (2005)
		numcoag = Beta2*num_partk*num_partj
		numcoagtot = Beta*num_partk*num_partj # all kernels, for loss calculation
		# effect change in number of particles
		Pn_gain = ((numcoag).sum().sum())
		Pn_lost = ((numcoagtot[sbi, :]).sum())+((numcoagtot[:, sbi]).sum())
		num_part[0, sbi] += Pn_gain-Pn_lost

		# fraction of number of k particles coagulating with j for production
		numfrack = (numcoag/num_partk)
		# sum over columns to get total fraction of k contributing to new particles
		numfrack = numfrack.sum(1)
		# fraction of number of j particles coagulating with k for production
		numfracj = (numcoag/num_partj)
		# sum over rows to get total fraction of j contributing to new particles
		numfracj = numfracj.sum(0)
		
		# now multiply each by number of molecules per size bin and sum over size bins
		# for total molecules going to new size bin
		molec_gain = (numfrack*num_molec+numfracj*num_molec).sum(1)
			
		# now, for molecules lost from k
		# first get number of particles lost from k and divide by total number
		numfrack = ((numcoagtot[sbi, :].sum()+numcoagtot[:, sbi].sum()))/num_partk[sbi, 0]
		
		molec_lost = numfrack*num_molec[:, sbi]
		molec_change = molec_gain-molec_lost
		num_molec2[:, sbi] += molec_change
			

	MV = (M[:, 0]/(rho)).reshape(num_comp, 1) # molar volume (cc/mol)
	# new volume of single particle per size bin (um3)
	ish = num_part[0, :]>1.0e-20 # only use size bins where particles reside
	Vnew = np.zeros((sbrn))
	Vnew[ish] = np.sum(((num_molec2[:, ish]/(6.0221409e+23*num_part[0, ish]))*MV*1.0e12), 
				0)

	# new radius per size bin (um)
	ish = num_part[0, :]<=1.0e-20 # only use size bins where particles reside
	rad = ((3.0*Vnew)/(4.0*np.pi))**(1.0/3.0)
	rad[ish] = rad0[ish]

	# return number of particles/cc(air) per size bin (columns) and
	# number of molecules (molecules/cc(air)) flattened into species followed by size bins
	return(np.squeeze(num_part), num_molec2.flatten(order='F'), rad, Gi, eta_ai, Vnew)