'''module to track particle number size distribution using moving centre size structure (p. 416 of Jacobson 2000)'''
# note, no limits on number of size bins a particle can move in
# this version of moving centre

import numpy as np
import scipy.constants as si

def mov_cen_main(n0, Vbou, Cn, sbn, nc, Vol0, t, tinc_count, MV):



	# input:---------------------------------------------------------
	# n0 - particle number concentration per size bin before time step
	# (# particle/cc (air)) (excluding wall)
	# Vbou - volume bounds per size bin (um3) (number of size bins +1)
	# Cn - particle phase concentration per component per size bin
	# (molecules/cc (air)), with rows representing
	# components and columns representing size bins
	# sbn - number of size bins excluding wall (if present)
	# nc - number of components
	# Vol0 - original volume of size bins (um3) (excluding wall)
	# t - integration time (s)
	# tinc_count - count on number of time steps since time interval last required 
	# decreasing
	# MV - molar volume per component (um3/mol)
	# ---------------------------------------------------------------
	
	NA = si.Avogadro # Avogadro's number (molecules/mol)	

	Vnew = np.zeros((sbn))
	ish = n0[:, 0] > 0. # size bins containing particle number concentration
	nmolC = np.zeros((nc, sbn))
	
	# number of moles of each component in a single particle (mol/cc (air))
	nmolC[:, ish] = ((Cn[:, ish]/(NA*n0[ish, 0])))
	MVrep = np.repeat(MV, nmolC.shape[1], axis=1)
	
	
	# new volume of single particle per size bin (um3) including volume of water
	Vnew[:] = np.sum(nmolC*MV, axis=0)
	Vnew[n0[:, 0] <= 1.e-10] = Vol0[0::][n0[:, 0] <= 1.e-10]

	# array of new particle number concentration (# particle/cc (air))
	num_part_new = np.zeros((sbn, 1))
	# array of new molecular concentration (# molecules/cc (air))
	num_molec_new = np.zeros((nc, sbn))
	# loop through size bins and reallocate particle number concentration 
	# (# particle/cc (air)) and molecular concentration (# molecules/cc (air))
	for i in range(sbn):
	
		# return and send message if beyond uppermost boundary
		if Vnew[i]>=Vbou[-1]:
			print('particle in size bin ' +str(i) + ' exceed uppermost volume bound, therefore will try reducing ode solver time step')
			redt = 1
			tnew = t/2.0 # new time for integration on next step (s)
			# note, that when called from water initiator or coag, this condition will
			# terminate the model, so the returned values are meaningless
			return(0, 0, 0, 0, 0, 0, 0)
			
		Vindx = (Vnew[i] >= Vbou[0:-1]).sum()-1 # index of where these particle fit
		num_part_new[Vindx, 0] += n0[i] # allocate number concentration
		num_molec_new[:, Vindx] += Cn[:, i] # allocate molecular concentration
		
	# remove any particles and their constituents if number concentration is relatively very low
	rm_indx = num_part_new[:, 0] < (sum(num_part_new)*1.e-8)
	num_part_new[rm_indx, 0] = 0.
	num_molec_new[:, rm_indx] = 0.
	
	
	# need to find new volumes of single particles (um3) because possible that particles
	# originally in separate size bins now share a size bin
	# total volume of components 
	# ((um3 (all particles)/cc (air))/(particle number/cc (air))) 
	# calculation is:
	# divide number of molecules/cc (air) by Na to get moles/cc(air), then 
	# multiply by um3/mol (MV) to get um3 (of each component)/cc (air),
	# then sum volume of components per size bin to get um3 (all particles)/cc (air)
	MVrep = np.repeat(MV, num_molec_new.shape[1], axis=1)
	
	Vtot = (np.sum((num_molec_new/NA)*(MVrep), 0)) # (um3)

	# then divide by particle number (#/cc (air)) to get volume of single particles
	# per size bin (um3)
	Vsing = np.zeros(len(Vtot))
	isb = num_part_new[:, 0]>0.
	Vsing[isb] = Vtot[isb]/num_part_new[isb, 0]
	
	# new radius of single particles per size bin (um)
	rad = ((3.0*Vsing)/(4.0*np.pi))**(1.0/3.0)
	
	isb = num_part_new[:, 0] == 0.
	
	# fill volume array elements for bins without particles with central volume (um3)
	Vsing[isb] = Vol0[isb]
	rad[isb] = ((3.*Vsing[isb])/(4.*np.pi))**(1./3.) # new radius (um)
	
	# flag to show no reduction in time step needed
	redt = 0
	tnew = t # maintain same time step (s)

	return(num_part_new, Vsing, np.ravel(np.transpose(num_molec_new)), rad, 
			redt, t, tnew)
