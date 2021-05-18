'''plots a replication of mass spectrum as reported by a chemical ionisation mass spectrometer (CIMS)'''
# simulation results are represented graphically

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import retr_out
import numpy as np
import scipy.constants as si

def plotter_CIMS(dir_path, mcr_res, tn, iont):

	# inputs: -----------------
	# dir_path - path to folder containing results files to plot
	# mcr_res - resolution of molar mass to charge ratio (g/mol/charge)
	# tn - time through experiment to plot at (s)
	# iont - type of ionisation
	# ---------------------------

	# retrieve results, note that num_sb (number of size bins)
	# includes wall if wall turned on
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_MW, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, H2Oi, _, _, _, RO2i) = retr_out.retr_out(dir_path)
	
	# convert to 2D numpy array
	y_MW = np.array((y_MW)).reshape(-1, 1)
	
	# get index of time wanted
	ti = (np.where(np.abs(timehr-tn/3600.) == np.min(np.abs(timehr-tn/3600.))))[0][0]
	
	# convert yrec from 1D to 2D with times in rows, then select time wanted
	yrec = (yrec.reshape(len(timehr), num_comp*(num_sb+1)))[ti, :]
	
	# get gas-phase concentrations (ug/m3, note starting concentration is ppb)
	gp = yrec[0:num_comp]*Cfac[ti]/si.N_A*y_MW[:, 0]*1.e12
	
	# get particle-phase concentrations (molecules/cm3)
	pp = yrec[num_comp:num_comp*(num_sb+1-wall_on)]
	
	# sum each component over size bins (molecules/cm3)
	pp = np.sum(pp.reshape(num_sb-wall_on, num_comp), axis=0)
	
	# convert to ug/m3
	pp = (pp/si.N_A)*y_MW[:, 0]*1.e12
	
	# remove water
	gp = np.append(gp[0:H2Oi], gp[H2Oi+1::])
	pp = np.append(pp[0:H2Oi], pp[H2Oi+1::])
	y_MW = np.append(y_MW[0:H2Oi, 0], y_MW[H2Oi+1::, 0])

	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
	ax0.semilogy(y_MW, gp, '+m', label = str('gas-phase'))
	ax0.semilogy(y_MW, pp, 'xb', label = str('particle-phase'))
	
	ax0.set_title(str('Mass spectrum at ' + str(timehr[ti]) + ' hours'), fontsize = 14)
	ax0.set_xlabel(r'molar mass:charge ($\mathrm{g\,mol^{-1}}$)', fontsize = 14)
	ax0.set_ylabel(r'Concentration ($\mu\,\mathrm{g\,m^{-3}}$)', fontsize = 14)
	ax0.legend(fontsize = 14)
	plt.show()

	return()