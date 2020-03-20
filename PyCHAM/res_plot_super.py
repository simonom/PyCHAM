'''module to plot results from PyCHAM, called on when 'Plot Results' button selected'''

# used this guide: https://matplotlib.org/3.1.0/gallery/ticks_and_spines/multiple_yaxis_with_spines.html

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import sys
import scipy.constants as si
import ipdb
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import user_input as ui

def run(testf):
	
	# --------------------------------------------------------------
	# inputs:
	# testf - flag for whether in test mode (=0 no test, =1, testing from test_PyCHAM.py)
	
	# --------------------------------------------------------------
	if testf==1:
		print('"Plot Results" button works fine')
		return()
	
	# module to ask, receive and return required inputs
	[fname, resfname, y_indx_plot, Comp0] = ui.run(1,0)
	
	# path to saved files
	dir_path = os.getcwd()    
		
	output_root = 'PyCHAM/output'
	filename = os.path.basename(fname)
	filename = os.path.splitext(filename)[0]
	# one folder for one simulation
	output_by_sim = os.path.join(dir_path, output_root, filename, resfname)
	
	# name of file where experiment constants saved
	fname = str(output_by_sim+'/model_and_component_constants')

	const_in = open(fname)
	const = {} # prepare to create dictionary
	for line in const_in.readlines():
	
		# convert to python list
		dlist = []
		for i in line.split(',')[1::]:
			if str(line.split(',')[0]) == 'number_of_size_bins':
				dlist.append(int(i))
			if str(line.split(',')[0]) == 'number_of_components':
				dlist.append(int(i))
			if str(line.split(',')[0]) == 'molecular_weights_g/mol_corresponding_to_component_names' or  str(line.split(',')[0]) == 'molecular_volumes_cm3/mol':
				i = i.strip('\n')
				i = i.strip('[')
				i = i.strip(']')
				i = i.strip(' ')
				dlist.append(float(i))
			if str(line.split(',')[0]) == 'component_names':
				i = i.strip('\n')
				i = i.strip('[')
				i = i.strip(']')
				i = i.strip(' ')
				i = i.strip('\'')
				dlist.append(str(i))
		const[str(line.split(',')[0])] = dlist
	
	num_sb = (const['number_of_size_bins'])[0] # number of size bins
	num_speci = (const['number_of_components'])[0] # number of components
	y_mw = const['molecular_weights_g/mol_corresponding_to_component_names']
	y_MV = const['molecular_volumes_cm3/mol']
	PyCHAM_names = const['component_names']
	
	# name of file where concentration (molecules/cc (air)) results saved
	fname = str(output_by_sim+'/concentrations_all_components_all_times_gas_particle_wall')
	y = np.loadtxt(fname, delimiter=',', skiprows=2) 
	if y.ndim==1: # occurs if only one time step saved
		y = np.array(y.reshape(1, -1))
	
	# withdraw times
	fname = str(output_by_sim+'/time')
	t_array = np.loadtxt(fname, delimiter=',', skiprows=1) # skiprows=1 omits header)	
	if t_array.ndim==0: # occurs if only one time step saved
		print('Please note only results for one time step have been saved; number size distribution contours will not be plotted')
		t_array = np.array(t_array.reshape(-1, 1))	
	
	plt.ion() # show figures on screen

	if num_sb>1:
		# name of file where concentration (# particles/cc (air)) results saved
		fname = str(output_by_sim+'/particle_number_concentration')
		N = np.loadtxt(fname, delimiter=',', skiprows=2) # skiprows=1 omits header)
		if N.ndim==1: # occurs if only one time step saved
			N = np.array(N.reshape(1, num_sb-1))	
	
		# name of file where particle size results saved
		fname = str(output_by_sim+'/size_bin_radius')
		x = np.loadtxt(fname, delimiter=',', skiprows=2)
		if x.ndim==1: # occurs if only one time step saved
			x = np.array(x.reshape(1, num_sb-1))
	
		# calculate radius bounds (um) for the number size distribution contour plot 
		# (assumes a linear size distribution)
		xbound = x[0, :]+(x[0, 1]-x[0, 0])/2.0
		xbound = np.append(0.0, xbound) # bottom bound
	
		# size bin radius bounds (um) for calculating dN/dlog10D (/cc (air))
		# name of file where particle size results saved
		fname = str(output_by_sim+'/size_bin_bounds')
		sbb = np.loadtxt(fname, delimiter=',', skiprows=2)
	
	
		def make_patch_spines_invisible(ax):
			ax.set_frame_on(True)
			ax.patch.set_visible(False)
			for sp in ax.spines.values():
				sp.set_visible(False)
	
	
	
	
		# ------------------------
		# plotting number-size distribution
	
		# convert particle concentration to dN/d(log10(D))
		# the log 10 of size bin diameter bounds
		# undo change from Size_distributions.py where the upper bin boundary was increased by two
		# orders of magnitude to reduce possibility of particles growing beyond this
		sbb[-1] = sbb[-1]*1e-2
	
		# don't use the first boundary as it's zero, so will error when log10 taken
		log10D = np.log10(sbb[1::]*2.0)
		# note, can't append zero to start of log10D to cover first size bin as the log10 of the
		# non-zero boundaries give negative results due to the value being below 1, so instead
		# assume same log10 distance as the next pair
		log10D = np.append((log10D[0]-(log10D[1]-log10D[0])).reshape(1, 1), log10D.reshape(-1,1), axis=0)
		# radius distance covered by each size bin (log10(um))
		dlog10D = (log10D[1::]-log10D[0:-1]).reshape(1, -1)
		# repeat over times
		dlog10D = np.repeat(dlog10D, N.shape[0], axis=0)
	
		ish = (dlog10D != 0)
		# prepare number size distribution contours (/cc (air))
		dNdlog10D = np.zeros((N.shape[0], N.shape[1]))
	
	
		dNdlog10D[ish] = N[ish]/dlog10D[ish]
		# ------------------------------------------------------------------------------
		# plot to show number concentration as contours on a plot with time along
		# the x axis and radius along the y axis
	
	
		# smooth particle concentration results with a four-point moving average (# particles/cm3)
		# dNdlog10D_smooth = np.zeros((np.shape(dNdlog10D)))
		# dNdlog10D_smooth[:, 1:-2] = (dNdlog10D[:, 0:-3]+dNdlog10D[:, 1:-2]+dNdlog10D[:, 2:-1]+dNdlog10D[:, 3::])/4.0
		# dNdlog10D_smooth[:, 0] = dNdlog10D[:, 0]
		# dNdlog10D_smooth[:, -1] = dNdlog10D[:, -1]
	
		# smooth particle concentration results with a three-point moving average (# particles/cm3)
		# dNdlog10D_smooth = np.zeros((np.shape(dNdlog10D)))
		# dNdlog10D_smooth[:, 1:-1] = (dNdlog10D[:, 0:-2]+dNdlog10D[:, 1:-1]+dNdlog10D[:, 2::])/3.0
		# dNdlog10D_smooth[:, 0] = dNdlog10D[:, 0]
		# dNdlog10D_smooth[:, -1] = dNdlog10D[:, -1]
	
		# smooth particle concentration results with a two-point moving average (# particles/cm3)
		#dNdlog10D_smooth = np.zeros((np.shape(dNdlog10D)))
		#dNdlog10D_smooth[:, 1::] = (dNdlog10D[:, 0:-1]+dNdlog10D[:, 1::])/2.0
		#dNdlog10D_smooth[:, 0] = dNdlog10D[:, 0]
		#dNdlog10D_smooth[:, -1] = dNdlog10D[:, -1]
	
		dNdlog10D_smooth = np.zeros((np.shape(dNdlog10D)))
		dNdlog10D_smooth[:,:] = dNdlog10D[:,:] # don't smooth option
	
		fig, host = plt.subplots(figsize=(14, 7))
		fig.subplots_adjust(right=0.75)
	
		par1 = host.twinx() # first parasite axis
		par2 = host.twinx() # second parasite axis
	
		# Offset the right spine of par2.  The ticks and label have already been
		# placed on the right by twinx above.
		par2.spines["right"].set_position(("axes", 1.2))
		# Having been created by twinx, par2 has its frame off, so the line of its
		# detached spine is invisible.  First, activate the frame but make the patch
		# and spines invisible.
		make_patch_spines_invisible(par2)
		# Second, show the right spine.
		par2.spines["right"].set_visible(True)
	
		# transpose number concentration results, so time on x axis and diameter on y
		dNdlog10D_smooth = dNdlog10D_smooth.transpose()
	
		# mask the nan values so they're not plotted
		z = np.ma.masked_where(np.isnan(dNdlog10D_smooth), dNdlog10D_smooth)
	
		# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
		colors = [(0.60, 0.0, 0.70), (0, 0, 1), (0, 1.0, 1.0), (0, 1.0, 0.0), (1.0, 1.0, 0.0), (1.0, 0.0, 0.0)]  # R -> G -> B
		n_bin = 100  # Discretizes the colormap interpolation into bins
		cmap_name = 'my_list'
		# Create the colormap
		cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
	
		# set contour levels
		# levels = (MaxNLocator(nbins = 20).tick_values(np.min(z[~np.isnan(z)]), 
		# 			np.max(z[~np.isnan(z)])))
		levels = (MaxNLocator(nbins = 100).tick_values(np.min(z[~np.isnan(z)]), 
				np.max(z[~np.isnan(z)])))
	
		# associate colours and contour levels
		norm1 = BoundaryNorm(levels, ncolors = cm.N, clip=True)
		# contour plot with times along x axis and particle diameters along y axis
		#p1 = host.contour(t_array/3600.0, (x[:,0]*2*1e3), z[:, :], cmap=cm, norm=norm1)
		p1 = host.pcolormesh(t_array/3600.0, (xbound*2*1e3), z[:, :], cmap=cm, norm=norm1)
	
	
		host.set_xlabel('Time (hrs)', size=18)
		host.set_ylabel('Diameter (nm)', size=18)
		host.xaxis.set_tick_params(labelsize=18)
		host.yaxis.set_tick_params(labelsize=18)
	
		# function for doing colorbar tick labels in standard notation
		def fmt(x, pos):
			a, b = '{:.1e}'.format(x).split('e')
			b = int(b)
			return r'${} \times 10^{{{}}}$'.format(a, b)
		cb = plt.colorbar(p1, format=ticker.FuncFormatter(fmt), pad=0.25)
		cb.ax.tick_params(labelsize=18)   
		# colour bar label
		cb.set_label('dN/dlog10(D) $\mathrm{(cm^{-3})}$', size=18, rotation=270, labelpad=20)
	
	
		# ----------------------------------------------------------------------------------------
		# total particle number concentration #/cm3
	
		# include total number concentration (# particles/cc (air)) on contour plot
		# first identify size bins with radius exceeding 3nm
		Nvs_time = np.zeros((N.shape[0])) # empty array for holding total number of particles
		for i in range(num_sb-1):
			ish = x[:, i]>3.0e-3 # get the times when bin exceeds 3nm
			Nvs_time[ish] += N[ish, i] # sum number
	
		p3, = par1.plot(t_array/3600.0, Nvs_time, '+k', label = 'N (sim)')
	
		par1.set_ylabel('N (# $\mathrm{cm^{-3})}$', size=18, rotation=270, labelpad=20) # vertical axis label
		par1.ticklabel_format(style='sci', scilimits=(0,0)) # scientific notation used on vertical axis
		par1.yaxis.set_tick_params(labelsize=18)
	
	
	
		# ----------------------------------------------------------------------------------------
	
		# SOA mass concentration
		# array for SOA sum with time
		SOAvst = np.zeros((1, len(t_array)))
	
		for i in range(num_sb): # size bin loop, including wall
	
			# sum of organics in condensed-phase at end of simulation (ug/m3 (air))
			if i < num_sb-1:
			
				# to replicate the SMPS results when using 
				# low size bin resolution in the model, find the volume of particles, then
				# assume a density of 1.3 g/cm3
				SOAvst[0, :] += np.sum((y[:, ((i+1)*num_speci):((i+2)*num_speci-2)]/si.N_A*(y_MV[0:-2])*1.0e12), axis = 1)
			
	
		p5, = par2.plot(t_array/3600.0, SOAvst[0, :], 'xr', label = 'M (sim)')
		par2.set_ylabel('[SOA] ($\mathrm{\mu g\, m^{-3}})$', rotation=270, size=18, labelpad=25)
		# set label, tick font and [SOA] vertical axis to red to match scatter plot presentation
		par2.yaxis.label.set_color('red')
		par2.tick_params(axis='y', colors='red')
		par2.spines['right'].set_color('red')
		par2.yaxis.set_tick_params(labelsize=18)
		plt.legend(fontsize=18, handles=[p3, p5] ,loc=4)

		plt.savefig(str(output_by_sim+'/contours.png'), transparent=True)

	# -----------------------------------------------------------------------------------
	# gas-phase concentrations
	
	plt.figure(figsize=(8,7))
	for i in range(len(y_indx_plot)):
		plt.plot(t_array/3600, y[:, y_indx_plot[i]], '+',linewidth=4.0, 
					label=str(str(Comp0[i])+' (sim)'))
		
	plt.ylabel(r'Gas-phase concentration (ppb)', fontsize=18)
	plt.xlabel(r'Time (hours)', fontsize=18)
	plt.yticks(fontsize=18)
	plt.xticks(fontsize=18)
	
	
	plt.legend(fontsize=18)
	plt.savefig(str(output_by_sim+'/gas_ppb.png'), transparent=True)
	
	os.remove('PyCHAM/var_store.pkl') # remove pickle file
	
	return()