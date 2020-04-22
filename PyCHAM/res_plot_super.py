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
import user_input as ui

def run(testf):
	
	# --------------------------------------------------------------
	# inputs:
	# testf - flag for whether in test mode (=0 no test, =1 testing from test_PyCHAM.py,
	#			=2 testing from test_res_plot_super.py)
	
	# --------------------------------------------------------------
	if testf==1:
		print('"Plot Results" button works fine')
		return()
	if testf==2:
		print('in testing mode and called from test_res_plot_super.py, will try plotting results from example_run output, note that test_res_plot_super.py must be called from the PyCHAM home directory')
		cwd = os.getcwd() # address of current working directory
		fname = str(cwd+'/PyCHAM/output/Example_Run')
		resfname = 'Example_output'
		y_indx_plot = np.array(([1, 312, 2, 3]))
		Comp0 = np.array((['O3', 'APINENE', 'NO', 'NO2']))
	if testf==0:
		# module to ask, receive and return required inputs
		[fname, resfname, y_indx_plot, Comp0] = ui.run(1,0)
	
	# path to saved files
	dir_path = os.getcwd()    
		
	output_root = 'PyCHAM/output'
	filename = os.path.basename(fname)
	filename = os.path.splitext(filename)[0]
	
	# one folder for one simulation
	output_by_sim = os.path.join(dir_path, output_root, filename, resfname)
	
	# name of file where experiment constants saved (number of size bins and whether wall 
	# included)
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
			if str(line.split(',')[0]) == 'factor_for_multiplying_ppb_to_get_molec/cm3':
				dlist.append(float(i))
				
		const[str(line.split(',')[0])] = dlist
	
	num_sb = int((const['number_of_size_bins'])[0]) # number of size bins
	num_speci = int((const['number_of_components'])[0]) # number of species
	y_mw = const['molecular_weights_g/mol_corresponding_to_component_names']
	y_MV = const['molecular_volumes_cm3/mol']
	PyCHAM_names = const['component_names']
	# conversion factor to change gas-phase concentrations from molecules/cc 
	# (air) into ppb
	Cfactor = float((const['factor_for_multiplying_ppb_to_get_molec/cm3'])[0])

	
	# name of file where concentration (molecules/cc (air)) results saved
	fname = str(output_by_sim+'/concentrations_all_components_all_times_gas_particle_wall')
	y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
	if y.ndim==1: # occurs if only one time step saved
		y = np.array(y.reshape(1, -1))
	
	# withdraw times
	fname = str(output_by_sim+'/time')
	t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)	
	if t_array.ndim==0: # occurs if only one time step saved
		print('Please note only results for one time step have been saved; number size distribution contours will not be plotted')
		t_array = np.array(t_array.reshape(-1, 1))	
	
	if testf == 0:
		plt.ion() # show figures on screen

	if num_sb>1:
		# name of file where concentration (# particles/cc (air)) results saved
		fname = str(output_by_sim+'/particle_number_concentration_dry')
		N = np.loadtxt(fname, delimiter=',', skiprows=1) # skiprows=1 omits header
		if t_array.ndim==0: # occurs if only one time step saved
			N = np.array(N.reshape(1, num_sb-1))
		if num_sb==2: # just one particle size bin (wall included in num_sb)
			N = np.array(N.reshape(len(t_array), num_sb-1))
			
	
		# name of file where particle size results saved
		fname = str(output_by_sim+'/size_bin_radius')
		x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header
		if t_array.ndim==0: # occurs if only one time step saved
			x = np.array(x.reshape(1, num_sb-1))
		if num_sb==2: # just one particle size bin (wall included in num_sb)
			x = np.array(x.reshape(len(t_array), num_sb-1))
	
		# size bin radius bounds (um) (for the number size distribution contour plot)
		fname = str(output_by_sim+'/size_bin_bounds')
		sbb = np.loadtxt(fname, delimiter=',', skiprows=1) # skiprows=1 omits header)
		if t_array.ndim==0: # occurs if only one time step saved
			sbb = np.array(sbb.reshape(1, num_sb))
	
	
		def make_patch_spines_invisible(ax):
			ax.set_frame_on(True)
			ax.patch.set_visible(False)
			for sp in ax.spines.values():
				sp.set_visible(False)
	
	
		# ------------------------
		# plotting number-size distribution
	
		# don't use the first boundary as it's zero, so will error when log10 taken
		log10D = np.log10(sbb[1::]*2.0)
		if num_sb>2:
			# note, can't append zero to start of log10D to cover first size bin as the log10 of the
			# non-zero boundaries give negative results due to the value being below 1, so instead
			# assume same log10 distance as the next pair
			log10D = np.append((log10D[0]-(log10D[1]-log10D[0])).reshape(1, 1), log10D.reshape(-1,1), axis=0)
			# radius distance covered by each size bin (log10(um))
			dlog10D = (log10D[1::]-log10D[0:-1]).reshape(1, -1)
		if num_sb==2: # number of size bins includes wall
			# assume lower radius bound is ten times smaller than upper
			dlog10D = (log10D-np.log10((sbb[1::]/10.0)*2.0)).reshape(1, -1)
			
		
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
	
		fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(14, 7))
		fig.subplots_adjust(right=0.75)
	
		par1 = ax0.twinx() # first parasite axis
		par2 = ax0.twinx() # second parasite axis
	
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
		p1 = ax0.pcolormesh(t_array/3600.0, (sbb*2*1e3), z[:, :], cmap=cm, norm=norm1)
	
		ax0.set_ylabel('Diameter (nm)', size=14)
		ax0.xaxis.set_tick_params(labelsize=14)
		ax0.yaxis.set_tick_params(labelsize=14)
		ax0.text(x=t_array[0]/3600-(t_array[-1]/3600-t_array[0]/3600)/10.0, y=max(sbb*2*1e3)*1.05, s='a)', size=14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize=14)
		# function for doing colorbar tick labels in standard notation
		def fmt(x, pos):
			a, b = '{:.1e}'.format(x).split('e')
			b = int(b)
			return r'${} \times 10^{{{}}}$'.format(a, b)
		cb = plt.colorbar(p1, format=ticker.FuncFormatter(fmt), pad=0.25)
		cb.ax.tick_params(labelsize=14)   
		# colour bar label
		cb.set_label('dN/dlog10(D) $\mathrm{(cm^{-3})}$', size=14, rotation=270, labelpad=20)
	
	
		# ----------------------------------------------------------------------------------------
		# total particle number concentration #/cm3
	
		# include total number concentration (# particles/cc (air)) on contour plot
		# first identify size bins with radius exceeding 3nm
		# empty array for holding total number of particles
		Nvs_time = np.zeros((N.shape[0]))
		for i in range(num_sb-1):
			# get the times when bin exceeds 3nm - might be wanted to deal with particle counter detection limits
# 			ish = x[:, i]>3.0e-3
# 			Nvs_time[ish] += N[ish, i] # sum number
			Nvs_time[:] += N[:, i] # sum number
	
		p3, = par1.plot(t_array/3600.0, Nvs_time, '+k', label = 'N (sim)')
	
		par1.set_ylabel('N (# $\mathrm{cm^{-3})}$', size=14, rotation=270, labelpad=20) # vertical axis label
		par1.ticklabel_format(style='sci', scilimits=(0,0)) # scientific notation used on vertical axis
		par1.yaxis.set_tick_params(labelsize=14)
	
	
	
		# ----------------------------------------------------------------------------------------
	
		# SOA mass concentration
		# array for SOA sum with time
		SOAvst = np.zeros((1, len(t_array)))
		
		final_i = 0
		# check whether water and/or core is present
		if PyCHAM_names[-2] == 'H2O': # if both present
			final_i = 2
		if PyCHAM_names[-1] == 'H2O': # if just water
			final_i = 1
		# note that the seed component is only registered in init_conc_func if initial
		# particle concentration (pconc) exceeds zero, therefore particle-phase material 
		# must be present at start of experiment (row 0 in y)
		if final_i == 0 and y[0, num_speci:(num_speci*(num_sb+1))].sum()>1.0e-10: 
			final_i = -1
		
		for i in range(num_sb): # size bin loop, including wall
	
			# sum of organics in condensed-phase at end of simulation (ug/m3 (air))
			if i < num_sb-1:
			
				# to replicate the SMPS results when using 
				# low size bin resolution in the model, find the volume of particles, then
				# assume a density of 1.0 g/cm3
				SOAvst[0, :] += np.sum((y[:, ((i+1)*num_speci):((i+2)*num_speci-final_i)]/si.N_A*(y_MV[0:-final_i])*1.0e12), axis = 1)
		
		# log10 of maximum in SOA
		SOAmax = int(np.log10(max(SOAvst[0, :])))
		# transform SOA so no standard notation required
		SOAvst[0, :] = SOAvst[0, :]/(10**(SOAmax))
		
	
		p5, = par2.plot(t_array/3600.0, SOAvst[0, :], 'xr', label = 'M (sim)')
		par2.set_ylabel(str('[SOA]/ ' + str(10**(SOAmax)) + ' ($\mathrm{\mu g\, m^{-3}})$'), rotation=270, size=16, labelpad=25)
		# set label, tick font and [SOA] vertical axis to red to match scatter plot presentation
		par2.yaxis.label.set_color('red')
		par2.tick_params(axis='y', colors='red')
		par2.spines['right'].set_color('red')
		par2.yaxis.set_tick_params(labelsize=16)
		par2.text((t_array/3600.0)[0], max(SOAvst[0, :])/2.0, 'assumed particle density for [SOA] = 1.0 $\mathrm{g\, cm^{-3}}$')
		plt.legend(fontsize=14, handles=[p3, p5] ,loc=4)
		

	# -----------------------------------------------------------------------------------
	# gas-phase concentrations
	
	for i in range(len(y_indx_plot)):
		ax1.semilogy(t_array/3600, y[:, y_indx_plot[i]], '+',linewidth=4.0, 
					label=str(str(Comp0[i])+' (sim)'))
		
		ish = y[:, y_indx_plot[i]]>0.0 # prevent log10 of zero
		
		if i == 0:
			miny = min(np.log10(y[ish, y_indx_plot[i]]))
			maxy = max(np.log10(y[ish, y_indx_plot[i]]))
		if min(np.log10(y[ish, y_indx_plot[i]]))<miny:
			miny = min(np.log10(y[ish, y_indx_plot[i]]))
		if max(np.log10(y[ish, y_indx_plot[i]]))>maxy:
			maxy = max(np.log10(y[ish, y_indx_plot[i]]))
			
	ax1.set_ylabel(r'Gas-phase concentration (ppb)', fontsize=14)
	ax1.set_xlabel(r'Time through simulation (hours)', fontsize=14)
	ax1.yaxis.set_tick_params(labelsize=14)
	ax1.xaxis.set_tick_params(labelsize=14)
	ax1.legend(fontsize=14)
	
	ax1.text(x=t_array[0]/3600-(t_array[-1]/3600-t_array[0]/3600)/10.0, y=maxy, s='b)', size=14)
	
	if testf == 0:
		plt.savefig(str(output_by_sim+'/'+str(resfname)+'_output_plot.png'), transparent=True)
		os.remove('PyCHAM/var_store.pkl') # remove pickle file
		
	if testf == 2:
		plt.show()
	
	return()