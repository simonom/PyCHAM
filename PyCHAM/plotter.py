'''plots results'''
# simulation results are represented graphically

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import os
import user_input as ui
import retr_out
import numpy as np
import scipy.constants as si

def plotter(caller):
	
	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0) or tests (2) are the calling module
	# --------------------------------------------------------------------------

	# retrieve useful information from pickle file
	[sav_name, sch_name, indx_plot, Comp0] = ui.share(1)

	if (sav_name == 'default_res_name'):
		print('Default results name was used, therefore results not saved and nothing to plot')
		return()

	dir_path = os.getcwd() # current working directory
	# obtain just part of the path up to PyCHAM home directory
	for i in range(len(dir_path)):
		if dir_path[i:i+7] == 'PyCHAM':
			dir_path = dir_path[0:i+7]
			break
	# isolate the scheme name from path to scheme
	for i in range(len(sch_name)-1, 0, -1):
		if sch_name[i] == '/':
			sch_name = sch_name[i+1::]
			break
	# remove any file formats
	for i in range(len(sch_name)-1, 0, -1):
		if sch_name[i] == '.':
			sch_name = sch_name[0:i]
			break
	

	dir_path = str(dir_path+'/PyCHAM/output/'+sch_name+'/'+sav_name)

	# chamber condition ---------------------------------------------------------

	# retrieve results
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, comp_names, 
		_, _, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(dir_path)
	
	# number of actual particle size bins
	num_asb = num_sb-wall_on

	if caller == 0:
		plt.ion() # show results to screen and turn on interactive mode

	# prepare sub-plots depending on whether particles present
	if (num_asb) == 0: # no particle size bins
		if not (indx_plot): # check whether there are any gaseous components to plot
			print('Please note no initial gas-phase concentrations were received and no particle size bins were present, therefore there is nothing for the standard plot to show')
			return()
		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))

	if (num_asb > 0): 
		if not (indx_plot):
			print('Please note, no initial gas-phase concentrations were registered, therefore the gas-phase standard plot will not be shown')
			fig, (ax1) = plt.subplots(1, 1, figsize=(14, 7))
		else:
			fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(14, 7))

		par1 = ax1.twinx() # first parasite axis
		par2 = ax1.twinx() # second parasite axis
		
		# Offset the right spine of par2.  The ticks and label have already been
		# placed on the right by twinx above.
		par2.spines["right"].set_position(("axes", 1.2))
		# Having been created by twinx, par2 has its frame off, so the line of its
		# detached spine is invisible.  First, activate the frame but make the patch
		# and spines invisible.
		make_patch_spines_invisible(par2)
		# Second, show the right spine.
		par2.spines["right"].set_visible(True)	

	if (indx_plot):
		# gas-phase concentration sub-plot ---------------------------------------------	
		for i in range(len(indx_plot)):
		
			ax0.semilogy(timehr, yrec[:, indx_plot[i]], '+',linewidth=4.0, 
						label=str(str(Comp0[i])))

		ax0.set_ylabel(r'Gas-phase concentration (ppb)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize=14)

		# find maximum and minimum of plotted concentrations for sub-plot label		
		maxy = max(yrec[:, indx_plot].flatten())
		miny = min(yrec[:, indx_plot].flatten())	
			
		ax0.text(x = timehr[0]-(timehr[-1]-timehr[0])/10., y = maxy+((maxy-miny)/10.), s='a)', size=14)

		# end of gas-phase concentration sub-plot ---------------------------------------
	
	# particle properties sub-plot --------------------------------------------------
	if (num_asb > 0): # if particles present		

		if timehr.ndim == 0: # occurs if only one time step saved
			Ndry = np.array(Ndry.reshape(1, num_asb))
		if num_asb == 1: # just one particle size bin (wall included in num_sb)
			Ndry = np.array(Ndry.reshape(len(timehr), num_asb))

		if timehr.ndim == 0: # occurs if only one time step saved
			x = np.array(x.reshape(1, num_asb))
		if num_asb == 1: # just one particle size bin (wall included in num_sb)
			x = np.array(x.reshape(len(timehr), num_asb))

		if timehr.ndim==0: # occurs if only one time step saved
			rbou_rec = np.array(rbou_rec.reshape(1, num_sb))

		
		# plotting number size distribution --------------------------------------
	
		# don't use the first boundary as it's zero, so will error when log10 taken
		log10D = np.log10(rbou_rec[:, 1::]*2.)
		if (num_asb > 1) :
			# note, can't append zero to start of log10D to cover first size bin as the log10 of the
			# non-zero boundaries give negative results due to the value being below 1, so instead
			# assume same log10 distance as the next pair
			log10D = np.append((log10D[:, 0]-(log10D[:, 1]-log10D[:, 0])).reshape(-1, 1), log10D, axis=1)
			# radius distance covered by each size bin (log10(um))
			dlog10D = (log10D[:, 1::]-log10D[:, 0:-1]).reshape(log10D.shape[0], log10D.shape[1]-1)
		if (num_asb == 1): # single particle size bin
			# assume lower radius bound is ten times smaller than upper
			dlog10D = (log10D[:, 0]-np.log10((rbou_rec[:, 1]/10.)*2.)).reshape(log10D.shape[0], 1)
			
	
		# number size distribution contours (/cc (air))
		dNdlog10D = np.zeros((Ndry.shape[0], Ndry.shape[1]))
		dNdlog10D[:, :] = Ndry[:, :]/dlog10D[:, :]
		# transpose ready for contour plot
		dNdlog10D = np.transpose(dNdlog10D)
		
		# mask the nan values so they're not plotted
		z = np.ma.masked_where(np.isnan(dNdlog10D), dNdlog10D)
	
		# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
		colors = [(0.60, 0.0, 0.70), (0, 0, 1), (0, 1.0, 1.0), (0, 1.0, 0.0), (1.0, 1.0, 0.0), (1.0, 0.0, 0.0)]  # R -> G -> B
		n_bin = 100  # Discretizes the colormap interpolation into bins
		cmap_name = 'my_list'
		# Create the colormap
		cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
	
		# set contour levels
		levels = (MaxNLocator(nbins = 100).tick_values(np.min(z[~np.isnan(z)]), 
				np.max(z[~np.isnan(z)])))
	
		# associate colours and contour levels
		norm1 = BoundaryNorm(levels, ncolors = cm.N, clip=True)
		
		# contour plot with times (hours) along x axis and 
		# particle diameters (nm) along y axis
		for ti in range(len(timehr)-1): # loop through times
			p1 = ax1.pcolormesh(timehr[ti:ti+2], (rbou_rec[ti, :]*2*1e3), z[:, ti].reshape(-1, 1), cmap=cm, norm=norm1)
	
		# if logarithmic spacing of size bins specified, plot vertical axis 
		# logarithmically
		if space_mode == 'log':
			ax1.set_yscale("log")
		ax1.set_ylabel('Diameter (nm)', size = 14)
		ax1.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax1.yaxis.set_tick_params(labelsize = 14, direction = 'in')

		# label according to whether gas-phase plot also displayed		
		if (indx_plot):
			ax1.text(x=timehr[0]-(timehr[-1]-timehr[0])/11., y = np.amax(rbou_rec*2*1e3)*1.05, s='b)', size=14)
		else:
			ax1.text(x=timehr[0]-(timehr[-1]-timehr[0])/11., y = np.amax(rbou_rec*2*1e3)*1.3, s='a)', size=14)
		ax1.set_xlabel(r'Time through simulation (hours)', fontsize=14)
		
		cb = plt.colorbar(p1, format=ticker.FuncFormatter(fmt), pad=0.25)
		cb.ax.tick_params(labelsize=14)   
		# colour bar label
		cb.set_label('dN/dlog10(D) $\mathrm{(cm^{-3})}$', size=14, rotation=270, labelpad=20)

		# ----------------------------------------------------------------------------------------
		# total particle number concentration #/cm3
	
		# include total number concentration (# particles/cc (air)) on contour plot
		# first identify size bins with radius exceeding 3nm
		# empty array for holding total number of particles
		Nvs_time = np.zeros((Ndry.shape[0]))
	
		for i in range(num_asb): # size bin loop
			# get the times when bin exceeds 3nm - might be wanted to deal with particle counter detection limits
# 			ish = x[:, i]>3.0e-3
# 			Nvs_time[ish] += Ndry[ish, i] # sum number
			Nvs_time[:] += Ndry[:, i] # sum number
	
		p3, = par1.plot(timehr, Nvs_time, '+k', label = 'N')
	
		par1.set_ylabel('N (# $\mathrm{cm^{-3})}$', size=14, rotation=270, labelpad=20) # vertical axis label
		par1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e')) # set tick format for vertical axis
		par1.yaxis.set_tick_params(labelsize=14)

		# SOA mass concentration ---------------------------------------------------------------
		# array for SOA sum with time
		SOAvst = np.zeros((1, len(timehr)))
		
		final_i = 0
		# check whether water and/or core is present
		if comp_names[-2] == 'H2O': # if both present
			final_i = 2
		if comp_names[-1] == 'H2O': # if just water
			final_i = 1
		# note that the seed component is only registered in init_conc_func if initial
		# particle concentration (pconc) exceeds zero, therefore particle-phase material 
		# must be present at start of experiment (row 0 in y)
		if final_i == 0 and y[0, num_speci:(num_speci*(num_asb+1))].sum()>1.0e-10: 
			final_i = 1
		
		for i in range(num_asb): # particle size bin loop
	
			# sum of organics in condensed-phase at end of simulation (ug/m3 (air))
			
			# to replicate the SMPS results, find the volume of particles then
			# assume a density of 1.0 g/cm3
			SOAvst[0, :] += np.sum((yrec[:, ((i+1)*num_comp):((i+2)*num_comp-final_i)]/si.N_A*(y_MV[0:-final_i])*1.0e12), axis = 1)

		# log10 of maximum in SOA
		if (max(SOAvst[0, :]) > 0):
			SOAmax = int(np.log10(max(SOAvst[0, :])))
		else:
			SOAmax = 0.
		# transform SOA so no standard notation required
		SOAvst[0, :] = SOAvst[0, :]
		
	
		p5, = par2.plot(timehr, SOAvst[0, :], 'xk', label = '[secondary]')
		par2.set_ylabel(str('[secondary] ($\mathrm{\mu g\, m^{-3}})$'), rotation=270, size=16, labelpad=25)
		# set label, tick font and [SOA] vertical axis to red to match scatter plot presentation
		par2.yaxis.label.set_color('black')
		par2.tick_params(axis='y', colors='black')
		par2.spines['right'].set_color('black')
		par2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e')) # set tick format for vertical axis
		par2.yaxis.set_tick_params(labelsize=16)
		par2.text((timehr)[0], max(SOAvst[0, :])/2.0, 'assumed particle density = 1.0 $\mathrm{g\, cm^{-3}}$')
		plt.legend(fontsize=14, handles=[p3, p5] ,loc=4)	

	# end of particle properties sub-plot -----------------------------------

	# save and display
	plt.savefig(str(dir_path+'/'+sav_name+'_output_plot.png'))
	if caller == 2:
		plt.show()	

	return()

# function to makes spines invisible
def make_patch_spines_invisible(ax):
	ax.set_frame_on(True)
	ax.patch.set_visible(False)
	for sp in ax.spines.values():
		sp.set_visible(False)

# function for doing colorbar tick labels in standard notation
def fmt(x, pos):
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	return r'${} \times 10^{{{}}}$'.format(a, b)
