'''plots a replication of number size distribution as reported by a Scanning Mobility Particle Spectrometer (SMPS)'''
# simulation results are represented graphically

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import os
import retr_out
import numpy as np
import scipy.constants as si
from scipy import interpolate
import math
import importlib

def plotter(caller, dir_path, self, dryf, cdt, sdt, min_size, max_size, csbn, p_rho):
	
	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0) or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# self - reference to GUI
	# dryf - whether particles dried (0) or not (1)
	# cdt - particle number concentration detection limit (particles/cm3)
	# sdt - particle size at 50 % detection efficiency (nm), 
	# 	width factor for detection efficiency dependence on particle size
	# min_size - minimum size measure by counter (nm)
	# max_size - maximum size measure by counter (nm)
	# csbn - number of size bins for counter
	# p_rho - assumed density of particles (g/cm3)
	# --------------------------------------------------------------------------

	# chamber condition ---------------------------------------------------------
	# retrieve results
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_mw, Nwet, _, y_MV, _, wall_on, space_mode, indx_plot, 
		comp0, _, PsatPa, OC, H2Oi, _, siz_str, _) = retr_out.retr_out(dir_path)
	
	# number of actual particle size bins
	num_asb = (num_sb-wall_on)

	if (caller == 0):
		plt.ion() # show results to screen and turn on interactive mode

	
	# prepare figure -------------------------------------------
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
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
	
	# ---------------------------------------------------------------
		
	# affect whether particle number concentration is for dry or wet inlet to instrument (# particles/cm3)
	if (dryf == 0):
		Nuse = Ndry
	else:
		Nuse = Nwet
	
	# if just one size bin, ensure two dimensional
	if (num_sb-wall_on) == 1:
		Nuse = Nuse.reshape(-1, 1)
		x = x.reshape(-1, 1)
	
	# account for minimum particle size detectable by instrument (um)
	Nuse[x<(min_size*1.e-3)] = 0.
	
	# if moving centre used rather than full moving then 
	# change Nuse, x and rbou_rec to a two-point moving average
	if (siz_str[0] == 0):
		# two point moving average number concentration
		Nuse = (Nuse[:, 0:-1]+Nuse[:, 1::])/2.
	
		if (space_mode == 'log'):
			# two point moving average size (radius) bin centres (um)
			x = 10.**(np.log10(x[:, 0:-1])+(np.log10(x[:, 1::])-np.log10(x[:, 0:-1]))/2.)
			# two point moving average size (radius) bin bounds (um)
			rbou_rec = 10.**(np.log10(rbou_rec[:, 0:-1])+(np.log10(rbou_rec[:, 1::])-np.log10(rbou_rec[:, 0:-1]))/2.)
			# fixed point centre (radius) of size bins (um)
			xf = 10.**(np.log10(rbou_rec[:, 0:-1])+(np.log10(rbou_rec[:, 1::])-np.log10(rbou_rec[:, 0:-1]))/2.)
	
		if (space_mode == 'lin'):
			# two point moving average size (radius) bin centres (um)
			x = x[:, 0:-1]+(x[:, 1::]-x[:, 0:-1])/2.
			# two point moving average size (radius) bin bounds (um)
			rbou_rec = rbou_rec[:, 0:-1]+(rbou_rec[:, 1::]-rbou_rec[:, 0:-1])/2.
			# fixed point centre (radius) of size bins (um)
			xf = rbou_rec[:, 0:-1]+(rbou_rec[:, 1::]-rbou_rec[:, 0:-1])/2.
	
	else: # if using full-moving then set the size bin centre radius (um) as the known
		xf = x 
	
	# difference in the log10 of size (diameter) bin widths of model (um) 
	# (could vary with time depending on size structure)
	sbwm = (np.log10(rbou_rec[:, 1::]*2.))-(np.log10(rbou_rec[:, 0:-1]*2.))

	# normalise number concentration by size (diameter) bin width (# particles/cm3/um)
	Nuse = Nuse/sbwm
	
	# interpolate particle number concentration to the counter size bins --------------------
	
	Nint = np.zeros((len(timehr), csbn)) # empty array to hold interpolated results (#/m3)
	
	# size bin bounds of SMPS (diameter) (um)
	csbb = (np.logspace(np.log10(min_size), np.log10(max_size), csbn+1, base = 10.))*1.e-3

	# size bin centres of instrument (diameter) (um)
	csbc = (10.**(np.log10(csbb[0:-1])+(np.log10(csbb[1::])-np.log10(csbb[0:-1]))/2.))
	
	# difference in log10 of size bins (diameter) bin widths of instrument (um)
	sbwc = np.log10(csbb[1::])-np.log10(csbb[0:-1])

	# interpolation based on fixed size bin centres - note that this is preferred over changeable
	# size bin centres when calculating number size distribution, total particle number 
	# concentration and total mass concentration, since the area under the
	# number vs size line is comparable between times with changeable size bin centres, note when
	# this is not the case, even when a particle grows when using moving centre size structure
	# the area beneath the curve may decrease because it becomes more jagged, thereby 
	# unrealistically affecting particle number concentration during interpolation
	for ti in range(len(timehr)):
		Nint[ti, :] = np.interp(csbc, xf[ti, :]*2., Nuse[ti, :]) # (# particles/cm3/difference in log 10(size(um)))
		Nint[ti, :] = Nint[ti, :]*sbwc # correct for size bin width (# particles/cm3)
	
	
	# account for minimum detectable particle concentration (# particles/cm3)
	Nint[Nint<cdt] = 0.
	
	# get detection efficiency as a function of particle size (nm)
	[Dp, ce] = count_eff_plot(3, 0, self, sdt)
	
	# interpolate detection efficiency (fraction) to instrument size bin centres
	ce = np.interp(csbc, Dp, ce)
	Nint = Nint*ce # correct for detection efficiency

	# plotting number size distribution --------------------------------------
	
	# take log10 of instrument size (diameter) bin boundaries
	log10D = np.log10(csbb)
	# take difference of log 10 of diameter at size bin boundaries
	dlog10D = log10D[1::]-log10D[0:-1]
	dlog10D = np.tile(dlog10D, [len(timehr), 1]) # tile over times
			
	# number size distribution contours ((# particle/cm3))
	dNdlog10D = np.zeros((Nint.shape[0], Nint.shape[1]))
	dNdlog10D[:, :] = Nint/dlog10D
	# transpose ready for contour plot
	dNdlog10D = np.transpose(dNdlog10D)
	
	# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
	colors = [(0.6, 0., 0.7), (0, 0, 1), (0, 1., 1.), (0, 1., 0.), (1., 1., 0.), (1., 0., 0.)]  # R -> G -> B
	n_bin = 100  # discretizes the colormap interpolation into bins
	cmap_name = 'my_list'
	# create the colormap
	cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
	
	# set contour levels
	levels = (MaxNLocator(nbins = 100).tick_values(np.min(dNdlog10D), np.max(dNdlog10D)))
	
	# associate colours and contour levels
	norm1 = BoundaryNorm(levels, ncolors=cm.N, clip=True)
		
	# contour plot with times (hours) along x axis and 
	# particle diameters (nm) along y axis
	for ti in range(len(timehr)-1): # loop through times
		p1 = ax0.pcolormesh(timehr[ti:ti+2], (csbb*1e3), dNdlog10D[:, ti].reshape(-1, 1), cmap=cm, norm=norm1)
	
	# plot vertical axis logarithmically
	ax0.set_yscale("log")
	
	# set tick format for vertical axis
	ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
	ax0.set_ylabel('Diameter (nm)', size = 14)
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which= 'both')
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which= 'both')

		
	ax0.set_xlabel(r'Time through simulation (hours)', fontsize=14)
	
	cb = plt.colorbar(p1, format=ticker.FuncFormatter(fmt), pad=0.25)
	cb.ax.tick_params(labelsize=14)   
	# colour bar label
	cb.set_label('dN (#$\,$$\mathrm{cm^{-3}}$)/d$\,$log$_{10}$(D ($\mathrm{\mu m}$))', size=14, rotation=270, labelpad=20)

	# total particle number concentration (# particles/cm3) -------------------------
	
	# include total number concentration (# particles/cm3 (air)) on contour plot
	Nvs_time = Nint.sum(axis = 1)
	
	p3, = par1.plot(timehr, Nvs_time, '-+k', label = 'N')
	
	par1.set_ylabel('N (#$\,$ $\mathrm{cm^{-3})}$', size=14, rotation=270, labelpad=20) # vertical axis label
	par1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e')) # set tick format for vertical axis
	par1.yaxis.set_tick_params(labelsize=14)

	# mass concentration of particles (ug/m3) ---------------------------------------------------------------
	
	MCvst = np.zeros((len(timehr))) # empty array for total mass concentration (ug/m3)
	
	# assumed volumes of particles per size bin (um3)
	Vn = ((4./3.)*np.pi)*((csbc/2.)**3.)
	Vn = Vn*1.e-12 # convert to cm3
	Vn = np.tile(Vn, [len(timehr), 1]) # tile over times
	# convert number concentration to volume concentration (cm3 (particle)/cm3 (air))
	Vc = Nint*Vn
	# convert volume concentration to total mass concentration (ug/m3)
	MCvst = ((Vc*p_rho)*1.e12).sum(axis=1)
		
	# log10 of maximum in mass concentration
	if (max(MCvst[:]) > 0):
		MCmax = int(np.log10(max(MCvst[:])))
	else:
		MCmax = 0.
	
	p5, = par2.plot(timehr, MCvst[:], '-xk', label = 'Total Particle Mass Concentration')
	par2.set_ylabel(str('Mass Concentration ($\mathrm{\mu g\, m^{-3}})$'), rotation=270, size=16, labelpad=25)
	# set colour of label, tick font and corresponding vertical axis to match scatter plot presentation
	par2.yaxis.label.set_color('black')
	par2.tick_params(axis='y', colors='black')
	par2.spines['right'].set_color('black')
	par2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e')) # set tick format for vertical axis
	par2.yaxis.set_tick_params(labelsize=16)
	plt.legend(fontsize=14, handles=[p3, p5] ,loc=4)

	if (caller == 2): # display when in test mode
		plt.show()	


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

# condensation particle counter (CPC)
def cpc_plotter(caller, dir_path, self, dryf, cdt, max_dt, sdt, max_size, uncert, 
		delays, wfuncs, Hz, loss_func_str, losst):

	import rad_resp_hum
	import inlet_loss
	
	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0) or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# self - reference to GUI
	# dryf - relative humidity of aerosol at entrance to condensing unit of CPC (fraction 0-1)
	# cdt - false background counts (# particles/cm3)
	# max_dt - maximum detectable concentration (# particles/cm3)
	# sdt - particle size at 50 % detection efficiency (nm), 
	# 	width factor for detection efficiency dependence on particle size
	# max_size - maximum size measure by counter (nm)
	# uncert - uncertainty (%) around counts by counter
	# delays - the significant response times for counter
	# wfuncs - the weighting as a function of time for particles of different age
	# Hz - temporal frequency of output
	# loss_func_str - string stating loss rate (fraction/s) as a 
	#			function of particle size (um)
	# losst - time of passage through inlet (s)
	# --------------------------------------------------------------------------

	# required outputs ---------------------------------------------------------
	# retrieve results
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_mw, Nwet, _, y_MV, _, wall_on, space_mode, indx_plot, 
		comp0, _, PsatPa, OC, H2Oi, _, siz_str, _) = retr_out.retr_out(dir_path)
	# ------------------------------------------------------------------------------
	
	# condition wet particles assuming equilibrium with relative humidity at 
	# entrance to condensing unit of CPC.  Get new radius at size bin centre (um)
	[xn, yrec[:, num_comp:(num_sb-wall_on+1)*(num_comp)]]  = rad_resp_hum.rad_resp_hum(yrec[:, num_comp:(num_sb-wall_on+1)*(num_comp)], x, dryf, H2Oi, num_comp, (num_sb-wall_on), Nwet, y_MV)
	
	# remove particles lost during transit through inlet (# particles/cm3)
	[Nwet,  yrec[:, num_comp:(num_sb-wall_on+1)*(num_comp)]]= inlet_loss.inlet_loss(Nwet, xn, yrec[:, num_comp:(num_sb-wall_on+1)*(num_comp)], loss_func_str, losst, num_comp)
	
	# all CPC output times, assuming first report is at 0 s through experiment
	times = np.arange(0, timehr[-1]*3600., 1./Hz)
	
	# empty array for holding corrected concentrations (# particles/cm3)
	Nwetn = np.zeros((len(times), Nwet.shape[1]))
	xnn = np.zeros((len(times), Nwet.shape[1]))
	
	# interpolate simulation output to instrument output frequency (# particles/cm3)
	# loop through size bins
	for sbi in range(num_sb-wall_on):
		Nwetn[:, sbi] = np.interp(times, timehr*3600., Nwet[:, sbi])
		xnn[:, sbi] = np.interp(times, timehr*3600., xn[:, sbi])
	
	Nwet = Nwetn # rename Nwet (# particles/cm3)
	xn = xnn # rename xn (um)
	
	# number of simulation outputs within the instrument response time
	rt_num = (times[1]-times[0])/delays[2]
	
	# if more than one output within response time, then loop through times to correct
	# for response time and any mixing of ages of particle
	if (rt_num > 1):
	
		# account for response time and mixing of particles of different ages
		[weight, weightt] = resp_time_func(3, delays, wfuncs)
	
		# empty array for holding corrected concentrations (# particles/cm3)
		Nwetn = np.zeros((Nwet.shape[0], Nwet.shape[1]))
		xnn = np.zeros((Nwet.shape[0], Nwet.shape[1]))
	
		for it in range(1, len(times)): # loop through times
			# number of time points to consider
			trel = (times >= (times[it]-delays[2]))*(times <= times[it])
			tsim = times[trel] # extract relevant time points (s)
			tsim = np.abs(tsim-tsim[-1]) # time difference with present (s)
			Nsim = Nwet[trel, :] # extract relevant number concentrations (# particles/cm3)
			xsim = xn[trel, :]
			# interpolate weights, use flip to align times
			weightn = np.flip(np.interp(np.flip(tsim), weightt, weight))
			# tile across size bins
			weightn = np.tile(weightn.reshape(-1, 1), [1, num_sb-wall_on])
			# corrected concentration
			Nwetn[it, :] = np.sum(Nsim*weightn, axis=0)
			xnn[it, :] = np.sum(xsim*weightn, axis=0)
		
		Nwet = Nwetn # rename Nwet
		xn = xnn # size bin radii (um)
	
	
	# account for size dependent detection efficiency below one
	# get detection efficiency as a function of particle size (nm)
	[Dp, ce] = count_eff_plot(3, 0, self, sdt)
	
	# empty array to hold detection efficiencies across times and simulation size bins
	# Dp is in um
	ce_t = np.zeros((len(times), xn.shape[1]))
	
	
	# loop through times
	for it in range(len(times)):
		# interpolate detection efficiency (fraction) to simulation size bin centres
		# Dp is in um
		ce_t[it, :] = np.interp(xn[it, :]*2., Dp, ce)
		
		# correct for upper size range of instrument, note conversion of
		# upper size from nm to um
		size_indx = (xn[it, :]*2. > max_size*1.e-3)
		Nwet[it, size_indx] = 0.
	
	Nwet = Nwet*ce_t # correct for detection efficiency
	
	
	# ------------
	# sum particle number concentration across size bins
	Nwet = Nwet.sum(axis=1)
	
	plt.ion() # show results to screen and turn on interactive mode
	

	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	ax0.plot(times/3600., Nwet, label = 'uncertainty mid-point')
	
	# plot vertical axis logarithmically
	ax0.set_yscale("log")
	
	# include uncertainty region, note conversion of uncertainty from percentage to fraction
	ax0.fill_between(times/3600., Nwet-Nwet*uncert/100., Nwet+Nwet*uncert/100., alpha=0.3, label = 'uncertainty bounds')
	
	# set tick format for vertical axis
	ax0.set_xlabel(r'Time through simulation (hours)', fontsize=14)
	ax0.set_ylabel('Total Number Concentration (#$\mathrm{particles\, cm^{-3}}$)', size = 14)
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which= 'both')
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which= 'both')
	ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
	ax0.set_title('Simulated total particle concentration convolved to represent \ncondensation particle counter (CPC) measurements')
	ax0.legend()
	
	return()
	# ------------
	
	
	
	Nwet = Nwet.sum(axis=1) # sum particle concentrations (# particles/cm3)
	
	
	
	# account for false background counts 
	# (minimum detectable particle concentration)  (# particles/cm3)
	Nwet[Nwet<cdt] = cdt
	
	# account for maximum particle concentration (# particles/cm3)
	Nwet[Nwet>max_dt] = max_dt
	
	if (caller == 0): # when called from gui
		plt.ion() # show results to screen and turn on interactive mode
	
	# plot temporal profile of total particle number concentration (# particles/cm3)
	# prepare figure -------------------------------------------
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
	ax0.plot(times/3600.0, Nwet, label = 'uncertainty mid-point')
	
	# plot vertical axis logarithmically
	ax0.set_yscale("log")
	
	# include uncertainty region, note conversion of uncertainty from percentage to fraction
	ax0.fill_between(times/3600., Nwet-Nwet*uncert/100., Nwet+Nwet*uncert/100., alpha=0.3, label = 'uncertainty bounds')
	
	# set tick format for vertical axis
	ax0.set_xlabel(r'Time through simulation (hours)', fontsize=14)
	ax0.set_ylabel('Total Number Concentration (#$\mathrm{particles\, cm^{-3}}$)', size = 14)
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which= 'both')
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which= 'both')
	ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
	ax0.set_title('Simulated total particle concentration convolved to represent \ncondensation particle counter (CPC) measurements')
	ax0.legend()
	if (caller == 2): # display when in test mode
		plt.show()
	
	return()

# function to plot detection efficiency as a function of particle size (% vs. nm)
def count_eff_plot(caller, dir_path, self, sdt):

	# note that this function assumes a logarithmic (in diameter space) sigmoid
	# for detection efficiency dependence on particle diameter, which is
	# observed by Niida et al. (1988) doi.org/10.1016/0021-8502(88)90188-7
	# and is also presented in TSI's documents (which references Niida et al. (1988): 
	# https://www.tsi.com/getmedia/952fc58b-ec23-4f18-86e4-58d23636e56d/SMPS-003appnote?ext=.pdf

	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0) or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# self - reference to GUI
	# sdt - particle diameter at 50 % detection efficiency (nm), 
	# 	width factor for detection efficiency dependence on particle size
	# --------------------------------------------------------------------------
	
	# prepare figure -------------------------------------------
	if (caller == 0): # if calling from gui
		plt.ion()
	
	if (caller != 3): # if called from something that does need a plot
		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	# --------------------------------------------------------------
	
	# particle size at 50 % detection efficiency (nm)
	b = sdt[0]
	# width factor for detection efficiency dependence on particle size
	c = sdt[1]
	
	# diameters (without transformation)
	Dp = np.linspace(-10, 10, 100)
	# sigmoid function, converting from fraction to percent
	ce = ((np.exp(Dp)/(np.exp(Dp)+1)))*1.e2
	
	# convert diameters so that sigmoid begins at the origin and ends at 1
	Dp = Dp/20.+0.5
		
	# apply stretch factor
	Dp = Dp*c
	
	# find diameter at 50 % point
	Dpcen = Dp[np.abs(ce-50.)==np.min(np.abs(ce-50.))]
	
	# apply the shift factor so that the 50 % point is at the specified value
	Dp += (np.log10(b)-Dpcen)
	
	# now map onto logarithmic space (nm)
	Dp = 10**Dp
	
	if (caller != 3): # if called from something that does not need a plot
		ax0.semilogx(Dp, ce, '-x') # state plotting points
		# set axis formats
		ax0.set_ylabel('Detection efficiency (%)', size = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.set_xlabel(r'Particle diameter (nm)', fontsize=14)
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
		ax0.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
	
	if (caller == 2): # if calling from test
		plt.show()
	
	# before returning convert diameter to um and ce to fraction
	Dp = Dp*1.e-3
	ce = ce/100.
	
	return(Dp, ce)

# function for plotting weighting of particles by age due to instrument response time
def resp_time_func(caller, delays, wfuncs):

	import datetime
	import cpc_response_eqs

	# inputs: --------------------------
	# caller - flag for the calling function
	# delays - the delays to when particles pass through counter
	# wfuncs - weighting functions
	# ------------------------------------

	# create new  file - will contain module response time weighting function
	f = open('PyCHAM/cpc_response_eqs.py', mode='w')

	f.write('\'\'\'solving the weighting of particles of different ages over response time of instrument\'\'\'\n')
	f.write('# module to estimate the weighting of particles of different ages during the response time of instrument to represent the mixing of particles of different ages due to differential flow prior to counter \n')
	f.write('# File Created at %s\n' %(datetime.datetime.now()))	
	f.write('\n')
	f.write('import numpy as np\n')
	f.write('\n')
	f.write('# function for weighting\n')
	f.write('def cpc_response(delays, wfuncs):\n')
	f.write('	\n')
	f.write('	# inputs: -----------------\n')
	f.write('	# ---------------------------\n')
	f.write('	\n')
	f.write('	# remember all times (s) \n')
	f.write('	t_all = np.zeros((100)) \n')
	f.write('	sht = %s # shortest delay\n' %(delays[0]))
	f.write('	resp_timepeak = %s # delay at peak weighting\n' %(delays[1]))
	f.write('	# time range for increasing weighting (s)\n')
	f.write('	t = np.arange(sht, (resp_timepeak), (resp_timepeak-sht)/50.)\n')
	f.write('	wpre = %s \n' %(wfuncs[0]))
	f.write('	t_all[0:50] = t[:] \n')
	f.write('	\n')
	f.write('	lot = %s # longest delay (s)\n' %(delays[2]))
	f.write('	# time range for decreasing weighting (s)\n')
	f.write('	t = np.arange((resp_timepeak+(lot-resp_timepeak)/50.), (lot), (lot-resp_timepeak)/51.)\n')
	f.write('	wpro = %s \n' %(wfuncs[1]))
	f.write('	t_all[50::] = t[:] \n')
	f.write('	\n')
	f.write('	# join weighting\n')
	f.write('	w = np.append(wpre, wpro)\n')
	f.write('	\n')
	f.write('	# integrate weight curve\n')
	f.write('	area = np.trapz(w, t_all)\n')
	f.write('	# normalise so that integral is one\n')
	f.write('	w = w/area\n')
	f.write('	\n')
	f.write('	return(w, t_all)')
	f.close() # close file
	
	# weighting function over this time
	importlib.reload(cpc_response_eqs)
	[w, t] = cpc_response_eqs.cpc_response(delays, wfuncs)
	
	
	# prepare figure -------------------------------------------
	if (caller == 0): # if calling from gui
		plt.ion()
	
	if (caller != 3): # if called from something that does need a plot
		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	# --------------------------------------------------------------
	
		# plot weighting dependency against response time
		ax0.plot(t, w)
		# plot details
		ax0.set_title('Weighting of simulated particle ages with instrument response time')
		ax0.set_ylabel('Weighting (fraction (0-1))', size = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.set_xlabel('Response time, with 0 s representing counting of particles \nsimultaneous with presence in atmosphere (s)', fontsize=14)
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')

	return(w, t)