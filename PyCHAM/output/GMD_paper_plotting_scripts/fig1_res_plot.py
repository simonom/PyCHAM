'''Code to plot results for limonene oxidation example (Fig. 1 of GMD paper)'''
# aim is to exemplify the coupled integration of gas-phase chemistry and partitioning to
# particles and wall, call from the PyCHAM home directory if inside PyCHAM, or, if 
# inside the GMD paper folder, call from the Results folder

# For the chemical scheme file use: limonene MCM PRAM sheme (fig01_scheme.txt)
# For model variable inputs use: GMD_paper/Results/fig01_mod_var.txt
# Which, for the NOx present case should include:
# C0 = 10.0, 0.0, 0.0, 22.0, 0.0, 500.0e3
# Comp0 = LIMONENE, N2O5, NO3, NO2, O3, CO
# injectt = 5400.0, 14400.0
# Compt = LIMONENE, N2O5, NO3, NO2, O3
# Ct = 0.0, 10.0; 0.0, 0.0; 0.0, 0.0; 8.0, 19.0; 38.0, 45.0
# Use the xml file in the inputs folder as this has been updated for PRAM components

# import required modules
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as si
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import os
from retrieve_PyCHAM_outputs import retrieve_outputs as retr

# testing plot ----------------------------------

# import results

# empty dictionaries to hold results
num_sb_dict = {}
num_speci_dict = {}
Cfac_dict = {}
y_dict = {}
N_dict = {}
sbb_dict = {}
x_dict = {}
thr_dict = {}
comp_names = {}
mw_dict = {}


cwd = os.getcwd() # get current working directory
try: # in case calling from inside the GMD paper folder results section
	(num_sb_dict['num_sb0'], num_speci_dict['num_speci0'], Cfac_dict['Cfac0'], y_dict['y0'], 		N_dict['N0'], sbb_dict['sbb0'], x_dict['x0'], thr_dict['thr0'], comp_names['cn0'], mw_dict['mw0'], 		_, _, _, _) = retr(str(cwd + '/fig01_data/'))
except: # in case calling from PyCHAM home folder
	(num_sb_dict['num_sb0'], num_speci_dict['num_speci0'], Cfac_dict['Cfac0'], y_dict['y0'], 		N_dict['N0'], sbb_dict['sbb0'], x_dict['x0'], thr_dict['thr0'], comp_names['cn0'], mw_dict['mw0'], 		_, _, _, _) = retr(str(cwd + 'PyCHAM/output/GMD_paper_plotting_scripts/fig01_data/'))

# ----------------------------------------------------------------------------------------
# rename variables for use below
num_speci = num_speci_dict['num_speci0']
num_sb = num_sb_dict['num_sb0']

# ----------------------------------------------------------------------------------------
# get indices of components with a nitrate group
org_nitr_ind = np.empty(0, dtype='int')
inorg_nitr_ind = np.empty(0, dtype='int')
PyCHAM_names = comp_names['cn0']
for i in range(len(PyCHAM_names)):
	if len(PyCHAM_names[i])<3:
		continue
	# don't include inorganic nitrate
	if PyCHAM_names[i] == 'HNO3' or PyCHAM_names[i] == 'N2O5':
		inorg_nitr_ind = np.append(inorg_nitr_ind, i)
		continue
	else:
		for i2 in range(len(PyCHAM_names[i])):
			if PyCHAM_names[i][i2:i2+3] == 'NO3':
				org_nitr_ind = np.append(org_nitr_ind, i)

# ----------------------------------------------------------------------------------------
# prepare plot
# note, we have to use the gridspec approach because the parasite axis on the lower plot
# means that sharex won't bring the x axis of both plots in line
fig = plt.figure(figsize = (13,7))
gs = fig.add_gridspec(2, 15)

ax0 = fig.add_subplot(gs[0, :9])
ax1 = fig.add_subplot(gs[1, :])

# parasite axis for particle-phase plot
par1 = ax1.twinx() # first parasite axis
par2 = ax1.twinx() # second parasite axis

def make_patch_spines_invisible(ax):
	ax.set_frame_on(True)
	ax.patch.set_visible(False)
	for sp in ax.spines.values():
		sp.set_visible(False)

# Offset the right spine of par2.  The ticks and label have already been
# placed on the right by twinx above.
par2.spines["right"].set_position(("axes", 1.2))
# Having been created by twinx, par2 has its frame off, so the line of its
# detached spine is invisible.  First, activate the frame but make the patch
# and spines invisible.
make_patch_spines_invisible(par2)
# Second, show the right spine.
par2.spines["right"].set_visible(True)
# ----------------------------------------------------------------------------------------
# gas-phase concentrations

# get index of components of interest
targ_comps = ['O3', 'NO2', 'LIMONENE', 'N2O5', 'NO3']
comp_index = []
for i in range(len(targ_comps)):
	index = PyCHAM_names.index(targ_comps[i])
	comp_index.append(index)

y = y_dict['y0']

# plot gas-phase concentrations of O3 and NO2
ax0.semilogy(thr_dict['thr0'], y[:, comp_index[0]], label = 'O3')
ax0.semilogy(thr_dict['thr0'], y[:, comp_index[1]], label = 'NO2')
ax0.semilogy(thr_dict['thr0'], y[:, comp_index[2]], label = 'LIMONENE')
# ax0.semilogy(thr_dict['thr0'], y[:, comp_index[3]], label = 'N2O5')
# ax0.semilogy(thr_dict['thr0'], y[:, comp_index[4]], label = 'NO3')
ax0.set_ylim(1.0e-1, 1.0e2)
ax0.set_ylabel(r'Gas-phase concentration (ppb)', fontsize=12)
ax0.yaxis.set_tick_params(labelsize=14, direction = 'in', which = 'both')
#ax0.set_xlabel(r'Time through simulation (hours)', fontsize=12)
ax0.xaxis.set_tick_params(direction = 'in', which = 'both')
ax0.xaxis.set_ticklabels([])


ax0.plot([1.5, 1.5], [0.0, 500.0], color='black', linewidth=1, linestyle='dashed')
ax0.text(x=1.2, y=1, s='injection 1', size=12, rotation=90, color='black')
ax0.plot([4.0, 4.0], [0.0, 500.0], color='black', linewidth=1, linestyle='dashed')
ax0.text(x=3.7, y=1, s='injection 2', size=12, rotation=90, color='black')

ax0.legend(fontsize=10)
ax0.text(x=(thr_dict['thr0'])[0]-0.2, y=1.6e2, s='(a)', size=12)

# ----------------------------------------------------------------------------------------
# particle-phase properties

# number size distribution (dN/dlog10(Dp))
sbb = sbb_dict['sbb0']
N = N_dict['N0']

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

# repeat size bin diameter widths over times
dlog10D = np.repeat(dlog10D, N.shape[0], axis=0)

# prepare number size distribution contours (/cc (air))
dNdlog10D = np.zeros((N.shape[0], N.shape[1]))
dNdlog10D = N/dlog10D

# two-size bin moving average
dNdlog10D[:, 1::] = (dNdlog10D[:, 0:-1]+dNdlog10D[:, 1::])/2.0

# transpose number concentration results, so time on x axis and diameter on y
dNdlog10D = dNdlog10D.transpose()

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

# contour plot with times along x axis and particle diameters along y axis
p1 = ax1.pcolormesh(thr_dict['thr0'], (sbb*2*1e3), z[:, :], cmap=cm, norm=norm1)


ax1.set_ylabel('Diameter ($D_p$, nm)', size=12)
ax1.xaxis.set_tick_params(labelsize=14, direction = 'in', which = 'both')
ax1.yaxis.set_tick_params(labelsize=14, direction = 'in', which = 'both')
ax1.set_xlabel(r'Time through simulation (hours)', fontsize=12)
ax1.set_yscale('log')
ax1.set_ylim((sbb[1]*2*1e3/0.7), (sbb[-1]*2*1e3))
# set y axis to standard notation
#ax1.ticklabel_format(axis='y', style='sci')

# function for doing colorbar tick labels in standard notation
def fmt(x, pos):
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	return r'${} \times 10^{{{}}}$'.format(a, b)

cb = plt.colorbar(p1, format=ticker.FuncFormatter(fmt), pad=0.25)
cb.ax.tick_params(labelsize=14) 

# colour bar label
cb.set_label('d$N$ ($\#\,\mathrm{cm^{-3}}$)/dlog$\mathrm{_{10}}$($D_p$)', size=12, rotation=270, labelpad=20)

ax1.text(x=(thr_dict['thr0'])[0]-0.1, y=(sbb[-1]*2*1e3)*1.05, s='(b)', size=12)

# set line and label denoting injections
ax1.plot([1.5, 1.5], [0.0, 500.0], color='white', linewidth=1, linestyle='dashed')
ax1.text(x=1.2, y=40, s='injection 1', size=12, rotation=90, color='white')
ax1.plot([4.0, 4.0], [0.0, 500.0], color='white', linewidth=1, linestyle='dashed')
ax1.text(x=3.7, y=150, s='injection 2', size=12, rotation=90, color='white')

# ----------------------------------------------------------------------------------------
# organic aerosol mass concentration in particles (ug/m3)
final_i = 0

y_mw = mw_dict['mw0']

# check whether water and/or core is present
if PyCHAM_names[-2] == 'H2O': # if both present
	final_i = 2
if PyCHAM_names[-1] == 'H2O': # if just water
	final_i = 1

SOAvst = np.zeros((y.shape[0], 1))
NO3vst = np.zeros((y.shape[0], 1))
inorg_NO3vst = np.zeros((y.shape[0], 1))

for i in range(1, num_sb):

	# isolate concentrations in this size bin (molecules/cc (air))
	y_now = y[:, num_speci*i:num_speci*(i+1)-final_i]
	
	# calculate SOA (*1.0E-12 to convert from g/cc (air) to ug/m3 (air) assuming a density
	# of 1.0 g/cc)
	SOAvst[:, 0] += (((y_now/si.N_A)*y_mw[0:-final_i]*1.0e12).sum(axis=1))
	# sum of inorganic nitrates
	for i2 in inorg_nitr_ind:
		inorg_NO3vst[:, 0] += ((y_now[:, i2]/si.N_A)*y_mw[i2]*1.0e12) 	
	# get just the sum concentrations of organic nitrate components in this size bin 
	for i2 in org_nitr_ind:
		NO3vst[:, 0] += ((y_now[:, i2]/si.N_A)*y_mw[i2]*1.0e12)

SOAvst[:, 0] -=  inorg_NO3vst[:, 0] # remove inorganic nitrate contribution from SOA

# mass concentration of organics on wall
# (*1.0E-12 to convert from g/cc (air) to ug/m3 (air) assuming a density
# of 1.0 g/cc)
Owall = (((y[:, num_speci*num_sb:num_speci*(num_sb+1)-final_i]/si.N_A)*y_mw[0:-final_i]*1.0e12).sum(axis=1))

# 
# # log10 of maximum in SOA
# SOAmax = int(np.log10(max(SOAvst)))
# Owallmax = int(np.log10(max(Owall)))
# SOAmax = max(SOAmax, Owallmax)
# # transform SOA so no standard notation required
# SOAvst = SOAvst/(10**(SOAmax))
# # transform Owall so no standard notation required
# Owall = Owall/(10**(SOAmax))
# # transform NO3vst so no standard notation required
# NO3vst = NO3vst/(10**(SOAmax))
# # transform inorganic nitrates so no standard notation required
# inorg_NO3vst = inorg_NO3vst/(10**(SOAmax))

p5, = par2.semilogy(thr_dict['thr0'], SOAvst, '--k', label = 'SOPM')
p6, = par2.semilogy(thr_dict['thr0'], Owall, '-.k', label = 'wall deposit')
p7, = par2.semilogy(thr_dict['thr0'], NO3vst, ':k', label = 'particle organic nitrate')
p8, = par2.semilogy(thr_dict['thr0'], inorg_NO3vst, '.k', label = 'particle inorganic nitrate')
par2.set_ylim(5.0e-1, 1.0e2)
# par2.set_ylabel(str('[organics]x ' + str(10**(SOAmax)) + ' ($\mathrm{\mu g\, m^{-3}})$'), rotation=270, size=12, labelpad=25)
par2.set_ylabel(str('Particulate & wall concentration' + ' ($\mathrm{\mu g\, m^{-3}})$'), rotation=270, size=12, labelpad=25)

par2.yaxis.set_tick_params(labelsize=12, direction = 'in', which = 'both')

# ----------------------------------------------------------------------------------------
# total particle number concentration

p3, = par1.plot(thr_dict['thr0'], N.sum(axis=1), '-k', label = '$N$')
# set tick format for vertical axis
par1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e'))
par1.set_ylabel('$N$ ($\#\,\mathrm{cm^{-3}}$)', rotation=270, size=12, labelpad=25)
par1.yaxis.set_tick_params(labelsize=12, direction = 'in', which = 'both')

plt.legend(fontsize=12, handles=[p3, p5, p6, p7, p8] ,loc=2, framealpha=0.25)
fig.savefig('fig01.png')
plt.show()

