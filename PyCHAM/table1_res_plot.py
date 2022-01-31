'''script to plot results for the comparison of modelled nucleation against observations'''
# using the output produced by fig11_mod_var.txt and fig11_simp_scheme.txt, 
# assumed calling from the PyCHAM home folder - note that 24 size bin used for
# moving-centre and 12 size bins used for full-moving, with the save file name
# and size strucure inputs changed accordingly also

# import required modules
import  matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import numpy as np
import scipy.constants as si
import sys
import os
# ensure modules can be seen 
# (assumes calling from the home folder, but also works from the GMD paper Results 
# folder)
sys.path.append(str(os.getcwd() + '/PyCHAM'))
import retr_out

# ----------------------------------------------------------------------------------------
# function for doing colorbar tick labels in standard notation
def fmt(x, pos):
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	return r'${} \times 10^{{{}}}$'.format(a, b)

# ----------------------------------------------------------------------------------------
# create figure to plot results

fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(12, 6.5), sharex='all')
fig.subplots_adjust(hspace = 0.05)


# ------------------------------------------------------------------------------
# results
cwd = os.getcwd() # get current working directory

# no wall losses (wall_on = 0) # nuc_vsobs_output_mc_24sb_walloff
# particle loss to wall, no gas loss to wall
# inflectDp = 1.e-6
#Grad_pre_inflect = 1.
#Grad_post_inflect = 1.
#Rate_at_inflect = 1.e-4

# better particle loss to wall, no gas loss to wall
# inflectDp = 1.e-6
#Grad_pre_inflect = 1.
#Grad_post_inflect = 1.
#Rate_at_inflect = 6.e-6

# high gas loss to wall, no particle loss
# Kw = 1.e-2
# Cw = 1.0
# inflectDp = 1.e-6
#Grad_pre_inflect = 0.
#Grad_post_inflect = 0.
#Rate_at_inflect = 0.

# better gas loss to wall, no particle loss
#mass_trans_coeff = 1.e-4
#eff_abs_wall_massC = 1.e0
#inflectDp = 1.e-6
#Grad_pre_inflect = 0.
#Grad_post_inflect = 0.
#Rate_at_inflect = 0.

# best particle and gas loss to wall
#mass_trans_coeff = 1.e-6
#eff_abs_wall_massC = 1.e0
#inflectDp = 1.e-6
#Grad_pre_inflect = 1.
#Grad_post_inflect = 1.
#Rate_at_inflect = 6.e-6

output_by_sim = str(cwd + '/PyCHAM/output/fig11_full_scheme/nuc_vsobs_output_mc_24sb_gpm')
# required outputs
(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, xfm, time_array, comp_names, 
	_, N, _, y_MV, _, wall_on, space_mode, indx_plot, comp0, 
	yrec_p2w, PsatPa, OC, H2Oi, seedi, siz_str, cham_env, group_indx, 
	tot_in_res) = retr_out.retr_out(output_by_sim)

	
dlog10D = np.log10(rbou_rec[:, 1::]*2.)-np.log10(rbou_rec[:, 0:-1]*2.)
dNdD = Ndry/dlog10D # normalised number size distribution (# particles/cm3 (air))

# observation part ---------------------------------------------------------------
import openpyxl # for opening xlsx file

# if calling from the PyCHAM home folder
wb = openpyxl.load_workbook(str(cwd+'/PyCHAM/output/GMD_paper_plotting_scripts/obs_fig11.xlsx'))
wb = wb['20190628_Dark_limonene_LowNOx'] # take just the required sheet


sr = 14 # starting row for desired time since O3 injection

obst = np.zeros((46-sr, 1)) # empty array for observation times (s)
obsN = np.zeros((46-sr, 23)) # empty array for particle number concentration (# particles/cm3 (air)) 
obsDp= np.zeros((1, 23)) # empty array for diameters (nm)

# get data
rcn = 0 # row count
rc = 0 # row count
for r in wb.iter_rows(): # row loop
	if (rcn > sr): # desired times (s)
		
		obstn = r[1].value
		obstn = str(obstn) # ensure time is string
		if (isinstance(obstn, str)): # convert string type times to fraction of day
			splt = obstn.split(':')
			obstn = float(splt[0])*(1/24)+float(splt[1])*(1/(24*60))+float(splt[2])*(1/(24*60*60))
		
		obst[rc, 0] = obstn*(3600.*24.) # convert to seconds

	cc = 0 # column count
	for c in r[2::]: # column loop
		if (rcn == 0):
			obsDp[0, cc] = c.value
			
		if (rcn > sr):
			obsN[rc, cc] = c.value
		cc += 1 # column count

	if (rcn > sr):
		rc += 1 # row count
	rcn += 1

# normalised number concentration - note values in excel already normalised by d(log10(Dp))
dNdDp = obsN

# consider only observations within the time period of interest, note that
# for Table 1 of the GMD paper, the first hour of the experiment was considered
t_indx_min = int(sum((obst/3600.) < 2.))
t_indx_max = int(sum((obst/3600.) <= 5.))

obst = obst[t_indx_min:t_indx_max, :]
dNdDp = dNdDp[t_indx_min:t_indx_max, :]

# make simulation and observation results comparable by mapping simulation results onto 
# observation times and sizes (the latter by interpolation)
dNdD_interp = np.zeros((dNdDp.shape[0], dNdDp.shape[1]))

# observation diameters (um)
obss = obsDp*1.e-3

for it in range(len(obst)-1): # loop through observation time steps
	# observation time (hours)
	obstn = obst[it]/3600.

	# index of simulation result at this time
	simi = np.where(time_array==((time_array[time_array>=obstn])[0]))[0]
	
	# fill model matrix - convert simulation radii to diameter
	dNdD_interp[it, :] = np.interp(obss[0, :], np.squeeze(xfm[simi, :]*2.), np.squeeze(dNdD[simi, :]))


# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
colors = [(0.60, 0.0, 0.70), (0, 0, 1), (0, 1.0, 1.0), (0, 1.0, 0.0), (1.0, 1.0, 0.0), (1.0, 0.0, 0.0)]  # R -> G -> B
n_bin = 100  # discretizes the colormap interpolation into bins
cmap_name = 'my_list'
# Create the colormap
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N = n_bin)
	
# set contour levels
levels = (MaxNLocator(nbins = 100).tick_values(np.min(dNdDp), np.max(dNdDp)))
	
# associate colours and contour levels
norm1 = BoundaryNorm(levels, ncolors = cm.N, clip=True)
p1 = ax0.contour(obst[:, 0]/3600., obsDp[0, :], dNdD_interp.transpose(), cmap=cm, norm=norm1)
ax0.set_yscale("log")
ax0.set_ylim(1.e1, 4.e2) # set vertical axis limits
p2 = ax1.contour(obst[:, 0]/3600., obsDp[0, :], dNdDp.transpose(), cmap=cm, norm=norm1)
ax1.set_yscale("log")
ax1.set_ylim(1.e1, 4.e2) # set vertical axis limits

# colorbar
cax = plt.axes([0.30, 0.92, 0.40, 0.01])
cb = plt.colorbar(p1, cax = cax, ticks=[0., 6.0e4, 1.2e5, 1.8e5], format = ticker.FuncFormatter(fmt), orientation = 'horizontal')

plt.show()

res = np.zeros((obst.shape[0], obsDp.shape[1]))
# loop through times and size bins to get residuals
for it in range(obst.shape[0]):
	for isb in range(obsDp.shape[1]):
		if ((dNdD_interp[it, isb] >= 0.) and (dNdDp[it, isb] >= 0.)):
			res[it, isb] = (np.abs(dNdD_interp[it, isb]-dNdDp[it, isb]))

# denominator is the sum of number concentrations from all size bins at all times for each simulation
# results to use
dNdD_interp_sp = dNdD_interp[0:obst.shape[0], :]
den_cnt = (sum(dNdD_interp_sp[dNdD_interp_sp>0.]))
# results to use
dNdDp_sp = dNdDp[0:obst.shape[0], :]
den_cnt += (sum(dNdDp_sp[dNdDp_sp>0.]))

# sum over times and size bins to get percentage error
resT = ((sum(sum(res)))/(den_cnt))*100.

print(resT)