'''module to test the moving-centre method of PyCHAM against benchmarks'''
# benchmark simulations are Jacobson (2005) and Zhang et al (1999)
# can be called from the PyCHAM home folder or GMD paper Results folder
print('The script to create plots for assessment of size structures; note needs to be called from either the PyCHAM home folder or GMD paper Results folder')

import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import os
import sys
# ensure modules can be seen 
sys.path.append(str(os.getcwd() + '/PyCHAM'))
import retr_out
from scipy import stats # Import the scipy.stats module
import numpy as np

# associated files --------------------------
# Use the fig06_scheme for the chemical scheme and the model variables provided below,
# possibly input to the fig06_mod_var.txt file, along with the example xml file.

# ----------------------------------------------------------------------------------------
# when comparing against Fig. 13.8 of Jacobson, note that preferential growth of larger
# particles over smaller ones is a function of the Kelvin effect, for which the density 
# and molar weight of the seed material is important, and of the supersaturation of water
# since it needs to be only slightly supersaturated for larger particles to reduce to
# subsaturated due to condensation, therefore the specified vapour pressure and injection
# concentration of water must be accurate and precise   

# model variables input file for Jacobson, note the size_structure will change:

#res_file_name = test_Jacobs_16sb_mc
#total_model_time = 600.
#update_step = 1.
#recording_time_step = 60.
#size_structure = 0
#number_size_bins = 16
#lower_part_size = 4.e-2
#upper_part_size = 5.e1
#space_mode = log
#coag_on = 1
#mass_trans_coeff = 0.
#eff_abs_wall_massC = 0.
#temperature = 300.00
#tempt = 0.
#p_init = 98000
#rh = 0.
#Ct = 3.6154e7, 3.6154e7, 3.6154e7, 3.6154e7, 3.6154e7, 3.6154e7, 3.6154e7, 3.6154e7, 3.6154e7, 3.6154e7
#Compt = O
#injectt = 0., 60., 120., 180., 240., 300., 360., 420., 480., 540.
#pconct = 0.
#pconc = 5.e3 : 9.5e2 : 4.0
#mean_rad = 0.070 : 0.2 : 1.5
#std = 1.24 : 1.40 : 1.70
#vol_Comp = O
#volP = 3536.01
#seed_name = core
#seed_mw = 350.
#seed_dens = 0.9
#core_diss = 1.
#chem_scheme_markers = %, RO2, +, , , ;, +, ;, , %, :, ;

# when testing against Fig. 3 of Zhang et al. (1999): doi.org/10.1080/027868299304039
# remember that the density and molecular weight of sulphuric acid are automatically set
# by PyCHAM and that these determine the resulting volumetric condensation rate

# inputs for Zhang (note siz_structure will vary)

#res_file_name = test_Zhang_100sb_mc
#total_model_time = 43200.
#update_step = 120.
#recording_time_step = 120.
#size_structure = 0
#number_size_bins = 100
#lower_part_size = 5.e-4
#upper_part_size = 1.e1
#space_mode = log
#coag_on = 1
#mass_trans_coeff = 0.
#eff_abs_wall_massC = 0.
#temperature = 273.15
#tempt = 0.0
#p_init = 98000
#rh = 0.000
#const_comp = 
#const_infl = SA
#const_infl_t = 0.
#Cinfl = 2.65e-5
#vol_Comp = SA
#volP = 0.
#pconct = 0.
#pconc = 2350. : 3600. : 3.95
#mean_rad = 0.0205 : 0.0611: 0.87
#std = 1.2 : 1.8: 2.2
#seed_name = core
#seed_mw = 200.
#seed_dens = 1.
#core_diss = 1.
#light_time = 0., 43200. 
#light_status = 0, 0
#chem_scheme_markers = %, RO2, +, , , ;, +, ;, , %, :, ;



# ----------------------------------------------------------------------------------------
# create figure to plot results
fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(10, 6.5))
fig.subplots_adjust(right=0.9, hspace=0.3)

par0 = ax0.twinx() # first parasite axis
par1 = par0.twiny() # first parasite axis

# ----------------------------------------------------------------------------------------
# Zhang et al. (1999) plot

# open full-moving
cwd = os.getcwd() # get current working directory
try: # if calling from the PyCHAM home directory
	output_by_sim = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig06_data/test_Zhang_100sb_fm')
	# required outputs from full-moving
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, xfm, time_array, comp_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)
except:
	output_by_sim = str(cwd + '/fig06_data/test_Zhang_100sb_fm')
	# required outputs from full-moving
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, xfm, time_array, comp_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)


import scipy.constants as si
Dpbb = rbou_rec*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[:, 1::]-(np.log10(Dpbb))[:, 0:-1] # difference in the log10 of bin bounds

# volume concentration of particles per size bin per time step (um3/cc (air))
dV = ((4.0/3.0)*np.pi*xfm**3.0)*N


# volumes normalised by size bin width (um3/cm)
dv0fm = dV[0, :]/dDp[0, :] # at initial time
dv1fm = dV[-1, :]/dDp[-1, :] # at final time step
dn0fm = N[0, :]/dDp[0, :] # initial time, number distribution
dn1fm = N[-1, :]/dDp[-1, :] # final time, number distribution

# open saved files
try:
	output_by_sim = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig06_data/test_Zhang_100sb_mc')
	# required outputs from full-moving
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, time_array, comp_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)
except:
	output_by_sim = str(cwd + '/fig06_data/test_Zhang_100sb_mc')
	# required outputs from full-moving
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, time_array, comp_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)


Dpbb = rbou_rec*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[:, 1::]-(np.log10(Dpbb))[:, 0:-1] # difference in the log10 of bin bounds

# volume concentration of particles per size bin per time step (um3/cc (air))
dV = ((4.0/3.0)*np.pi*x**3.0)*N

# volumes normalised by size bin width (um3/cm)
dv0mc = dV[0, :]/dDp[0, :] # at initial time
dv1mc = dV[-1, :]/dDp[-1, :] # at final time step
dn0mc = N[0, :]/dDp[0, :] # initial time, number distribution
dn1mc = N[-1, :]/dDp[-1, :] # final time, number distribution

n_ma = 3 # number of points to combine for moving average
dvn = np.zeros((len(dv1mc)-(n_ma-1)))
dvnx = np.zeros((len(dv1mc)-(n_ma-1)))
dnn = np.zeros((len(dn1mc)-(n_ma-1)))
dnnx = np.zeros((len(dn1mc)-(n_ma-1)))
for i in range(len(dvn)):
	dvn[i] = sum(dv1mc[i:i+n_ma])/n_ma
	dvnx[i] = sum(x[-1, i:i+n_ma])/n_ma
	dnn[i] = sum(dn1mc[i:i+n_ma])/n_ma
	dnnx[i] = sum(x[-1, i:i+n_ma])/n_ma

# plotting part ----------------------------------------------------------------------


# volume distribution
ax0.semilogx(xfm[0, :]*2., dv0fm, '-k', label='initial') # initial
ax0.semilogx(xfm[-1, :]*2., dv1fm, '--k', label='full-moving')# full-moving
ax0.semilogx(dvnx*2., dvn, ':k', label='moving-centre') # moving-centre

# number distribution
p3, = par1.semilogx(xfm[0, xfm[0, :]>5.e-3]*2., dn0fm[xfm[0, :]>5.e-3], '-b') # initial
p4, = par1.semilogx(xfm[-1, xfm[-1, :]>5.e-3]*2, dn1fm[xfm[-1, :]>5.e-3], '--b') # full-moving
p4, = par1.semilogx(dnnx[dnnx>5.e-3]*2., dnn[dnnx>5.e-3], ':b') # moving-centre

# standard axis formatting
ax0.set_ylabel(r'd$\nu$ ($\rm{\mu m^{3}\, cm^{-3}})/dlog_{10}{\mathit{D}\rm_{p}}$', fontsize=14)
ax0.set_xlabel(r'${\mathit{D}\rm_{p}}$ ($\mu$m)', fontsize=14)

ax0.set_xlim([1.0e-3, 1.0e1])
ax0.set_xscale('log')
ax0.text(x=6.0e-4, y=34, s='(a)', size = 14)
ax0.yaxis.set_tick_params(labelsize = 14, direction='in', which='both')
ax0.xaxis.set_tick_params(labelsize = 14, direction='in', which='both')
ax0.legend(fontsize = 12)

# parasite axis formatting
par0.set_ylabel(r'dN ($\rm{cm^{-3}})/dlog_{10}{\mathit{D}\rm_{p}}$', size=14, rotation=270, 
labelpad=20, color='b') # vertical axis label
par1.set_xlabel(r'${\mathit{D}\rm_{p}}$ ($\mu$m)', fontsize=14, color='b')
par1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0e')) # set tick format
par0.yaxis.set_tick_params(labelsize=14, colors='b', direction='in', which='both')
par1.spines['right'].set_color('b')
par1.xaxis.set_tick_params(labelsize=14, colors='b', direction='in', which='both')
par1.spines['top'].set_color('b')
par0.set_ylim([-1.e3, 3.8e4])

# Jacobson 2005 part ----------------------------------------------------------------------------
# open full-moving
cwd = os.getcwd() # get current working directory
try:
	output_by_sim = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig06_data/test_Jacobs_16sb_fm')
	# required outputs from full-moving
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, xfm, time_array, comp_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)
except:
	output_by_sim = str(cwd + '/fig06_data/test_Jacobs_16sb_fm')
	# required outputs from full-moving
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, xfm, time_array, comp_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)


import scipy.constants as si
Dpbb = rbou_rec*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[:, 1::]-(np.log10(Dpbb))[:, 0:-1] # difference in the log10 of bin bounds

# volume concentration of particles per size bin per time step (um3/cc (air))
dV = ((4.0/3.0)*np.pi*xfm**3.0)*N

# volumes normalised by size bin width (um3/cm3)
dv0fm = dV[0, :]/dDp[0, :] # at initial time
dv1fm = dV[-1, :]/dDp[-1, :] # at final time step

# number concentrations normalised by size bin width (/cm3)
dn0fm = N[0, :]/dDp[0, :] # at initial time


# open moving-centre
cwd = os.getcwd() # get current working directory
try:
	output_by_sim = str(cwd + '/PyCHAM/output/GMD_paper_plotting_scripts/fig06_data/test_Jacobs_16sb_mc')
	# required outputs from full-moving
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, xmc, time_array, comp_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)
except:
	output_by_sim = str(cwd + '/fig06_data/test_Jacobs_16sb_mc')
	# required outputs from full-moving
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, xmc, time_array, comp_names, 
		_, N, _, y_MV, _, wall_on, space_mode) = retr_out.retr_out(output_by_sim)


Dpbb = rbou_rec*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[:, 1::]-(np.log10(Dpbb))[:, 0:-1] # difference in the log10 of bin bounds

# volume concentration of particles per size bin per time step (um3/cc (air))
dV = ((4.0/3.0)*np.pi*xmc**3.0)*N

# volumes normalised by size bin width (um3/cm3)
dv0mc = dV[0, :]/dDp[0, :] # at initial time
dv1mc = dV[-1, :]/dDp[-1, :] # at final time step

# number concentrations normalised by size bin width (/cm3)
dn0mc = N[0, :]/dDp[0, :] # at initial time
dn1mc = N[-1, :]/dDp[-1, :] # at initial time

n_ma = 2 # number of points to combine for moving average
dvn = np.zeros((len(dv1mc)-(n_ma-1)))
dvnx = np.zeros((len(dv1mc)-(n_ma-1)))
dnn = np.zeros((len(dn1mc)-(n_ma-1)))
dnnx = np.zeros((len(dn1mc)-(n_ma-1)))
for i in range(len(dvn)):
	dvn[i] = sum(dv1mc[i:i+n_ma])/n_ma
	dvnx[i] = sum(xmc[-1, i:i+n_ma])/n_ma
	dnn[i] = sum(dn1mc[i:i+n_ma])/n_ma
	dnnx[i] = sum(xmc[-1, i:i+n_ma])/n_ma

# Jacobson (2005) result
# diameters (um)
Jacob_x = np.array((0.1, 0.11, 0.2, 0.28, 0.4, 0.41, 16.0, 16.1, 20.0, 21.0, 22.0, 26.0, 30.0, 70.0, 100.0))
Jacob_y = np.array((0.4, 0.4, 1.25, 110.0, 300.0, 1.0e-10, 1.0e-10, 3.0e7, 2.0e6, 3.0e5, 4.0e5, 1.1e4, 1.1e3, 2.0, 0.3))


ax1.loglog(xfm[0, :]*2., dv0fm, '-+k', label = 'initial')
ax1.loglog(xfm[-1, :]*2., dv1fm, 'xk', label = 'full-moving')
ax1.loglog(dvnx*2., dvn, '*k', label = 'moving-centre')
ax1.loglog(Jacob_x, Jacob_y, '--k', label='Jacobson (2005) Fig. 13.4')
ax1.legend(fontsize = 12)
ax1.text(x=0.04, y=4.e7, s='(b)', size=14)
ax1.set_ylim([1.e-1, 3.e7])
ax1.yaxis.set_tick_params(labelsize = 14, direction='in', which='both')
ax1.xaxis.set_tick_params(labelsize = 14, direction='in', which='both')

ax1.set_ylabel(r'd$\nu$ ($\rm{\mu m^{3}\, cm^{-3}})/dlog_{10}{\mathit{D}\rm_{p}}$', fontsize=14)
ax1.set_xlabel(r'${\mathit{D}\rm_{p}}$ ($\mu$m)', fontsize=14)


fig.savefig('fig06.png')

plt.show()
