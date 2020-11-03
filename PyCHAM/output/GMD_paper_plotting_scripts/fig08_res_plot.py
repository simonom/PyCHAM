'''code to plot results from spatial and temporal tests of coagulation sensitivity for GMD paper'''
# for creating the coagulation sensitivity plot, currently can only be called from the GMD 
# paper Results folder

import numpy as np
import matplotlib.pyplot as plt
import os

# Figure 13.4 Jacobson (2005)

# use fig08_scheme.txt as the input chemical scheme from Results folder in GMD_paper, setting the one 
# reaction rate coefficient to zero for the no reaction cases and to 5.6e-17 for the with reaction case

# total time for simulation: 21600 s
# seed_name = core
# seed_mw = 132.00
# seed_dens = 1.77

# for all number of size bins we use the radius bounds (um) input:
# lower_part_size = 0.005
# upper_part_size = 6
# space_mode = log
# 
# # in the case of 8 size bins we use the following input:
# pconc = 45000.0, 13000.0, 2500.0, 3200.0, 200.0, 3.5e-1, 5.0e-1, 6.0e-2
# # in the case of 32 size bins we use the following input:
# pconc = 20000.0, 15000.0, 11000.0, 8000.0, 5000.0, 3000.0, 1700.0, 1100.0, 700.0, 550.0, 700.0, 950.0, 1090.0, 960.0, 830.0, 300.0, 80.0, 120.0, 160.0, 80.0, 20.0, 2.0, 1.0e-1, 1.5e-1, 1.5e-1, 1.4e-1, 1.3e-1, 9.0e-2, 5.0e-2, 2.0e-2, 8.0e-3, 2.0e-3
# # for 128 size bins, used the following input:
# pconc = 2800.0, 2500.0, 2200.0, 1950.0, 1700.0, 1500.0, 1300.0, 1150.0, 1000.0, 900.0, 800.0, 720.0, 640.0, 580.0, 520.0, 460.0, 430.0, 400.0, 370.0, 345.0, 320.0, 295.0, 275.0, 255.0, 235.0, 220.0, 205.0, 190.0, 180.0, 170.0, 164.0, 162.0, 158.0, 154.0, 150.0, 158.0, 166.0, 174.0, 182.0, 190.0, 198.0, 206.0, 214.0, 220.0, 228.0, 235.0, 226.0, 217.0, 209.0, 201.0, 194.0, 187.0, 170.0, 146.0, 124.0, 104.0, 86.0, 70.0, 56.0, 43.0, 31.0, 20.0, 22.0, 24.0, 26.0, 28.0, 32.0, 28.0, 24.0, 20.0, 17.0, 14.0, 11.0, 9.0, 7.0, 5.0, 4.0, 3.0e0, 1.5e0, 7.25e-1, 3.6e-1, 1.0e-1, 1.53e-2, 1.62e-2, 1.70e-2, 1.78e-2, 1.85e-2, 1.92e-2, 1.98e-2, 2.04e-2, 2.09e-2, 2.13e-2, 2.17e-2, 2.20e-2, 2.22e-2, 2.23e-2, 2.22e-2, 2.21e-2, 2.20e-2, 2.19e-2, 2.18e-2, 2.07e-2, 1.97e-2, 1.87e-2, 1.78e-2, 1.69e-2, 1.60e-2, 1.51e-2, 1.43e-2, 1.35e-2, 1.27e-2, 1.20e-2, 1.13e-2, 1.05e-2, 9.7e-3, 8.9e-3, 8.1e-3, 7.3e-3, 6.7e-3, 6.1e-3, 5.5e-3, 4.9e-3, 4.3e-3, 3.7e-3, 3.4e-3, 3.1e-3, 2.8e-3, 2.5e-3

# prepare number size distributions
fig, ((ax0, ax1, ax2), (ax3, ax4, ax5)) = plt.subplots(2, 3, sharey=True, sharex=True, figsize=(12,7))

# function to calculate percentage change in non-volatile particle-phase concentration
# between start and end of simulation
def nv_change_func(y):
	# inputs:
	# y - concentrations of components (columns) against time (rows), molecules/cc (air)
	
	# assuming the coag_simple_chem.txt scheme used, there are five components presents
	# (APINENE, OZONE, APINOOA, WATER & CORE), so the first index of particle-phase
	# core component is 9 and then it occurs every 5 elements
	nv_change = np.abs(y[0, 9::5].sum()-y[-1, 9::5].sum())/y[0, 9::5].sum()
	
	return(nv_change)

# open results and plot number size distributions

cwd = os.getcwd() # get current working directory

# 8 size bins, 60 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd + '/fig08_data/coag_resol_8sb_60s')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname, delimiter=',', skiprows=1)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname, delimiter=',', skiprows=1)

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax0.loglog(x[0,:]*2.0, N[0,:]/dDp, '-x', markersize=10, label='Initial')
ax0.loglog(x[-1,:]*2.0, N[-1,:]/dDp, 'o', markersize=8, label=r'$\Delta$ t = 60 s')
ax0.text(0.1,3e5, '$\Delta nv_{60 s} = $' + str(nv_change) + '$\,\%$', size=12)


# 8 size bins, 600 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd + '/fig08_data/coag_resol_8sb_600s')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

ax0.loglog(x[-1,:]*2.0, N[-1,:]/dDp, '^', markersize=6, label=r'$\Delta$ t = 600 s')
ax0.text(0.1,1e5, '$\Delta nv_{600 s} = $' + str(nv_change) + '$\,\%$', size=12)

# 8 size bins, 6000 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd + '/fig08_data/coag_resol_8sb_6000s')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax0.loglog(x[-1,:]*2.0, N[-1,:]/dDp, 'x', markersize=4, label=r'$\Delta$ t = 6000 s')
ax0.set_ylim([1.0e-1, 1.0e6])
ax0.text(0.1,3.7e4, '$\Delta nv_{6000 s} = $' + str(nv_change) + '$\,\%$', size=12)

ax0.set_title('8 size bins', size=14)
ax0.legend(fontsize=12)
ax0.yaxis.set_tick_params(labelsize=14)
ax0.set_ylabel('d$N$ $\mathrm{(\#\, cm^{-3})}$/dlog$\mathrm{_{10}}(D_p)$', size=14)

# 32 size bins, 60 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd + '/fig08_data/coag_resol_32sb_60s')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax1.loglog(x[0,:]*2.0, N[0,:]/dDp, '-x', markersize=10, label='Initial 32 size bins')
ax1.loglog(x[-1,:]*2.0, N[-1,:]/dDp, 'o', markersize=8, label=r'$\Delta$ t = 60 s, 32 size bins')
ax1.text(0.1,3.0e5, '$\Delta nv_{60 s} = $' + str(nv_change) + '$\,\%$', size=12)
ax1.set_title('32 size bins', size = 14)

# 32 size bins, 600 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd + '/fig08_data/coag_resol_32sb_600s')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax1.loglog(x[-1,:]*2.0, N[-1,:]/dDp, '^', markersize=6, label=r'$\Delta$ t = 600 s, 32 size bins')
ax1.text(0.1,1.0e5, '$\Delta nv_{600 s} = $' + str(nv_change) + '$\,\%$', size=12)

# 32 size bins, 6000 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd +  '/fig08_data/coag_resol_32sb_6000s')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax1.loglog(x[-1,:]*2.0, N[-1,:]/dDp, 'x', markersize=4, label=r'$\Delta$ t = 6000 s, 32 size bins')
ax1.text(0.1,3.7e4, '$\Delta nv_{6000 s} = $' + str(nv_change) + '$\,\%$', size=12)

# 128 size bins, 60 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd +  '/fig08_data/coag_resol_128sb_60s')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax2.loglog(x[0,:]*2.0, N[0,:]/dDp, '-x', markersize=10, label='Initial 128 size bins')
ax2.loglog(x[-1,:]*2.0, N[-1,:]/dDp, 'o', markersize=8, label=r'$\Delta$ t = 60 s, 128 size bins')
ax2.text(0.1,3.0e5, '$\Delta nv_{60 s} = $' + str(nv_change) + '$\,\%$', size=12)
ax2.set_title('128 size bins', size=14)

# 128 size bins, 600 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd +  '/fig08_data/coag_resol_128sb_600s')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax2.loglog(x[-1,:]*2.0, N[-1,:]/dDp, '^', markersize=6, label=r'$\Delta$ t = 600 s, 128 size bins')
ax2.text(0.1,1.0e5, '$\Delta nv_{600 s} = $' + str(nv_change) + '$\,\%$', size=12)

# 128 size bins, 6000 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd +  '/fig08_data/coag_resol_128sb_6000s')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax2.loglog(x[-1,:]*2.0, N[-1,:]/dDp, 'x', markersize=4, label=r'$\Delta$ t = 6000 s, 128 size bins')
ax2.text(0.1,3.7e4, '$\Delta nv_{6000 s} = $' + str(nv_change) + '$\,\%$', size=12)

# ----------------------------------------------------------------------------------------
# second row is coagulation with SOA formation, using the simple chemical mechanism:

#% 8.05D-16*EXP(-640/TEMP)*0.60 : APINENE + O3 = APINOOA ;

# and the following relevant inputs:
#C0 = 30.0, 30.0
#Comp0 = APINENE, O3
#voli = 2
#volP = 1.0E-6

# 8 size bins, 60 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd +  '/fig08_data/coag_resol_8sb_60s_SOA')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax3.loglog(x[0,:]*2.0, N[0,:]/dDp, '-x', markersize=10, label='Initial')
ax3.loglog(x[-1,:]*2.0, N[-1,:]/dDp, 'o', markersize=8, label=r'$\Delta$ t = 60 s')
ax3.text(0.4,3.0e5, '$\Delta nv_{60 s} = $' + str(nv_change) + '$\,\%$', size=12)
ax3.yaxis.set_tick_params(labelsize=14)
ax3.xaxis.set_tick_params(labelsize=14)

# 8 size bins, 600 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd +  '/fig08_data/coag_resol_8sb_600s_SOA')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax3.loglog(x[-1,:]*2.0, N[-1,:]/dDp, '^', markersize=6, label=r'$\Delta$ t = 600 s')
ax3.text(0.4,1.0e5, '$\Delta nv_{600 s} = $' + str(nv_change) + '$\,\%$', size=12)

# 8 size bins, 6000 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd +  '/fig08_data/coag_resol_8sb_6000s_SOA')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax3.loglog(x[-1,:]*2.0, N[-1,:]/dDp, 'x', markersize=4, label=r'$\Delta$ t = 6000 s')
ax3.text(0.4,3.7e4, '$\Delta nv_{6000 s} = $' + str(nv_change) + '$\,\%$', size=12)
ax3.legend(fontsize=12)
ax3.set_xlabel('$D_p\, \mathrm{(\mu m)}$', size=14)
ax3.set_ylabel('d$N$ $\mathrm{(\#\, cm^{-3})}$/dlog$\mathrm{_{10}}(D_p)$', size=14)

# 32 size bins, 60 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd +  '/fig08_data/coag_resol_32sb_60s_SOA')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax4.loglog(x[0,:]*2.0, N[0,:]/dDp, '-x', markersize=10, label='Initial')
ax4.loglog(x[-1,:]*2.0, N[-1,:]/dDp, 'o', markersize=8, label=r'$\Delta$ t = 60 s')
ax4.text(0.4,3.0e5, '$\Delta nv_{60 s} = $' + str(nv_change) + '$\,\%$', size=12)
ax4.xaxis.set_tick_params(labelsize=14)

# 32 size bins, 600 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd +  '/fig08_data/coag_resol_32sb_600s_SOA')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax4.loglog(x[-1,:]*2.0, N[-1,:]/dDp, '^', markersize=8, label=r'$\Delta$ t = 600 s')
ax4.text(0.4,1.0e5, '$\Delta nv_{600 s} = $' + str(nv_change) + '$\,\%$', size=12)
# 32 size bins, 6000 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd +  '/fig08_data/coag_resol_32sb_6000s_SOA')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax4.loglog(x[-1,:]*2.0, N[-1,:]/dDp, 'x', markersize=6, label=r'$\Delta$ t = 6000 s')
ax4.text(0.4,3.7e4, '$\Delta nv_{6000 s} = $' + str(nv_change) + '$\,\%$', size=12)

ax4.set_xlabel('$D_p\, \mathrm{(\mu m)}$', size=14)


# 128 size bins, 60 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd +  '/fig08_data/coag_resol_128sb_60s_SOA')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax5.loglog(x[0,:]*2.0, N[0,:]/dDp, '-x', markersize=10, label='Initial')
ax5.loglog(x[-1,:]*2.0, N[-1,:]/dDp, 'o', markersize=8, label=r'$\Delta$ t = 60 s')
ax5.text(0.5,3.0e5, '$\Delta nv_{60 s} = $' + str(nv_change) + '$\,\%$', size=12)
ax5.xaxis.set_tick_params(labelsize=14)

# 128 size bins, 600 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd +  '/fig08_data/coag_resol_128sb_600s_SOA')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax5.loglog(x[-1,:]*2.0, N[-1,:]/dDp, '^', markersize=6, label=r'$\Delta$ t = 600 s')
ax5.text(0.5,1.0e5, '$\Delta nv_{600 s} = $' + str(nv_change) + '$\,\%$', size=12)

# 128 size bins, 6000 s ----------------------------------------------------------------------
# open saved files
output_by_sim = str(cwd +  '/fig08_data/coag_resol_128sb_6000s_SOA')

# withdraw times (s)
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw number-size distributions (# particles/cc (air))
fname = str(output_by_sim+'/N_dry')
N = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw radii at size bin centre (um)
fname = str(output_by_sim+'/x')
x = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw size bin bounds, represented by radii (um)
fname = str(output_by_sim+'/sbb')
sbb = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# withdraw molecular concentrations (molecules/cc)
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
nv_change = nv_change_func(y)
nv_change = '{:.2e}'.format(float(nv_change))

Dpbb = sbb*2.0 # diameters at size bin boundaries (um)
dDp = (np.log10(Dpbb))[1::]-(np.log10(Dpbb))[0:-1] # difference in the log10 of bin bounds

ax5.loglog(x[-1,:]*2.0, N[-1,:]/dDp, 'x', markersize=4, label=r'$\Delta$ t = 6000 s')
ax5.text(0.5,3.7e4, '$\Delta nv_{6000 s} = $' + str(nv_change) + '$\,\%$', size=12)
ax5.set_xlabel('$D_p\, \mathrm{(\mu m)}$', size=14)

plt.show()
fig.savefig('fig08.png')