'''code to plot results from temporal resolution tests, please call from inside the GMD_paper/Results folder of PyCHAM'''
# important that we can exemplify the sensitivity of model estimates to the resolution of 
# the maximum integration time step and the operator-split time step, therefore this code
# is responsible for recording inputs to tests and plotting corresponding results

import numpy as np
import matplotlib.pyplot as plt
import os
from retrieve_PyCHAM_outputs import retrieve_outputs as retr
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import scipy.constants as si

# for moving-centre, the maximum integration time step is adapted depending on the 
# structure's tolerance for size bin change

# use the tr_tests_chem.txt scheme for these tests (saved in GMD_paper/Results)

# input file is saved as: tr_tests_input.txt

# parameters to assess sensitivity over:
# 128 size bins, Riverside (polydisperse) distribution, with pconc 
# given in coag_resol_test_res_plot.py
# 128 size bins pconc input in model variables input file: pconc = 2800.0, 2500.0, 2200.0, 1950.0, 1700.0, 1500.0, 1300.0, 1150.0, 1000.0, 900.0, 800.0, 720.0, 640.0, 580.0, 520.0, 460.0, 430.0, 400.0, 370.0, 345.0, 320.0, 295.0, 275.0, 255.0, 235.0, 220.0, 205.0, 190.0, 180.0, 170.0, 164.0, 162.0, 158.0, 154.0, 150.0, 158.0, 166.0, 174.0, 182.0, 190.0, 198.0, 206.0, 214.0, 220.0, 228.0, 235.0, 226.0, 217.0, 209.0, 201.0, 194.0, 187.0, 170.0, 146.0, 124.0, 104.0, 86.0, 70.0, 56.0, 43.0, 31.0, 20.0, 22.0, 24.0, 26.0, 28.0, 32.0, 28.0, 24.0, 20.0, 17.0, 14.0, 11.0, 9.0, 7.0, 5.0, 4.0, 3.0e0, 1.5e0, 7.25e-1, 3.6e-1, 1.0e-1, 1.53e-2, 1.62e-2, 1.70e-2, 1.78e-2, 1.85e-2, 1.92e-2, 1.98e-2, 2.04e-2, 2.09e-2, 2.13e-2, 2.17e-2, 2.20e-2, 2.22e-2, 2.23e-2, 2.22e-2, 2.21e-2, 2.20e-2, 2.19e-2, 2.18e-2, 2.07e-2, 1.97e-2, 1.87e-2, 1.78e-2, 1.69e-2, 1.60e-2, 1.51e-2, 1.43e-2, 1.35e-2, 1.27e-2, 1.20e-2, 1.13e-2, 1.05e-2, 9.7e-3, 8.9e-3, 8.1e-3, 7.3e-3, 6.7e-3, 6.1e-3, 5.5e-3, 4.9e-3, 4.3e-3, 3.7e-3, 3.4e-3, 3.1e-3, 2.8e-3, 2.5e-3
# pconct input in model variables input file: pconct = 0
# # in the case of 8 size bins we use the following input:
# pconc = 15000.0, 3600.0, 3150.0, 1200.0, 140.0, 3.0e-1, 2.7e-1, 6.8e-2
# in the case of 32 size bins we use the following input:
# pconc = 10000.0, 6100.0, 3700.0, 2300.0, 1650.0, 1180.0, 800.0, 600.0, 600.0, 640.0, 750.0, 950.0, 850.0, 650.0, 250.0, 80.0, 110.0, 80.0, 30.0, 10.0, 0.2, 7.3e-2, 8.0e-2, 9.0e-2, 9.0e-2, 8.0e-2, 7.0e-2, 5.0e-2, 4.0e-2, 2.6e-2, 1.7e-2, 1.2e-2


# for experiments with partitioning:
# const_infl = two-methylglyceric_acid
# const_infl_t = 0.0, 86400.0
# Cinfl = 0.01, 0.01

# for the nucleation experiment:
# nucv1 = 1.0e4
# nucv2 = -10
# nucv3 = 100
# nuc_comp = ELVOC_o3
# C0 = 0.0, 1.0
# Comp0 = two-methylglyceric_acid, ELVOC_o3
# vol_Comp = ELVOC_o3
# volP = 1.0E-15

# set up plots ---------------------------------------------------------
fig, ((ax0, ax1, ax2, ax3), (ax4, ax5, ax6, ax7)) = plt.subplots(2, 4, figsize=(14, 7), sharey=True)
# ensure sufficient spacing between subplots
fig.subplots_adjust(top=0.90, bottom = 0.10, hspace=0.5, wspace=0.2)

# 128 size bin, polydisperse -------------------------------------------------------------

# outputs from 128 size bin, no partitioning with variable operator-split time step and
# only coagulation affecting particles


# empty dictionary of results from each simulation
num_sb_dict = {}
num_speci_dict = {}
Cfac_dict = {}
y_dict = {}
N_dict = {}
sbb_dict = {}
x_dict = {}
thr_dict = {}
spd_dict = {}

# get current working directory - assumes module called inside the GMD_paper/Results 
# directory of PyCHAM
cwd = os.getcwd()

# number of points in moving average
np_mvav = 1

# coagulation only, no partitioning ------------------------------------------------------
(num_sb_dict['num_sb0'], num_speci_dict['num_speci0'], Cfac_dict['Cfac0'], y_dict['y0'], N_dict['N0'], sbb_dict['sbb0'], x_dict['x0'], thr_dict['thr0'], _, _, _, _, _, spd_dict['spd0']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_43200s_opsplt_128sb_tempconst_seeded'))
(num_sb_dict['num_sb1'], num_speci_dict['num_speci1'], Cfac_dict['Cfac1'], y_dict['y1'], N_dict['N1'], sbb_dict['sbb1'], x_dict['x1'], thr_dict['thr1'], _, _, _, _, _, spd_dict['spd1']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_10800s_opsplt_128sb_tempconst_seeded'))
(num_sb_dict['num_sb2'], num_speci_dict['num_speci2'], Cfac_dict['Cfac2'], y_dict['y2'], N_dict['N2'], sbb_dict['sbb2'], x_dict['x2'], thr_dict['thr2'], _, _, _, _, _, spd_dict['spd2']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_3600s_opsplt_128sb_tempconst_seeded'))
(num_sb_dict['num_sb3'], num_speci_dict['num_speci3'], Cfac_dict['Cfac3'], y_dict['y3'], N_dict['N3'], sbb_dict['sbb3'], x_dict['x3'], thr_dict['thr3'], _, _, _, _, _, spd_dict['spd3']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_1800s_opsplt_128sb_tempconst_seeded'))
(num_sb_dict['num_sb4'], num_speci_dict['num_speci4'], Cfac_dict['Cfac4'], y_dict['y4'], N_dict['N4'], sbb_dict['sbb4'], x_dict['x4'], thr_dict['thr4'], _, _, _, _, _, spd_dict['spd4']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_720s_opsplt_128sb_tempconst_seeded'))
(num_sb_dict['num_sb5'], num_speci_dict['num_speci5'], Cfac_dict['Cfac5'], y_dict['y5'], N_dict['N5'], sbb_dict['sbb5'], x_dict['x5'], thr_dict['thr5'], _, _, _, _, _, spd_dict['spd5']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_360s_opsplt_128sb_tempconst_seeded'))
(num_sb_dict['num_sb6'], num_speci_dict['num_speci6'], Cfac_dict['Cfac6'], y_dict['y6'], N_dict['N6'], sbb_dict['sbb6'], x_dict['x6'], thr_dict['thr6'], _, _, _, _, _, spd_dict['spd6']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_180s_opsplt_128sb_tempconst_seeded'))
(num_sb_dict['num_sb7'], num_speci_dict['num_speci7'], Cfac_dict['Cfac7'], y_dict['y7'], N_dict['N7'], sbb_dict['sbb7'], x_dict['x7'], thr_dict['thr7'], _, _, _, _, _, spd_dict['spd7']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_120s_opsplt_128sb_tempconst_seeded'))
(num_sb_dict['num_sb8'], num_speci_dict['num_speci8'], Cfac_dict['Cfac8'], y_dict['y8'], N_dict['N8'], sbb_dict['sbb8'], x_dict['x8'], thr_dict['thr8'], _, _, _, _, _, spd_dict['spd8']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_60s_opsplt_128sb_tempconst_seeded'))
(num_sb_dict['num_sb9'], num_speci_dict['num_speci9'], Cfac_dict['Cfac9'], y_dict['y9'], N_dict['N9'], sbb_dict['sbb9'], x_dict['x9'], thr_dict['thr9'], _, _, _, _, _, spd_dict['spd9']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_20s_opsplt_128sb_tempconst_seeded'))
(num_sb_dict['num_sb10'], num_speci_dict['num_speci10'], Cfac_dict['Cfac10'], y_dict['y10'], N_dict['N10'], sbb_dict['sbb10'], x_dict['x10'], thr_dict['thr10'], _, _, _, _, _, spd_dict['spd10']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_6s_opsplt_128sb_tempconst_seeded'))

per_err = np.zeros((len(num_sb_dict), 3)) # empty store for percentage error results
per_err[:, 0] = (43200, 10800, 3600, 1800, 720, 360, 180, 120, 60, 20, 6)
speed = np.zeros((1, len(num_sb_dict))) # empty array for simulation times (s)

for resi in np.linspace(len(num_sb_dict)-1, 0, len(num_sb_dict)):
	resi = int(resi)# ensure integer not float
	
	speed[0, resi] = spd_dict[str('spd' + str(resi))] # simulation times (s)
	N = N_dict[str('N' + str(resi))] # particle number concentration
	sbb = sbb_dict[str('sbb' + str(resi))] # size bin bounds
	num_sb = num_sb_dict[str('num_sb' + str(resi))]-1 # number of size bins
	x = x_dict[str('x' + str(resi))] # radii at size bin centre
	t_array = thr_dict[str('thr' + str(resi))]*3600.0 # times results saved at (s)
	t_step = t_array[1]-t_array[0] # operator-split time step
	if resi == len(num_sb_dict)-1:
		t_array_bench = t_array # times to interpolate to (s)
		t_step_bench = t_step
		t_bench_len = int(len(t_array_bench)) # number of times for reference
	
	# empty matrix for holding interpolated values
	Nint = np.zeros((int(86400.0/t_step_bench+1), num_sb))
	t_count = 0 # initial count on rows of matrix containing interpolation results



	for ti in range(int(len(t_array)-1)): # loop through operator-split time steps
		
		# start and finish times for interpolation (hours)
		t_interp = np.zeros((2))
		t_interp[0] = t_array[ti]
		t_interp[1] = t_array[ti+1]
		# points to interpolate
		N2int = np.zeros((2, N.shape[1]))
		N2int[0, :] = N[ti, :]
		N2int[1, :] = N[ti+1, :]
		
		# index of reference time array
		ref_indx = (np.where((t_interp[1]-t_array_bench)>=0))[0][-1]
		# number of reference time points covered
		num_ref_points = int(ref_indx-t_count)
		
		for ti2 in range(num_ref_points): # loop through times to interpolate to

			# loop through size bins
			for sbi in range(num_sb):
				Nint[t_count, sbi] = np.interp(t_array_bench[t_count], t_interp, N2int[:, sbi]) 
			t_count += 1
	

	
	# new array for holding number concentrations		
	N = np.zeros((Nint.shape[0], Nint.shape[1]))
	N[:, :] = Nint[:, :]
	# new array for holding dN/dlog10(Dp)
	dN = np.zeros((N.shape[0], N.shape[1]))
	# loop through times to normalise number concentrations by size bin width
	for ti in range(N.shape[0]):
		dN[ti, :] = N[ti, :]/(np.log10(sbb[1::]*2.0)-np.log10(sbb[0:-1]*2.0))

	dN_av = np.zeros((dN.shape[0], num_sb-(np_mvav-1))) # dN/dlog10(Dp) moving average
	
	for i in range(np_mvav): # loop through number concentration points to get average
	
		# moving average number concentration, times in rows, size bins in columns
		dN_av[:, :] += (dN[:, i:(num_sb)-(np_mvav-(i+1))])/np_mvav
	
	
	if resi == len(num_sb_dict)-1:
		# remember benchmark values
		dN_av_bench = np.zeros((dN_av.shape[0], dN_av.shape[1]))
		dN_av_bench[:, :] = dN_av[:, :]
		dN_av_bench[dN_av_bench==0.0] = 1.0e-40
		N_bench = N.sum(axis=1)
		N_bench[N_bench==0.0] = 1.0e-40 # set zeros to very low number
		continue # no need to plot root-mean square deviation of benchmark against itself

	# get total number concentration per time step (# particles/cc (air))
	N = N.sum(axis=1)
	N[N==0.0] = 1.0e-40 # set zeros to very low number	
	dN_av[dN_av==0.0] = 1.0e-40 # set zeros to very low number
	
	# loop through time steps of the low temporal resolution case to calculate divergence
	# prepare denominator matrices for divergence estimation
	den_matrdN = np.zeros((dN_av.shape[1]))
	den_matrN = np.zeros((1))

	num_times = 0 # number of times to average over

	for ti in range(t_bench_len): # loop through time steps
		if ti == 0:
			continue # skip results at simulation start

		# denominator matrices to limit error to 100 %
		den_matrdN = np.zeros((dN_av_bench.shape[1]))
		den_matrdN[:] = dN_av_bench[ti, :]
		den_matrN = np.zeros((1))
		den_matrN[0] = N_bench[ti]

		# index of size bins where reference value is below comparator value
		lt_indx = dN_av_bench[ti, :]<dN_av[ti, :]
		den_matrdN[lt_indx] = dN_av[ti, lt_indx]
		if N_bench[ti]<N[ti]:
			den_matrN[0] = N[ti]
		
		# index of size bins where both reference and test cases have particles present
		sb_indx = (dN_av_bench[ti, :]>1.0e-10)+(dN_av[ti, :]>1.0e-10)
		# absolute percentage error in number size distribution
		per_error_nsd_ti = ((np.abs(dN_av_bench[ti, sb_indx]-dN_av[ti, sb_indx]))/den_matrdN[sb_indx])*100.0
		
		# absolute percentage error in total particle number concentration
		per_error_N_ti = ((np.abs(N_bench[ti]-N[ti]))/den_matrN)*100.0

		
		if sum(sb_indx) > 0:
			num_times += 1 # count on number of times to average over
			# for number size distribution sum with previous time steps average across size bins
			per_err[resi, 1] += (sum(per_error_nsd_ti)/(sum(sb_indx)))
			# for total number average sum with previous time steps
			per_err[resi, 2] += (per_error_N_ti)
	
	# average over time steps
	per_err[resi, 1] = per_err[resi, 1]/num_times
	per_err[resi, 2] = per_err[resi, 2]/num_times

	

# plot of absolute percentage error, summed across size and times, excluding the reference
# (highest temporal resolution) case
p1, = ax0.loglog(per_err[0:-1, 0], per_err[0:-1, 1], '-k', linewidth=2, label='NSD, coag.')
p2, = ax0.loglog(per_err[0:-1, 0], per_err[0:-1, 2], '-b', linewidth=2, label='$N$, coag.')
ax0.set_title('Seeded, no partitioning', size = 12)

# add contour plot to show computer processing time --------------------------------------
cont_y = np.array((1.5e-2, 1.1e2))
cont_x = np.append(per_err[:, 0], 2.0)
cont_x[0] = 5.0e4

# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
colors = [(0.00, 0.5, 0.0), (0, 1, 0)]#, (0, 1.0, 1.0), (0, 1.0, 0.0), (1.0, 1.0, 0.0), (1.0, 0.0, 0.0)]  # R -> G -> B
n_bin = 100  # Discretizes the colormap interpolation into bins
cmap_name = 'my_list'
# Create the colormap
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
# set contour levels
levels = (MaxNLocator(nbins = n_bin).tick_values(0.85, 4.65))
# associate colours and contour levels
norm1 = BoundaryNorm(levels, ncolors = cm.N, clip=True)
# contour plot with times along x axis and particle diameters along y axis
pc1 = ax0.pcolormesh(cont_x, cont_y, np.log10(speed[:, :]), cmap=cm, norm=norm1)


# end of contour plot section ------------------------------------------------------------

# particle coagulation and loss to wall, no partitioning ---------------------------------
num_sb_dict = {}
num_speci_dict = {}
Cfac_dict = {}
y_dict = {}
N_dict = {}
sbb_dict = {}
x_dict = {}
thr_dict = {}

(num_sb_dict['num_sb0'], num_speci_dict['num_speci0'], Cfac_dict['Cfac0'], y_dict['y0'], N_dict['N0'], sbb_dict['sbb0'], x_dict['x0'], thr_dict['thr0'], _, _, _, _, _, spd_dict['spd0']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_43200s_opspltwlcoag_128sb_tempconst_seeded'))
(num_sb_dict['num_sb1'], num_speci_dict['num_speci1'], Cfac_dict['Cfac1'], y_dict['y1'], N_dict['N1'], sbb_dict['sbb1'], x_dict['x1'], thr_dict['thr1'], _, _, _, _, _, spd_dict['spd1']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_15360s_opspltwlcoag_128sb_tempconst_seeded'))
(num_sb_dict['num_sb2'], num_speci_dict['num_speci2'], Cfac_dict['Cfac2'], y_dict['y2'], N_dict['N2'], sbb_dict['sbb2'], x_dict['x2'], thr_dict['thr2'], _, _, _, _, _, spd_dict['spd2']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_7680s_opspltwlcoag_128sb_tempconst_seeded'))
(num_sb_dict['num_sb3'], num_speci_dict['num_speci3'], Cfac_dict['Cfac3'], y_dict['y3'], N_dict['N3'], sbb_dict['sbb3'], x_dict['x3'], thr_dict['thr3'], _, _, _, _, _, spd_dict['spd3']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_3840s_opspltwlcoag_128sb_tempconst_seeded'))
(num_sb_dict['num_sb4'], num_speci_dict['num_speci4'], Cfac_dict['Cfac4'], y_dict['y4'], N_dict['N4'], sbb_dict['sbb4'], x_dict['x4'], thr_dict['thr4'], _, _, _, _, _, spd_dict['spd4']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_1920s_opspltwlcoag_128sb_tempconst_seeded'))
(num_sb_dict['num_sb5'], num_speci_dict['num_speci5'], Cfac_dict['Cfac5'], y_dict['y5'], N_dict['N5'], sbb_dict['sbb5'], x_dict['x5'], thr_dict['thr5'], _, _, _, _, _, spd_dict['spd5']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_960s_opspltwlcoag_128sb_tempconst_seeded'))
(num_sb_dict['num_sb6'], num_speci_dict['num_speci6'], Cfac_dict['Cfac6'], y_dict['y6'], N_dict['N6'], sbb_dict['sbb6'], x_dict['x6'], thr_dict['thr6'], _, _, _, _, _, spd_dict['spd6']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_480s_opspltwlcoag_128sb_tempconst_seeded'))
(num_sb_dict['num_sb7'], num_speci_dict['num_speci7'], Cfac_dict['Cfac7'], y_dict['y7'], N_dict['N7'], sbb_dict['sbb7'], x_dict['x7'], thr_dict['thr7'], _, _, _, _, _, spd_dict['spd7']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_240s_opspltwlcoag_128sb_tempconst_seeded'))
(num_sb_dict['num_sb8'], num_speci_dict['num_speci8'], Cfac_dict['Cfac8'], y_dict['y8'], N_dict['N8'], sbb_dict['sbb8'], x_dict['x8'], thr_dict['thr8'], _, _, _, _, _, spd_dict['spd8']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_120s_opspltwlcoag_128sb_tempconst_seeded'))
(num_sb_dict['num_sb9'], num_speci_dict['num_speci9'], Cfac_dict['Cfac9'], y_dict['y9'], N_dict['N9'], sbb_dict['sbb9'], x_dict['x9'], thr_dict['thr9'], _, _, _, _, _, spd_dict['spd9']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_60s_opspltwlcoag_128sb_tempconst_seeded'))
(num_sb_dict['num_sb10'], num_speci_dict['num_speci10'], Cfac_dict['Cfac10'], y_dict['y10'], N_dict['N10'], sbb_dict['sbb10'], x_dict['x10'], thr_dict['thr10'], _, _, _, _, _, spd_dict['spd10']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_48s_opspltwlcoag_128sb_tempconst_seeded'))
(num_sb_dict['num_sb11'], num_speci_dict['num_speci11'], Cfac_dict['Cfac11'], y_dict['y11'], N_dict['N11'], sbb_dict['sbb11'], x_dict['x11'], thr_dict['thr11'], _, _, _, _, _, spd_dict['spd11']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_24s_opspltwlcoag_128sb_tempconst_seeded'))
(num_sb_dict['num_sb12'], num_speci_dict['num_speci12'], Cfac_dict['Cfac12'], y_dict['y12'], N_dict['N12'], sbb_dict['sbb12'], x_dict['x12'], thr_dict['thr12'], _, _, _, _, _, spd_dict['spd12']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_12s_opspltwlcoag_128sb_tempconst_seeded'))
(num_sb_dict['num_sb13'], num_speci_dict['num_speci13'], Cfac_dict['Cfac13'], y_dict['y13'], N_dict['N13'], sbb_dict['sbb13'], x_dict['x13'], thr_dict['thr13'], _, _, _, _, _, spd_dict['spd13']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_6s_opspltwlcoag_128sb_tempconst_seeded'))
per_err = np.zeros((len(num_sb_dict), 3)) # empty store for percentage error results
per_err[:, 0] = (43200, 15360, 7680, 3840, 1920, 960, 480, 240, 120, 60, 48, 24, 12, 6)

for resi in np.linspace(len(num_sb_dict)-1, 0, len(num_sb_dict)):
	resi = int(resi)# ensure integer not float
	
	N = N_dict[str('N' + str(resi))] # particle number concentration
	sbb = sbb_dict[str('sbb' + str(resi))] # size bin bounds
	num_sb = num_sb_dict[str('num_sb' + str(resi))]-1 # number of size bins
	x = x_dict[str('x' + str(resi))] # radii at size bin centre
	t_array = thr_dict[str('thr' + str(resi))]*3600.0 # times results saved at (s)
	t_step = t_array[1]-t_array[0] # operator-split time step
	if resi == len(num_sb_dict)-1:
		t_array_bench = t_array # times to interpolate to (s)
		t_step_bench = t_step
		t_bench_len = int(len(t_array_bench)) # number of times for reference
	
	# empty matrix for holding interpolated values
	Nint = np.zeros((int(86400.0/t_step_bench+1), num_sb))
	t_count = 0 # initial count on rows of matrix containing interpolation results



	for ti in range(int(len(t_array)-1)): # loop through operator-split time steps
		
		# start and finish times for interpolation (hours)
		t_interp = np.zeros((2))
		t_interp[0] = t_array[ti]
		t_interp[1] = t_array[ti+1]
		# points to interpolate
		N2int = np.zeros((2, N.shape[1]))
		N2int[0, :] = N[ti, :]
		N2int[1, :] = N[ti+1, :]
		
		# index of reference time array
		ref_indx = (np.where((t_interp[1]-t_array_bench)>=0))[0][-1]
		# number of reference time points covered
		num_ref_points = int(ref_indx-t_count)
		
		for ti2 in range(num_ref_points): # loop through times to interpolate to

			# loop through size bins
			for sbi in range(num_sb):
				Nint[t_count, sbi] = np.interp(t_array_bench[t_count], t_interp, N2int[:, sbi]) 
			t_count += 1
	

	
	# new array for holding number concentrations		
	N = np.zeros((Nint.shape[0], Nint.shape[1]))
	N[:, :] = Nint[:, :]
	# new array for holding dN/dlog10(Dp)
	dN = np.zeros((N.shape[0], N.shape[1]))
	# loop through times to normalise number concentrations by size bin width
	for ti in range(N.shape[0]):
		dN[ti, :] = N[ti, :]/(np.log10(sbb[1::]*2.0)-np.log10(sbb[0:-1]*2.0))

	dN_av = np.zeros((dN.shape[0], num_sb-(np_mvav-1))) # dN/dlog10(Dp) moving average
	
	for i in range(np_mvav): # loop through number concentration points to get average
	
		# moving average number concentration, times in rows, size bins in columns
		dN_av[:, :] += (dN[:, i:(num_sb)-(np_mvav-(i+1))])/np_mvav

	
	if resi == len(num_sb_dict)-1:
		# remember benchmark values
		dN_av_bench = np.zeros((dN_av.shape[0], dN_av.shape[1]))
		dN_av_bench[:, :] = dN_av[:, :]
		dN_av_bench[dN_av_bench==0.0] = 1.0e-40
		N_bench = N.sum(axis=1)
		N_bench[N_bench==0.0] = 1.0e-40 # set zeros to very low number
		continue # no need to plot root-mean square deviation of benchmark against itself

	# get total number concentration per time step (# particles/cc (air))
	N = N.sum(axis=1)
	N[N==0.0] = 1.0e-40 # set zeros to very low number	
	dN_av[dN_av==0.0] = 1.0e-40 # set zeros to very low number
	
	# loop through time steps of the low temporal resolution case to calculate divergence
	# prepare denominator matrices for divergence estimation
	den_matrdN = np.zeros((dN_av.shape[1]))
	den_matrN = np.zeros((1))

	num_times = 0 # number of times to average over

	for ti in range(t_bench_len): # loop through time steps
		if ti == 0:
			continue # skip results at simulation start

		# denominator matrices to limit error to 100 %
		den_matrdN = np.zeros((dN_av_bench.shape[1]))
		den_matrdN[:] = dN_av_bench[ti, :]
		den_matrN = np.zeros((1))
		den_matrN[0] = N_bench[ti]

		# index of size bins where reference value is below comparator value
		lt_indx = dN_av_bench[ti, :]<dN_av[ti, :]
		den_matrdN[lt_indx] = dN_av[ti, lt_indx]
		if N_bench[ti]<N[ti]:
			den_matrN[0] = N[ti]
		
		# index of size bins where both reference and test cases have particles present
		sb_indx = (dN_av_bench[ti, :]>1.0e-10)+(dN_av[ti, :]>1.0e-10)
		# absolute percentage error in number size distribution
		per_error_nsd_ti = ((np.abs(dN_av_bench[ti, sb_indx]-dN_av[ti, sb_indx]))/den_matrdN[sb_indx])*100.0
		
		# absolute percentage error in total particle number concentration
		per_error_N_ti = ((np.abs(N_bench[ti]-N[ti]))/den_matrN)*100.0

		
		if sum(sb_indx) > 0:
			num_times += 1 # count on number of times to average over
			# for number size distribution sum with previous time steps average across size bins
			per_err[resi, 1] += (sum(per_error_nsd_ti)/(sum(sb_indx)))
			# for total number average sum with previous time steps
			per_err[resi, 2] += (per_error_N_ti)
	
	# average over time steps
	per_err[resi, 1] = per_err[resi, 1]/num_times
	per_err[resi, 2] = per_err[resi, 2]/num_times
	
# plot of absolute percentage error, summed across size and times
p5, = ax0.loglog(per_err[0:-1, 0], per_err[0:-1, 1], '-.k', linewidth=2, label='NSD, coag. & wall')
p6, = ax0.loglog(per_err[0:-1, 0], per_err[0:-1, 2], '-.b', linewidth=2, label='$N$, coag. & wall')

lines = [p1, p2, p5, p6]
## ax0.legend(lines, [l.get_label() for l in lines])
# ax0.set_ylim(4.5e-2, 1.05e2)

# ax0.set_xlabel('operator-split time step (s)', size=16) # vertical axis label
ax0.set_ylabel('Deviation (%)\n(temporal resolution)', size=18) # vertical axis label
ax0.set_xlim(2.0e0, 5.0e4)
ax0.set_ylim(1.5e-2, 1.1e2)
ax0.xaxis.set_tick_params(labelsize=16)
ax0.yaxis.set_tick_params(labelsize=16)
ax0.text(x=1.0e0, y=3.0e2, s='(a)', size=14) # plot label
		


# end of no partitioning plot, onto partitioning for coagulation and wall loss effects ---

# coagulation only with partitioning -----------------------------------------------------
num_sb_dict = {}
num_speci_dict = {}
Cfac_dict = {}
y_dict = {}
N_dict = {}
sbb_dict = {}
x_dict = {}
thr_dict = {}
MV_dict = {}
spd_dict = {}

(num_sb_dict['num_sb0'], num_speci_dict['num_speci0'], Cfac_dict['Cfac0'], y_dict['y0'], N_dict['N0'], sbb_dict['sbb0'], x_dict['x0'], thr_dict['thr0'], _, _, _, _, MV_dict['MV0'], spd_dict['spd0']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_43200s_opspltcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb1'], num_speci_dict['num_speci1'], Cfac_dict['Cfac1'], y_dict['y1'], N_dict['N1'], sbb_dict['sbb1'], x_dict['x1'], thr_dict['thr1'], _, _, _, _, MV_dict['MV1'], spd_dict['spd1']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_18000s_opspltcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb2'], num_speci_dict['num_speci2'], Cfac_dict['Cfac2'], y_dict['y2'], N_dict['N2'], sbb_dict['sbb2'], x_dict['x2'], thr_dict['thr2'], _, _, _, _, MV_dict['MV2'], spd_dict['spd2']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_6000s_opspltcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb3'], num_speci_dict['num_speci3'], Cfac_dict['Cfac3'], y_dict['y3'], N_dict['N3'], sbb_dict['sbb3'], x_dict['x3'], thr_dict['thr3'], _, _, _, _, MV_dict['MV3'], spd_dict['spd3']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_1800s_opspltcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb4'], num_speci_dict['num_speci4'], Cfac_dict['Cfac4'], y_dict['y4'], N_dict['N4'], sbb_dict['sbb4'], x_dict['x4'], thr_dict['thr4'], _, _, _, _, MV_dict['MV4'], spd_dict['spd4']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_600s_opspltcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb5'], num_speci_dict['num_speci5'], Cfac_dict['Cfac5'], y_dict['y5'], N_dict['N5'], sbb_dict['sbb5'], x_dict['x5'], thr_dict['thr5'], _, _, _, _, MV_dict['MV5'], spd_dict['spd5']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_180s_opspltcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb6'], num_speci_dict['num_speci6'], Cfac_dict['Cfac6'], y_dict['y6'], N_dict['N6'], sbb_dict['sbb6'], x_dict['x6'], thr_dict['thr6'], _, _, _, _, MV_dict['MV6'], spd_dict['spd6']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_60s_opspltCOAG_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb7'], num_speci_dict['num_speci7'], Cfac_dict['Cfac7'], y_dict['y7'], N_dict['N7'], sbb_dict['sbb7'], x_dict['x7'], thr_dict['thr7'], _, _, _, _, MV_dict['MV7'], spd_dict['spd7']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_30s_opspltcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb8'], num_speci_dict['num_speci8'], Cfac_dict['Cfac8'], y_dict['y8'], N_dict['N8'], sbb_dict['sbb8'], x_dict['x8'], thr_dict['thr8'], _, _, _, _, MV_dict['MV8'], spd_dict['spd8']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_12s_opspltcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb9'], num_speci_dict['num_speci9'], Cfac_dict['Cfac9'], y_dict['y9'], N_dict['N9'], sbb_dict['sbb9'], x_dict['x9'], thr_dict['thr9'], _, _, _, _, MV_dict['MV9'], spd_dict['spd9']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_6s_opspltcoag_128sb_tempconst_seeded_sv'))

per_err = np.zeros((len(num_sb_dict), 4)) # empty store for percentage error results
per_err[:, 0] = (43200.0, 18000.0, 6000.0, 1800.0, 600.0, 180.0, 60.0, 30.0, 12.0, 6.0)

for resi in np.linspace(len(num_sb_dict)-1, 0, len(num_sb_dict)):
	resi = int(resi)# ensure integer not float
	
	num_sb = num_sb_dict[str('num_sb' + str(resi))]-1 # number of size bins
	ncomp = num_speci_dict[str('num_speci' + str(resi))] # number of components
	N = N_dict[str('N' + str(resi))] # particle number concentration
	y = y_dict[str('y' + str(resi))] # component molecular concentrations (molecules/cm3)
	MV = MV_dict[str('MV' + str(resi))] # molar volume (cm3/mol)
	MV = np.array((MV))
	# prepare for multiplying with molar concentrations
	MV = (MV.reshape(1,-1).repeat(num_sb, axis=1)).repeat(y.shape[0], axis=0)
	# particle mass concentrations (ug/m3), note, indexing of y means only particle-phase
	# is considered, then multiply by 1.0e6 to convert from molecules/cm3 to molecules/m3,
	# then divide by Avogadro's constant to get mol/m3, then multiply by molar volume to
	# get cm3/m3, then assume a density of 1.0e6 ug/cm3 and multiply by this to get ug/m3,
	# this is for each component per size bin (columns) per time step (rows)
	y = y[:, ncomp:-ncomp]*1.0e6/si.N_A*MV*1.0e6
	# obtain just the mass concentration of the semi-volatile
	y = y[:, 2::3]
	sbb = sbb_dict[str('sbb' + str(resi))] # size bin bounds
	
	x = x_dict[str('x' + str(resi))] # radii at size bin centre
	t_array = thr_dict[str('thr' + str(resi))]*3600.0 # times results saved at (s)
	t_step = t_array[1]-t_array[0] # operator-split time step
	if resi == len(num_sb_dict)-1:
		t_array_bench = t_array # times to interpolate to (s)
		t_step_bench = t_step
		t_bench_len = int(len(t_array_bench)) # number of times for reference
	
	# empty matrix for holding interpolated values of number concentration
	Nint = np.zeros((int(86400.0/t_step_bench+1), num_sb))
	# empty matrix for holding interpolated values of SOA mass concentration
	SOAint = np.zeros((int(86400.0/t_step_bench+1), num_sb))
	t_count = 0 # initial count on rows of matrix containing interpolation results

	for ti in range(int(len(t_array)-1)): # loop through operator-split time steps
		
		# start and finish times for interpolation (hours)
		t_interp = np.zeros((2))
		t_interp[0] = t_array[ti]
		t_interp[1] = t_array[ti+1]
		# points to interpolate
		N2int = np.zeros((2, N.shape[1]))
		N2int[0, :] = N[ti, :]
		N2int[1, :] = N[ti+1, :]
		SOA2int = np.zeros((2, N.shape[1]))
		SOA2int[0, :] = y[ti, :]
		SOA2int[1, :] = y[ti+1, :]
		
		# index of reference time array
		ref_indx = (np.where((t_interp[1]-t_array_bench)>=0))[0][-1]
		# number of reference time points covered
		num_ref_points = int(ref_indx-t_count)
		
		for ti2 in range(num_ref_points): # loop through times to interpolate to

			# loop through size bins
			for sbi in range(num_sb):
				Nint[t_count, sbi] = np.interp(t_array_bench[t_count], t_interp, N2int[:, sbi])
				SOAint[t_count, sbi] = np.interp(t_array_bench[t_count], t_interp, SOA2int[:, sbi])
			t_count += 1
	

	
	# new array for holding number concentrations		
	N = np.zeros((Nint.shape[0], Nint.shape[1]))
	N[:, :] = Nint[:, :]
	# new array for holding dN/dlog10(Dp)
	dN = np.zeros((N.shape[0], N.shape[1]))
	# new array for holding SOA mass concentration per size bin (columns) against time
	# (rows)
	SOA = np.zeros((SOAint.shape[0], SOAint.shape[1]))
	SOA[:, :] = SOAint[:, :]
	
	# loop through times to normalise number concentrations by size bin width
	for ti in range(N.shape[0]):
		dN[ti, :] = N[ti, :]/(np.log10(sbb[1::]*2.0)-np.log10(sbb[0:-1]*2.0))

	dN_av = np.zeros((dN.shape[0], num_sb-(np_mvav-1))) # dN/dlog10(Dp) moving average
	
	for i in range(np_mvav): # loop through number concentration points to get average
	
		# moving average number concentration, times in rows, size bins in columns
		dN_av[:, :] += (dN[:, i:(num_sb)-(np_mvav-(i+1))])/np_mvav

	
	if resi == len(num_sb_dict)-1:
		# remember benchmark values
		dN_av_bench = np.zeros((dN_av.shape[0], dN_av.shape[1]))
		dN_av_bench[:, :] = dN_av[:, :]
		dN_av_bench[dN_av_bench==0.0] = 1.0e-40
		N_bench = N.sum(axis=1)
		N_bench[N_bench==0.0] = 1.0e-40 # set zeros to very low number
		SOA_bench = np.zeros((SOA.shape[0], 1))
		SOA_bench[:, 0] = SOA.sum(axis=1)
		continue # no need to plot root-mean square deviation of benchmark against itself

	# get total number concentration per time step (# particles/cc (air))
	N = N.sum(axis=1)
	# get total SOA mass concentration per time step (ug/m3 (air))
	SOA = SOA.sum(axis=1)
	N[N==0.0] = 1.0e-40 # set zeros to very low number	
	dN_av[dN_av==0.0] = 1.0e-40 # set zeros to very low number
	SOA[SOA==0.0] = 1.0e-40 # set zeros to very low number
	
	# loop through time steps of the low temporal resolution case to calculate divergence
	# prepare denominator matrices for divergence estimation
	den_matrdN = np.zeros((dN_av.shape[1]))
	den_matrN = np.zeros((1))

	num_times = 0 # number of times to average over

	for ti in range(t_bench_len): # loop through time steps
		if ti == 0:
			continue # skip results at simulation start

		# denominator matrices to limit error to 100 %
		den_matrdN = np.zeros((dN_av_bench.shape[1]))
		den_matrdN[:] = dN_av_bench[ti, :]
		den_matrN = np.zeros((1))
		den_matrN[0] = N_bench[ti]
		den_matrSOA = np.zeros((1))
		den_matrSOA[0] = SOA_bench[ti]
		# index of size bins where reference value is below comparator value
		lt_indx = dN_av_bench[ti, :]<dN_av[ti, :]
		den_matrdN[lt_indx] = dN_av[ti, lt_indx]
		if N_bench[ti]<N[ti]:
			den_matrN[0] = N[ti]
		if SOA_bench[ti]<SOA[ti]:
			den_matrSOA[0] = SOA[ti]
		
		# index of size bins where both reference and test cases have particles present
		sb_indx = (dN_av_bench[ti, :]>1.0e-10)+(dN_av[ti, :]>1.0e-10)
		# absolute percentage error in number size distribution
		per_error_nsd_ti = ((np.abs(dN_av_bench[ti, sb_indx]-dN_av[ti, sb_indx]))/den_matrdN[sb_indx])*100.0
		
		# absolute percentage error in total particle number concentration
		per_error_N_ti = ((np.abs(N_bench[ti]-N[ti]))/den_matrN)*100.0
		# absolute percentage error in SOA
		per_error_SOA_ti = ((np.abs(SOA_bench[ti]-SOA[ti]))/den_matrSOA)*100.0
		
		if sum(sb_indx) > 0:
			num_times += 1 # count on number of times to average over
			# for number size distribution sum with previous time steps average across size bins
			per_err[resi, 1] += (sum(per_error_nsd_ti)/(sum(sb_indx)))
			# for total number average sum with previous time steps
			per_err[resi, 2] += (per_error_N_ti)
			# for SOA average across times and store
			per_err[resi, 3] += (per_error_SOA_ti)
	
	# average over time steps
	per_err[resi, 1] = per_err[resi, 1]/num_times
	per_err[resi, 2] = per_err[resi, 2]/num_times
	per_err[resi, 3] = per_err[resi, 3]/num_times
	
# plot of absolute percentage error, summed across size and times
p7, = ax1.loglog(per_err[0:-1, 0], per_err[0:-1, 1], '-k', linewidth=2, label='NSD, coag.')
p8, = ax1.loglog(per_err[0:-1, 0], per_err[0:-1, 2], '-b', linewidth=2, label='$N$, coag.')
p9, = ax1.loglog(per_err[0:-1, 0], per_err[0:-1, 3], '-r', linewidth=2, label='[SOA], coag.')
ax1.set_title('Seeded, with partitioning', size=12)



# coagulation & wall loss with partitioning ----------------------------------------------
num_sb_dict = {}
num_speci_dict = {}
Cfac_dict = {}
y_dict = {}
N_dict = {}
sbb_dict = {}
x_dict = {}
thr_dict = {}

(num_sb_dict['num_sb0'], num_speci_dict['num_speci0'], Cfac_dict['Cfac0'], y_dict['y0'], N_dict['N0'], sbb_dict['sbb0'], x_dict['x0'], thr_dict['thr0'], _, _, _, _, _, spd_dict['spd0']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_43200s_opspltwlcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb1'], num_speci_dict['num_speci1'], Cfac_dict['Cfac1'], y_dict['y1'], N_dict['N1'], sbb_dict['sbb1'], x_dict['x1'], thr_dict['thr1'], _, _, _, _, _, spd_dict['spd1']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_21600s_opspltwlcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb2'], num_speci_dict['num_speci2'], Cfac_dict['Cfac2'], y_dict['y2'], N_dict['N2'], sbb_dict['sbb2'], x_dict['x2'], thr_dict['thr2'], _, _, _, _, _, spd_dict['spd2']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_10800s_opspltwlcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb3'], num_speci_dict['num_speci3'], Cfac_dict['Cfac3'], y_dict['y3'], N_dict['N3'], sbb_dict['sbb3'], x_dict['x3'], thr_dict['thr3'], _, _, _, _, _, spd_dict['spd3']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_3000s_opspltwlcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb4'], num_speci_dict['num_speci4'], Cfac_dict['Cfac4'], y_dict['y4'], N_dict['N4'], sbb_dict['sbb4'], x_dict['x4'], thr_dict['thr4'], _, _, _, _, _, spd_dict['spd4']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_1200s_opspltwlcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb5'], num_speci_dict['num_speci5'], Cfac_dict['Cfac5'], y_dict['y5'], N_dict['N5'], sbb_dict['sbb5'], x_dict['x5'], thr_dict['thr5'], _, _, _, _, _, spd_dict['spd5']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_600s_opspltwlcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb6'], num_speci_dict['num_speci6'], Cfac_dict['Cfac6'], y_dict['y6'], N_dict['N6'], sbb_dict['sbb6'], x_dict['x6'], thr_dict['thr6'], _, _, _, _, _, spd_dict['spd6']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_480s_opspltwlcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb7'], num_speci_dict['num_speci7'], Cfac_dict['Cfac7'], y_dict['y7'], N_dict['N7'], sbb_dict['sbb7'], x_dict['x7'], thr_dict['thr7'], _, _, _, _, _, spd_dict['spd7']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_180s_opspltwlcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb8'], num_speci_dict['num_speci8'], Cfac_dict['Cfac8'], y_dict['y8'], N_dict['N8'], sbb_dict['sbb8'], x_dict['x8'], thr_dict['thr8'], _, _, _, _, _, spd_dict['spd8']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_20s_opspltwlcoag_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb9'], num_speci_dict['num_speci9'], Cfac_dict['Cfac9'], y_dict['y9'], N_dict['N9'], sbb_dict['sbb9'], x_dict['x9'], thr_dict['thr9'], _, _, _, _, _, spd_dict['spd9']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_6s_opspltwlcoag_128sb_tempconst_seeded_sv'))

per_err = np.zeros((len(num_sb_dict), 4)) # empty store for percentage error results
per_err[:, 0] = (43200.0, 21600, 10800.0, 3000.0, 1200.0, 600.0, 480.0, 180.0, 18.0, 6.0)
speed = np.zeros((1, len(num_sb_dict))) # empty array for simulation times (s)
	
	
for resi in np.linspace(len(num_sb_dict)-1, 0, len(num_sb_dict)):
	resi = int(resi)# ensure integer not float
	
	speed[0, resi] = spd_dict[str('spd' + str(resi))] # simulation times (s)
	num_sb = num_sb_dict[str('num_sb' + str(resi))]-1 # number of size bins
	ncomp = num_speci_dict[str('num_speci' + str(resi))] # number of components
	N = N_dict[str('N' + str(resi))] # particle number concentration
	y = y_dict[str('y' + str(resi))] # component molecular concentrations (molecules/cm3)
	MV = MV_dict[str('MV' + str(resi))] # molar volume (cm3/mol)
	MV = np.array((MV))
	# prepare for multiplying with molar concentrations
	MV = (MV.reshape(1,-1).repeat(num_sb, axis=1)).repeat(y.shape[0], axis=0)
	# particle mass concentrations (ug/m3), note, indexing of y means only particle-phase
	# is considered, then multiply by 1.0e6 to convert from molecules/cm3 to molecules/m3,
	# then divide by Avogadro's constant to get mol/m3, then multiply by molar volume to
	# get cm3/m3, then assume a density of 1.0e6 ug/cm3 and multiply by this to get ug/m3,
	# this is for each component per size bin (columns) per time step (rows)
	y = y[:, ncomp:-ncomp]*1.0e6/si.N_A*MV*1.0e6
	# obtain just the mass concentration of the semi-volatile
	y = y[:, 2::3]
	sbb = sbb_dict[str('sbb' + str(resi))] # size bin bounds
	
	x = x_dict[str('x' + str(resi))] # radii at size bin centre
	t_array = thr_dict[str('thr' + str(resi))]*3600.0 # times results saved at (s)
	t_step = t_array[1]-t_array[0] # operator-split time step
	if resi == len(num_sb_dict)-1:
		t_array_bench = t_array # times to interpolate to (s)
		t_step_bench = t_step
		t_bench_len = int(len(t_array_bench)) # number of times for reference
	
	# empty matrix for holding interpolated values of number concentration
	Nint = np.zeros((int(86400.0/t_step_bench+1), num_sb))
	# empty matrix for holding interpolated values of SOA mass concentration
	SOAint = np.zeros((int(86400.0/t_step_bench+1), num_sb))
	t_count = 0 # initial count on rows of matrix containing interpolation results

	for ti in range(int(len(t_array)-1)): # loop through operator-split time steps
		
		# start and finish times for interpolation (hours)
		t_interp = np.zeros((2))
		t_interp[0] = t_array[ti]
		t_interp[1] = t_array[ti+1]
		# points to interpolate
		N2int = np.zeros((2, N.shape[1]))
		N2int[0, :] = N[ti, :]
		N2int[1, :] = N[ti+1, :]
		SOA2int = np.zeros((2, N.shape[1]))
		SOA2int[0, :] = y[ti, :]
		SOA2int[1, :] = y[ti+1, :]
		
		# index of reference time array
		ref_indx = (np.where((t_interp[1]-t_array_bench)>=0))[0][-1]
		# number of reference time points covered
		num_ref_points = int(ref_indx-t_count)
		
		for ti2 in range(num_ref_points): # loop through times to interpolate to

			# loop through size bins
			for sbi in range(num_sb):
				Nint[t_count, sbi] = np.interp(t_array_bench[t_count], t_interp, N2int[:, sbi])
				SOAint[t_count, sbi] = np.interp(t_array_bench[t_count], t_interp, SOA2int[:, sbi])
			t_count += 1
	

	
	# new array for holding number concentrations		
	N = np.zeros((Nint.shape[0], Nint.shape[1]))
	N[:, :] = Nint[:, :]
	# new array for holding dN/dlog10(Dp)
	dN = np.zeros((N.shape[0], N.shape[1]))
	# new array for holding SOA mass concentration per size bin (columns) against time
	# (rows)
	SOA = np.zeros((SOAint.shape[0], SOAint.shape[1]))
	SOA[:, :] = SOAint[:, :]
	
	# loop through times to normalise number concentrations by size bin width
	for ti in range(N.shape[0]):
		dN[ti, :] = N[ti, :]/(np.log10(sbb[1::]*2.0)-np.log10(sbb[0:-1]*2.0))

	dN_av = np.zeros((dN.shape[0], num_sb-(np_mvav-1))) # dN/dlog10(Dp) moving average
	
	for i in range(np_mvav): # loop through number concentration points to get average
	
		# moving average number concentration, times in rows, size bins in columns
		dN_av[:, :] += (dN[:, i:(num_sb)-(np_mvav-(i+1))])/np_mvav

	
	if resi == len(num_sb_dict)-1:
		# remember benchmark values
		dN_av_bench = np.zeros((dN_av.shape[0], dN_av.shape[1]))
		dN_av_bench[:, :] = dN_av[:, :]
		dN_av_bench[dN_av_bench==0.0] = 1.0e-40
		N_bench = N.sum(axis=1)
		N_bench[N_bench==0.0] = 1.0e-40 # set zeros to very low number
		SOA_bench = np.zeros((SOA.shape[0], 1))
		SOA_bench[:, 0] = SOA.sum(axis=1)
		continue # no need to plot root-mean square deviation of benchmark against itself

	# get total number concentration per time step (# particles/cc (air))
	N = N.sum(axis=1)
	# get total SOA mass concentration per time step (ug/m3 (air))
	SOA = SOA.sum(axis=1)
	N[N==0.0] = 1.0e-40 # set zeros to very low number	
	dN_av[dN_av==0.0] = 1.0e-40 # set zeros to very low number
	SOA[SOA==0.0] = 1.0e-40 # set zeros to very low number
	
	# loop through time steps of the low temporal resolution case to calculate divergence
	# prepare denominator matrices for divergence estimation
	den_matrdN = np.zeros((dN_av.shape[1]))
	den_matrN = np.zeros((1))

	num_times = 0 # number of times to average over

	for ti in range(t_bench_len): # loop through time steps
		if ti == 0:
			continue # skip results at simulation start

		# denominator matrices to limit error to 100 %
		den_matrdN = np.zeros((dN_av_bench.shape[1]))
		den_matrdN[:] = dN_av_bench[ti, :]
		den_matrN = np.zeros((1))
		den_matrN[0] = N_bench[ti]
		den_matrSOA = np.zeros((1))
		den_matrSOA[0] = SOA_bench[ti]
		# index of size bins where reference value is below comparator value
		lt_indx = dN_av_bench[ti, :]<dN_av[ti, :]
		den_matrdN[lt_indx] = dN_av[ti, lt_indx]
		if N_bench[ti]<N[ti]:
			den_matrN[0] = N[ti]
		if SOA_bench[ti]<SOA[ti]:
			den_matrSOA[0] = SOA[ti]
		
		# index of size bins where both reference and test cases have particles present
		sb_indx = (dN_av_bench[ti, :]>1.0e-10)+(dN_av[ti, :]>1.0e-10)
		# absolute percentage error in number size distribution
		per_error_nsd_ti = ((np.abs(dN_av_bench[ti, sb_indx]-dN_av[ti, sb_indx]))/den_matrdN[sb_indx])*100.0
		
		# absolute percentage error in total particle number concentration
		per_error_N_ti = ((np.abs(N_bench[ti]-N[ti]))/den_matrN)*100.0
		# absolute percentage error in SOA
		per_error_SOA_ti = ((np.abs(SOA_bench[ti]-SOA[ti]))/den_matrSOA)*100.0
		
		if sum(sb_indx) > 0:
			num_times += 1 # count on number of times to average over
			# for number size distribution sum with previous time steps average across size bins
			per_err[resi, 1] += (sum(per_error_nsd_ti)/(sum(sb_indx)))
			# for total number average sum with previous time steps
			per_err[resi, 2] += (per_error_N_ti)
			# for SOA average across times and store
			per_err[resi, 3] += (per_error_SOA_ti)
	
	# average over time steps
	per_err[resi, 1] = per_err[resi, 1]/num_times
	per_err[resi, 2] = per_err[resi, 2]/num_times
	per_err[resi, 3] = per_err[resi, 3]/num_times
	
# plot of absolute percentage error, summed across size and times
p11, = ax1.loglog(per_err[0:-1, 0], per_err[0:-1, 1], '-.k', linewidth=2, label='NSD, coag. & wall')
p12, = ax1.loglog(per_err[0:-1, 0], per_err[0:-1, 2], '-.b', linewidth=2, label='$N$, coag. & wall')
p13, = ax1.loglog(per_err[0:-1, 0], per_err[0:-1, 3], '-.r', linewidth=2, label='[SOA], coag. & wall')


ax1.set_xlabel('update time interval (s)', size=16) # horizontal axis label
# ax1.set_yticklabels([]) # turn off labels for number size distribution
ax1.xaxis.set_tick_params(labelsize=16)
# par2.set_ylabel('$|\%\, \Delta|$ N', rotation=270, size=18, color='b') # vertical axis label
# par2.yaxis.set_tick_params(labelsize=16)
# par2.yaxis.label.set_color('blue')
# par2.tick_params(axis='y', colors='blue')
# par2.spines['right'].set_color('blue')

ax1.set_ylim(1.5e-2, 1.1e2)
ax1.set_xlim(2.0e0, 5.0e4)
ax1.text(x=1.0e0, y=3.0e2, s='(b)', size=14)

lines = [p7, p8, p9, p11, p12, p13]
## ax1.legend(lines, [l.get_label() for l in lines])
# ax1.set_ylim(4.5e-2, 1.05e2)
# add contour plot to show computer processing time --------------------------------------
cont_y = np.array((1.5e-2, 1.1e2))
cont_x = np.append(per_err[:, 0], 2.0)
cont_x[0] = 5.0e4

# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
n_bin = 100  # Discretizes the colormap interpolation into bins
cmap_name = 'my_list'
# Create the colormap
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
# set contour levels
levels = (MaxNLocator(nbins = n_bin).tick_values(0.85, 4.65))
# associate colours and contour levels
norm1 = BoundaryNorm(levels, ncolors = cm.N, clip=True)
# contour plot with times along x axis and particle diameters along y axis
pc1 = ax1.pcolormesh(cont_x, cont_y, np.log10(speed[:, :]), cmap=cm, norm=norm1)



# end of contour plot section ------------------------------------------------------------

# partitioning with nucleation section ---------------------------------------------------

num_sb_dict = {}
num_speci_dict = {}
Cfac_dict = {}
y_dict = {}
N_dict = {}
sbb_dict = {}
x_dict = {}
thr_dict = {}
MV_dict = {}
spd_dict = {}

(num_sb_dict['num_sb0'], num_speci_dict['num_speci0'], Cfac_dict['Cfac0'], y_dict['y0'], N_dict['N0'], sbb_dict['sbb0'], x_dict['x0'], thr_dict['thr0'], _, _, _, _, MV_dict['MV0'], spd_dict['spd0']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_43200s_opspltcoagnuc_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb1'], num_speci_dict['num_speci1'], Cfac_dict['Cfac1'], y_dict['y1'], N_dict['N1'], sbb_dict['sbb1'], x_dict['x1'], thr_dict['thr1'], _, _, _, _, MV_dict['MV1'], spd_dict['spd1']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_21600s_opspltcoagnuc_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb2'], num_speci_dict['num_speci2'], Cfac_dict['Cfac2'], y_dict['y2'], N_dict['N2'], sbb_dict['sbb2'], x_dict['x2'], thr_dict['thr2'], _, _, _, _, MV_dict['MV2'], spd_dict['spd2']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_10800s_opspltcoagnuc_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb3'], num_speci_dict['num_speci3'], Cfac_dict['Cfac3'], y_dict['y3'], N_dict['N3'], sbb_dict['sbb3'], x_dict['x3'], thr_dict['thr3'], _, _, _, _, MV_dict['MV3'], spd_dict['spd3']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_3600s_opspltcoagnuc_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb4'], num_speci_dict['num_speci4'], Cfac_dict['Cfac4'], y_dict['y4'], N_dict['N4'], sbb_dict['sbb4'], x_dict['x4'], thr_dict['thr4'], _, _, _, _, MV_dict['MV4'], spd_dict['spd4']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_900s_opspltcoagnuc_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb5'], num_speci_dict['num_speci5'], Cfac_dict['Cfac5'], y_dict['y5'], N_dict['N5'], sbb_dict['sbb5'], x_dict['x5'], thr_dict['thr5'], _, _, _, _, MV_dict['MV5'], spd_dict['spd5']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_360s_opspltcoagnuc_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb6'], num_speci_dict['num_speci6'], Cfac_dict['Cfac6'], y_dict['y6'], N_dict['N6'], sbb_dict['sbb6'], x_dict['x6'], thr_dict['thr6'], _, _, _, _, MV_dict['MV6'], spd_dict['spd6']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_60s_opspltcoagnuc_128sb_tempconst_sv'))
(num_sb_dict['num_sb7'], num_speci_dict['num_speci7'], Cfac_dict['Cfac7'], y_dict['y7'], N_dict['N7'], sbb_dict['sbb7'], x_dict['x7'], thr_dict['thr7'], _, _, _, _, MV_dict['MV7'], spd_dict['spd7']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_20s_opspltcoagnuc_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb8'], num_speci_dict['num_speci8'], Cfac_dict['Cfac8'], y_dict['y8'], N_dict['N8'], sbb_dict['sbb8'], x_dict['x8'], thr_dict['thr8'], _, _, _, _, MV_dict['MV8'], spd_dict['spd8']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_6s_opspltcoagnuc_128sb_tempconst_seeded_sv'))

per_err = np.zeros((len(num_sb_dict), 4)) # empty store for percentage error results
per_err[:, 0] = (43200, 21600, 10800, 3600, 900, 360, 60, 20, 6)
speed = np.zeros((1, len(num_sb_dict))) # empty array for simulation times (s)
	
	
for resi in np.linspace(len(num_sb_dict)-1, 0, len(num_sb_dict)):
	resi = int(resi)# ensure integer not float
	
	speed[0, resi] = spd_dict[str('spd' + str(resi))] # simulation times (s)
	num_sb = num_sb_dict[str('num_sb' + str(resi))]-1 # number of size bins
	ncomp = num_speci_dict[str('num_speci' + str(resi))] # number of components
	N = N_dict[str('N' + str(resi))] # particle number concentration
	y = y_dict[str('y' + str(resi))] # component molecular concentrations (molecules/cm3)
	MV = MV_dict[str('MV' + str(resi))] # molar volume (cm3/mol)
	MV = np.array((MV))
	# prepare for multiplying with molar concentrations
	MV = (MV.reshape(1,-1).repeat(num_sb, axis=1)).repeat(y.shape[0], axis=0)
	# particle mass concentrations (ug/m3), note, indexing of y means only particle-phase
	# is considered, then multiply by 1.0e6 to convert from molecules/cm3 to molecules/m3,
	# then divide by Avogadro's constant to get mol/m3, then multiply by molar volume to
	# get cm3/m3, then assume a density of 1.0e6 ug/cm3 and multiply by this to get ug/m3,
	# this is for each component per size bin (columns) per time step (rows)
	y = y[:, ncomp:-ncomp]*1.0e6/si.N_A*MV*1.0e6
	# obtain just the mass concentration of the semi-volatile
	y = y[:, 2::3]
	sbb = sbb_dict[str('sbb' + str(resi))] # size bin bounds
	
	x = x_dict[str('x' + str(resi))] # radii at size bin centre
	t_array = thr_dict[str('thr' + str(resi))]*3600.0 # times results saved at (s)
	t_step = t_array[1]-t_array[0] # operator-split time step
	if resi == len(num_sb_dict)-1:
		t_array_bench = t_array # times to interpolate to (s)
		t_step_bench = 60.0
		t_bench_len = int(len(t_array_bench)) # number of times for reference
	
	# empty matrix for holding interpolated values of number concentration
	Nint = np.zeros((int(86400.0/t_step_bench+1), num_sb))
	# empty matrix for holding interpolated values of SOA mass concentration
	SOAint = np.zeros((int(86400.0/t_step_bench+1), num_sb))
	t_count = 0 # initial count on rows of matrix containing interpolation results

	for ti in range(int(len(t_array)-1)): # loop through operator-split time steps
		
		# start and finish times for interpolation (hours)
		t_interp = np.zeros((2))
		t_interp[0] = t_array[ti]
		t_interp[1] = t_array[ti+1]
		# points to interpolate
		N2int = np.zeros((2, N.shape[1]))
		N2int[0, :] = N[ti, :]
		N2int[1, :] = N[ti+1, :]
		SOA2int = np.zeros((2, N.shape[1]))
		SOA2int[0, :] = y[ti, :]
		SOA2int[1, :] = y[ti+1, :]
		
		# index of reference time array
		ref_indx = (np.where((t_interp[1]-t_array_bench)>=0))[0][-1]
		# number of reference time points covered
		num_ref_points = int(ref_indx-t_count)
		
		for ti2 in range(num_ref_points): # loop through times to interpolate to

			# loop through size bins
			for sbi in range(num_sb):
				Nint[t_count, sbi] = np.interp(t_array_bench[t_count], t_interp, N2int[:, sbi])
				SOAint[t_count, sbi] = np.interp(t_array_bench[t_count], t_interp, SOA2int[:, sbi])
			t_count += 1
	

	
	# new array for holding number concentrations		
	N = np.zeros((Nint.shape[0], Nint.shape[1]))
	N[:, :] = Nint[:, :]
	# new array for holding dN/dlog10(Dp)
	dN = np.zeros((N.shape[0], N.shape[1]))
	# new array for holding SOA mass concentration per size bin (columns) against time
	# (rows)
	SOA = np.zeros((SOAint.shape[0], SOAint.shape[1]))
	SOA[:, :] = SOAint[:, :]
	
	# loop through times to normalise number concentrations by size bin width
	for ti in range(N.shape[0]):
		dN[ti, :] = N[ti, :]/(np.log10(sbb[1::]*2.0)-np.log10(sbb[0:-1]*2.0))

	dN_av = np.zeros((dN.shape[0], num_sb-(np_mvav-1))) # dN/dlog10(Dp) moving average
	
	for i in range(np_mvav): # loop through number concentration points to get average
	
		# moving average number concentration, times in rows, size bins in columns
		dN_av[:, :] += (dN[:, i:(num_sb)-(np_mvav-(i+1))])/np_mvav

	
	if resi == len(num_sb_dict)-1:
		# remember benchmark values
		dN_av_bench = np.zeros((dN_av.shape[0], dN_av.shape[1]))
		dN_av_bench[:, :] = dN_av[:, :]
		dN_av_bench[dN_av_bench==0.0] = 1.0e-40
		N_bench = N.sum(axis=1)
		N_bench[N_bench==0.0] = 1.0e-40 # set zeros to very low number
		SOA_bench = np.zeros((SOA.shape[0], 1))
		SOA_bench[:, 0] = SOA.sum(axis=1)
		continue # no need to plot root-mean square deviation of benchmark against itself

	# get total number concentration per time step (# particles/cc (air))
	N = N.sum(axis=1)
	# get total SOA mass concentration per time step (ug/m3 (air))
	SOA = SOA.sum(axis=1)
	N[N==0.0] = 1.0e-40 # set zeros to very low number	
	dN_av[dN_av==0.0] = 1.0e-40 # set zeros to very low number
	SOA[SOA==0.0] = 1.0e-40 # set zeros to very low number
	
	# loop through time steps of the low temporal resolution case to calculate divergence
	# prepare denominator matrices for divergence estimation
	den_matrdN = np.zeros((dN_av.shape[1]))
	den_matrN = np.zeros((1))
	num_times = 0 # number of times to average over

	for ti in range(t_bench_len): # loop through time steps
		if ti == 0:
			continue # skip results at simulation start

		# denominator matrices to limit error to 100 %
		den_matrdN = np.zeros((dN_av_bench.shape[1]))
		den_matrdN[:] = dN_av_bench[ti, :]
		den_matrN = np.zeros((1))
		den_matrN[0] = N_bench[ti]
		den_matrSOA = np.zeros((1))
		den_matrSOA[0] = SOA_bench[ti]
		# index of size bins where reference value is below comparator value
		lt_indx = dN_av_bench[ti, :]<dN_av[ti, :]
		den_matrdN[lt_indx] = dN_av[ti, lt_indx]
		if N_bench[ti]<N[ti]:
			den_matrN[0] = N[ti]
		if SOA_bench[ti]<SOA[ti]:
			den_matrSOA[0] = SOA[ti]
		
		# index of size bins where both reference and test cases have particles present
		sb_indx = (dN_av_bench[ti, :]>1.0e-10)+(dN_av[ti, :]>1.0e-10)
		# absolute percentage error in number size distribution
		per_error_nsd_ti = ((np.abs(dN_av_bench[ti, sb_indx]-dN_av[ti, sb_indx]))/den_matrdN[sb_indx])*100.0
		
		# absolute percentage error in total particle number concentration
		per_error_N_ti = ((np.abs(N_bench[ti]-N[ti]))/den_matrN)*100.0
		# absolute percentage error in SOA
		per_error_SOA_ti = ((np.abs(SOA_bench[ti]-SOA[ti]))/den_matrSOA)*100.0
		
		if sum(sb_indx) > 0:
			num_times += 1 # count on number of times to average over
			# for number size distribution sum with previous time steps average across size bins
			per_err[resi, 1] += (sum(per_error_nsd_ti)/(sum(sb_indx)))
			# for total number average sum with previous time steps
			per_err[resi, 2] += (per_error_N_ti)
			# for SOA average across times and store
			per_err[resi, 3] += (per_error_SOA_ti)
	
	# average over time steps
	per_err[resi, 1] = per_err[resi, 1]/num_times
	per_err[resi, 2] = per_err[resi, 2]/num_times
	per_err[resi, 3] = per_err[resi, 3]/num_times
	
# plot of absolute percentage error, summed across size and times
p14, = ax2.loglog(per_err[0:-1, 0], per_err[0:-1, 1], '-k', linewidth=2, label='NSD, coag.')
p15, = ax2.loglog(per_err[0:-1, 0], per_err[0:-1, 2], '-b', linewidth=2, label='$N$, coag.')
p16, = ax2.loglog(per_err[0:-1, 0], per_err[0:-1, 3], '-r', linewidth=2, label='[SOA], coag.')


# ax2.set_xlabel('operator-split time step (s)', size=16) # vertical axis label
# ax1.set_yticklabels([]) # turn off labels for number size distribution
ax2.xaxis.set_tick_params(labelsize=16)
ax2.text(x=1.0e0, y=3.0e2, s='(c)', size=14)
ax2.set_title('Nucleation, with partitioning', size=12)
# par2.set_ylabel('$|\%\, \Delta|$ N', rotation=270, size=18, color='b') # vertical axis label
# par2.yaxis.set_tick_params(labelsize=16)
# par2.yaxis.label.set_color('blue')
# par2.tick_params(axis='y', colors='blue')
# par2.spines['right'].set_color('blue')


# add contour plot to show computer processing time --------------------------------------
cont_y = np.array((1.5e-2, 1.1e2))
cont_x = np.append(per_err[:, 0], 2.0)
cont_x[0] = 5.0e4

# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
n_bin = 100  # Discretizes the colormap interpolation into bins
cmap_name = 'my_list'
# Create the colormap
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
# set contour levels
levels = (MaxNLocator(nbins = n_bin).tick_values(0.85, 4.65))
# associate colours and contour levels
norm1 = BoundaryNorm(levels, ncolors = cm.N, clip=True)
# contour plot with times along x axis and particle diameters along y axis
pc1 = ax2.pcolormesh(cont_x, cont_y, np.log10(speed[:, :]), cmap=cm, norm=norm1)

# end of contour plot section ------------------------------------------------------------

# partitioning with nucleation, coagulation and wall loss section ------------------------

num_sb_dict = {}
num_speci_dict = {}
Cfac_dict = {}
y_dict = {}
N_dict = {}
sbb_dict = {}
x_dict = {}
thr_dict = {}
MV_dict = {}
spd_dict = {}

(num_sb_dict['num_sb0'], num_speci_dict['num_speci0'], Cfac_dict['Cfac0'], y_dict['y0'], N_dict['N0'], sbb_dict['sbb0'], x_dict['x0'], thr_dict['thr0'], _, _, _, _, MV_dict['MV0'], spd_dict['spd0']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_43200s_opspltcoagnucwl_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb1'], num_speci_dict['num_speci1'], Cfac_dict['Cfac1'], y_dict['y1'], N_dict['N1'], sbb_dict['sbb1'], x_dict['x1'], thr_dict['thr1'], _, _, _, _, MV_dict['MV1'], spd_dict['spd1']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_21600s_opspltcoagnucwl_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb2'], num_speci_dict['num_speci2'], Cfac_dict['Cfac2'], y_dict['y2'], N_dict['N2'], sbb_dict['sbb2'], x_dict['x2'], thr_dict['thr2'], _, _, _, _, MV_dict['MV2'], spd_dict['spd2']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_10800s_opspltcoagnucwl_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb3'], num_speci_dict['num_speci3'], Cfac_dict['Cfac3'], y_dict['y3'], N_dict['N3'], sbb_dict['sbb3'], x_dict['x3'], thr_dict['thr3'], _, _, _, _, MV_dict['MV3'], spd_dict['spd3']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_3600s_opspltcoagnucwl_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb4'], num_speci_dict['num_speci4'], Cfac_dict['Cfac4'], y_dict['y4'], N_dict['N4'], sbb_dict['sbb4'], x_dict['x4'], thr_dict['thr4'], _, _, _, _, MV_dict['MV4'], spd_dict['spd4']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_900s_opspltcoagnucwl_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb5'], num_speci_dict['num_speci5'], Cfac_dict['Cfac5'], y_dict['y5'], N_dict['N5'], sbb_dict['sbb5'], x_dict['x5'], thr_dict['thr5'], _, _, _, _, MV_dict['MV5'], spd_dict['spd5']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_360s_opspltcoagnucwl_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb6'], num_speci_dict['num_speci6'], Cfac_dict['Cfac6'], y_dict['y6'], N_dict['N6'], sbb_dict['sbb6'], x_dict['x6'], thr_dict['thr6'], _, _, _, _, MV_dict['MV6'], spd_dict['spd6']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_60s_opspltcoagnucwl_128sb_tempconst_sv'))
(num_sb_dict['num_sb7'], num_speci_dict['num_speci7'], Cfac_dict['Cfac7'], y_dict['y7'], N_dict['N7'], sbb_dict['sbb7'], x_dict['x7'], thr_dict['thr7'], _, _, _, _, MV_dict['MV7'], spd_dict['spd7']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_20s_opspltcoagnucwl_128sb_tempconst_seeded_sv'))
(num_sb_dict['num_sb8'], num_speci_dict['num_speci8'], Cfac_dict['Cfac8'], y_dict['y8'], N_dict['N8'], sbb_dict['sbb8'], x_dict['x8'], thr_dict['thr8'], _, _, _, _, MV_dict['MV8'], spd_dict['spd8']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_6s_opspltcoagnucwl_128sb_tempconst_seeded_sv'))

per_err = np.zeros((len(num_sb_dict), 4)) # empty store for percentage error results
per_err[:, 0] = (43200, 21600, 10800, 3600, 900, 360, 60, 20, 6)
speed = np.zeros((1, len(num_sb_dict))) # empty array for simulation times (s)
	
	
for resi in np.linspace(len(num_sb_dict)-1, 0, len(num_sb_dict)):
	resi = int(resi)# ensure integer not float
	
	speed[0, resi] = spd_dict[str('spd' + str(resi))] # simulation times (s)
	num_sb = num_sb_dict[str('num_sb' + str(resi))]-1 # number of size bins
	ncomp = num_speci_dict[str('num_speci' + str(resi))] # number of components
	N = N_dict[str('N' + str(resi))] # particle number concentration
	y = y_dict[str('y' + str(resi))] # component molecular concentrations (molecules/cm3)
	MV = MV_dict[str('MV' + str(resi))] # molar volume (cm3/mol)
	MV = np.array((MV))
	# prepare for multiplying with molar concentrations
	MV = (MV.reshape(1,-1).repeat(num_sb, axis=1)).repeat(y.shape[0], axis=0)
	# particle mass concentrations (ug/m3), note, indexing of y means only particle-phase
	# is considered, then multiply by 1.0e6 to convert from molecules/cm3 to molecules/m3,
	# then divide by Avogadro's constant to get mol/m3, then multiply by molar volume to
	# get cm3/m3, then assume a density of 1.0e6 ug/cm3 and multiply by this to get ug/m3,
	# this is for each component per size bin (columns) per time step (rows)
	y = y[:, ncomp:-ncomp]*1.0e6/si.N_A*MV*1.0e6
	# obtain just the mass concentration of the semi-volatile
	y = y[:, 2::3]
	sbb = sbb_dict[str('sbb' + str(resi))] # size bin bounds
	
	x = x_dict[str('x' + str(resi))] # radii at size bin centre
	t_array = thr_dict[str('thr' + str(resi))]*3600.0 # times results saved at (s)
	t_step = t_array[1]-t_array[0] # operator-split time step
	if resi == len(num_sb_dict)-1:
		t_array_bench = t_array # times to interpolate to (s)
		t_step_bench = 60.0
		t_bench_len = int(len(t_array_bench)) # number of times for reference
	
	# empty matrix for holding interpolated values of number concentration
	Nint = np.zeros((int(86400.0/t_step_bench+1), num_sb))
	# empty matrix for holding interpolated values of SOA mass concentration
	SOAint = np.zeros((int(86400.0/t_step_bench+1), num_sb))
	t_count = 0 # initial count on rows of matrix containing interpolation results

	for ti in range(int(len(t_array)-1)): # loop through operator-split time steps
		
		# start and finish times for interpolation (hours)
		t_interp = np.zeros((2))
		t_interp[0] = t_array[ti]
		t_interp[1] = t_array[ti+1]
		# points to interpolate
		N2int = np.zeros((2, N.shape[1]))
		N2int[0, :] = N[ti, :]
		N2int[1, :] = N[ti+1, :]
		SOA2int = np.zeros((2, N.shape[1]))
		SOA2int[0, :] = y[ti, :]
		SOA2int[1, :] = y[ti+1, :]
		
		# index of reference time array
		ref_indx = (np.where((t_interp[1]-t_array_bench)>=0))[0][-1]
		# number of reference time points covered
		num_ref_points = int(ref_indx-t_count)
		
		for ti2 in range(num_ref_points): # loop through times to interpolate to

			# loop through size bins
			for sbi in range(num_sb):
				Nint[t_count, sbi] = np.interp(t_array_bench[t_count], t_interp, N2int[:, sbi])
				SOAint[t_count, sbi] = np.interp(t_array_bench[t_count], t_interp, SOA2int[:, sbi])
			t_count += 1
	

	
	# new array for holding number concentrations		
	N = np.zeros((Nint.shape[0], Nint.shape[1]))
	N[:, :] = Nint[:, :]
	# new array for holding dN/dlog10(Dp)
	dN = np.zeros((N.shape[0], N.shape[1]))
	# new array for holding SOA mass concentration per size bin (columns) against time
	# (rows)
	SOA = np.zeros((SOAint.shape[0], SOAint.shape[1]))
	SOA[:, :] = SOAint[:, :]
	
	# loop through times to normalise number concentrations by size bin width
	for ti in range(N.shape[0]):
		dN[ti, :] = N[ti, :]/(np.log10(sbb[1::]*2.0)-np.log10(sbb[0:-1]*2.0))

	dN_av = np.zeros((dN.shape[0], num_sb-(np_mvav-1))) # dN/dlog10(Dp) moving average
	
	for i in range(np_mvav): # loop through number concentration points to get average
	
		# moving average number concentration, times in rows, size bins in columns
		dN_av[:, :] += (dN[:, i:(num_sb)-(np_mvav-(i+1))])/np_mvav

	
	if resi == len(num_sb_dict)-1:
		# remember benchmark values
		dN_av_bench = np.zeros((dN_av.shape[0], dN_av.shape[1]))
		dN_av_bench[:, :] = dN_av[:, :]
		dN_av_bench[dN_av_bench==0.0] = 1.0e-40
		N_bench = N.sum(axis=1)
		N_bench[N_bench==0.0] = 1.0e-40 # set zeros to very low number
		SOA_bench = np.zeros((SOA.shape[0], 1))
		SOA_bench[:, 0] = SOA.sum(axis=1)
		continue # no need to plot root-mean square deviation of benchmark against itself

	# get total number concentration per time step (# particles/cc (air))
	N = N.sum(axis=1)
	# get total SOA mass concentration per time step (ug/m3 (air))
	SOA = SOA.sum(axis=1)
	N[N==0.0] = 1.0e-40 # set zeros to very low number	
	dN_av[dN_av==0.0] = 1.0e-40 # set zeros to very low number
	SOA[SOA==0.0] = 1.0e-40 # set zeros to very low number
	
	# loop through time steps of the low temporal resolution case to calculate divergence
	# prepare denominator matrices for divergence estimation
	den_matrdN = np.zeros((dN_av.shape[1]))
	den_matrN = np.zeros((1))
	num_times = 0 # number of times to average over

	for ti in range(t_bench_len): # loop through time steps
		if ti == 0:
			continue # skip results at simulation start

		# denominator matrices to limit error to 100 %
		den_matrdN = np.zeros((dN_av_bench.shape[1]))
		den_matrdN[:] = dN_av_bench[ti, :]
		den_matrN = np.zeros((1))
		den_matrN[0] = N_bench[ti]
		den_matrSOA = np.zeros((1))
		den_matrSOA[0] = SOA_bench[ti]
		# index of size bins where reference value is below comparator value
		lt_indx = dN_av_bench[ti, :]<dN_av[ti, :]
		den_matrdN[lt_indx] = dN_av[ti, lt_indx]
		if N_bench[ti]<N[ti]:
			den_matrN[0] = N[ti]
		if SOA_bench[ti]<SOA[ti]:
			den_matrSOA[0] = SOA[ti]
		
		# index of size bins where both reference and test cases have particles present
		sb_indx = (dN_av_bench[ti, :]>1.0e-10)+(dN_av[ti, :]>1.0e-10)
		# absolute percentage error in number size distribution
		per_error_nsd_ti = ((np.abs(dN_av_bench[ti, sb_indx]-dN_av[ti, sb_indx]))/den_matrdN[sb_indx])*100.0
		
		# absolute percentage error in total particle number concentration
		per_error_N_ti = ((np.abs(N_bench[ti]-N[ti]))/den_matrN)*100.0
		# absolute percentage error in SOA
		per_error_SOA_ti = ((np.abs(SOA_bench[ti]-SOA[ti]))/den_matrSOA)*100.0
		
		if sum(sb_indx) > 0:
			num_times += 1 # count on number of times to average over
			# for number size distribution sum with previous time steps average across size bins
			per_err[resi, 1] += (sum(per_error_nsd_ti)/(sum(sb_indx)))
			# for total number average sum with previous time steps
			per_err[resi, 2] += (per_error_N_ti)
			# for SOA average across times and store
			per_err[resi, 3] += (per_error_SOA_ti)
	
	# average over time steps
	per_err[resi, 1] = per_err[resi, 1]/num_times
	per_err[resi, 2] = per_err[resi, 2]/num_times
	per_err[resi, 3] = per_err[resi, 3]/num_times
	
# plot of absolute percentage error, summed across size and times
p17, = ax2.loglog(per_err[0:-1, 0], per_err[0:-1, 1], '-.k', linewidth=2, label='NSD, coag. & wall')
p18, = ax2.loglog(per_err[0:-1, 0], per_err[0:-1, 2], '-.b', linewidth=2, label='$N$, coag. & wall')
p19, = ax2.loglog(per_err[0:-1, 0], per_err[0:-1, 3], '-.r', linewidth=2, label='[SOA], coag. & wall')


# par2.set_ylim(1.0e-2, 1.05e2)
ax2.set_xlim(2.0e0, 5.0e4)
# 
lines = [p14, p15, p16, p17, p18, p19]
# ax2.legend(lines, [l.get_label() for l in lines])
ax2.legend(fontsize=12, loc = [1.55, 0.4])

# set fourth plot to be invisible
ax3.set_visible(False)



# sensitivity to spatial resolution, coagulation only no partitioning --------------------
num_sb_dict = {}
num_speci_dict = {}
Cfac_dict = {}
y_dict = {}
N_dict = {}
sbb_dict = {}
x_dict = {}
thr_dict = {}
MV_dict = {}
spd_dict = {}
nam_dict = {}

(num_sb_dict['num_sb0'], num_speci_dict['num_speci0'], Cfac_dict['Cfac0'], y_dict['y0'], N_dict['N0'], sbb_dict['sbb0'], x_dict['x0'], thr_dict['thr0'], nam_dict['nam0'], _, _, _, MV_dict['MV0'], spd_dict['spd0']) = retr(str(cwd + '/tr_tests_data/sens2sr_60s_opspltcoag_8sb_seeded'))
(num_sb_dict['num_sb1'], num_speci_dict['num_speci1'], Cfac_dict['Cfac1'], y_dict['y1'], N_dict['N1'], sbb_dict['sbb1'], x_dict['x1'], thr_dict['thr1'], nam_dict['nam1'], _, _, _, MV_dict['MV1'], spd_dict['spd1']) = retr(str(cwd + '/tr_tests_data/sens2sr_60s_opspltcoag_32sb_seeded'))
(num_sb_dict['num_sb2'], num_speci_dict['num_speci2'], Cfac_dict['Cfac2'], y_dict['y2'], N_dict['N2'], sbb_dict['sbb2'], x_dict['x2'], thr_dict['thr2'], nam_dict['nam2'], _, _, _, MV_dict['MV2'], spd_dict['spd2']) = retr(str(cwd + '/tr_tests_data/sens2sr_60s_opspltcoag_64sb_seeded'))
(num_sb_dict['num_sb3'], num_speci_dict['num_speci3'], Cfac_dict['Cfac3'], y_dict['y3'], N_dict['N3'], sbb_dict['sbb3'], x_dict['x3'], thr_dict['thr3'], nam_dict['nam3'], _, _, _, MV_dict['MV3'], spd_dict['spd3']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_60s_opsplt_128sb_tempconst_seeded'))


per_err = np.zeros((len(num_sb_dict), 4)) # empty store for percentage error results
per_err[:, 0] = (8, 32, 64, 128)
speed = np.zeros((1, len(num_sb_dict))) # empty array for simulation times (s)

for resi in np.linspace(len(num_sb_dict)-1, 0, len(num_sb_dict)):
	resi = int(resi)# ensure integer not float
	
	speed[0, resi] = spd_dict[str('spd' + str(resi))] # simulation times (s)
	N = N_dict[str('N' + str(resi))] # particle number concentration (#/cm3 (air))
	sbb = sbb_dict[str('sbb' + str(resi))] # size bin bounds
	num_sb = num_sb_dict[str('num_sb' + str(resi))]-1 # number of size bins
	if resi==len(num_sb_dict)-1:
		num_sb_bench = num_sb
	ncomp = num_speci_dict[str('num_speci' + str(resi))] # number of components
	comp_names = nam_dict[str('nam' + str(resi))]
	x = x_dict[str('x' + str(resi))] # radii at size bin centre
	t_array = thr_dict[str('thr' + str(resi))]*3600.0 # times results saved at (s)
	t_step = t_array[1]-t_array[0] # operator-split time step
	
	y = y_dict[str('y' + str(resi))] # component molecular concentrations (molecules/cm3)
	MV = MV_dict[str('MV' + str(resi))] # molar volume (cm3/mol)
	MV = np.array((MV))
	y = y[:, ncomp:-ncomp] # get just particle-phase
	# just the semi-volatile
	for comp_indx in range(len(comp_names)):
		if comp_names[comp_indx]=='two-methylglyceric_acid':
			ystart = comp_indx
			break
	y = y[:, ystart::ncomp]
	# SOA mass concentrations (ug/m3)
	# multiply by 1.0e6 to convert from molecules/cm3 to molecules/m3,
	# then divide by Avogadro's constant to get mol/m3, then multiply by molar volume to
	# get cm3/m3, then assume a density of 1.0e6 ug/cm3 and multiply by this to get ug/m3
	y = y*1.0e6/si.N_A*MV[ystart]*1.0e6
	
	
	if resi == len(num_sb_dict)-1:
		sb_array_bench = x[0, :] # size bins to interpolate to (um)
		sb_bench_len = int(len(sb_array_bench)) # number of size bins for reference
	
	
	# empty matrix for holding interpolated values
	dNint = np.zeros((int(86400/60+1), sb_bench_len))
	# empty matrix for holding interpolated values of SOA mass concentration
	SOAint = np.zeros((int(86400/60+1), sb_bench_len))
	int_count = 0 # count on rows of matrix containing interpolation results

	# new array for holding dN/dlog10(Dp)
	dN = np.zeros((N.shape[0], N.shape[1]))
	# loop through times to normalise number concentrations by size bin width
	for ti in range(N.shape[0]):
		dN[ti, :] = N[ti, :]/(np.log10(sbb[1::]*2.0)-np.log10(sbb[0:-1]*2.0))


	for sbi in range(int(num_sb-1)): # loop through operator-split size bins
		
		# start and finish size bins for interpolation (um)
		sb_interp = np.zeros((2))
		sb_interp[0] = x[0, sbi]
		sb_interp[1] = x[0, sbi+1]
		# points to interpolate
		dN2int = np.zeros((N.shape[0], 2))
		dN2int[:, 0] = dN[:, sbi]
		dN2int[:, 1] = dN[:, sbi+1]
		SOA2int = np.zeros((N.shape[0], 2))
		SOA2int[:, 0] = y[:, sbi]
		SOA2int[:, 1] = y[:, sbi+1]
		
		# index of reference array
		ref_indx = (np.where((sb_interp[1]-sb_array_bench)>=0))[0][-1]
		# plus one to account for python indexing
		ref_indx += 1
		# number of reference points covered
		num_ref_points = int(ref_indx-int_count)

		
		for sbi2 in range(num_ref_points): # loop through size bins to interpolate to

			# loop through times
			for ti in range(len(t_array)):
				dNint[ti, int_count] = np.interp(sb_array_bench[int_count], sb_interp, dN2int[ti, :])
				SOAint[ti, int_count] = np.interp(sb_array_bench[int_count], sb_interp, SOA2int[ti, :])
			int_count += 1 # count on number of size bins considered

	dN_av = np.zeros((dNint.shape[0], num_sb_bench-(np_mvav-1))) # dN/dlog10(Dp) moving average
	
	for i in range(np_mvav): # loop through number concentration points to get average
	
		# moving average number concentration, times in rows, size bins in columns
		dN_av[:, :] += (dNint[:, i:(num_sb_bench)-(np_mvav-(i+1))])/np_mvav
	
	
	if resi == len(num_sb_dict)-1:
		# remember benchmark values
		dN_av_bench = np.zeros((dN_av.shape[0], dN_av.shape[1]))
		dN_av_bench[:, :] = dN_av[:, :]
		dN_av_bench[dN_av_bench==0.0] = 1.0e-40
		dNtot_bench = dNint.sum(axis=1)
		SOA_bench = np.zeros((SOAint.shape[0], 1))
		SOA_bench[:, 0] = SOAint.sum(axis=1)
		t_bench_len = len(t_array)
		continue # no need to plot root-mean square deviation of benchmark against itself

	# get total number concentration per time step (# particles/cc (air))
	dNtot = dNint.sum(axis=1)
	# get total SOA mass concentration per time step (ug/m3 (air))
	SOAtot = SOAint.sum(axis=1)
	dNtot[dNtot==0.0] = 1.0e-40 # set zeros to very low number	
	dN_av[dN_av==0.0] = 1.0e-40 # set zeros to very low number
	SOAtot[SOAtot==0.0] = 1.0e-40 # set zeros to very low number
	
	# count on number of times to average divergence results over
	num_times = 0

	for ti in range(t_bench_len): # loop through time steps
		if ti == 0:
			continue # skip results at simulation start

		# denominator matrices to limit error to 100 %
		den_matrdN = np.zeros((dN_av_bench.shape[1]))
		den_matrdN[:] = dN_av_bench[ti, :]
		den_matrNtot = np.zeros((1))
		den_matrNtot[0] = dNtot_bench[ti]
		den_matrSOAtot = np.zeros((1))
		den_matrSOAtot[0] = SOA_bench[ti]
		
		# index of size bins where reference value is below comparator value
		lt_indx = dN_av_bench[ti, :]<dN_av[ti, :]
		den_matrdN[lt_indx] = dN_av[ti, lt_indx]
		if dNtot_bench[ti]<dNtot[ti]:
			den_matrNtot[0] = dNtot[ti]
		if SOA_bench[ti]<SOAtot[ti]:
			den_matrSOAtot[0] = SOAtot[ti]
		
		# index of size bins where both reference and test cases have particles present
		sb_indx = (dN_av_bench[ti, :]>1.0e-10)+(dN_av[ti, :]>1.0e-10)
		# absolute percentage error in number size distribution
		per_error_nsd_ti = ((np.abs(dN_av_bench[ti, sb_indx]-dN_av[ti, sb_indx]))/den_matrdN[sb_indx])*100.0
		
			
		# absolute percentage error in total particle number concentration
		per_error_N_ti = ((np.abs(dNtot[ti]-dNtot_bench[ti]))/den_matrNtot)*100.0
		# absolute percentage error in SOA
		per_error_SOA_ti = ((np.abs(SOA_bench[ti]-SOAtot[ti]))/den_matrSOAtot)*100.0

		if sum(sb_indx) > 0:
			num_times += 1 # count on number of times to average over
			# for number size distribution sum with previous time steps average across size bins
			per_err[resi, 1] += (sum(per_error_nsd_ti)/(sum(sb_indx)))
			# for total number average sum with previous time steps
			per_err[resi, 2] += (per_error_N_ti)
			# for SOA average across times and store
			per_err[resi, 3] += (per_error_SOA_ti)
			
	# average over time steps
	per_err[resi, 1] = per_err[resi, 1]/num_times
	per_err[resi, 2] = per_err[resi, 2]/num_times
	per_err[resi, 3] = per_err[resi, 3]/num_times
	
# plot of absolute percentage error, summed across size and times, excluding the reference
# (highest temporal resolution) case
p20, = ax4.loglog(per_err[0:-1, 0], per_err[0:-1, 1], '-k', linewidth=2, label='NSD, coag.')
p21, = ax4.loglog(per_err[0:-1, 0], per_err[0:-1, 2], '-b', linewidth=2, label='$N$, coag.')

ax4.set_ylabel('Deviation (%)\n(spatial resolution)', size=18) # vertical axis label
ax4.xaxis.set_tick_params(labelsize=16)
ax4.yaxis.set_tick_params(labelsize=16)
ax4.text(x=1.0e0, y=2.0e2, s='(d)', size=14)


# add contour plot to show computer processing time --------------------------------------
cont_y = np.array((1.5e-2, 1.1e2))
cont_x = np.append(2.0, per_err[:, 0])

# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
n_bin = 100  # Discretizes the colormap interpolation into bins
cmap_name = 'my_list'
# Create the colormap
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
# set contour levels
levels = (MaxNLocator(nbins = n_bin).tick_values(0.85, 4.65))
# associate colours and contour levels
norm1 = BoundaryNorm(levels, ncolors = cm.N, clip=True)
# contour plot with times along x axis and particle diameters along y axis
pc1 = ax4.pcolormesh(cont_x, cont_y, np.log10(speed[:, :]), cmap=cm, norm=norm1)

cbaxes = fig.add_axes([0.72, 0.1, 0.03, 0.8]) # axes for colorbar
cb = plt.colorbar(pc1, cax = cbaxes)
cb.ax.tick_params(labelsize=14)   
# colour bar label
cb.set_label('$\mathrm{log_{10}}$ (simulation time (s))', size=14, rotation=270, labelpad=20)

# end of contour plot section ------------------------------------------------------------

# spatial resolution sensitivity section, coagulation and wall loss no partitioning ------

num_sb_dict = {}
num_speci_dict = {}
Cfac_dict = {}
y_dict = {}
N_dict = {}
sbb_dict = {}
x_dict = {}
thr_dict = {}
MV_dict = {}
spd_dict = {}
nam_dict = {}

(num_sb_dict['num_sb0'], num_speci_dict['num_speci0'], Cfac_dict['Cfac0'], y_dict['y0'], N_dict['N0'], sbb_dict['sbb0'], x_dict['x0'], thr_dict['thr0'], nam_dict['nam0'], _, _, _, MV_dict['MV0'], spd_dict['spd0']) = retr(str(cwd + '/tr_tests_data/sens2sr_60s_opspltcoagwl_8sb_seeded'))
(num_sb_dict['num_sb1'], num_speci_dict['num_speci1'], Cfac_dict['Cfac1'], y_dict['y1'], N_dict['N1'], sbb_dict['sbb1'], x_dict['x1'], thr_dict['thr1'], nam_dict['nam1'], _, _, _, MV_dict['MV1'], spd_dict['spd1']) = retr(str(cwd + '/tr_tests_data/sens2sr_60s_opspltcoagwl_32sb_seeded'))
(num_sb_dict['num_sb2'], num_speci_dict['num_speci2'], Cfac_dict['Cfac2'], y_dict['y2'], N_dict['N2'], sbb_dict['sbb2'], x_dict['x2'], thr_dict['thr2'], nam_dict['nam2'], _, _, _, MV_dict['MV2'], spd_dict['spd2']) = retr(str(cwd + '/tr_tests_data/sens2sr_60s_opspltcoagwl_64sb_seeded'))
(num_sb_dict['num_sb3'], num_speci_dict['num_speci3'], Cfac_dict['Cfac3'], y_dict['y3'], N_dict['N3'], sbb_dict['sbb3'], x_dict['x3'], thr_dict['thr3'], nam_dict['nam3'], _, _, _, MV_dict['MV3'], spd_dict['spd3']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_60s_opspltwl_128sb_tempconst_seeded'))


per_err = np.zeros((len(num_sb_dict), 4)) # empty store for percentage error results
per_err[:, 0] = (8, 32, 64, 128)
speed = np.zeros((1, len(num_sb_dict))) # empty array for simulation times (s)

for resi in np.linspace(len(num_sb_dict)-1, 0, len(num_sb_dict)):
	resi = int(resi)# ensure integer not float
	
	speed[0, resi] = spd_dict[str('spd' + str(resi))] # simulation times (s)
	N = N_dict[str('N' + str(resi))] # particle number concentration (#/cm3 (air))
	sbb = sbb_dict[str('sbb' + str(resi))] # size bin bounds
	num_sb = num_sb_dict[str('num_sb' + str(resi))]-1 # number of size bins
	if resi==len(num_sb_dict)-1:
		num_sb_bench = num_sb
	ncomp = num_speci_dict[str('num_speci' + str(resi))] # number of components
	comp_names = nam_dict[str('nam' + str(resi))]
	x = x_dict[str('x' + str(resi))] # radii at size bin centre
	t_array = thr_dict[str('thr' + str(resi))]*3600.0 # times results saved at (s)
	t_step = t_array[1]-t_array[0] # operator-split time step
	
	y = y_dict[str('y' + str(resi))] # component molecular concentrations (molecules/cm3)
	MV = MV_dict[str('MV' + str(resi))] # molar volume (cm3/mol)
	MV = np.array((MV))
	y = y[:, ncomp:-ncomp] # get just particle-phase
	# just the semi-volatile
	for comp_indx in range(len(comp_names)):
		if comp_names[comp_indx]=='two-methylglyceric_acid':
			ystart = comp_indx
			break
	y = y[:, ystart::ncomp]
	# SOA mass concentrations (ug/m3)
	# multiply by 1.0e6 to convert from molecules/cm3 to molecules/m3,
	# then divide by Avogadro's constant to get mol/m3, then multiply by molar volume to
	# get cm3/m3, then assume a density of 1.0e6 ug/cm3 and multiply by this to get ug/m3
	y = y*1.0e6/si.N_A*MV[ystart]*1.0e6
	
	
	if resi == len(num_sb_dict)-1:
		sb_array_bench = x[0, :] # size bins to interpolate to (um)
		sb_bench_len = int(len(sb_array_bench)) # number of size bins for reference
	
	
	# empty matrix for holding interpolated values
	dNint = np.zeros((int(86400/60+1), sb_bench_len))
	# empty matrix for holding interpolated values of SOA mass concentration
	SOAint = np.zeros((int(86400/60+1), sb_bench_len))
	int_count = 0 # count on rows of matrix containing interpolation results

	# new array for holding dN/dlog10(Dp)
	dN = np.zeros((N.shape[0], N.shape[1]))
	# loop through times to normalise number concentrations by size bin width
	for ti in range(N.shape[0]):
		dN[ti, :] = N[ti, :]/(np.log10(sbb[1::]*2.0)-np.log10(sbb[0:-1]*2.0))


	for sbi in range(int(num_sb-1)): # loop through operator-split size bins
		
		# start and finish size bins for interpolation (um)
		sb_interp = np.zeros((2))
		sb_interp[0] = x[0, sbi]
		sb_interp[1] = x[0, sbi+1]
		# points to interpolate
		dN2int = np.zeros((N.shape[0], 2))
		dN2int[:, 0] = dN[:, sbi]
		dN2int[:, 1] = dN[:, sbi+1]
		SOA2int = np.zeros((N.shape[0], 2))
		SOA2int[:, 0] = y[:, sbi]
		SOA2int[:, 1] = y[:, sbi+1]
		
		# index of reference array
		ref_indx = (np.where((sb_interp[1]-sb_array_bench)>=0))[0][-1]
		# plus one to account for python indexing
		ref_indx += 1
		# number of reference points covered
		num_ref_points = int(ref_indx-int_count)

		
		for sbi2 in range(num_ref_points): # loop through size bins to interpolate to

			# loop through times
			for ti in range(len(t_array)):
				dNint[ti, int_count] = np.interp(sb_array_bench[int_count], sb_interp, dN2int[ti, :])
				SOAint[ti, int_count] = np.interp(sb_array_bench[int_count], sb_interp, SOA2int[ti, :])
			int_count += 1 # count on number of size bins considered

	dN_av = np.zeros((dNint.shape[0], num_sb_bench-(np_mvav-1))) # dN/dlog10(Dp) moving average
	
	for i in range(np_mvav): # loop through number concentration points to get average
	
		# moving average number concentration, times in rows, size bins in columns
		dN_av[:, :] += (dNint[:, i:(num_sb_bench)-(np_mvav-(i+1))])/np_mvav
	
	
	if resi == len(num_sb_dict)-1:
		# remember benchmark values
		dN_av_bench = np.zeros((dN_av.shape[0], dN_av.shape[1]))
		dN_av_bench[:, :] = dN_av[:, :]
		dN_av_bench[dN_av_bench==0.0] = 1.0e-40
		dNtot_bench = dNint.sum(axis=1)
		SOA_bench = np.zeros((SOAint.shape[0], 1))
		SOA_bench[:, 0] = SOAint.sum(axis=1)
		t_bench_len = len(t_array)
		continue # no need to plot root-mean square deviation of benchmark against itself

	# get total number concentration per time step (# particles/cc (air))
	dNtot = dNint.sum(axis=1)
	# get total SOA mass concentration per time step (ug/m3 (air))
	SOAtot = SOAint.sum(axis=1)
	dNtot[dNtot==0.0] = 1.0e-40 # set zeros to very low number	
	dN_av[dN_av==0.0] = 1.0e-40 # set zeros to very low number
	SOAtot[SOAtot==0.0] = 1.0e-40 # set zeros to very low number
	
	# count on number of times to average divergence results over
	num_times = 0

	for ti in range(t_bench_len): # loop through time steps
		if ti == 0:
			continue # skip results at simulation start

		# denominator matrices to limit error to 100 %
		den_matrdN = np.zeros((dN_av_bench.shape[1]))
		den_matrdN[:] = dN_av_bench[ti, :]
		den_matrNtot = np.zeros((1))
		den_matrNtot[0] = dNtot_bench[ti]
		den_matrSOAtot = np.zeros((1))
		den_matrSOAtot[0] = SOA_bench[ti]
		
		# index of size bins where reference value is below comparator value
		lt_indx = dN_av_bench[ti, :]<dN_av[ti, :]
		den_matrdN[lt_indx] = dN_av[ti, lt_indx]
		if dNtot_bench[ti]<dNtot[ti]:
			den_matrNtot[0] = dNtot[ti]
		if SOA_bench[ti]<SOAtot[ti]:
			den_matrSOAtot[0] = SOAtot[ti]
		
		# index of size bins where both reference and test cases have particles present
		sb_indx = (dN_av_bench[ti, :]>1.0e-10)+(dN_av[ti, :]>1.0e-10)
		# absolute percentage error in number size distribution
		per_error_nsd_ti = ((np.abs(dN_av_bench[ti, sb_indx]-dN_av[ti, sb_indx]))/den_matrdN[sb_indx])*100.0
		
			
		# absolute percentage error in total particle number concentration
		per_error_N_ti = ((np.abs(dNtot[ti]-dNtot_bench[ti]))/den_matrNtot)*100.0
		# absolute percentage error in SOA
		per_error_SOA_ti = ((np.abs(SOA_bench[ti]-SOAtot[ti]))/den_matrSOAtot)*100.0

		if sum(sb_indx) > 0:
			num_times += 1 # count on number of times to average over
			# for number size distribution sum with previous time steps average across size bins
			per_err[resi, 1] += (sum(per_error_nsd_ti)/(sum(sb_indx)))
			# for total number average sum with previous time steps
			per_err[resi, 2] += (per_error_N_ti)
			# for SOA average across times and store
			per_err[resi, 3] += (per_error_SOA_ti)
			
	# average over time steps
	per_err[resi, 1] = per_err[resi, 1]/num_times
	per_err[resi, 2] = per_err[resi, 2]/num_times
	per_err[resi, 3] = per_err[resi, 3]/num_times
	
# plot of absolute percentage error, summed across size and times, excluding the reference
# (highest temporal resolution) case
p22, = ax4.loglog(per_err[0:-1, 0], per_err[0:-1, 1], '-.k', linewidth=2, label='NSD, coag. & wall')
p23, = ax4.loglog(per_err[0:-1, 0], per_err[0:-1, 2], '-.b', linewidth=2, label='$N$, coag.  & wall')

# lines = [p20, p21, p22, p23]
# ax4.legend(lines, [l.get_label() for l in lines])
# ----------------------------------------------------------------------------------------




# sensitivity to spatial resolution, coagulation only with partitioning ------------------
num_sb_dict = {}
num_speci_dict = {}
Cfac_dict = {}
y_dict = {}
N_dict = {}
sbb_dict = {}
x_dict = {}
thr_dict = {}
MV_dict = {}
spd_dict = {}
nam_dict = {}

(num_sb_dict['num_sb0'], num_speci_dict['num_speci0'], Cfac_dict['Cfac0'], y_dict['y0'], N_dict['N0'], sbb_dict['sbb0'], x_dict['x0'], thr_dict['thr0'], nam_dict['nam0'], _, _, _, MV_dict['MV0'], spd_dict['spd0']) = retr(str(cwd + '/tr_tests_data/sens2sr_60s_opspltcoag_8sb_seeded_sv'))
(num_sb_dict['num_sb1'], num_speci_dict['num_speci1'], Cfac_dict['Cfac1'], y_dict['y1'], N_dict['N1'], sbb_dict['sbb1'], x_dict['x1'], thr_dict['thr1'], nam_dict['nam1'], _, _, _, MV_dict['MV1'], spd_dict['spd1']) = retr(str(cwd + '/tr_tests_data/sens2sr_60s_opspltcoag_32sb_seeded_sv'))
(num_sb_dict['num_sb2'], num_speci_dict['num_speci2'], Cfac_dict['Cfac2'], y_dict['y2'], N_dict['N2'], sbb_dict['sbb2'], x_dict['x2'], thr_dict['thr2'], nam_dict['nam2'], _, _, _, MV_dict['MV2'], spd_dict['spd2']) = retr(str(cwd + '/tr_tests_data/sens2sr_60s_opspltcoag_64sb_seeded_sv'))
(num_sb_dict['num_sb3'], num_speci_dict['num_speci3'], Cfac_dict['Cfac3'], y_dict['y3'], N_dict['N3'], sbb_dict['sbb3'], x_dict['x3'], thr_dict['thr3'], nam_dict['nam3'], _, _, _, MV_dict['MV3'], spd_dict['spd3']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_60s_opspltcoag_128sb_tempconst_seeded_sv'))


per_err = np.zeros((len(num_sb_dict), 4)) # empty store for percentage error results
per_err[:, 0] = (8, 32, 64, 128)
speed = np.zeros((1, len(num_sb_dict))) # empty array for simulation times (s)

for resi in np.linspace(len(num_sb_dict)-1, 0, len(num_sb_dict)):
	resi = int(resi)# ensure integer not float
	
	speed[0, resi] = spd_dict[str('spd' + str(resi))] # simulation times (s)
	N = N_dict[str('N' + str(resi))] # particle number concentration (#/cm3 (air))
	sbb = sbb_dict[str('sbb' + str(resi))] # size bin bounds
	num_sb = num_sb_dict[str('num_sb' + str(resi))]-1 # number of size bins
	if resi==len(num_sb_dict)-1:
		num_sb_bench = num_sb
	ncomp = num_speci_dict[str('num_speci' + str(resi))] # number of components
	comp_names = nam_dict[str('nam' + str(resi))]
	x = x_dict[str('x' + str(resi))] # radii at size bin centre
	t_array = thr_dict[str('thr' + str(resi))]*3600.0 # times results saved at (s)
	t_step = t_array[1]-t_array[0] # operator-split time step
	
	y = y_dict[str('y' + str(resi))] # component molecular concentrations (molecules/cm3)
	MV = MV_dict[str('MV' + str(resi))] # molar volume (cm3/mol)
	MV = np.array((MV))
	y = y[:, ncomp:-ncomp] # get just particle-phase
	# just the semi-volatile
	for comp_indx in range(len(comp_names)):
		if comp_names[comp_indx]=='two-methylglyceric_acid':
			ystart = comp_indx
			break
	y = y[:, ystart::ncomp]
	# SOA mass concentrations (ug/m3)
	# multiply by 1.0e6 to convert from molecules/cm3 to molecules/m3,
	# then divide by Avogadro's constant to get mol/m3, then multiply by molar volume to
	# get cm3/m3, then assume a density of 1.0e6 ug/cm3 and multiply by this to get ug/m3
	y = y*1.0e6/si.N_A*MV[ystart]*1.0e6
	
	
	if resi == len(num_sb_dict)-1:
		sb_array_bench = x[0, :] # size bins to interpolate to (um)
		sb_bench_len = int(len(sb_array_bench)) # number of size bins for reference
	
	
	# empty matrix for holding interpolated values
	dNint = np.zeros((int(86400/60+1), sb_bench_len))
	int_count = 0 # count on rows of matrix containing interpolation results

	# new array for holding dN/dlog10(Dp)
	dN = np.zeros((N.shape[0], N.shape[1]))
	# loop through times to normalise number concentrations by size bin width
	for ti in range(N.shape[0]):
		dN[ti, :] = N[ti, :]/(np.log10(sbb[1::]*2.0)-np.log10(sbb[0:-1]*2.0))


	for sbi in range(int(num_sb-1)): # loop through operator-split size bins
		
		# start and finish size bins for interpolation (um)
		sb_interp = np.zeros((2))
		sb_interp[0] = x[0, sbi]
		sb_interp[1] = x[0, sbi+1]
		# points to interpolate
		dN2int = np.zeros((N.shape[0], 2))
		dN2int[:, 0] = dN[:, sbi]
		dN2int[:, 1] = dN[:, sbi+1]
		
		# index of reference array
		ref_indx = (np.where((sb_interp[1]-sb_array_bench)>=0))[0][-1]
		# plus one to account for python indexing
		ref_indx += 1
		# number of reference points covered
		num_ref_points = int(ref_indx-int_count)

		
		for sbi2 in range(num_ref_points): # loop through size bins to interpolate to

			# loop through times
			for ti in range(len(t_array)):
				dNint[ti, int_count] = np.interp(sb_array_bench[int_count], sb_interp, dN2int[ti, :])
			int_count += 1 # count on number of size bins considered

	dN_av = np.zeros((dNint.shape[0], num_sb_bench-(np_mvav-1))) # dN/dlog10(Dp) moving average
	
	for i in range(np_mvav): # loop through number concentration points to get average
	
		# moving average number concentration, times in rows, size bins in columns
		dN_av[:, :] += (dNint[:, i:(num_sb_bench)-(np_mvav-(i+1))])/np_mvav
	
	
	if resi == len(num_sb_dict)-1:
		# remember benchmark values
		dN_av_bench = np.zeros((dN_av.shape[0], dN_av.shape[1]))
		dN_av_bench[:, :] = dN_av[:, :]
		dN_av_bench[dN_av_bench==0.0] = 1.0e-40
		dNtot_bench = N.sum(axis=1)
		SOA_bench = np.zeros((y.shape[0], 1))
		SOA_bench[:, 0] = y.sum(axis=1)
		t_bench_len = len(t_array)
		continue # no need to plot root-mean square deviation of benchmark against itself

	# get total number concentration per time step (# particles/cc (air))
	dNtot = N.sum(axis=1)
	# get total SOA mass concentration per time step (ug/m3 (air))
	SOAtot = y.sum(axis=1)
	dNtot[dNtot==0.0] = 1.0e-40 # set zeros to very low number	
	dN_av[dN_av==0.0] = 1.0e-40 # set zeros to very low number
	SOAtot[SOAtot==0.0] = 1.0e-40 # set zeros to very low number
	
	# count on number of times to average divergence results over
	num_times = 0

	for ti in range(t_bench_len): # loop through time steps
		if ti == 0:
			continue # skip results at simulation start

		# denominator matrices to limit error to 100 %
		den_matrdN = np.zeros((dN_av_bench.shape[1]))
		den_matrdN[:] = dN_av_bench[ti, :]
		den_matrNtot = np.zeros((1))
		den_matrNtot[0] = dNtot_bench[ti]
		den_matrSOAtot = np.zeros((1))
		den_matrSOAtot[0] = SOA_bench[ti]
		
		# index of size bins where reference value is below comparator value
		lt_indx = dN_av_bench[ti, :]<dN_av[ti, :]
		den_matrdN[lt_indx] = dN_av[ti, lt_indx]
		if dNtot_bench[ti]<dNtot[ti]:
			den_matrNtot[0] = dNtot[ti]
		if SOA_bench[ti]<SOAtot[ti]:
			den_matrSOAtot[0] = SOAtot[ti]
		
		# index of size bins where both reference and test cases have particles present
		sb_indx = (dN_av_bench[ti, :]>1.0e-10)+(dN_av[ti, :]>1.0e-10)
		# absolute percentage error in number size distribution
		per_error_nsd_ti = ((np.abs(dN_av_bench[ti, sb_indx]-dN_av[ti, sb_indx]))/den_matrdN[sb_indx])*100.0
		
			
		# absolute percentage error in total particle number concentration
		per_error_N_ti = ((np.abs(dNtot[ti]-dNtot_bench[ti]))/den_matrNtot)*100.0
		# absolute percentage error in SOA
		per_error_SOA_ti = ((np.abs(SOA_bench[ti]-SOAtot[ti]))/den_matrSOAtot)*100.0

		if sum(sb_indx) > 0:
			num_times += 1 # count on number of times to average over
			# for number size distribution sum with previous time steps average across size bins
			per_err[resi, 1] += (sum(per_error_nsd_ti)/(sum(sb_indx)))
			# for total number average sum with previous time steps
			per_err[resi, 2] += (per_error_N_ti)
			# for SOA average across times and store
			per_err[resi, 3] += (per_error_SOA_ti)
			
	# average over time steps
	per_err[resi, 1] = per_err[resi, 1]/num_times
	per_err[resi, 2] = per_err[resi, 2]/num_times
	per_err[resi, 3] = per_err[resi, 3]/num_times
	
# plot of absolute percentage error, summed across size and times, excluding the reference
# (highest temporal resolution) case
p24, = ax5.loglog(per_err[0:-1, 0], per_err[0:-1, 1], '-k', linewidth=2, label='NSD, coag.')
p25, = ax5.loglog(per_err[0:-1, 0], per_err[0:-1, 2], '-b', linewidth=2, label='$N$, coag.')
p26, = ax5.loglog(per_err[0:-1, 0], per_err[0:-1, 3], '-r', linewidth=2, label='[SOA], coag.')

ax5.xaxis.set_tick_params(labelsize=16)
ax5.text(x=1.0e0, y=2.0e2, s='(e)', size=14)



# add contour plot to show computer processing time --------------------------------------
cont_y = np.array((1.5e-2, 1.1e2))
cont_x = np.append(2.0, per_err[:, 0])

# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
n_bin = 100  # Discretizes the colormap interpolation into bins
cmap_name = 'my_list'
# Create the colormap
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
# set contour levels
levels = (MaxNLocator(nbins = n_bin).tick_values(0.85, 4.65))
# associate colours and contour levels
norm1 = BoundaryNorm(levels, ncolors = cm.N, clip=True)
# contour plot with times along x axis and particle diameters along y axis
pc1 = ax5.pcolormesh(cont_x, cont_y, np.log10(speed[:, :]), cmap=cm, norm=norm1)
  
# colour bar label
cb.set_label('$\mathrm{log_{10}}$ (simulation time (s))', size=14, rotation=270, labelpad=20)

# end of contour plot section ------------------------------------------------------------

# spatial resolution sensitivity section, coagulation and wall loss no partitioning ------

num_sb_dict = {}
num_speci_dict = {}
Cfac_dict = {}
y_dict = {}
N_dict = {}
sbb_dict = {}
x_dict = {}
thr_dict = {}
MV_dict = {}
spd_dict = {}
nam_dict = {}

(num_sb_dict['num_sb0'], num_speci_dict['num_speci0'], Cfac_dict['Cfac0'], y_dict['y0'], N_dict['N0'], sbb_dict['sbb0'], x_dict['x0'], thr_dict['thr0'], nam_dict['nam0'], _, _, _, MV_dict['MV0'], spd_dict['spd0']) = retr(str(cwd + '/tr_tests_data/sens2sr_60s_opspltcoagwl_8sb_seeded_sv'))
(num_sb_dict['num_sb1'], num_speci_dict['num_speci1'], Cfac_dict['Cfac1'], y_dict['y1'], N_dict['N1'], sbb_dict['sbb1'], x_dict['x1'], thr_dict['thr1'], nam_dict['nam1'], _, _, _, MV_dict['MV1'], spd_dict['spd1']) = retr(str(cwd + '/tr_tests_data/sens2sr_60s_opspltcoagwl_32sb_seeded_sv'))
(num_sb_dict['num_sb2'], num_speci_dict['num_speci2'], Cfac_dict['Cfac2'], y_dict['y2'], N_dict['N2'], sbb_dict['sbb2'], x_dict['x2'], thr_dict['thr2'], nam_dict['nam2'], _, _, _, MV_dict['MV2'], spd_dict['spd2']) = retr(str(cwd + '/tr_tests_data/sens2sr_60s_opspltcoagwl_64sb_seeded_sv'))
(num_sb_dict['num_sb3'], num_speci_dict['num_speci3'], Cfac_dict['Cfac3'], y_dict['y3'], N_dict['N3'], sbb_dict['sbb3'], x_dict['x3'], thr_dict['thr3'], nam_dict['nam3'], _, _, _, MV_dict['MV3'], spd_dict['spd3']) = retr(str(cwd + '/tr_tests_data/sens2tr_60s_opspltcoagwl_128sb_seeded_sv'))


per_err = np.zeros((len(num_sb_dict), 4)) # empty store for percentage error results
per_err[:, 0] = (8, 32, 64, 128)
speed = np.zeros((1, len(num_sb_dict))) # empty array for simulation times (s)

for resi in np.linspace(len(num_sb_dict)-1, 0, len(num_sb_dict)):
	resi = int(resi)# ensure integer not float
	
	speed[0, resi] = spd_dict[str('spd' + str(resi))] # simulation times (s)
	N = N_dict[str('N' + str(resi))] # particle number concentration (#/cm3 (air))
	sbb = sbb_dict[str('sbb' + str(resi))] # size bin bounds
	num_sb = num_sb_dict[str('num_sb' + str(resi))]-1 # number of size bins
	if resi==len(num_sb_dict)-1:
		num_sb_bench = num_sb
	ncomp = num_speci_dict[str('num_speci' + str(resi))] # number of components
	comp_names = nam_dict[str('nam' + str(resi))]
	x = x_dict[str('x' + str(resi))] # radii at size bin centre
	t_array = thr_dict[str('thr' + str(resi))]*3600.0 # times results saved at (s)
	t_step = t_array[1]-t_array[0] # operator-split time step
	
	y = y_dict[str('y' + str(resi))] # component molecular concentrations (molecules/cm3)
	MV = MV_dict[str('MV' + str(resi))] # molar volume (cm3/mol)
	MV = np.array((MV))
	y = y[:, ncomp:-ncomp] # get just particle-phase
	# just the semi-volatile
	for comp_indx in range(len(comp_names)):
		if comp_names[comp_indx]=='two-methylglyceric_acid':
			ystart = comp_indx
			break
	y = y[:, ystart::ncomp]
	# SOA mass concentrations (ug/m3)
	# multiply by 1.0e6 to convert from molecules/cm3 to molecules/m3,
	# then divide by Avogadro's constant to get mol/m3, then multiply by molar volume to
	# get cm3/m3, then assume a density of 1.0e6 ug/cm3 and multiply by this to get ug/m3
	y = y*1.0e6/si.N_A*MV[ystart]*1.0e6
	
	
	if resi == len(num_sb_dict)-1:
		sb_array_bench = x[0, :] # size bins to interpolate to (um)
		sb_bench_len = int(len(sb_array_bench)) # number of size bins for reference
	
	
	# empty matrix for holding interpolated values
	dNint = np.zeros((int(86400/60+1), sb_bench_len))
	# empty matrix for holding interpolated values of SOA mass concentration
	int_count = 0 # count on rows of matrix containing interpolation results

	# new array for holding dN/dlog10(Dp)
	dN = np.zeros((N.shape[0], N.shape[1]))
	# loop through times to normalise number concentrations by size bin width
	for ti in range(N.shape[0]):
		dN[ti, :] = N[ti, :]/(np.log10(sbb[1::]*2.0)-np.log10(sbb[0:-1]*2.0))


	for sbi in range(int(num_sb-1)): # loop through operator-split size bins
		
		# start and finish size bins for interpolation (um)
		sb_interp = np.zeros((2))
		sb_interp[0] = x[0, sbi]
		sb_interp[1] = x[0, sbi+1]
		# points to interpolate
		dN2int = np.zeros((N.shape[0], 2))
		dN2int[:, 0] = dN[:, sbi]
		dN2int[:, 1] = dN[:, sbi+1]

		
		# index of reference array
		ref_indx = (np.where((sb_interp[1]-sb_array_bench)>=0))[0][-1]
		# plus one to account for python indexing
		ref_indx += 1
		# number of reference points covered
		num_ref_points = int(ref_indx-int_count)

		
		for sbi2 in range(num_ref_points): # loop through size bins to interpolate to

			# loop through times
			for ti in range(len(t_array)):
				dNint[ti, int_count] = np.interp(sb_array_bench[int_count], sb_interp, dN2int[ti, :])
			int_count += 1 # count on number of size bins considered

	dN_av = np.zeros((dNint.shape[0], num_sb_bench-(np_mvav-1))) # dN/dlog10(Dp) moving average
	
	for i in range(np_mvav): # loop through number concentration points to get average
	
		# moving average number concentration, times in rows, size bins in columns
		dN_av[:, :] += (dNint[:, i:(num_sb_bench)-(np_mvav-(i+1))])/np_mvav
	
	
	if resi == len(num_sb_dict)-1:
		# remember benchmark values
		dN_av_bench = np.zeros((dN_av.shape[0], dN_av.shape[1]))
		dN_av_bench[:, :] = dN_av[:, :]
		dN_av_bench[dN_av_bench==0.0] = 1.0e-40
		dNtot_bench = N.sum(axis=1)
		SOA_bench = np.zeros((y.shape[0], 1))
		SOA_bench[:, 0] = y.sum(axis=1)
		t_bench_len = len(t_array)
		continue # no need to plot root-mean square deviation of benchmark against itself

	# get total number concentration per time step (# particles/cc (air))
	dNtot = N.sum(axis=1)
	# get total SOA mass concentration per time step (ug/m3 (air))
	SOAtot = y.sum(axis=1)
	dNtot[dNtot==0.0] = 1.0e-40 # set zeros to very low number	
	dN_av[dN_av==0.0] = 1.0e-40 # set zeros to very low number
	SOAtot[SOAtot==0.0] = 1.0e-40 # set zeros to very low number
	
	# count on number of times to average divergence results over
	num_times = 0

	for ti in range(t_bench_len): # loop through time steps
		if ti == 0:
			continue # skip results at simulation start

		# denominator matrices to limit error to 100 %
		den_matrdN = np.zeros((dN_av_bench.shape[1]))
		den_matrdN[:] = dN_av_bench[ti, :]
		den_matrNtot = np.zeros((1))
		den_matrNtot[0] = dNtot_bench[ti]
		den_matrSOAtot = np.zeros((1))
		den_matrSOAtot[0] = SOA_bench[ti]
		
		# index of size bins where reference value is below comparator value
		lt_indx = dN_av_bench[ti, :]<dN_av[ti, :]
		den_matrdN[lt_indx] = dN_av[ti, lt_indx]
		if dNtot_bench[ti]<dNtot[ti]:
			den_matrNtot[0] = dNtot[ti]
		if SOA_bench[ti]<SOAtot[ti]:
			den_matrSOAtot[0] = SOAtot[ti]
		
		# index of size bins where both reference and test cases have particles present
		sb_indx = (dN_av_bench[ti, :]>1.0e-10)+(dN_av[ti, :]>1.0e-10)
		# absolute percentage error in number size distribution
		per_error_nsd_ti = ((np.abs(dN_av_bench[ti, sb_indx]-dN_av[ti, sb_indx]))/den_matrdN[sb_indx])*100.0
		
			
		# absolute percentage error in total particle number concentration
		per_error_N_ti = ((np.abs(dNtot[ti]-dNtot_bench[ti]))/den_matrNtot)*100.0
		# absolute percentage error in SOA
		per_error_SOA_ti = ((np.abs(SOA_bench[ti]-SOAtot[ti]))/den_matrSOAtot)*100.0

		if sum(sb_indx) > 0:
			num_times += 1 # count on number of times to average over
			# for number size distribution sum with previous time steps average across size bins
			per_err[resi, 1] += (sum(per_error_nsd_ti)/(sum(sb_indx)))
			# for total number average sum with previous time steps
			per_err[resi, 2] += (per_error_N_ti)
			# for SOA average across times and store
			per_err[resi, 3] += (per_error_SOA_ti)
			
	# average over time steps
	per_err[resi, 1] = per_err[resi, 1]/num_times
	per_err[resi, 2] = per_err[resi, 2]/num_times
	per_err[resi, 3] = per_err[resi, 3]/num_times
	
# plot of absolute percentage error, summed across size and times, excluding the reference
# (highest temporal resolution) case
p27, = ax5.loglog(per_err[0:-1, 0], per_err[0:-1, 1], '-.k', linewidth=2, label='NSD, coag. & wall')
p28, = ax5.loglog(per_err[0:-1, 0], per_err[0:-1, 2], '-.b', linewidth=2, label='$N$, coag.  & wall')
p29, = ax5.loglog(per_err[0:-1, 0], per_err[0:-1, 3], '-.r', linewidth=2, label='[SOA], coag.  & wall')

# lines = [p24, p25, p26, p27, p28, p29]
# ax5.legend(lines, [l.get_label() for l in lines])
ax5.set_xlabel('number of size bins', size=18) # horizontal axis label

# ----------------------------------------------------------------------------------------




# sensitivity to spatial resolution, coagulation only with partitioning and nucleation ---
num_sb_dict = {}
num_speci_dict = {}
Cfac_dict = {}
y_dict = {}
N_dict = {}
sbb_dict = {}
x_dict = {}
thr_dict = {}
MV_dict = {}
spd_dict = {}
nam_dict = {}

(num_sb_dict['num_sb0'], num_speci_dict['num_speci0'], Cfac_dict['Cfac0'], y_dict['y0'], N_dict['N0'], sbb_dict['sbb0'], x_dict['x0'], thr_dict['thr0'], nam_dict['nam0'], _, _, _, MV_dict['MV0'], spd_dict['spd0']) = retr(str(cwd + '/tr_tests_data/sens2sr_60s_opspltcoagnuc_8sb_sv'))
(num_sb_dict['num_sb1'], num_speci_dict['num_speci1'], Cfac_dict['Cfac1'], y_dict['y1'], N_dict['N1'], sbb_dict['sbb1'], x_dict['x1'], thr_dict['thr1'], nam_dict['nam1'], _, _, _, MV_dict['MV1'], spd_dict['spd1']) = retr(str(cwd + '/tr_tests_data/sens2tr_60s_opspltcoagnuc_32sb_sv'))
(num_sb_dict['num_sb2'], num_speci_dict['num_speci2'], Cfac_dict['Cfac2'], y_dict['y2'], N_dict['N2'], sbb_dict['sbb2'], x_dict['x2'], thr_dict['thr2'], nam_dict['nam2'], _, _, _, MV_dict['MV2'], spd_dict['spd2']) = retr(str(cwd + '/tr_tests_data/sens2tr_60s_opspltcoagnuc_64sb_sv'))
(num_sb_dict['num_sb3'], num_speci_dict['num_speci3'], Cfac_dict['Cfac3'], y_dict['y3'], N_dict['N3'], sbb_dict['sbb3'], x_dict['x3'], thr_dict['thr3'], nam_dict['nam3'], _, _, _, MV_dict['MV3'], spd_dict['spd3']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_60s_opspltcoagnuc_128sb_tempconst_sv'))


per_err = np.zeros((len(num_sb_dict), 4)) # empty store for percentage error results
per_err[:, 0] = (8, 32, 64, 128)
speed = np.zeros((1, len(num_sb_dict))) # empty array for simulation times (s)

for resi in np.linspace(len(num_sb_dict)-1, 0, len(num_sb_dict)):
	resi = int(resi)# ensure integer not float
	
	speed[0, resi] = spd_dict[str('spd' + str(resi))] # simulation times (s)
	N = N_dict[str('N' + str(resi))] # particle number concentration (#/cm3 (air))
	sbb = sbb_dict[str('sbb' + str(resi))] # size bin bounds
	num_sb = num_sb_dict[str('num_sb' + str(resi))]-1 # number of size bins
	if resi==len(num_sb_dict)-1:
		num_sb_bench = num_sb
	ncomp = num_speci_dict[str('num_speci' + str(resi))] # number of components
	comp_names = nam_dict[str('nam' + str(resi))]
	x = x_dict[str('x' + str(resi))] # radii at size bin centre
	t_array = thr_dict[str('thr' + str(resi))]*3600.0 # times results saved at (s)
	t_step = t_array[1]-t_array[0] # operator-split time step
	
	y = y_dict[str('y' + str(resi))] # component molecular concentrations (molecules/cm3)
	MV = MV_dict[str('MV' + str(resi))] # molar volume (cm3/mol)
	MV = np.array((MV))
	y = y[:, ncomp:-ncomp] # get just particle-phase
	# just the semi-volatile
	for comp_indx in range(len(comp_names)):
		if comp_names[comp_indx]=='two-methylglyceric_acid':
			ystart = comp_indx
			break
	y = y[:, ystart::ncomp]
	# SOA mass concentrations (ug/m3)
	# multiply by 1.0e6 to convert from molecules/cm3 to molecules/m3,
	# then divide by Avogadro's constant to get mol/m3, then multiply by molar volume to
	# get cm3/m3, then assume a density of 1.0e6 ug/cm3 and multiply by this to get ug/m3
	y = y*1.0e6/si.N_A*MV[ystart]*1.0e6
	
	
	if resi == len(num_sb_dict)-1:
		sb_array_bench = x[0, :] # size bins to interpolate to (um)
		sb_bench_len = int(len(sb_array_bench)) # number of size bins for reference
	
	
	# empty matrix for holding interpolated values
	dNint = np.zeros((int(86400/60+1), sb_bench_len))
	int_count = 0 # count on rows of matrix containing interpolation results

	# new array for holding dN/dlog10(Dp)
	dN = np.zeros((N.shape[0], N.shape[1]))
	# loop through times to normalise number concentrations by size bin width
	for ti in range(N.shape[0]):
		dN[ti, :] = N[ti, :]/(np.log10(sbb[1::]*2.0)-np.log10(sbb[0:-1]*2.0))


	for sbi in range(int(num_sb-1)): # loop through operator-split size bins
		
		# start and finish size bins for interpolation (um)
		sb_interp = np.zeros((2))
		sb_interp[0] = x[0, sbi]
		sb_interp[1] = x[0, sbi+1]
		# points to interpolate
		dN2int = np.zeros((N.shape[0], 2))
		dN2int[:, 0] = dN[:, sbi]
		dN2int[:, 1] = dN[:, sbi+1]

		
		# index of reference array
		ref_indx = (np.where((sb_interp[1]-sb_array_bench)>=0))[0][-1]
		# plus one to account for python indexing
		ref_indx += 1
		# number of reference points covered
		num_ref_points = int(ref_indx-int_count)

		
		for sbi2 in range(num_ref_points): # loop through size bins to interpolate to

			# loop through times
			for ti in range(len(t_array)):
				dNint[ti, int_count] = np.interp(sb_array_bench[int_count], sb_interp, dN2int[ti, :])
			int_count += 1 # count on number of size bins considered

	dN_av = np.zeros((dNint.shape[0], num_sb_bench-(np_mvav-1))) # dN/dlog10(Dp) moving average
	
	for i in range(np_mvav): # loop through number concentration points to get average
	
		# moving average number concentration, times in rows, size bins in columns
		dN_av[:, :] += (dNint[:, i:(num_sb_bench)-(np_mvav-(i+1))])/np_mvav
	
	
	if resi == len(num_sb_dict)-1:
		# remember benchmark values
		dN_av_bench = np.zeros((dN_av.shape[0], dN_av.shape[1]))
		dN_av_bench[:, :] = dN_av[:, :]
		dN_av_bench[dN_av_bench==0.0] = 1.0e-40
		dNtot_bench = N.sum(axis=1)
		SOA_bench = np.zeros((y.shape[0], 1))
		SOA_bench[:, 0] = y.sum(axis=1)
		t_bench_len = len(t_array)
		continue # no need to plot root-mean square deviation of benchmark against itself

	# get total number concentration per time step (# particles/cc (air))
	dNtot = N.sum(axis=1)
	# get total SOA mass concentration per time step (ug/m3 (air))
	SOAtot = y.sum(axis=1)
	dNtot[dNtot==0.0] = 1.0e-40 # set zeros to very low number	
	dN_av[dN_av==0.0] = 1.0e-40 # set zeros to very low number
	SOAtot[SOAtot==0.0] = 1.0e-40 # set zeros to very low number
	
	# count on number of times to average divergence results over
	num_times = 0

	for ti in range(t_bench_len): # loop through time steps
		if ti == 0:
			continue # skip results at simulation start

		# denominator matrices to limit error to 100 %
		den_matrdN = np.zeros((dN_av_bench.shape[1]))
		den_matrdN[:] = dN_av_bench[ti, :]
		den_matrNtot = np.zeros((1))
		den_matrNtot[0] = dNtot_bench[ti]
		den_matrSOAtot = np.zeros((1))
		den_matrSOAtot[0] = SOA_bench[ti]
		
		# index of size bins where reference value is below comparator value
		lt_indx = dN_av_bench[ti, :]<dN_av[ti, :]
		den_matrdN[lt_indx] = dN_av[ti, lt_indx]
		if dNtot_bench[ti]<dNtot[ti]:
			den_matrNtot[0] = dNtot[ti]
		if SOA_bench[ti]<SOAtot[ti]:
			den_matrSOAtot[0] = SOAtot[ti]
		
		# index of size bins where both reference and test cases have particles present
		sb_indx = (dN_av_bench[ti, :]>1.0e-10)+(dN_av[ti, :]>1.0e-10)
		# absolute percentage error in number size distribution
		per_error_nsd_ti = ((np.abs(dN_av_bench[ti, sb_indx]-dN_av[ti, sb_indx]))/den_matrdN[sb_indx])*100.0
		
			
		# absolute percentage error in total particle number concentration
		per_error_N_ti = ((np.abs(dNtot[ti]-dNtot_bench[ti]))/den_matrNtot)*100.0
		# absolute percentage error in SOA
		per_error_SOA_ti = ((np.abs(SOA_bench[ti]-SOAtot[ti]))/den_matrSOAtot)*100.0

		if sum(sb_indx) > 0:
			num_times += 1 # count on number of times to average over
			# for number size distribution sum with previous time steps average across size bins
			per_err[resi, 1] += (sum(per_error_nsd_ti)/(sum(sb_indx)))
			# for total number average sum with previous time steps
			per_err[resi, 2] += (per_error_N_ti)
			# for SOA average across times and store
			per_err[resi, 3] += (per_error_SOA_ti)
			
	# average over time steps
	per_err[resi, 1] = per_err[resi, 1]/num_times
	per_err[resi, 2] = per_err[resi, 2]/num_times
	per_err[resi, 3] = per_err[resi, 3]/num_times
	
# plot of absolute percentage error, summed across size and times, excluding the reference
# (highest temporal resolution) case
p24, = ax6.loglog(per_err[0:-1, 0], per_err[0:-1, 1], '-k', linewidth=2, label='NSD, coag.')
p25, = ax6.loglog(per_err[0:-1, 0], per_err[0:-1, 2], '-b', linewidth=2, label='$N$, coag.')
p26, = ax6.loglog(per_err[0:-1, 0], per_err[0:-1, 3], '-r', linewidth=2, label='[SOA], coag.')

ax6.xaxis.set_tick_params(labelsize=16)
ax6.text(x=1.0e0, y=2.0e2, s='(f)', size=14)


# spatial resolution sensitivity section, coagulation and wall loss with partitioning and 
# nucleation------------------------------------------------------------------------------

num_sb_dict = {}
num_speci_dict = {}
Cfac_dict = {}
y_dict = {}
N_dict = {}
sbb_dict = {}
x_dict = {}
thr_dict = {}
MV_dict = {}
spd_dict = {}
nam_dict = {}

(num_sb_dict['num_sb0'], num_speci_dict['num_speci0'], Cfac_dict['Cfac0'], y_dict['y0'], N_dict['N0'], sbb_dict['sbb0'], x_dict['x0'], thr_dict['thr0'], nam_dict['nam0'], _, _, _, MV_dict['MV0'], spd_dict['spd0']) = retr(str(cwd + '/tr_tests_data/sens2sr_60s_opspltcoagwlnuc_8sb_sv'))
(num_sb_dict['num_sb1'], num_speci_dict['num_speci1'], Cfac_dict['Cfac1'], y_dict['y1'], N_dict['N1'], sbb_dict['sbb1'], x_dict['x1'], thr_dict['thr1'], nam_dict['nam1'], _, _, _, MV_dict['MV1'], spd_dict['spd1']) = retr(str(cwd + '/tr_tests_data/sens2tr_60s_opspltcoagwlnuc_32sb_sv'))
(num_sb_dict['num_sb2'], num_speci_dict['num_speci2'], Cfac_dict['Cfac2'], y_dict['y2'], N_dict['N2'], sbb_dict['sbb2'], x_dict['x2'], thr_dict['thr2'], nam_dict['nam2'], _, _, _, MV_dict['MV2'], spd_dict['spd2']) = retr(str(cwd + '/tr_tests_data/sens2tr_60s_opspltcoagwlnuc_64sb_sv'))
(num_sb_dict['num_sb3'], num_speci_dict['num_speci3'], Cfac_dict['Cfac3'], y_dict['y3'], N_dict['N3'], sbb_dict['sbb3'], x_dict['x3'], thr_dict['thr3'], nam_dict['nam3'], _, _, _, MV_dict['MV3'], spd_dict['spd3']) = retr(str(cwd + '/tr_tests_data/mov_cen_sens2tr_60s_opspltcoagnucwl_128sb_tempconst_sv'))


per_err = np.zeros((len(num_sb_dict), 4)) # empty store for percentage error results
per_err[:, 0] = (8, 32, 64, 128)
speed = np.zeros((1, len(num_sb_dict))) # empty array for simulation times (s)

for resi in np.linspace(len(num_sb_dict)-1, 0, len(num_sb_dict)):
	resi = int(resi)# ensure integer not float
	
	speed[0, resi] = spd_dict[str('spd' + str(resi))] # simulation times (s)
	N = N_dict[str('N' + str(resi))] # particle number concentration (#/cm3 (air))
	sbb = sbb_dict[str('sbb' + str(resi))] # size bin bounds
	num_sb = num_sb_dict[str('num_sb' + str(resi))]-1 # number of size bins
	if resi==len(num_sb_dict)-1:
		num_sb_bench = num_sb
	ncomp = num_speci_dict[str('num_speci' + str(resi))] # number of components
	comp_names = nam_dict[str('nam' + str(resi))]
	x = x_dict[str('x' + str(resi))] # radii at size bin centre
	t_array = thr_dict[str('thr' + str(resi))]*3600.0 # times results saved at (s)
	t_step = t_array[1]-t_array[0] # operator-split time step
	
	y = y_dict[str('y' + str(resi))] # component molecular concentrations (molecules/cm3)
	MV = MV_dict[str('MV' + str(resi))] # molar volume (cm3/mol)
	MV = np.array((MV))
	y = y[:, ncomp:-ncomp] # get just particle-phase
	# just the semi-volatile
	for comp_indx in range(len(comp_names)):
		if comp_names[comp_indx]=='two-methylglyceric_acid':
			ystart = comp_indx
			break
	y = y[:, ystart::ncomp]
	# SOA mass concentrations (ug/m3)
	# multiply by 1.0e6 to convert from molecules/cm3 to molecules/m3,
	# then divide by Avogadro's constant to get mol/m3, then multiply by molar volume to
	# get cm3/m3, then assume a density of 1.0e6 ug/cm3 and multiply by this to get ug/m3
	y = y*1.0e6/si.N_A*MV[ystart]*1.0e6
	
	
	if resi == len(num_sb_dict)-1:
		sb_array_bench = x[0, :] # size bins to interpolate to (um)
		sb_bench_len = int(len(sb_array_bench)) # number of size bins for reference
	
	
	# empty matrix for holding interpolated values
	dNint = np.zeros((int(86400/60+1), sb_bench_len))
	int_count = 0 # count on rows of matrix containing interpolation results

	# new array for holding dN/dlog10(Dp)
	dN = np.zeros((N.shape[0], N.shape[1]))
	# loop through times to normalise number concentrations by size bin width
	for ti in range(N.shape[0]):
		dN[ti, :] = N[ti, :]/(np.log10(sbb[1::]*2.0)-np.log10(sbb[0:-1]*2.0))


	for sbi in range(int(num_sb-1)): # loop through operator-split size bins
		
		# start and finish size bins for interpolation (um)
		sb_interp = np.zeros((2))
		sb_interp[0] = x[0, sbi]
		sb_interp[1] = x[0, sbi+1]
		# points to interpolate
		dN2int = np.zeros((N.shape[0], 2))
		dN2int[:, 0] = dN[:, sbi]
		dN2int[:, 1] = dN[:, sbi+1]

		
		# index of reference array
		ref_indx = (np.where((sb_interp[1]-sb_array_bench)>=0))[0][-1]
		# plus one to account for python indexing
		ref_indx += 1
		# number of reference points covered
		num_ref_points = int(ref_indx-int_count)

		
		for sbi2 in range(num_ref_points): # loop through size bins to interpolate to

			# loop through times
			for ti in range(len(t_array)):
				dNint[ti, int_count] = np.interp(sb_array_bench[int_count], sb_interp, dN2int[ti, :])
				
			int_count += 1 # count on number of size bins considered

	dN_av = np.zeros((dNint.shape[0], num_sb_bench-(np_mvav-1))) # dN/dlog10(Dp) moving average
	
	for i in range(np_mvav): # loop through number concentration points to get average
	
		# moving average number concentration, times in rows, size bins in columns
		dN_av[:, :] += (dNint[:, i:(num_sb_bench)-(np_mvav-(i+1))])/np_mvav
	
	
	if resi == len(num_sb_dict)-1:
		# remember benchmark values
		dN_av_bench = np.zeros((dN_av.shape[0], dN_av.shape[1]))
		dN_av_bench[:, :] = dN_av[:, :]
		dN_av_bench[dN_av_bench==0.0] = 1.0e-40
		dNtot_bench = N.sum(axis=1)
		SOA_bench = np.zeros((y.shape[0], 1))
		SOA_bench[:, 0] = y.sum(axis=1)
		t_bench_len = len(t_array)
		continue # no need to plot root-mean square deviation of benchmark against itself

	# get total number concentration per time step (# particles/cc (air))
	dNtot = N.sum(axis=1)
	# get total SOA mass concentration per time step (ug/m3 (air))
	SOAtot = y.sum(axis=1)
	dNtot[dNtot==0.0] = 1.0e-40 # set zeros to very low number	
	dN_av[dN_av==0.0] = 1.0e-40 # set zeros to very low number
	SOAtot[SOAtot==0.0] = 1.0e-40 # set zeros to very low number
	
	# count on number of times to average divergence results over
	num_times = 0

	for ti in range(t_bench_len): # loop through time steps
		if ti == 0:
			continue # skip results at simulation start

		# denominator matrices to limit error to 100 %
		den_matrdN = np.zeros((dN_av_bench.shape[1]))
		den_matrdN[:] = dN_av_bench[ti, :]
		den_matrNtot = np.zeros((1))
		den_matrNtot[0] = dNtot_bench[ti]
		den_matrSOAtot = np.zeros((1))
		den_matrSOAtot[0] = SOA_bench[ti]
		
		# index of size bins where reference value is below comparator value
		lt_indx = dN_av_bench[ti, :]<dN_av[ti, :]
		den_matrdN[lt_indx] = dN_av[ti, lt_indx]
		if dNtot_bench[ti]<dNtot[ti]:
			den_matrNtot[0] = dNtot[ti]
		if SOA_bench[ti]<SOAtot[ti]:
			den_matrSOAtot[0] = SOAtot[ti]
		
		# index of size bins where both reference and test cases have particles present
		sb_indx = (dN_av_bench[ti, :]>1.0e-10)+(dN_av[ti, :]>1.0e-10)
		# absolute percentage error in number size distribution
		per_error_nsd_ti = ((np.abs(dN_av_bench[ti, sb_indx]-dN_av[ti, sb_indx]))/den_matrdN[sb_indx])*100.0
		
			
		# absolute percentage error in total particle number concentration
		per_error_N_ti = ((np.abs(dNtot[ti]-dNtot_bench[ti]))/den_matrNtot)*100.0
		# absolute percentage error in SOA
		per_error_SOA_ti = ((np.abs(SOA_bench[ti]-SOAtot[ti]))/den_matrSOAtot)*100.0

		if sum(sb_indx) > 0:
			num_times += 1 # count on number of times to average over
			# for number size distribution sum with previous time steps average across size bins
			per_err[resi, 1] += (sum(per_error_nsd_ti)/(sum(sb_indx)))
			# for total number average sum with previous time steps
			per_err[resi, 2] += (per_error_N_ti)
			# for SOA average across times and store
			per_err[resi, 3] += (per_error_SOA_ti)
			
	# average over time steps
	per_err[resi, 1] = per_err[resi, 1]/num_times
	per_err[resi, 2] = per_err[resi, 2]/num_times
	per_err[resi, 3] = per_err[resi, 3]/num_times
	
# plot of absolute percentage error, summed across size and times, excluding the reference
# (highest temporal resolution) case
p27, = ax6.loglog(per_err[0:-1, 0], per_err[0:-1, 1], '-.k', label='NSD, coag. & wall')
p28, = ax6.loglog(per_err[0:-1, 0], per_err[0:-1, 2], '-.b', label='$N$, coag.  & wall')
p29, = ax6.loglog(per_err[0:-1, 0], per_err[0:-1, 3], '-.r', label='[SOA], coag.  & wall')

# add contour plot to show computer processing time --------------------------------------
cont_y = np.array((1.5e-2, 1.1e2))
cont_x = np.append(2.0, per_err[:, 0])

# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
n_bin = 100  # Discretizes the colormap interpolation into bins
cmap_name = 'my_list'
# Create the colormap
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
# set contour levels
levels = (MaxNLocator(nbins = n_bin).tick_values(0.85, 4.65))
# associate colours and contour levels
norm1 = BoundaryNorm(levels, ncolors = cm.N, clip=True)
# contour plot with times along x axis and particle diameters along y axis
pc1 = ax6.pcolormesh(cont_x, cont_y, np.log10(speed[:, :]), cmap=cm, norm=norm1)
  
# colour bar label
cb.set_label('$\mathrm{log_{10}}$ (simulation time (s))', size=14, rotation=270, labelpad=20)

# end of contour plot section ------------------------------------------------------------
# ----------------------------------------------------------------------------------------


# lines = [p24, p25, p26, p27, p28, p29]
# ax6.legend(lines, [l.get_label() for l in lines])

ax7.set_visible(False)

plt.savefig('fig11.png', transparent=True)
plt.show()