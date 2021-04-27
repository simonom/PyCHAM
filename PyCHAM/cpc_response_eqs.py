'''solving the weighting of particles of different ages over response time of instrument'''
# module to estimate the weighting of particles of different ages during the response time of instrument to represent the mixing of particles of different ages due to differential flow prior to counter 
# File Created at 2021-03-16 15:58:36.852870

import numpy as np

# function for weighting
def cpc_response(delays, wfuncs):
	
	# inputs: -----------------
	# ---------------------------
	
	# remember all times (s) 
	t_all = np.zeros((100)) 
	sht = 1.0 # shortest delay
	resp_timepeak = 1.0 # delay at peak weighting
	# time range for increasing weighting (s)
	if ((resp_timepeak-sht) > 0):
		t = np.arange(sht, (resp_timepeak), (resp_timepeak-sht)/50.)
	if (resp_timepeak == sht):
		t = np.ones((50))*resp_timepeak
	wpre = 1.*t 
	t_all[0:50] = t[:] 
	
	lot = 1.0 # longest delay (s)
	# time range for decreasing weighting (s)
	if ((lot-resp_timepeak) > 0):
		t = np.arange((resp_timepeak+(lot-resp_timepeak)/50.), (lot), (lot-resp_timepeak)/51.)
	if (resp_timepeak == lot):
		t = np.ones((50))*resp_timepeak
	wpro = 1.*t 
	t_all[50::] = t[:] 
	
	# join weighting
	w = np.append(wpre, wpro)
	
	# integrate weight curve
	area = np.trapz(w, t_all)
	# normalise so that integral is one
	if (area > 0):
		w = w/area
	
	return(w, t_all)