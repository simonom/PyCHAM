'''solving the weighting of particles of different ages over response time of instrument'''
# module to estimate the weighting of particles of different ages during the response time of instrument to represent the mixing of particles of different ages due to differential flow prior to counter 
# File Created at 2021-03-10 17:31:31.470657

import numpy as np

# function for weighting
def cpc_response(delays, wfuncs):
	
	# inputs: -----------------
	# ---------------------------
	
	# remember all times (s) 
	t_all = np.zeros((100)) 
	sht = 0.1 # shortest delay
	resp_timepeak = 1.3 # delay at peak weighting
	# time range for increasing weighting (s)
	t = np.arange(sht, (resp_timepeak), (resp_timepeak-sht)/50.)
	wpre = np.exp(t) 
	t_all[0:50] = t[:] 
	
	lot = 2.5 # longest delay (s)
	# time range for decreasing weighting (s)
	t = np.arange((resp_timepeak+(lot-resp_timepeak)/50.), (lot), (lot-resp_timepeak)/51.)
	wpro = np.exp(np.flip(t-1.3)) 
	t_all[50::] = t[:] 
	
	# join weighting
	w = np.append(wpre, wpro)
	
	# integrate weight curve
	area = np.trapz(w, t_all)
	# normalise so that integral is one
	w = w/area
	
	return(w, t_all)