import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as si
import pandas as pd


# open saved files
output_by_sim = '/Users/Simon_OMeara/Documents/Manchester/postdoc_stuff/box-model/PyCHAM_Gitw/PyCHAM/output/Example_Run_2016/2016-10-20_V3p3_nowall_Jall_MCM_0p1AF'

# name of file where experiment constants saved (number of size bins and whether wall 
# included)
fname = str(output_by_sim+'/constants')
const = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

num_sb = int(const[0,0]) # number of size bins
num_speci = int(const[0,1]) # number of species
y_mw = (const[1,:])
y_MV = (const[2,:])

# name of file where concentration (molecules/cc (air)) results saved
fname = str(output_by_sim+'/y')
y = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)


# withdraw times
fname = str(output_by_sim+'/t')
t_array = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)

# open measurements
fnamepre = '/Users/Simon_OMeara/Documents/Manchester/postdoc_stuff/box-model/PyCHAM_Gitw/PyCHAM/output/MACobs2016-10-20'
# plot simulated and measured O3

# measured time
fname = str(fnamepre+'/gastime.csv')
time = np.loadtxt(fname,delimiter=',',skiprows=0)

# measured O3
fname = str(fnamepre+'/O3ppb.csv')
O3ppb = np.loadtxt(fname,delimiter=',',skiprows=0)

# measured NO2
fname = str(fnamepre+'/NO2ppb.csv')
NO2ppb = np.loadtxt(fname,delimiter=',',skiprows=0)

# measured NO
fname = str(fnamepre+'/NOppb.csv')
NOppb = np.loadtxt(fname,delimiter=',',skiprows=0)

plt.plot(t_array/3600.0, y[:,1], label='simO3')
plt.plot(t_array/3600.0, y[:,3], label='simNO2')
plt.plot(t_array/3600.0, y[:,2], label='simNO')
# start at 862:2762 for gas-phase concentration arrays for 2016-10-19, [0:1900] for time
# start at 813:2610 for 2016-10-20, [0:1797] for time
plt.plot(time[0:1797]/3600.0, O3ppb[813:2610], label='obsO3')
plt.plot(time[0:1797]/3600.0, NO2ppb[813:2610], label='obsNO2')
plt.plot(time[0:1797]/3600.0, NOppb[813:2610], label='obsNO')
plt.title('MAC 20/10/2016, a-Pinene ozonolysis')
plt.ylabel(r'Gas-phase concentration (ppb)', fontsize=18)
plt.xlabel(r'Time (hours)', fontsize=18)	
plt.legend()
plt.show()

# open tendency to change
# fname = str(output_by_sim+'/NO2_dydt')
# dydt = np.loadtxt(fname,delimiter=',',skiprows=1) # skiprows=1 omits header)
# for i in range(0, len(dydt[0,:])):
# 	if int(i)<(len(dydt[0,:])-2):
# 		plt.plot(t_array/3600.0, dydt[1::, int(i)], label='sim_tend_'+str(int(dydt[0,int(i)])))
# 	elif int(i)==(len(dydt[0,:])-2): # particle partitioning
# 		plt.plot(t_array/3600.0, dydt[1::, int(i)], '--', label='sim_tend_part')
# 	elif int(i)==(len(dydt[0,:])-1): # particle partitioning
# 		plt.plot(t_array/3600.0, dydt[1::, int(i)], '--', label='sim_tend_wall')
# plt.title('MAC 20/10/2016, a-Pinene ozonolysis')
# plt.ylabel(r'Gas-phase concentration tendency (%/hour)', fontsize=18)
# plt.xlabel(r'Time (hours)', fontsize=18)	
# plt.legend(loc="best")
# plt.show()