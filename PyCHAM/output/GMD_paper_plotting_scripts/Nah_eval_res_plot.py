'''For analysing results relating to the computer time for PyCHAM when solving the alpha-pinene ozonolysis scheme, which also provides an evaluation of PyCHAM against Nah et al. 2016 results'''
# Nah et al. 2016: https://www.atmos-chem-phys.net/16/9361/2016/acp-16-9361-2016.pdf
# supplement: https://www.atmos-chem-phys.net/16/9361/2016/acp-16-9361-2016-supplement.pdf

# whilst fitting nucleation parameters need to check the temporal profile of total
# particle volume in addition to total particle number

import matplotlib.pyplot as plt
import os
from retrieve_PyCHAM_outputs import retrieve_outputs as retr
import numpy as np

# set up plots ---------------------------------------------------------
fig, ((ax0)) = plt.subplots(1, 1, figsize=(13, 6))
# ensure sufficient spacing between subplots
# fig.subplots_adjust(top=0.90, bottom = 0.10, hspace=0.5, wspace=0.2)
# ----------------------------------------------------------------------

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

# retrieve results
(num_sb_dict['num_sb0'], num_speci_dict['num_speci0'], Cfac_dict['Cfac0'], y_dict['y0'], N_dict['N0'], sbb_dict['sbb0'], x_dict['x0'], thr_dict['thr0'], _, _, _, _, _, spd_dict['spd0']) = retr(str(cwd + '/Nah_eval_output/Nah_eval_60s_opsl_2sb_medwallloss_4pcELVOC'))

# temporal profile of number concentration of dry particles (#/cm3)
N = N_dict['N0']
# temporal profile of average radius of particle per size bin (um)
r = x_dict['x0']
# convert to temporal profile of total volume of particles (um/cm3)
Vtot = ((4.0/3.0*np.pi*r**3.0)*N).sum(axis=1)
ax0.plot(thr_dict['thr0']*60.0, Vtot)
ax0.set_xlabel('Time (min)')
ax0.set_ylabel('Volume (um/cm3)')
plt.show()

# ----------------------------------------------------------------------------
# check that particle loss to wall dependence on particle diameter matches that in 
# Fig. 1 of Nah et al. 2016
# inflectDp = 2.0e-7
# pwl_xpre = 2.0
# pwl_xpro = 3.0
# inflectk = 2.0e-6
# Dp = np.logspace(-8, -6, 100) # diameter (m)
# 
# Beta = np.zeros((Dp.shape[0], 1))
# Beta[Dp<inflectDp, 0] = 10**((np.log10(inflectDp)-
# 						np.log10(Dp[Dp<inflectDp]))*pwl_xpre+np.log10(inflectk))
# Beta[Dp>=inflectDp, 0] = 10**((np.log10(Dp[Dp>=inflectDp])-
# 						np.log10(inflectDp))*pwl_xpro+np.log10(inflectk))
# ax0.loglog(Dp, Beta)
# plt.show()
# ----------------------------------------------------------------------------