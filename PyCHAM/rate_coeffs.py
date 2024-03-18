##########################################################################################
#                                                                                        #
#    Copyright (C) 2018-2024 Simon O'Meara : simon.omeara@manchester.ac.uk               #
#                                                                                        #
#    All Rights Reserved.                                                                #
#    This file is part of PyCHAM                                                         #
#                                                                                        #
#    PyCHAM is free software: you can redistribute it and/or modify it under             #
#    the terms of the GNU General Public License as published by the Free Software       #
#    Foundation, either version 3 of the License, or (at your option) any later          #
#    version.                                                                            #
#                                                                                        #
#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT               #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              #
#    details.                                                                            #
#                                                                                        #
#    You should have received a copy of the GNU General Public License along with        #
#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                #
#                                                                                        #
##########################################################################################
'''module for calculating reaction rate coefficients (automatically generated)'''
# module to hold expressions for calculating rate coefficients # 
# created at 2024-03-08 15:35:02.839816

import numpy
import photolysisRates

def evaluate_rates(RO2, H2O, TEMP, time, M, N2, O2, Jlen, NO, HO2, NO3, sumt, self):

	# inputs: ------------------------------------------------------------------
	# RO2 - total concentration of alkyl peroxy radicals (# molecules/cm3) 
	# M - third body concentration (# molecules/cm3 (air))
	# N2 - nitrogen concentration (# molecules/cm3 (air))
	# O2 - oxygen concentration (# molecules/cm3 (air))
	# H2O, TEMP: given by the user
	# self.light_stat_now: given by the user and is 0 for lights off and >1 for on
	# reaction rate coefficients and their names parsed in eqn_parser.py 
	# Jlen - number of photolysis reactions
	# self.tf - sunlight transmission factor
	# NO - NO concentration (# molecules/cm3 (air))
	# HO2 - HO2 concentration (# molecules/cm3 (air))
	# NO3 - NO3 concentration (# molecules/cm3 (air))
	# self.tf_UVC - transmission factor for 254 nm wavelength light (0-1) 
	# ------------------------------------------------------------------------

	erf = 0; err_mess = '' # begin assuming no errors

	# calculate any generic reaction rate coefficients given by chemical scheme

	try:
		gprn=0
		gprn += 1 # keep count on reaction number
		KRO2NO=2.7e-12*numpy.exp(360/TEMP) 
		gprn += 1 # keep count on reaction number
		KRO2HO2=2.91e-13*numpy.exp(1300/TEMP) 
		gprn += 1 # keep count on reaction number
		KAPHO2=5.2e-13*numpy.exp(980/TEMP) 
		gprn += 1 # keep count on reaction number
		KAPNO=7.5e-12*numpy.exp(290/TEMP) 
		gprn += 1 # keep count on reaction number
		KRO2NO3=2.3e-12 
		gprn += 1 # keep count on reaction number
		KNO3AL=1.44e-12*numpy.exp(-1862/TEMP) 
		gprn += 1 # keep count on reaction number
		KDEC=1.00e+06 
		gprn += 1 # keep count on reaction number
		KROPRIM=2.50e-14*numpy.exp(-300/TEMP) 
		gprn += 1 # keep count on reaction number
		KROSEC=2.50e-14*numpy.exp(-300/TEMP) 
		gprn += 1 # keep count on reaction number
		KCH3O2=1.03e-13*numpy.exp(365/TEMP) 
		gprn += 1 # keep count on reaction number
		K298CH3O2=3.5e-13 
		gprn += 1 # keep count on reaction number
		K14ISOM1=3.00e7*numpy.exp(-5300/TEMP) 
		gprn += 1 # keep count on reaction number
		ASA=1. 
		gprn += 1 # keep count on reaction number
		KD0=1.10e-05*M*numpy.exp(-10100/TEMP) 
		gprn += 1 # keep count on reaction number
		KDI=1.90e17*numpy.exp(-14100/TEMP) 
		gprn += 1 # keep count on reaction number
		KRD=KD0/KDI 
		gprn += 1 # keep count on reaction number
		FCD=0.30 
		gprn += 1 # keep count on reaction number
		NCD=0.75-1.27*(numpy.log10(FCD)) 
		gprn += 1 # keep count on reaction number
		FD=10**(numpy.log10(FCD)/(1+(numpy.log10(KRD)/NCD)**2)) 
		gprn += 1 # keep count on reaction number
		KBPAN=(KD0*KDI)*FD/(KD0+KDI) 
		gprn += 1 # keep count on reaction number
		KC0=3.28e-28*M*(TEMP/300)**-6.87 
		gprn += 1 # keep count on reaction number
		KCI=1.125e-11*(TEMP/300)**-1.105 
		gprn += 1 # keep count on reaction number
		KRC=KC0/KCI 
		gprn += 1 # keep count on reaction number
		FCC=0.30 
		gprn += 1 # keep count on reaction number
		NC=0.75-1.27*(numpy.log10(FCC)) 
		gprn += 1 # keep count on reaction number
		FC=10**(numpy.log10(FCC)/(1+(numpy.log10(KRC)/NC)**2)) 
		gprn += 1 # keep count on reaction number
		KFPAN=(KC0*KCI)*FC/(KC0+KCI) 
		gprn += 1 # keep count on reaction number
		K10=1.0e-31*M*(TEMP/300)**-1.6 
		gprn += 1 # keep count on reaction number
		K1I=5.0e-11*(TEMP/300)**-0.3 
		gprn += 1 # keep count on reaction number
		KR1=K10/K1I 
		gprn += 1 # keep count on reaction number
		FC1=0.85 
		gprn += 1 # keep count on reaction number
		NC1=0.75-1.27*(numpy.log10(FC1)) 
		gprn += 1 # keep count on reaction number
		F1=10**(numpy.log10(FC1)/(1+(numpy.log10(KR1)/NC1)**2)) 
		gprn += 1 # keep count on reaction number
		KMT01=(K10*K1I)*F1/(K10+K1I) 
		gprn += 1 # keep count on reaction number
		K20=1.3e-31*M*(TEMP/300)**-1.5 
		gprn += 1 # keep count on reaction number
		K2I=2.3e-11*(TEMP/300)**0.24 
		gprn += 1 # keep count on reaction number
		KR2=K20/K2I 
		gprn += 1 # keep count on reaction number
		FC2=0.6 
		gprn += 1 # keep count on reaction number
		NC2=0.75-1.27*(numpy.log10(FC2)) 
		gprn += 1 # keep count on reaction number
		F2=10**(numpy.log10(FC2)/(1+(numpy.log10(KR2)/NC2)**2)) 
		gprn += 1 # keep count on reaction number
		KMT02=(K20*K2I)*F2/(K20+K2I) 
		gprn += 1 # keep count on reaction number
		K30=3.6e-30*M*(TEMP/300)**-4.1 
		gprn += 1 # keep count on reaction number
		K3I=1.9e-12*(TEMP/300)**0.2 
		gprn += 1 # keep count on reaction number
		KR3=K30/K3I 
		gprn += 1 # keep count on reaction number
		FC3=0.35 
		gprn += 1 # keep count on reaction number
		NC3=0.75-1.27*(numpy.log10(FC3)) 
		gprn += 1 # keep count on reaction number
		F3=10**(numpy.log10(FC3)/(1+(numpy.log10(KR3)/NC3)**2)) 
		gprn += 1 # keep count on reaction number
		KMT03=(K30*K3I)*F3/(K30+K3I) 
		gprn += 1 # keep count on reaction number
		K40=1.3e-3*M*(TEMP/300)**-3.5*numpy.exp(-11000/TEMP) 
		gprn += 1 # keep count on reaction number
		K4I=9.7e+14*(TEMP/300)**0.1*numpy.exp(-11080/TEMP) 
		gprn += 1 # keep count on reaction number
		KR4=K40/K4I 
		gprn += 1 # keep count on reaction number
		FC4=0.35 
		gprn += 1 # keep count on reaction number
		NC4=0.75-1.27*(numpy.log10(FC4)) 
		gprn += 1 # keep count on reaction number
		F4=10**(numpy.log10(FC4)/(1+(numpy.log10(KR4)/NC4)**2)) 
		gprn += 1 # keep count on reaction number
		KMT04=(K40*K4I)*F4/(K40+K4I) 
		gprn += 1 # keep count on reaction number
		KMT05=1.44e-13*(1+(M/4.2e+19)) 
		gprn += 1 # keep count on reaction number
		KMT06=1+(1.40e-21*numpy.exp(2200/TEMP)*H2O) 
		gprn += 1 # keep count on reaction number
		K70=7.4e-31*M*(TEMP/300)**-2.4 
		gprn += 1 # keep count on reaction number
		K7I=3.3e-11*(TEMP/300)**-0.3 
		gprn += 1 # keep count on reaction number
		KR7=K70/K7I 
		gprn += 1 # keep count on reaction number
		FC7=0.81 
		gprn += 1 # keep count on reaction number
		NC7=0.75-1.27*(numpy.log10(FC7)) 
		gprn += 1 # keep count on reaction number
		F7=10**(numpy.log10(FC7)/(1+(numpy.log10(KR7)/NC7)**2)) 
		gprn += 1 # keep count on reaction number
		KMT07=(K70*K7I)*F7/(K70+K7I) 
		gprn += 1 # keep count on reaction number
		K80=3.2e-30*M*(TEMP/300)**-4.5 
		gprn += 1 # keep count on reaction number
		K8I=3.0e-11 
		gprn += 1 # keep count on reaction number
		KR8=K80/K8I 
		gprn += 1 # keep count on reaction number
		FC8=0.41 
		gprn += 1 # keep count on reaction number
		NC8=0.75-1.27*(numpy.log10(FC8)) 
		gprn += 1 # keep count on reaction number
		F8=10**(numpy.log10(FC8)/(1+(numpy.log10(KR8)/NC8)**2)) 
		gprn += 1 # keep count on reaction number
		KMT08=(K80*K8I)*F8/(K80+K8I) 
		gprn += 1 # keep count on reaction number
		K90=1.4e-31*M*(TEMP/300)**-3.1 
		gprn += 1 # keep count on reaction number
		K9I=4.0e-12 
		gprn += 1 # keep count on reaction number
		KR9=K90/K9I 
		gprn += 1 # keep count on reaction number
		FC9=0.4 
		gprn += 1 # keep count on reaction number
		NC9=0.75-1.27*(numpy.log10(FC9)) 
		gprn += 1 # keep count on reaction number
		F9=10**(numpy.log10(FC9)/(1+(numpy.log10(KR9)/NC9)**2)) 
		gprn += 1 # keep count on reaction number
		KMT09=(K90*K9I)*F9/(K90+K9I) 
		gprn += 1 # keep count on reaction number
		K100=4.10e-05*M*numpy.exp(-10650/TEMP) 
		gprn += 1 # keep count on reaction number
		K10I=6.0e+15*numpy.exp(-11170/TEMP) 
		gprn += 1 # keep count on reaction number
		KR10=K100/K10I 
		gprn += 1 # keep count on reaction number
		FC10=0.4 
		gprn += 1 # keep count on reaction number
		NC10=0.75-1.27*(numpy.log10(FC10)) 
		gprn += 1 # keep count on reaction number
		F10=10**(numpy.log10(FC10)/(1+(numpy.log10(KR10)/NC10)**2)) 
		gprn += 1 # keep count on reaction number
		KMT10=(K100*K10I)*F10/(K100+K10I) 
		gprn += 1 # keep count on reaction number
		K1=2.40e-14*numpy.exp(460/TEMP) 
		gprn += 1 # keep count on reaction number
		K3=6.50e-34*numpy.exp(1335/TEMP) 
		gprn += 1 # keep count on reaction number
		K4=2.70e-17*numpy.exp(2199/TEMP) 
		gprn += 1 # keep count on reaction number
		K2=(K3*M)/(1+(K3*M/K4)) 
		gprn += 1 # keep count on reaction number
		KMT11=K1+K2 
		gprn += 1 # keep count on reaction number
		K120=2.5e-31*M*(TEMP/300)**-2.6 
		gprn += 1 # keep count on reaction number
		K12I=2.0e-12 
		gprn += 1 # keep count on reaction number
		KR12=K120/K12I 
		gprn += 1 # keep count on reaction number
		FC12=0.53 
		gprn += 1 # keep count on reaction number
		NC12=0.75-1.27*(numpy.log10(FC12)) 
		gprn += 1 # keep count on reaction number
		F12=10**(numpy.log10(FC12)/(1.0+(numpy.log10(KR12)/NC12)**2)) 
		gprn += 1 # keep count on reaction number
		KMT12=(K120*K12I*F12)/(K120+K12I) 
		gprn += 1 # keep count on reaction number
		K130=2.5e-30*M*(TEMP/300)**-5.5 
		gprn += 1 # keep count on reaction number
		K13I=1.8e-11 
		gprn += 1 # keep count on reaction number
		KR13=K130/K13I 
		gprn += 1 # keep count on reaction number
		FC13=0.36 
		gprn += 1 # keep count on reaction number
		NC13=0.75-1.27*(numpy.log10(FC13)) 
		gprn += 1 # keep count on reaction number
		F13=10**(numpy.log10(FC13)/(1+(numpy.log10(KR13)/NC13)**2)) 
		gprn += 1 # keep count on reaction number
		KMT13=(K130*K13I)*F13/(K130+K13I) 
		gprn += 1 # keep count on reaction number
		K140=9.0e-5*numpy.exp(-9690/TEMP)*M 
		gprn += 1 # keep count on reaction number
		K14I=1.1e+16*numpy.exp(-10560/TEMP) 
		gprn += 1 # keep count on reaction number
		KR14=K140/K14I 
		gprn += 1 # keep count on reaction number
		FC14=0.36 
		gprn += 1 # keep count on reaction number
		NC14=0.75-1.27*(numpy.log10(FC14)) 
		gprn += 1 # keep count on reaction number
		F14=10**(numpy.log10(FC14)/(1+(numpy.log10(KR14)/NC14)**2)) 
		gprn += 1 # keep count on reaction number
		KMT14=(K140*K14I)*F14/(K140+K14I) 
		gprn += 1 # keep count on reaction number
		K150=8.6e-29*M*(TEMP/300)**-3.1 
		gprn += 1 # keep count on reaction number
		K15I=9.0e-12*(TEMP/300)**-0.85 
		gprn += 1 # keep count on reaction number
		KR15=K150/K15I 
		gprn += 1 # keep count on reaction number
		FC15=0.48 
		gprn += 1 # keep count on reaction number
		NC15=0.75-1.27*(numpy.log10(FC15)) 
		gprn += 1 # keep count on reaction number
		F15=10**(numpy.log10(FC15)/(1+(numpy.log10(KR15)/NC15)**2)) 
		gprn += 1 # keep count on reaction number
		KMT15=(K150*K15I)*F15/(K150+K15I) 
		gprn += 1 # keep count on reaction number
		K160=8e-27*M*(TEMP/300)**-3.5 
		gprn += 1 # keep count on reaction number
		K16I=3.0e-11*(TEMP/300)**-1 
		gprn += 1 # keep count on reaction number
		KR16=K160/K16I 
		gprn += 1 # keep count on reaction number
		FC16=0.5 
		gprn += 1 # keep count on reaction number
		NC16=0.75-1.27*(numpy.log10(FC16)) 
		gprn += 1 # keep count on reaction number
		F16=10**(numpy.log10(FC16)/(1+(numpy.log10(KR16)/NC16)**2)) 
		gprn += 1 # keep count on reaction number
		KMT16=(K160*K16I)*F16/(K160+K16I) 
		gprn += 1 # keep count on reaction number
		K170=5.0e-30*M*(TEMP/300)**-1.5 
		gprn += 1 # keep count on reaction number
		K17I=1.0e-12 
		gprn += 1 # keep count on reaction number
		KR17=K170/K17I 
		gprn += 1 # keep count on reaction number
		FC17=0.17*numpy.exp(-51/TEMP)+numpy.exp(-TEMP/204) 
		gprn += 1 # keep count on reaction number
		NC17=0.75-1.27*(numpy.log10(FC17)) 
		gprn += 1 # keep count on reaction number
		F17=10**(numpy.log10(FC17)/(1.0+(numpy.log10(KR17)/NC17)**2)) 
		gprn += 1 # keep count on reaction number
		KMT17=(K170*K17I*F17)/(K170+K17I) 
		gprn += 1 # keep count on reaction number
		KMT18=9.5e-39*O2*numpy.exp(5270/TEMP)/(1+7.5e-29*O2*numpy.exp(5610/TEMP)) 
		gprn += 1 # keep count on reaction number
		KPPN0=1.7e-03*numpy.exp(-11280/TEMP)*M 
		gprn += 1 # keep count on reaction number
		KPPNI=8.3e+16*numpy.exp(-13940/TEMP) 
		gprn += 1 # keep count on reaction number
		KRPPN=KPPN0/KPPNI 
		gprn += 1 # keep count on reaction number
		FCPPN=0.36 
		gprn += 1 # keep count on reaction number
		NCPPN=0.75-1.27*(numpy.log10(FCPPN)) 
		gprn += 1 # keep count on reaction number
		FPPN=10**(numpy.log10(FCPPN)/(1+(numpy.log10(KRPPN)/NCPPN)**2)) 
		gprn += 1 # keep count on reaction number
		KBPPN=(KPPN0*KPPNI)*FPPN/(KPPN0+KPPNI) 

	except:
		erf = 1 # flag error
		err_mess = str('Error: generic reaction rates failed to be calculated inside rate_coeffs.py at number ' + str(gprn) + ', please check chemical scheme and associated chemical scheme markers, which are stated in the model variables input file') # error message
		return([], erf, err_mess)
	# estimate and append photolysis rates
	J = photolysisRates.PhotolysisCalculation(TEMP, Jlen, sumt, self)

	if (self.light_stat_now == 0):
		J = [0]*len(J)
	rate_values = numpy.zeros((159))
	
	# if reactions have been found in the chemical scheme
	# gas-phase reactions
	gprn = 0 # keep count on reaction number
	try:
		gprn += 1 # keep count on reaction number
		rate_values[0] = 4.50e-11
		gprn += 1 # keep count on reaction number
		rate_values[1] = 5.6e-34*N2*(TEMP/300)**-2.6*O2
		gprn += 1 # keep count on reaction number
		rate_values[2] = 6.0e-34*O2*(TEMP/300)**-2.6*O2
		gprn += 1 # keep count on reaction number
		rate_values[3] = 8.0e-12*numpy.exp(-2060/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[4] = KMT01
		gprn += 1 # keep count on reaction number
		rate_values[5] = 5.5e-12*numpy.exp(188/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[6] = KMT02
		gprn += 1 # keep count on reaction number
		rate_values[7] = 3.2e-11*numpy.exp(67/TEMP)*O2
		gprn += 1 # keep count on reaction number
		rate_values[8] = 2.0e-11*numpy.exp(130/TEMP)*N2
		gprn += 1 # keep count on reaction number
		rate_values[9] = 1.4e-12*numpy.exp(-1310/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[10] = 1.4e-13*numpy.exp(-2470/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[11] = 3.3e-39*numpy.exp(530/TEMP)*O2
		gprn += 1 # keep count on reaction number
		rate_values[12] = 1.8e-11*numpy.exp(110/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[13] = 4.50e-14*numpy.exp(-1260/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[14] = KMT03
		gprn += 1 # keep count on reaction number
		rate_values[15] = 2.14e-10*H2O
		gprn += 1 # keep count on reaction number
		rate_values[16] = 1.70e-12*numpy.exp(-940/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[17] = 7.7e-12*numpy.exp(-2100/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[18] = KMT05
		gprn += 1 # keep count on reaction number
		rate_values[19] = 2.9e-12*numpy.exp(-160/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[20] = 2.03e-16*(TEMP/300)**4.57*numpy.exp(693/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[21] = 4.8e-11*numpy.exp(250/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[22] = 2.20e-13*KMT06*numpy.exp(600/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[23] = 1.90e-33*M*KMT06*numpy.exp(980/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[24] = KMT07
		gprn += 1 # keep count on reaction number
		rate_values[25] = KMT08
		gprn += 1 # keep count on reaction number
		rate_values[26] = 2.0e-11
		gprn += 1 # keep count on reaction number
		rate_values[27] = 3.45e-12*numpy.exp(270/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[28] = KMT09
		gprn += 1 # keep count on reaction number
		rate_values[29] = 3.2e-13*numpy.exp(690/TEMP)*1.0
		gprn += 1 # keep count on reaction number
		rate_values[30] = 4.0e-12
		gprn += 1 # keep count on reaction number
		rate_values[31] = 2.5e-12*numpy.exp(260/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[32] = KMT11
		gprn += 1 # keep count on reaction number
		rate_values[33] = 4.0e-32*numpy.exp(-1000/TEMP)*M
		gprn += 1 # keep count on reaction number
		rate_values[34] = KMT12
		gprn += 1 # keep count on reaction number
		rate_values[35] = 1.3e-12*numpy.exp(-330/TEMP)*O2
		gprn += 1 # keep count on reaction number
		rate_values[36] = 6.00e-06
		gprn += 1 # keep count on reaction number
		rate_values[37] = 4.00e-04
		gprn += 1 # keep count on reaction number
		rate_values[38] = 1.20e-15*H2O
		gprn += 1 # keep count on reaction number
		rate_values[39] = J[1]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[40] = J[2]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[41] = J[3]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[42] = J[4]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[43] = J[5]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[44] = J[6]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[45] = J[7]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[46] = J[8]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[47] = KMT04
		gprn += 1 # keep count on reaction number
		rate_values[48] = KMT10
		gprn += 1 # keep count on reaction number
		rate_values[49] = 1.9e-13*numpy.exp(580/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[50] = 1.2e-11*numpy.exp(-280/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[51] = 9.4e-11*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[52] = 1.6e-17
		gprn += 1 # keep count on reaction number
		rate_values[53] = KMT18
		gprn += 1 # keep count on reaction number
		rate_values[54] = 2.74e+7*numpy.exp(-5950/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[55] = KRO2HO2*0.387
		gprn += 1 # keep count on reaction number
		rate_values[56] = 4.9e-12*numpy.exp(260/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[57] = KRO2NO3
		gprn += 1 # keep count on reaction number
		rate_values[58] = KRO2NO
		gprn += 1 # keep count on reaction number
		rate_values[59] = 8.90e+10*numpy.exp(-6040/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[60] = 7.03e-11
		gprn += 1 # keep count on reaction number
		rate_values[61] = J[41]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[62] = KDEC
		gprn += 1 # keep count on reaction number
		rate_values[63] = 2.78e-11
		gprn += 1 # keep count on reaction number
		rate_values[64] = 1.11e-11
		gprn += 1 # keep count on reaction number
		rate_values[65] = J[15]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[66] = 4.2e+7*numpy.exp(-5390/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[67] = 4.9e-12*numpy.exp(260/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[68] = 1.13e-13*numpy.exp(1300/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[69] = J[41]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[70] = 1.36e-11
		gprn += 1 # keep count on reaction number
		rate_values[71] = 1.4e-12*(5.6e-3*TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[72] = 1.4e-12*(2.24e-3*TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[73] = J[41]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[74] = J[14]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[75] = 1.15e-12*numpy.exp(430/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[76] = 6.00e-11*numpy.exp(240/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[77] = 4.0e-13
		gprn += 1 # keep count on reaction number
		rate_values[78] = 1.2e-11
		gprn += 1 # keep count on reaction number
		rate_values[79] = 4.40e-14
		gprn += 1 # keep count on reaction number
		rate_values[80] = 6.10e-12*numpy.exp(800/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[81] = 6.00e-11*numpy.exp(240/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[82] = 1.15e-12*numpy.exp(430/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[83] = 1.20e-16*numpy.exp(1580/TEMP)*O2
		gprn += 1 # keep count on reaction number
		rate_values[84] = J[11]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[85] = J[12]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[86] = 5.5e-16
		gprn += 1 # keep count on reaction number
		rate_values[87] = 5.4e-12*numpy.exp(135/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[88] = KRO2HO2*0.387
		gprn += 1 # keep count on reaction number
		rate_values[89] = KRO2NO
		gprn += 1 # keep count on reaction number
		rate_values[90] = KRO2NO3
		gprn += 1 # keep count on reaction number
		rate_values[91] = 2.00e-12*RO2*0.2
		gprn += 1 # keep count on reaction number
		rate_values[92] = 2.00e-12*RO2*0.6
		gprn += 1 # keep count on reaction number
		rate_values[93] = 2.00e-12*RO2*0.2
		gprn += 1 # keep count on reaction number
		rate_values[94] = 9.00e-11
		gprn += 1 # keep count on reaction number
		rate_values[95] = 3.8e-13*numpy.exp(780/TEMP)*(1-1/(1+498*numpy.exp(-1160/TEMP)))
		gprn += 1 # keep count on reaction number
		rate_values[96] = 2.3e-12*numpy.exp(360/TEMP)*0.001
		gprn += 1 # keep count on reaction number
		rate_values[97] = 2.3e-12*numpy.exp(360/TEMP)*0.999
		gprn += 1 # keep count on reaction number
		rate_values[98] = KMT13
		gprn += 1 # keep count on reaction number
		rate_values[99] = 1.2e-12
		gprn += 1 # keep count on reaction number
		rate_values[100] = 2*KCH3O2*RO2*7.18*numpy.exp(-885/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[101] = 2*KCH3O2*RO2*0.5*(1-7.18*numpy.exp(-885/TEMP))
		gprn += 1 # keep count on reaction number
		rate_values[102] = 2*KCH3O2*RO2*0.5*(1-7.18*numpy.exp(-885/TEMP))
		gprn += 1 # keep count on reaction number
		rate_values[103] = 1.20e-11*0.25
		gprn += 1 # keep count on reaction number
		rate_values[104] = 1.20e-11*0.75
		gprn += 1 # keep count on reaction number
		rate_values[105] = 4.00e-13
		gprn += 1 # keep count on reaction number
		rate_values[106] = 3.12e-16*numpy.exp(1580/TEMP)*O2
		gprn += 1 # keep count on reaction number
		rate_values[107] = 1.1e-11
		gprn += 1 # keep count on reaction number
		rate_values[108] = 2.2e-11
		gprn += 1 # keep count on reaction number
		rate_values[109] = 5.60e+16*numpy.exp(-10870/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[110] = 3.50e+10*numpy.exp(-3560/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[111] = 1.26e-12
		gprn += 1 # keep count on reaction number
		rate_values[112] = 3.60e-12
		gprn += 1 # keep count on reaction number
		rate_values[113] = J[41]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[114] = KDEC
		gprn += 1 # keep count on reaction number
		rate_values[115] = 1.78e-12
		gprn += 1 # keep count on reaction number
		rate_values[116] = J[15]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[117] = 5.23e-13
		gprn += 1 # keep count on reaction number
		rate_values[118] = 1.40e-13
		gprn += 1 # keep count on reaction number
		rate_values[119] = J[41]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[120] = 5.3e-12*numpy.exp(190/TEMP)*0.6
		gprn += 1 # keep count on reaction number
		rate_values[121] = 5.3e-12*numpy.exp(190/TEMP)*0.4
		gprn += 1 # keep count on reaction number
		rate_values[122] = J[51]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[123] = 4.0e-13*numpy.exp(-845/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[124] = 7.2e-14*numpy.exp(-1080/TEMP)*O2
		gprn += 1 # keep count on reaction number
		rate_values[125] = KMT14
		gprn += 1 # keep count on reaction number
		rate_values[126] = 2.85e-12*numpy.exp(-345/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[127] = 3.00e-13
		gprn += 1 # keep count on reaction number
		rate_values[128] = 5.00e+13*numpy.exp(-9673/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[129] = 1.03e-16*numpy.exp(1580/TEMP)*O2
		gprn += 1 # keep count on reaction number
		rate_values[130] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		rate_values[131] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		rate_values[132] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		rate_values[133] = 1.00e-11
		gprn += 1 # keep count on reaction number
		rate_values[134] = 1.20e-12*(TEMP/300)**-0.9
		gprn += 1 # keep count on reaction number
		rate_values[135] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		rate_values[136] = 9.10e+10*numpy.exp(-3560/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[137] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		rate_values[138] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		rate_values[139] = 5.00e-11
		gprn += 1 # keep count on reaction number
		rate_values[140] = 5.00e+13*numpy.exp(-9946/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[141] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		rate_values[142] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		rate_values[143] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		rate_values[144] = 1.00e-11
		gprn += 1 # keep count on reaction number
		rate_values[145] = 1.20e-12*(TEMP/300)**-0.9
		gprn += 1 # keep count on reaction number
		rate_values[146] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		rate_values[147] = 3.01e+10*numpy.exp(-3560/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[148] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		rate_values[149] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		rate_values[150] = 9.00e-11
		gprn += 1 # keep count on reaction number
		rate_values[151] = J[41]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[152] = 1.00e-11
		gprn += 1 # keep count on reaction number
		rate_values[153] = 5.40e+16*numpy.exp(-13112/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[154] = 2.24e-14
		gprn += 1 # keep count on reaction number
		rate_values[155] = 3.60e-12
		gprn += 1 # keep count on reaction number
		rate_values[156] = J[41]*ASA
		gprn += 1 # keep count on reaction number
		rate_values[157] = 3.60e-13
		gprn += 1 # keep count on reaction number
		rate_values[158] = 5.40e+16*numpy.exp(-13112/TEMP)
	except:
		erf = 1 # flag error
		err_mess = (str('Error: Could not calculate rate coefficient for equation number ' + str(gprn)))
	
	# aqueous-phase reactions
	
	# surface (e.g. wall) reactions
	
	return(rate_values, erf, err_mess)
