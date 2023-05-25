##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2023 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
#                                                                                       											 #
#    All Rights Reserved.                                                                									 #
#    This file is part of PyCHAM                                                         									 #
#                                                                                        											 #
#    PyCHAM is free software: you can redistribute it and/or modify it under              						 #
#    the terms of the GNU General Public License as published by the Free Software       					 #
#    Foundation, either version 3 of the License, or (at your option) any later          						 #
#    version.                                                                            										 #
#                                                                                        											 #
#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT                						 #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       			 #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              				 #
#    details.                                                                            										 #
#                                                                                        											 #
#    You should have received a copy of the GNU General Public License along with        					 #
#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                 							 #
#                                                                                        											 #
##########################################################################################
'''module for calculating reaction rate coefficients (automatically generated)'''
# module to hold expressions for calculating rate coefficients # 
# created at 2023-05-25 16:13:48.594201

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
		KRO2NO=2.7e-12*numpy.exp(360/TEMP) 
		KRO2HO2=2.91e-13*numpy.exp(1300/TEMP) 
		KAPHO2=5.2e-13*numpy.exp(980/TEMP) 
		KAPNO=7.5e-12*numpy.exp(290/TEMP) 
		KRO2NO3=2.3e-12 
		KNO3AL=1.4e-12*numpy.exp(-1860/TEMP) 
		KDEC=1.00e+06 
		KROPRIM=2.50e-14*numpy.exp(-300/TEMP) 
		KROSEC=2.50e-14*numpy.exp(-300/TEMP) 
		KCH3O2=1.03e-13*numpy.exp(365/TEMP) 
		K298CH3O2=3.5e-13 
		K14ISOM1=3.00e7*numpy.exp(-5300/TEMP) 
		KD0=1.10e-05*M*numpy.exp(-10100/TEMP) 
		KDI=1.90e17*numpy.exp(-14100/TEMP) 
		KRD=KD0/KDI 
		FCD=0.30 
		NCD=0.75-1.27*(numpy.log10(FCD)) 
		FD=10**(numpy.log10(FCD)/(1+(numpy.log10(KRD)/NCD)**2)) 
		KBPAN=(KD0*KDI)*FD/(KD0+KDI) 
		KC0=3.28e-28*M*(TEMP/300)**-6.87 
		KCI=1.125e-11*(TEMP/300)**-1.105 
		KRC=KC0/KCI 
		FCC=0.30 
		NC=0.75-1.27*(numpy.log10(FCC)) 
		FC=10**(numpy.log10(FCC)/(1+(numpy.log10(KRC)/NC)**2)) 
		KFPAN=(KC0*KCI)*FC/(KC0+KCI) 
		K10=1.0e-31*M*(TEMP/300)**-1.6 
		K1I=5.0e-11*(TEMP/300)**-0.3 
		KR1=K10/K1I 
		FC1=0.85 
		NC1=0.75-1.27*(numpy.log10(FC1)) 
		F1=10**(numpy.log10(FC1)/(1+(numpy.log10(KR1)/NC1)**2)) 
		KMT01=(K10*K1I)*F1/(K10+K1I) 
		K20=1.3e-31*M*(TEMP/300)**-1.5 
		K2I=2.3e-11*(TEMP/300)**0.24 
		KR2=K20/K2I 
		FC2=0.6 
		NC2=0.75-1.27*(numpy.log10(FC2)) 
		F2=10**(numpy.log10(FC2)/(1+(numpy.log10(KR2)/NC2)**2)) 
		KMT02=(K20*K2I)*F2/(K20+K2I) 
		K30=3.6e-30*M*(TEMP/300)**-4.1 
		K3I=1.9e-12*(TEMP/300)**0.2 
		KR3=K30/K3I 
		FC3=0.35 
		NC3=0.75-1.27*(numpy.log10(FC3)) 
		F3=10**(numpy.log10(FC3)/(1+(numpy.log10(KR3)/NC3)**2)) 
		KMT03=(K30*K3I)*F3/(K30+K3I) 
		K40=1.3e-3*M*(TEMP/300)**-3.5*numpy.exp(-11000/TEMP) 
		K4I=9.7e+14*(TEMP/300)**0.1*numpy.exp(-11080/TEMP) 
		KR4=K40/K4I 
		FC4=0.35 
		NC4=0.75-1.27*(numpy.log10(FC4)) 
		F4=10**(numpy.log10(FC4)/(1+(numpy.log10(KR4)/NC4)**2)) 
		KMT04=(K40*K4I)*F4/(K40+K4I) 
		KMT05=1.44e-13*(1+(M/4.2e+19)) 
		KMT06=1+(1.40e-21*numpy.exp(2200/TEMP)*H2O) 
		K70=7.4e-31*M*(TEMP/300)**-2.4 
		K7I=3.3e-11*(TEMP/300)**-0.3 
		KR7=K70/K7I 
		FC7=0.81 
		NC7=0.75-1.27*(numpy.log10(FC7)) 
		F7=10**(numpy.log10(FC7)/(1+(numpy.log10(KR7)/NC7)**2)) 
		KMT07=(K70*K7I)*F7/(K70+K7I) 
		K80=3.2e-30*M*(TEMP/300)**-4.5 
		K8I=3.0e-11 
		KR8=K80/K8I 
		FC8=0.41 
		NC8=0.75-1.27*(numpy.log10(FC8)) 
		F8=10**(numpy.log10(FC8)/(1+(numpy.log10(KR8)/NC8)**2)) 
		KMT08=(K80*K8I)*F8/(K80+K8I) 
		K90=1.4e-31*M*(TEMP/300)**-3.1 
		K9I=4.0e-12 
		KR9=K90/K9I 
		FC9=0.4 
		NC9=0.75-1.27*(numpy.log10(FC9)) 
		F9=10**(numpy.log10(FC9)/(1+(numpy.log10(KR9)/NC9)**2)) 
		KMT09=(K90*K9I)*F9/(K90+K9I) 
		K100=4.10e-05*M*numpy.exp(-10650/TEMP) 
		K10I=6.0e+15*numpy.exp(-11170/TEMP) 
		KR10=K100/K10I 
		FC10=0.4 
		NC10=0.75-1.27*(numpy.log10(FC10)) 
		F10=10**(numpy.log10(FC10)/(1+(numpy.log10(KR10)/NC10)**2)) 
		KMT10=(K100*K10I)*F10/(K100+K10I) 
		K1=2.40e-14*numpy.exp(460/TEMP) 
		K3=6.50e-34*numpy.exp(1335/TEMP) 
		K4=2.70e-17*numpy.exp(2199/TEMP) 
		K2=(K3*M)/(1+(K3*M/K4)) 
		KMT11=K1+K2 
		K120=2.5e-31*M*(TEMP/300)**-2.6 
		K12I=2.0e-12 
		KR12=K120/K12I 
		FC12=0.53 
		NC12=0.75-1.27*(numpy.log10(FC12)) 
		F12=10**(numpy.log10(FC12)/(1.0+(numpy.log10(KR12)/NC12)**2)) 
		KMT12=(K120*K12I*F12)/(K120+K12I) 
		K130=2.5e-30*M*(TEMP/300)**-5.5 
		K13I=1.8e-11 
		KR13=K130/K13I 
		FC13=0.36 
		NC13=0.75-1.27*(numpy.log10(FC13)) 
		F13=10**(numpy.log10(FC13)/(1+(numpy.log10(KR13)/NC13)**2)) 
		KMT13=(K130*K13I)*F13/(K130+K13I) 
		K140=9.0e-5*numpy.exp(-9690/TEMP)*M 
		K14I=1.1e+16*numpy.exp(-10560/TEMP) 
		KR14=K140/K14I 
		FC14=0.36 
		NC14=0.75-1.27*(numpy.log10(FC14)) 
		F14=10**(numpy.log10(FC14)/(1+(numpy.log10(KR14)/NC14)**2)) 
		KMT14=(K140*K14I)*F14/(K140+K14I) 
		K150=8.6e-29*M*(TEMP/300)**-3.1 
		K15I=9.0e-12*(TEMP/300)**-0.85 
		KR15=K150/K15I 
		FC15=0.48 
		NC15=0.75-1.27*(numpy.log10(FC15)) 
		F15=10**(numpy.log10(FC15)/(1+(numpy.log10(KR15)/NC15)**2)) 
		KMT15=(K150*K15I)*F15/(K150+K15I) 
		K160=8e-27*M*(TEMP/300)**-3.5 
		K16I=3.0e-11*(TEMP/300)**-1 
		KR16=K160/K16I 
		FC16=0.5 
		NC16=0.75-1.27*(numpy.log10(FC16)) 
		F16=10**(numpy.log10(FC16)/(1+(numpy.log10(KR16)/NC16)**2)) 
		KMT16=(K160*K16I)*F16/(K160+K16I) 
		K170=5.0e-30*M*(TEMP/300)**-1.5 
		K17I=1.0e-12 
		KR17=K170/K17I 
		FC17=0.17*numpy.exp(-51/TEMP)+numpy.exp(-TEMP/204) 
		NC17=0.75-1.27*(numpy.log10(FC17)) 
		F17=10**(numpy.log10(FC17)/(1.0+(numpy.log10(KR17)/NC17)**2)) 
		KMT17=(K170*K17I*F17)/(K170+K17I) 
		KMT18=9.5e-39*O2*numpy.exp(5270/TEMP)/(1+7.5e-29*O2*numpy.exp(5610/TEMP)) 
		KPPN0=1.7e-03*numpy.exp(-11280/TEMP)*M 
		KPPNI=8.3e+16*numpy.exp(-13940/TEMP) 
		KRPPN=KPPN0/KPPNI 
		FCPPN=0.36 
		NCPPN=0.75-1.27*(numpy.log10(FCPPN)) 
		FPPN=10**(numpy.log10(FCPPN)/(1+(numpy.log10(KRPPN)/NCPPN)**2)) 
		KBPPN=(KPPN0*KPPNI)*FCPPN/(KPPN0+KPPNI) 

	except:
		erf = 1 # flag error
		err_mess = 'Error: generic reaction rates failed to be calculated inside rate_coeffs.py, please check chemical scheme and associated chemical scheme markers, which are stated in the model variables input file' # error message
	# estimate and append photolysis rates
	J = photolysisRates.PhotolysisCalculation(TEMP, Jlen, sumt, self)

	if (self.light_stat_now == 0):
		J = [0]*len(J)
	rate_values = numpy.zeros((68))
	
	# if reactions have been found in the chemical scheme
	# gas-phase reactions
	gprn = 0 # keep count on reaction number
	try:
		gprn += 1 # keep count on reaction number
		rate_values[0] = 5.6e-34*N2*(TEMP/300)**-2.6*O2+6.0e-34*O2*(TEMP/300)**-2.6*O2
		gprn += 1 # keep count on reaction number
		rate_values[1] = 8.0e-12*numpy.exp(-2060/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[2] = KMT01
		gprn += 1 # keep count on reaction number
		rate_values[3] = 5.5e-12*numpy.exp(188/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[4] = KMT02
		gprn += 1 # keep count on reaction number
		rate_values[5] = 3.2e-11*numpy.exp(67/TEMP)*O2+2.0e-11*numpy.exp(130/TEMP)*N2
		gprn += 1 # keep count on reaction number
		rate_values[6] = 1.4e-12*numpy.exp(-1310/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[7] = 1.4e-13*numpy.exp(-2470/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[8] = 3.3e-39*numpy.exp(530/TEMP)*O2
		gprn += 1 # keep count on reaction number
		rate_values[9] = 1.8e-11*numpy.exp(110/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[10] = 4.50e-14*numpy.exp(-1260/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[11] = KMT03
		gprn += 1 # keep count on reaction number
		rate_values[12] = 2.14e-10*H2O
		gprn += 1 # keep count on reaction number
		rate_values[13] = 1.70e-12*numpy.exp(-940/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[14] = 7.7e-12*numpy.exp(-2100/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[15] = KMT05
		gprn += 1 # keep count on reaction number
		rate_values[16] = 2.9e-12*numpy.exp(-160/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[17] = 2.03e-16*(TEMP/300)**4.57*numpy.exp(693/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[18] = 4.8e-11*numpy.exp(250/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[19] = 2.20e-13*KMT06*numpy.exp(600/TEMP)+1.90e-33*M*KMT06*numpy.exp(980/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[20] = KMT07
		gprn += 1 # keep count on reaction number
		rate_values[21] = KMT08
		gprn += 1 # keep count on reaction number
		rate_values[22] = 2.0e-11
		gprn += 1 # keep count on reaction number
		rate_values[23] = 3.45e-12*numpy.exp(270/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[24] = KMT09
		gprn += 1 # keep count on reaction number
		rate_values[25] = 3.2e-13*numpy.exp(690/TEMP)*1.0
		gprn += 1 # keep count on reaction number
		rate_values[26] = 4.0e-12
		gprn += 1 # keep count on reaction number
		rate_values[27] = 2.5e-12*numpy.exp(260/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[28] = KMT11
		gprn += 1 # keep count on reaction number
		rate_values[29] = 4.0e-32*numpy.exp(-1000/TEMP)*M
		gprn += 1 # keep count on reaction number
		rate_values[30] = KMT12
		gprn += 1 # keep count on reaction number
		rate_values[31] = 1.3e-12*numpy.exp(-330/TEMP)*O2
		gprn += 1 # keep count on reaction number
		rate_values[32] = 6.00e-06
		gprn += 1 # keep count on reaction number
		rate_values[33] = 4.00e-04
		gprn += 1 # keep count on reaction number
		rate_values[34] = 1.20e-15*H2O
		gprn += 1 # keep count on reaction number
		rate_values[35] = J[1]
		gprn += 1 # keep count on reaction number
		rate_values[36] = J[2]
		gprn += 1 # keep count on reaction number
		rate_values[37] = J[3]
		gprn += 1 # keep count on reaction number
		rate_values[38] = J[4]
		gprn += 1 # keep count on reaction number
		rate_values[39] = J[5]
		gprn += 1 # keep count on reaction number
		rate_values[40] = J[6]
		gprn += 1 # keep count on reaction number
		rate_values[41] = J[7]
		gprn += 1 # keep count on reaction number
		rate_values[42] = J[8]
		gprn += 1 # keep count on reaction number
		rate_values[43] = KMT04
		gprn += 1 # keep count on reaction number
		rate_values[44] = KMT10
		gprn += 1 # keep count on reaction number
		rate_values[45] = 6.6e-12*numpy.exp(-1240/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[46] = 1.85e-12*numpy.exp(-1690/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[47] = 3.8e-13*numpy.exp(780/TEMP)*(1-1/(1+498*numpy.exp(-1160/TEMP)))
		gprn += 1 # keep count on reaction number
		rate_values[48] = 3.8e-13*numpy.exp(780/TEMP)*(1/(1+498*numpy.exp(-1160/TEMP)))
		gprn += 1 # keep count on reaction number
		rate_values[49] = 2.3e-12*numpy.exp(360/TEMP)*0.001
		gprn += 1 # keep count on reaction number
		rate_values[50] = 2.3e-12*numpy.exp(360/TEMP)*0.999
		gprn += 1 # keep count on reaction number
		rate_values[51] = KMT13
		gprn += 1 # keep count on reaction number
		rate_values[52] = 1.2e-12
		gprn += 1 # keep count on reaction number
		rate_values[53] = 2*KCH3O2*RO2*7.18*numpy.exp(-885/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[54] = 2*KCH3O2*RO2*0.5*(1-7.18*numpy.exp(-885/TEMP))
		gprn += 1 # keep count on reaction number
		rate_values[55] = 2*KCH3O2*RO2*0.5*(1-7.18*numpy.exp(-885/TEMP))
		gprn += 1 # keep count on reaction number
		rate_values[56] = J[41]
		gprn += 1 # keep count on reaction number
		rate_values[57] = 5.3e-12*numpy.exp(190/TEMP)*0.6
		gprn += 1 # keep count on reaction number
		rate_values[58] = 5.3e-12*numpy.exp(190/TEMP)*0.4
		gprn += 1 # keep count on reaction number
		rate_values[59] = J[11]
		gprn += 1 # keep count on reaction number
		rate_values[60] = J[12]
		gprn += 1 # keep count on reaction number
		rate_values[61] = 5.5e-16
		gprn += 1 # keep count on reaction number
		rate_values[62] = 5.4e-12*numpy.exp(135/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[63] = J[51]
		gprn += 1 # keep count on reaction number
		rate_values[64] = 4.0e-13*numpy.exp(-845/TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[65] = 7.2e-14*numpy.exp(-1080/TEMP)*O2
		gprn += 1 # keep count on reaction number
		rate_values[66] = KMT14
		gprn += 1 # keep count on reaction number
		rate_values[67] = 2.85e-12*numpy.exp(-345/TEMP)
	except:
		erf = 1 # flag error
		err_mess = str('Error: estimating reaction rate for reaction number ' + str(gprn) + ' failed, please check chemical scheme (including whether definitions for generic rate coefficients have been included), and associated chemical scheme markers, which are stated in the model variables input file') # error message
	
	# aqueous-phase reactions
	
	# surface (e.g. wall) reactions
	
	return(rate_values, erf, err_mess)
