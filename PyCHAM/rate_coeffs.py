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
# created at 2024-02-08 16:05:08.711852

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
		KBPPN=(KPPN0*KPPNI)*FCPPN/(KPPN0+KPPNI) 
		gprn += 1 # keep count on reaction number
		KNO=KRO2NO*NO 
		gprn += 1 # keep count on reaction number
		KHO2=KRO2HO2*HO2*0.706 
		gprn += 1 # keep count on reaction number
		KRO2=1.26e-12*RO2 
		gprn += 1 # keep count on reaction number
		KNO3=KRO2NO3*NO3 
		gprn += 1 # keep count on reaction number
		KTR=KNO+KHO2+KRO2+KNO3 
		gprn += 1 # keep count on reaction number
		K16ISOM=(KTR*5.18e-04*numpy.exp(1308/TEMP))+(2.76e+07*numpy.exp(-6759/TEMP)) 

	except:
		erf = 1 # flag error
		err_mess = str('Error: generic reaction rates failed to be calculated inside rate_coeffs.py at number ' + str(gprn) + ', please check chemical scheme and associated chemical scheme markers, which are stated in the model variables input file') # error message
		return([], erf, err_mess)
	# estimate and append photolysis rates
	J = photolysisRates.PhotolysisCalculation(TEMP, Jlen, sumt, self)

	if (self.light_stat_now == 0):
		J = [0]*len(J)
	rate_values = numpy.zeros((74))
	
	# if reactions have been found in the chemical scheme
	# gas-phase reactions
	gprn = 0 # keep count on reaction number
	try:
		gprn += 1 # keep count on reaction number
		rate_values[0] = 5.6E-34*N2*(TEMP/300.)**(-2.6)*O2+6.0E-34*O2*(TEMP/300.)**(-2.6)*O2
		gprn += 1 # keep count on reaction number
		rate_values[1] = 8.0E-12*numpy.exp(-2060./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[2] = KMT01
		gprn += 1 # keep count on reaction number
		rate_values[3] = 5.5E-12*numpy.exp(188./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[4] = KMT02
		gprn += 1 # keep count on reaction number
		rate_values[5] = 3.2E-11*numpy.exp(67./TEMP)*O2+2.0E-11*numpy.exp(130./TEMP)*N2
		gprn += 1 # keep count on reaction number
		rate_values[6] = 1.4E-12*numpy.exp(-1310./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[7] = 1.4E-13*numpy.exp(-2470./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[8] = 3.3E-39*numpy.exp(530./TEMP)*O2
		gprn += 1 # keep count on reaction number
		rate_values[9] = 1.8E-11*numpy.exp(110./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[10] = 4.50E-14*numpy.exp(-1260./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[11] = KMT03
		gprn += 1 # keep count on reaction number
		rate_values[12] = 2.14E-10*H2O
		gprn += 1 # keep count on reaction number
		rate_values[13] = 1.70E-12*numpy.exp(-940./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[14] = 7.7E-12*numpy.exp(-2100./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[15] = KMT05
		gprn += 1 # keep count on reaction number
		rate_values[16] = 2.9E-12*numpy.exp(-160./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[17] = 2.03E-16*(TEMP/300.)**(4.57)*numpy.exp(693./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[18] = 4.8E-11*numpy.exp(250./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[19] = 2.20E-13*KMT06*numpy.exp(600./TEMP)+1.90E-33*M*KMT06*numpy.exp(980./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[20] = KMT07
		gprn += 1 # keep count on reaction number
		rate_values[21] = KMT08
		gprn += 1 # keep count on reaction number
		rate_values[22] = 2.0E-11
		gprn += 1 # keep count on reaction number
		rate_values[23] = 3.45E-12*numpy.exp(270./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[24] = KMT09
		gprn += 1 # keep count on reaction number
		rate_values[25] = 3.2E-13*numpy.exp(690./TEMP)*1.0
		gprn += 1 # keep count on reaction number
		rate_values[26] = 4.0E-12
		gprn += 1 # keep count on reaction number
		rate_values[27] = 2.5E-12*numpy.exp(260./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[28] = KMT11
		gprn += 1 # keep count on reaction number
		rate_values[29] = 4.0E-32*numpy.exp(-1000./TEMP)*M
		gprn += 1 # keep count on reaction number
		rate_values[30] = KMT12
		gprn += 1 # keep count on reaction number
		rate_values[31] = 1.3E-12*numpy.exp(-330./TEMP)*O2
		gprn += 1 # keep count on reaction number
		rate_values[32] = 6.00E-06
		gprn += 1 # keep count on reaction number
		rate_values[33] = 4.00E-04
		gprn += 1 # keep count on reaction number
		rate_values[34] = 1.20E-15*H2O
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
		rate_values[45] = 1.85E-12*numpy.exp(-1690./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[46] = 2.3E-12*numpy.exp(360./TEMP)*0.999
		gprn += 1 # keep count on reaction number
		rate_values[47] = 2.3E-12*numpy.exp(360./TEMP)*0.001
		gprn += 1 # keep count on reaction number
		rate_values[48] = 7.2E-14*numpy.exp(-1080./TEMP)*O2
		gprn += 1 # keep count on reaction number
		rate_values[49] = KMT13
		gprn += 1 # keep count on reaction number
		rate_values[50] = KMT14
		gprn += 1 # keep count on reaction number
		rate_values[51] = 1.2E-12
		gprn += 1 # keep count on reaction number
		rate_values[52] = 3.8E-13*numpy.exp(780./TEMP)*(1.-1./(1.+498.*numpy.exp(-1160./TEMP)))
		gprn += 1 # keep count on reaction number
		rate_values[53] = 2.*KCH3O2*RO2*7.18*numpy.exp(-885./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[54] = 2.*KCH3O2*RO2*0.5*(1.-7.18*numpy.exp(-885./TEMP))
		gprn += 1 # keep count on reaction number
		rate_values[55] = 2.*KCH3O2*RO2*0.5*(1.-7.18*numpy.exp(-885./TEMP))
		gprn += 1 # keep count on reaction number
		rate_values[56] = 4.0E-13*numpy.exp(-845./TEMP)
		gprn += 1 # keep count on reaction number
		rate_values[57] = J[51]
		gprn += 1 # keep count on reaction number
		rate_values[58] = 5.3E-12*numpy.exp(190./TEMP)*0.6
		gprn += 1 # keep count on reaction number
		rate_values[59] = 5.3E-12*numpy.exp(190./TEMP)*0.4
		gprn += 1 # keep count on reaction number
		rate_values[60] = J[41]
		gprn += 1 # keep count on reaction number
		rate_values[61] = 1.2E-11*numpy.exp(440./TEMP)*0.572
		gprn += 1 # keep count on reaction number
		rate_values[62] = 1.2E-11*numpy.exp(440./TEMP)*0.353
		gprn += 1 # keep count on reaction number
		rate_values[63] = 1.2E-11*numpy.exp(440./TEMP)*0.075
		gprn += 1 # keep count on reaction number
		rate_values[64] = KRO2NO*0.770
		gprn += 1 # keep count on reaction number
		rate_values[65] = KRO2NO*0.230
		gprn += 1 # keep count on reaction number
		rate_values[66] = KDEC
		gprn += 1 # keep count on reaction number
		rate_values[67] = KRO2NO3
		gprn += 1 # keep count on reaction number
		rate_values[68] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		rate_values[69] = 9.20E-14*RO2*0.7
		gprn += 1 # keep count on reaction number
		rate_values[70] = 9.20E-14*RO2*0.3
		gprn += 1 # keep count on reaction number
		rate_values[71] = 5.50E-12
		gprn += 1 # keep count on reaction number
		rate_values[72] = J[41]
		gprn += 1 # keep count on reaction number
		rate_values[73] = 1.83E-11
	except:
		erf = 1 # flag error
		err_mess = (str('Error: Could not calculate rate coefficient for equation number ' + str(gprn)))
	
	# aqueous-phase reactions
	
	# surface (e.g. wall) reactions
	
	return(rate_values, erf, err_mess)
