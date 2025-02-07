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
# created at 2025-02-06 21:27:00.030175

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
		KRO2NO=2.7e-12*numpy.exp(360/TEMP); 
		gprn += 1 # keep count on reaction number
		KRO2HO2=2.91e-13*numpy.exp(1300/TEMP); 
		gprn += 1 # keep count on reaction number
		KAPHO2=5.2e-13*numpy.exp(980/TEMP); 
		gprn += 1 # keep count on reaction number
		KAPNO=7.5e-12*numpy.exp(290/TEMP); 
		gprn += 1 # keep count on reaction number
		KRO2NO3=2.3e-12; 
		gprn += 1 # keep count on reaction number
		KNO3AL=1.44e-12*numpy.exp(-1862/TEMP); 
		gprn += 1 # keep count on reaction number
		KDEC=1.00e+06; 
		gprn += 1 # keep count on reaction number
		KROPRIM=2.50e-14*numpy.exp(-300/TEMP); 
		gprn += 1 # keep count on reaction number
		KROSEC=2.50e-14*numpy.exp(-300/TEMP); 
		gprn += 1 # keep count on reaction number
		KCH3O2=1.03e-13*numpy.exp(365/TEMP); 
		gprn += 1 # keep count on reaction number
		K298CH3O2=3.5e-13; 
		gprn += 1 # keep count on reaction number
		K14ISOM1=3.00e7*numpy.exp(-5300/TEMP); 
		gprn += 1 # keep count on reaction number
		KD0=1.10e-05*M*numpy.exp(-10100/TEMP); 
		gprn += 1 # keep count on reaction number
		KDI=1.90e17*numpy.exp(-14100/TEMP); 
		gprn += 1 # keep count on reaction number
		KRD=KD0/KDI; 
		gprn += 1 # keep count on reaction number
		FCD=0.30; 
		gprn += 1 # keep count on reaction number
		NCD=0.75-1.27*(numpy.log10(FCD)); 
		gprn += 1 # keep count on reaction number
		FD=10**(numpy.log10(FCD)/(1+(numpy.log10(KRD)/NCD)**2)); 
		gprn += 1 # keep count on reaction number
		KBPAN=(KD0*KDI)*FD/(KD0+KDI); 
		gprn += 1 # keep count on reaction number
		KC0=3.28e-28*M*(TEMP/300)**-6.87; 
		gprn += 1 # keep count on reaction number
		KCI=1.125e-11*(TEMP/300)**-1.105; 
		gprn += 1 # keep count on reaction number
		KRC=KC0/KCI; 
		gprn += 1 # keep count on reaction number
		FCC=0.30; 
		gprn += 1 # keep count on reaction number
		NC=0.75-1.27*(numpy.log10(FCC)); 
		gprn += 1 # keep count on reaction number
		FC=10**(numpy.log10(FCC)/(1+(numpy.log10(KRC)/NC)**2)); 
		gprn += 1 # keep count on reaction number
		KFPAN=(KC0*KCI)*FC/(KC0+KCI); 
		gprn += 1 # keep count on reaction number
		K10=1.0e-31*M*(TEMP/300)**-1.6; 
		gprn += 1 # keep count on reaction number
		K1I=5.0e-11*(TEMP/300)**-0.3; 
		gprn += 1 # keep count on reaction number
		KR1=K10/K1I; 
		gprn += 1 # keep count on reaction number
		FC1=0.85; 
		gprn += 1 # keep count on reaction number
		NC1=0.75-1.27*(numpy.log10(FC1)); 
		gprn += 1 # keep count on reaction number
		F1=10**(numpy.log10(FC1)/(1+(numpy.log10(KR1)/NC1)**2)); 
		gprn += 1 # keep count on reaction number
		KMT01=(K10*K1I)*F1/(K10+K1I); 
		gprn += 1 # keep count on reaction number
		K20=1.3e-31*M*(TEMP/300)**-1.5; 
		gprn += 1 # keep count on reaction number
		K2I=2.3e-11*(TEMP/300)**0.24; 
		gprn += 1 # keep count on reaction number
		KR2=K20/K2I; 
		gprn += 1 # keep count on reaction number
		FC2=0.6; 
		gprn += 1 # keep count on reaction number
		NC2=0.75-1.27*(numpy.log10(FC2)); 
		gprn += 1 # keep count on reaction number
		F2=10**(numpy.log10(FC2)/(1+(numpy.log10(KR2)/NC2)**2)); 
		gprn += 1 # keep count on reaction number
		KMT02=(K20*K2I)*F2/(K20+K2I); 
		gprn += 1 # keep count on reaction number
		K30=3.6e-30*M*(TEMP/300)**-4.1; 
		gprn += 1 # keep count on reaction number
		K3I=1.9e-12*(TEMP/300)**0.2; 
		gprn += 1 # keep count on reaction number
		KR3=K30/K3I; 
		gprn += 1 # keep count on reaction number
		FC3=0.35; 
		gprn += 1 # keep count on reaction number
		NC3=0.75-1.27*(numpy.log10(FC3)); 
		gprn += 1 # keep count on reaction number
		F3=10**(numpy.log10(FC3)/(1+(numpy.log10(KR3)/NC3)**2)); 
		gprn += 1 # keep count on reaction number
		KMT03=(K30*K3I)*F3/(K30+K3I); 
		gprn += 1 # keep count on reaction number
		K40=1.3e-3*M*(TEMP/300)**-3.5*numpy.exp(-11000/TEMP); 
		gprn += 1 # keep count on reaction number
		K4I=9.7e+14*(TEMP/300)**0.1*numpy.exp(-11080/TEMP); 
		gprn += 1 # keep count on reaction number
		KR4=K40/K4I; 
		gprn += 1 # keep count on reaction number
		FC4=0.35; 
		gprn += 1 # keep count on reaction number
		NC4=0.75-1.27*(numpy.log10(FC4)); 
		gprn += 1 # keep count on reaction number
		F4=10**(numpy.log10(FC4)/(1+(numpy.log10(KR4)/NC4)**2)); 
		gprn += 1 # keep count on reaction number
		KMT04=(K40*K4I)*F4/(K40+K4I); 
		gprn += 1 # keep count on reaction number
		KMT05=1.44e-13*(1+(M/4.2e+19)); 
		gprn += 1 # keep count on reaction number
		KMT06=1+(1.40e-21*numpy.exp(2200/TEMP)*H2O); 
		gprn += 1 # keep count on reaction number
		K70=7.4e-31*M*(TEMP/300)**-2.4; 
		gprn += 1 # keep count on reaction number
		K7I=3.3e-11*(TEMP/300)**-0.3; 
		gprn += 1 # keep count on reaction number
		KR7=K70/K7I; 
		gprn += 1 # keep count on reaction number
		FC7=0.81; 
		gprn += 1 # keep count on reaction number
		NC7=0.75-1.27*(numpy.log10(FC7)); 
		gprn += 1 # keep count on reaction number
		F7=10**(numpy.log10(FC7)/(1+(numpy.log10(KR7)/NC7)**2)); 
		gprn += 1 # keep count on reaction number
		KMT07=(K70*K7I)*F7/(K70+K7I); 
		gprn += 1 # keep count on reaction number
		K80=3.2e-30*M*(TEMP/300)**-4.5; 
		gprn += 1 # keep count on reaction number
		K8I=3.0e-11; 
		gprn += 1 # keep count on reaction number
		KR8=K80/K8I; 
		gprn += 1 # keep count on reaction number
		FC8=0.41; 
		gprn += 1 # keep count on reaction number
		NC8=0.75-1.27*(numpy.log10(FC8)); 
		gprn += 1 # keep count on reaction number
		F8=10**(numpy.log10(FC8)/(1+(numpy.log10(KR8)/NC8)**2)); 
		gprn += 1 # keep count on reaction number
		KMT08=(K80*K8I)*F8/(K80+K8I); 
		gprn += 1 # keep count on reaction number
		K90=1.4e-31*M*(TEMP/300)**-3.1; 
		gprn += 1 # keep count on reaction number
		K9I=4.0e-12; 
		gprn += 1 # keep count on reaction number
		KR9=K90/K9I; 
		gprn += 1 # keep count on reaction number
		FC9=0.4; 
		gprn += 1 # keep count on reaction number
		NC9=0.75-1.27*(numpy.log10(FC9)); 
		gprn += 1 # keep count on reaction number
		F9=10**(numpy.log10(FC9)/(1+(numpy.log10(KR9)/NC9)**2)); 
		gprn += 1 # keep count on reaction number
		KMT09=(K90*K9I)*F9/(K90+K9I); 
		gprn += 1 # keep count on reaction number
		K100=4.10e-05*M*numpy.exp(-10650/TEMP); 
		gprn += 1 # keep count on reaction number
		K10I=6.0e+15*numpy.exp(-11170/TEMP); 
		gprn += 1 # keep count on reaction number
		KR10=K100/K10I; 
		gprn += 1 # keep count on reaction number
		FC10=0.4; 
		gprn += 1 # keep count on reaction number
		NC10=0.75-1.27*(numpy.log10(FC10)); 
		gprn += 1 # keep count on reaction number
		F10=10**(numpy.log10(FC10)/(1+(numpy.log10(KR10)/NC10)**2)); 
		gprn += 1 # keep count on reaction number
		KMT10=(K100*K10I)*F10/(K100+K10I); 
		gprn += 1 # keep count on reaction number
		K1=2.40e-14*numpy.exp(460/TEMP); 
		gprn += 1 # keep count on reaction number
		K3=6.50e-34*numpy.exp(1335/TEMP); 
		gprn += 1 # keep count on reaction number
		K4=2.70e-17*numpy.exp(2199/TEMP); 
		gprn += 1 # keep count on reaction number
		K2=(K3*M)/(1+(K3*M/K4)); 
		gprn += 1 # keep count on reaction number
		KMT11=K1+K2; 
		gprn += 1 # keep count on reaction number
		K120=2.5e-31*M*(TEMP/300)**-2.6; 
		gprn += 1 # keep count on reaction number
		K12I=2.0e-12; 
		gprn += 1 # keep count on reaction number
		KR12=K120/K12I; 
		gprn += 1 # keep count on reaction number
		FC12=0.53; 
		gprn += 1 # keep count on reaction number
		NC12=0.75-1.27*(numpy.log10(FC12)); 
		gprn += 1 # keep count on reaction number
		F12=10**(numpy.log10(FC12)/(1.0+(numpy.log10(KR12)/NC12)**2)); 
		gprn += 1 # keep count on reaction number
		KMT12=(K120*K12I*F12)/(K120+K12I); 
		gprn += 1 # keep count on reaction number
		K130=2.5e-30*M*(TEMP/300)**-5.5; 
		gprn += 1 # keep count on reaction number
		K13I=1.8e-11; 
		gprn += 1 # keep count on reaction number
		KR13=K130/K13I; 
		gprn += 1 # keep count on reaction number
		FC13=0.36; 
		gprn += 1 # keep count on reaction number
		NC13=0.75-1.27*(numpy.log10(FC13)); 
		gprn += 1 # keep count on reaction number
		F13=10**(numpy.log10(FC13)/(1+(numpy.log10(KR13)/NC13)**2)); 
		gprn += 1 # keep count on reaction number
		KMT13=(K130*K13I)*F13/(K130+K13I); 
		gprn += 1 # keep count on reaction number
		K140=9.0e-5*numpy.exp(-9690/TEMP)*M; 
		gprn += 1 # keep count on reaction number
		K14I=1.1e+16*numpy.exp(-10560/TEMP); 
		gprn += 1 # keep count on reaction number
		KR14=K140/K14I; 
		gprn += 1 # keep count on reaction number
		FC14=0.36; 
		gprn += 1 # keep count on reaction number
		NC14=0.75-1.27*(numpy.log10(FC14)); 
		gprn += 1 # keep count on reaction number
		F14=10**(numpy.log10(FC14)/(1+(numpy.log10(KR14)/NC14)**2)); 
		gprn += 1 # keep count on reaction number
		KMT14=(K140*K14I)*F14/(K140+K14I); 
		gprn += 1 # keep count on reaction number
		K150=8.6e-29*M*(TEMP/300)**-3.1; 
		gprn += 1 # keep count on reaction number
		K15I=9.0e-12*(TEMP/300)**-0.85; 
		gprn += 1 # keep count on reaction number
		KR15=K150/K15I; 
		gprn += 1 # keep count on reaction number
		FC15=0.48; 
		gprn += 1 # keep count on reaction number
		NC15=0.75-1.27*(numpy.log10(FC15)); 
		gprn += 1 # keep count on reaction number
		F15=10**(numpy.log10(FC15)/(1+(numpy.log10(KR15)/NC15)**2)); 
		gprn += 1 # keep count on reaction number
		KMT15=(K150*K15I)*F15/(K150+K15I); 
		gprn += 1 # keep count on reaction number
		K160=8e-27*M*(TEMP/300)**-3.5; 
		gprn += 1 # keep count on reaction number
		K16I=3.0e-11*(TEMP/300)**-1; 
		gprn += 1 # keep count on reaction number
		KR16=K160/K16I; 
		gprn += 1 # keep count on reaction number
		FC16=0.5; 
		gprn += 1 # keep count on reaction number
		NC16=0.75-1.27*(numpy.log10(FC16)); 
		gprn += 1 # keep count on reaction number
		F16=10**(numpy.log10(FC16)/(1+(numpy.log10(KR16)/NC16)**2)); 
		gprn += 1 # keep count on reaction number
		KMT16=(K160*K16I)*F16/(K160+K16I); 
		gprn += 1 # keep count on reaction number
		K170=5.0e-30*M*(TEMP/300)**-1.5; 
		gprn += 1 # keep count on reaction number
		K17I=1.0e-12; 
		gprn += 1 # keep count on reaction number
		KR17=K170/K17I; 
		gprn += 1 # keep count on reaction number
		FC17=0.17*numpy.exp(-51/TEMP)+numpy.exp(-TEMP/204); 
		gprn += 1 # keep count on reaction number
		NC17=0.75-1.27*(numpy.log10(FC17)); 
		gprn += 1 # keep count on reaction number
		F17=10**(numpy.log10(FC17)/(1.0+(numpy.log10(KR17)/NC17)**2)); 
		gprn += 1 # keep count on reaction number
		KMT17=(K170*K17I*F17)/(K170+K17I); 
		gprn += 1 # keep count on reaction number
		KMT18=9.5e-39*O2*numpy.exp(5270/TEMP)/(1+7.5e-29*O2*numpy.exp(5610/TEMP)); 
		gprn += 1 # keep count on reaction number
		KPPN0=1.7e-03*numpy.exp(-11280/TEMP)*M; 
		gprn += 1 # keep count on reaction number
		KPPNI=8.3e+16*numpy.exp(-13940/TEMP); 
		gprn += 1 # keep count on reaction number
		KRPPN=KPPN0/KPPNI; 
		gprn += 1 # keep count on reaction number
		FCPPN=0.36; 
		gprn += 1 # keep count on reaction number
		NCPPN=0.75-1.27*(numpy.log10(FCPPN)); 
		gprn += 1 # keep count on reaction number
		FPPN=10**(numpy.log10(FCPPN)/(1+(numpy.log10(KRPPN)/NCPPN)**2)); 
		gprn += 1 # keep count on reaction number
		KBPPN=(KPPN0*KPPNI)*FPPN/(KPPN0+KPPNI); 
		gprn += 1 # keep count on reaction number
		KRO2=1.26e-12*RO2; 
		gprn += 1 # keep count on reaction number
		KNO3=KRO2NO3*NO3; 
		gprn += 1 # keep count on reaction number
		#ClosedshellformationfromRO:RO->R=O(#H->#H-1 
		gprn += 1 # keep count on reaction number
		#RO2reactingwithsumofRO2sformingR=ObyOHremoval: 
		gprn += 1 # keep count on reaction number
		#RO2"+sum(RO2)->R=O";CxHyOz->CxH(y-1)O(z-1) 
		gprn += 1 # keep count on reaction number
		#ClosedshellformationfromRO:RO->R=O(#H->#H-1 
		gprn += 1 # keep count on reaction number
		#RO2undergoingH-shiftandformingcarbonylradical(RC=O*) 
		gprn += 1 # keep count on reaction number
		#(RC=O*)undergoingCOscissionandaddingO2tofromC5RO2again: 
		gprn += 1 # keep count on reaction number
		#RO2radicalsabstractingHfromalphahydroxylcarbon->RC=O 

	except:
		erf = 1 # flag error
		err_mess = str('Error: generic reaction rates failed to be calculated inside rate_coeffs.py at number ' + str(gprn) + ', please check chemical scheme and associated chemical scheme markers, which are stated in the model variables input file') # error message
		return([], erf, err_mess)
	# estimate and append photolysis rates
	J = photolysisRates.PhotolysisCalculation(TEMP, Jlen, sumt, self)

	if (self.light_stat_now == 0):
		J = [0]*len(J)
	rate_values = numpy.zeros((5797))
	
	# if reactions have been found in the chemical scheme
	# gas-phase reactions
	gprn = 0 # keep count on reaction number
	try:
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.6e-34*N2*(TEMP/300)**-2.6*O2+6.0e-34*O2*(TEMP/300)**-2.6*O2' 
		rate_values[0] = 5.6e-34*N2*(TEMP/300)**-2.6*O2+6.0e-34*O2*(TEMP/300)**-2.6*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.0e-12*numpy.exp(-2060/TEMP)' 
		rate_values[1] = 8.0e-12*numpy.exp(-2060/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT01' 
		rate_values[2] = KMT01
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5e-12*numpy.exp(188/TEMP)' 
		rate_values[3] = 5.5e-12*numpy.exp(188/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT02' 
		rate_values[4] = KMT02
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.2e-11*numpy.exp(67/TEMP)*O2+2.0e-11*numpy.exp(130/TEMP)*N2' 
		rate_values[5] = 3.2e-11*numpy.exp(67/TEMP)*O2+2.0e-11*numpy.exp(130/TEMP)*N2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4e-12*numpy.exp(-1310/TEMP)' 
		rate_values[6] = 1.4e-12*numpy.exp(-1310/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4e-13*numpy.exp(-2470/TEMP)' 
		rate_values[7] = 1.4e-13*numpy.exp(-2470/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3e-39*numpy.exp(530/TEMP)*O2' 
		rate_values[8] = 3.3e-39*numpy.exp(530/TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8e-11*numpy.exp(110/TEMP)' 
		rate_values[9] = 1.8e-11*numpy.exp(110/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.50e-14*numpy.exp(-1260/TEMP)' 
		rate_values[10] = 4.50e-14*numpy.exp(-1260/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT03' 
		rate_values[11] = KMT03
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.14e-10*H2O' 
		rate_values[12] = 2.14e-10*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.70e-12*numpy.exp(-940/TEMP)' 
		rate_values[13] = 1.70e-12*numpy.exp(-940/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.7e-12*numpy.exp(-2100/TEMP)' 
		rate_values[14] = 7.7e-12*numpy.exp(-2100/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT05' 
		rate_values[15] = KMT05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.9e-12*numpy.exp(-160/TEMP)' 
		rate_values[16] = 2.9e-12*numpy.exp(-160/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.03e-16*(TEMP/300)**4.57*numpy.exp(693/TEMP)' 
		rate_values[17] = 2.03e-16*(TEMP/300)**4.57*numpy.exp(693/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.8e-11*numpy.exp(250/TEMP)' 
		rate_values[18] = 4.8e-11*numpy.exp(250/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.20e-13*KMT06*numpy.exp(600/TEMP)+1.90e-33*M*KMT06*numpy.exp(980/TEMP)' 
		rate_values[19] = 2.20e-13*KMT06*numpy.exp(600/TEMP)+1.90e-33*M*KMT06*numpy.exp(980/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT07' 
		rate_values[20] = KMT07
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT08' 
		rate_values[21] = KMT08
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0e-11' 
		rate_values[22] = 2.0e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.45e-12*numpy.exp(270/TEMP)' 
		rate_values[23] = 3.45e-12*numpy.exp(270/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT09' 
		rate_values[24] = KMT09
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.2e-13*numpy.exp(690/TEMP)*1.0' 
		rate_values[25] = 3.2e-13*numpy.exp(690/TEMP)*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0e-12' 
		rate_values[26] = 4.0e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.5e-12*numpy.exp(260/TEMP)' 
		rate_values[27] = 2.5e-12*numpy.exp(260/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT11' 
		rate_values[28] = KMT11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0e-32*numpy.exp(-1000/TEMP)*M' 
		rate_values[29] = 4.0e-32*numpy.exp(-1000/TEMP)*M
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT12' 
		rate_values[30] = KMT12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.3e-12*numpy.exp(-330/TEMP)*O2' 
		rate_values[31] = 1.3e-12*numpy.exp(-330/TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.00e-06' 
		rate_values[32] = 6.00e-06
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.00e-04' 
		rate_values[33] = 4.00e-04
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.20e-15*H2O' 
		rate_values[34] = 1.20e-15*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[1]' 
		rate_values[35] = J[1]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[2]' 
		rate_values[36] = J[2]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[3]' 
		rate_values[37] = J[3]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]' 
		rate_values[38] = J[4]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[5]' 
		rate_values[39] = J[5]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[6]' 
		rate_values[40] = J[6]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[7]' 
		rate_values[41] = J[7]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[8]' 
		rate_values[42] = J[8]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT04' 
		rate_values[43] = KMT04
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT10' 
		rate_values[44] = KMT10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-12*numpy.exp(490/TEMP)*0.65' 
		rate_values[45] = 1.2e-12*numpy.exp(490/TEMP)*0.65
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-12*numpy.exp(490/TEMP)*0.35' 
		rate_values[46] = 1.2e-12*numpy.exp(490/TEMP)*0.35
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.05e-16*numpy.exp(-640/TEMP)*0.6' 
		rate_values[47] = 8.05e-16*numpy.exp(-640/TEMP)*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.05e-16*numpy.exp(-640/TEMP)*0.4' 
		rate_values[48] = 8.05e-16*numpy.exp(-640/TEMP)*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-11*numpy.exp(440/TEMP)*0.558' 
		rate_values[49] = 1.2e-11*numpy.exp(440/TEMP)*0.558
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-11*numpy.exp(440/TEMP)*0.344' 
		rate_values[50] = 1.2e-11*numpy.exp(440/TEMP)*0.344
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-11*numpy.exp(440/TEMP)*0.073' 
		rate_values[51] = 1.2e-11*numpy.exp(440/TEMP)*0.073
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.6e-12*numpy.exp(-1240/TEMP)' 
		rate_values[52] = 6.6e-12*numpy.exp(-1240/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.85e-12*numpy.exp(-1690/TEMP)' 
		rate_values[53] = 1.85e-12*numpy.exp(-1690/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3e-12*numpy.exp(-190/TEMP)*0.352' 
		rate_values[54] = 2.3e-12*numpy.exp(-190/TEMP)*0.352
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3e-12*numpy.exp(-190/TEMP)*0.118' 
		rate_values[55] = 2.3e-12*numpy.exp(-190/TEMP)*0.118
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3e-12*numpy.exp(-190/TEMP)*0.53' 
		rate_values[56] = 2.3e-12*numpy.exp(-190/TEMP)*0.53
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[57] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[58] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[59] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.1*RO2' 
		rate_values[60] = 6.70e-15*0.1*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.9*RO2' 
		rate_values[61] = 6.70e-15*0.9*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[62] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[63] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[64] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*0.1*RO2' 
		rate_values[65] = 2.50e-13*0.1*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*0.8*RO2' 
		rate_values[66] = 2.50e-13*0.8*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*0.1*RO2' 
		rate_values[67] = 2.50e-13*0.1*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.50' 
		rate_values[68] = KDEC*0.50
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.50' 
		rate_values[69] = KDEC*0.50
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[70] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.230' 
		rate_values[71] = KRO2NO*0.230
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.770' 
		rate_values[72] = KRO2NO*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[73] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.7' 
		rate_values[74] = 9.20e-14*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.3' 
		rate_values[75] = 9.20e-14*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[76] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.230' 
		rate_values[77] = KRO2NO*0.230
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.770' 
		rate_values[78] = KRO2NO*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[79] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[80] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[81] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[82] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[83] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.125' 
		rate_values[84] = KRO2NO*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.875' 
		rate_values[85] = KRO2NO*0.875
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[86] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2*0.7' 
		rate_values[87] = 6.70e-15*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2*0.3' 
		rate_values[88] = 6.70e-15*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8e-13*numpy.exp(780/TEMP)*(1-1/(1+498*numpy.exp(-1160/TEMP)))' 
		rate_values[89] = 3.8e-13*numpy.exp(780/TEMP)*(1-1/(1+498*numpy.exp(-1160/TEMP)))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8e-13*numpy.exp(780/TEMP)*(1/(1+498*numpy.exp(-1160/TEMP)))' 
		rate_values[90] = 3.8e-13*numpy.exp(780/TEMP)*(1/(1+498*numpy.exp(-1160/TEMP)))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3e-12*numpy.exp(360/TEMP)*0.001' 
		rate_values[91] = 2.3e-12*numpy.exp(360/TEMP)*0.001
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3e-12*numpy.exp(360/TEMP)*0.999' 
		rate_values[92] = 2.3e-12*numpy.exp(360/TEMP)*0.999
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT13' 
		rate_values[93] = KMT13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-12' 
		rate_values[94] = 1.2e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KCH3O2*RO2*7.18*numpy.exp(-885/TEMP)' 
		rate_values[95] = 2*KCH3O2*RO2*7.18*numpy.exp(-885/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KCH3O2*RO2*0.5*(1-7.18*numpy.exp(-885/TEMP))' 
		rate_values[96] = 2*KCH3O2*RO2*0.5*(1-7.18*numpy.exp(-885/TEMP))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KCH3O2*RO2*0.5*(1-7.18*numpy.exp(-885/TEMP))' 
		rate_values[97] = 2*KCH3O2*RO2*0.5*(1-7.18*numpy.exp(-885/TEMP))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[98] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.082' 
		rate_values[99] = KRO2NO*0.082
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.918' 
		rate_values[100] = KRO2NO*0.918
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[101] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[102] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[103] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[104] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KNO3AL*2.75' 
		rate_values[105] = 2*KNO3AL*2.75
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-18' 
		rate_values[106] = 2.00e-18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.08e-11*0.31' 
		rate_values[107] = 6.08e-11*0.31
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.08e-11*0.69' 
		rate_values[108] = 6.08e-11*0.69
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.1*0.5' 
		rate_values[109] = J[4]*0.1*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.1*0.5' 
		rate_values[110] = J[4]*0.1*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8e-12*0.742' 
		rate_values[111] = 3.8e-12*0.742
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8e-12*0.258' 
		rate_values[112] = 3.8e-12*0.258
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.7e-13*numpy.exp(1220/TEMP)*0.06' 
		rate_values[113] = 4.7e-13*numpy.exp(1220/TEMP)*0.06
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.7e-13*numpy.exp(1220/TEMP)*0.8' 
		rate_values[114] = 4.7e-13*numpy.exp(1220/TEMP)*0.8
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.7e-13*numpy.exp(1220/TEMP)*0.14' 
		rate_values[115] = 4.7e-13*numpy.exp(1220/TEMP)*0.14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.87e-12' 
		rate_values[116] = 6.87e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[117] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[118] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.64e-12' 
		rate_values[119] = 3.64e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[120] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.23e-11' 
		rate_values[121] = 1.23e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[122] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROSEC*O2' 
		rate_values[123] = KROSEC*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.00e+05' 
		rate_values[124] = 4.00e+05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.50e-12' 
		rate_values[125] = 5.50e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.55e-12' 
		rate_values[126] = 5.55e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]' 
		rate_values[127] = J[55]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[128] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[129] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[130] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*0.7*RO2' 
		rate_values[131] = 9.20e-14*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*0.3*RO2' 
		rate_values[132] = 9.20e-14*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[133] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[134] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[135] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2*0.05' 
		rate_values[136] = 2.00e-12*RO2*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2*0.90' 
		rate_values[137] = 2.00e-12*RO2*0.90
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2*0.05' 
		rate_values[138] = 2.00e-12*RO2*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.20e-15' 
		rate_values[139] = 1.20e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-14' 
		rate_values[140] = 1.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-15' 
		rate_values[141] = 1.00e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00e-14' 
		rate_values[142] = 7.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.40e-17*H2O' 
		rate_values[143] = 1.40e-17*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-18*H2O' 
		rate_values[144] = 2.00e-18*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[145] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.157' 
		rate_values[146] = KRO2NO*0.157
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.843' 
		rate_values[147] = KRO2NO*0.843
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[148] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.6*RO2' 
		rate_values[149] = 1.30e-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.2*RO2' 
		rate_values[150] = 1.30e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.2*RO2' 
		rate_values[151] = 1.30e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.83e-11' 
		rate_values[152] = 1.83e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[153] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[154] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.49e-11' 
		rate_values[155] = 1.49e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.28e-11' 
		rate_values[156] = 3.28e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[157] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[158] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.18e-12' 
		rate_values[159] = 8.18e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.03e-10' 
		rate_values[160] = 1.03e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[161] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.87e-11' 
		rate_values[162] = 9.87e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]' 
		rate_values[163] = J[55]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[164] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.91e-11' 
		rate_values[165] = 9.91e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[166] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.3e-12*numpy.exp(190/TEMP)*0.6' 
		rate_values[167] = 5.3e-12*numpy.exp(190/TEMP)*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.3e-12*numpy.exp(190/TEMP)*0.4' 
		rate_values[168] = 5.3e-12*numpy.exp(190/TEMP)*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[169] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[12]' 
		rate_values[170] = J[12]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5e-16' 
		rate_values[171] = 5.5e-16
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.4e-12*numpy.exp(135/TEMP)' 
		rate_values[172] = 5.4e-12*numpy.exp(135/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[51]' 
		rate_values[173] = J[51]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0e-13*numpy.exp(-845/TEMP)' 
		rate_values[174] = 4.0e-13*numpy.exp(-845/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.2e-14*numpy.exp(-1080/TEMP)*O2' 
		rate_values[175] = 7.2e-14*numpy.exp(-1080/TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT14' 
		rate_values[176] = KMT14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.85e-12*numpy.exp(-345/TEMP)' 
		rate_values[177] = 2.85e-12*numpy.exp(-345/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.77e-11' 
		rate_values[178] = 9.77e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[179] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.30e-11' 
		rate_values[180] = 7.30e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[181] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.50' 
		rate_values[182] = KDEC*0.50
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.50' 
		rate_values[183] = KDEC*0.50
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.21e-10' 
		rate_values[184] = 1.21e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.16e-11' 
		rate_values[185] = 8.16e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[186] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[187] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[188] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[189] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[190] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[191] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[192] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[193] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[194] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KNO3AL*4.0' 
		rate_values[195] = 2*KNO3AL*4.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.32e-11' 
		rate_values[196] = 4.32e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]*2' 
		rate_values[197] = J[17]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.18' 
		rate_values[198] = KDEC*0.18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.125' 
		rate_values[199] = KDEC*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.125' 
		rate_values[200] = KDEC*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.57' 
		rate_values[201] = KDEC*0.57
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.77' 
		rate_values[202] = KRO2HO2*0.77
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.105' 
		rate_values[203] = KRO2NO*0.105
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.895' 
		rate_values[204] = KRO2NO*0.895
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[205] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[206] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.6*RO2' 
		rate_values[207] = 8.80e-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[208] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[209] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[210] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[211] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[212] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[213] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[214] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KNO3AL*2.0' 
		rate_values[215] = 2*KNO3AL*2.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-18' 
		rate_values[216] = 2.00e-18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.20e-11*0.83' 
		rate_values[217] = 5.20e-11*0.83
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.20e-11*0.17' 
		rate_values[218] = 5.20e-11*0.17
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.14*0.4' 
		rate_values[219] = J[4]*0.14*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.14*0.6' 
		rate_values[220] = J[4]*0.14*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[221] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[222] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[223] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[224] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[225] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.7' 
		rate_values[226] = 8.00e-13*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.3' 
		rate_values[227] = 8.00e-13*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.9e-11' 
		rate_values[228] = 9.9e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.2e-18' 
		rate_values[229] = 9.2e-18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0e-10' 
		rate_values[230] = 1.0e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[231] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[232] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[233] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.7' 
		rate_values[234] = 8.00e-13*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.3' 
		rate_values[235] = 8.00e-13*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0e-14' 
		rate_values[236] = 2.0e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2e-12*numpy.exp(600/TEMP)*0.772' 
		rate_values[237] = 5.2e-12*numpy.exp(600/TEMP)*0.772
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2e-12*numpy.exp(600/TEMP)*0.228' 
		rate_values[238] = 5.2e-12*numpy.exp(600/TEMP)*0.228
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[239] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[240] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[241] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[242] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2' 
		rate_values[243] = 6.70e-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[244] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[245] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[246] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[247] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[248] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[249] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.7*RO2' 
		rate_values[250] = 1.00e-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.3*RO2' 
		rate_values[251] = 1.00e-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.01e-11' 
		rate_values[252] = 3.01e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]' 
		rate_values[253] = J[41]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[254] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.66e-11' 
		rate_values[255] = 2.66e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[256] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.47e-11' 
		rate_values[257] = 5.47e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]' 
		rate_values[258] = J[41]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[259] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.80' 
		rate_values[260] = KDEC*0.80
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.20' 
		rate_values[261] = KDEC*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.47e-11' 
		rate_values[262] = 5.47e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]+J[15]' 
		rate_values[263] = J[34]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.45e-11' 
		rate_values[264] = 4.45e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[265] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[266] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.65e-12' 
		rate_values[267] = 6.65e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[268] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[269] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-11' 
		rate_values[270] = 1.30e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[22]' 
		rate_values[271] = J[41]+J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.88e-12' 
		rate_values[272] = 2.88e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[53]+J[22]' 
		rate_values[273] = J[53]+J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.20e+10*numpy.exp(-3523/TEMP)' 
		rate_values[274] = 4.20e+10*numpy.exp(-3523/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.67e-12' 
		rate_values[275] = 7.67e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[276] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*8.5' 
		rate_values[277] = KNO3AL*8.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.64e-11' 
		rate_values[278] = 2.64e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[279] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.8e-12*numpy.exp(-1320/TEMP)+1.7e-14*numpy.exp(423/TEMP)' 
		rate_values[280] = 8.8e-12*numpy.exp(-1320/TEMP)+1.7e-14*numpy.exp(423/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[21]' 
		rate_values[281] = J[21]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.19e-10' 
		rate_values[282] = 1.19e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[283] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.278' 
		rate_values[284] = KRO2NO*0.278
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.722' 
		rate_values[285] = KRO2NO*0.722
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[286] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2*0.6' 
		rate_values[287] = 2.50e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2*0.2' 
		rate_values[288] = 2.50e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2*0.2' 
		rate_values[289] = 2.50e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[31]' 
		rate_values[290] = J[31]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[33]' 
		rate_values[291] = J[33]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[32]' 
		rate_values[292] = J[32]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL' 
		rate_values[293] = KNO3AL
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.1e-12*numpy.exp(340/TEMP)' 
		rate_values[294] = 3.1e-12*numpy.exp(340/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.00e-13' 
		rate_values[295] = 3.00e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.20e-19' 
		rate_values[296] = 2.20e-19
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.45e-11' 
		rate_values[297] = 4.45e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[298] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[299] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[300] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[301] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[302] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[303] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.30*RO2' 
		rate_values[304] = 1.00e-11*0.30*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.70*RO2' 
		rate_values[305] = 1.00e-11*0.70*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.06e-11' 
		rate_values[306] = 4.06e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]+J[19]' 
		rate_values[307] = J[18]+J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.37e-11' 
		rate_values[308] = 4.37e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[18]+J[19]' 
		rate_values[309] = J[41]+J[18]+J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.05e-11' 
		rate_values[310] = 4.05e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[311] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[312] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[313] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[314] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[315] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[316] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[317] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[318] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[319] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.520' 
		rate_values[320] = KRO2HO2*0.520
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[321] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[322] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[323] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[324] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[325] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-15' 
		rate_values[326] = 1.2e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0e-14' 
		rate_values[327] = 1.0e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0e-15' 
		rate_values[328] = 1.0e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.0e-14' 
		rate_values[329] = 7.0e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0e-18*H2O' 
		rate_values[330] = 6.0e-18*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0e-17*H2O' 
		rate_values[331] = 1.0e-17*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.31e-10' 
		rate_values[332] = 1.31e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]*2' 
		rate_values[333] = J[41]+J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.38e-11' 
		rate_values[334] = 4.38e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[335] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[336] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[337] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-11' 
		rate_values[338] = 9.20e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2+J[22]' 
		rate_values[339] = J[15]*2+J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.23e-11' 
		rate_values[340] = 8.23e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[341] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.52e-11' 
		rate_values[342] = 7.52e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[343] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]+J[19]' 
		rate_values[344] = J[18]+J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[345] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.90e-11' 
		rate_values[346] = 4.90e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]+J[18]+J[19]' 
		rate_values[347] = J[34]+J[18]+J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.75e-11' 
		rate_values[348] = 7.75e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]+J[19]' 
		rate_values[349] = J[18]+J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[350] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[351] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[352] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.20*RO2' 
		rate_values[353] = 8.80e-13*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.20*RO2' 
		rate_values[354] = 8.80e-13*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.60*RO2' 
		rate_values[355] = 8.80e-13*0.60*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.00e-14' 
		rate_values[356] = 9.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.00e-13' 
		rate_values[357] = 9.00e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[358] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[359] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[360] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[361] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.07e-10' 
		rate_values[362] = 1.07e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[363] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[364] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[365] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.04e-10' 
		rate_values[366] = 1.04e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[367] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[368] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[369] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[370] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.16e-10' 
		rate_values[371] = 1.16e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[372] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.71' 
		rate_values[373] = KDEC*0.71
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.29' 
		rate_values[374] = KDEC*0.29
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.13e-10' 
		rate_values[375] = 1.13e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[376] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.050' 
		rate_values[377] = KRO2NO*0.050
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.950' 
		rate_values[378] = KRO2NO*0.950
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[379] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.7*RO2' 
		rate_values[380] = 6.70e-15*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.3*RO2' 
		rate_values[381] = 6.70e-15*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.94e-12' 
		rate_values[382] = 5.94e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[383] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[384] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.73e-12' 
		rate_values[385] = 9.73e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[22]' 
		rate_values[386] = J[41]+J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.66e-12' 
		rate_values[387] = 3.66e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[388] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[389] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.125' 
		rate_values[390] = KRO2NO*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.875' 
		rate_values[391] = KRO2NO*0.875
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[392] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.7*RO2' 
		rate_values[393] = 6.70e-15*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.3*RO2' 
		rate_values[394] = 6.70e-15*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[395] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[396] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[397] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[398] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[399] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[400] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[401] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[402] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[403] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[404] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[405] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[406] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[407] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[408] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[409] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[410] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[411] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[412] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[413] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*RO2' 
		rate_values[414] = 1.30e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[415] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[416] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[417] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.7*RO2' 
		rate_values[418] = 6.70e-15*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.3*RO2' 
		rate_values[419] = 6.70e-15*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[420] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[421] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[422] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[423] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[424] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[425] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[426] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[427] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[428] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2' 
		rate_values[429] = 6.70e-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.36e-13*numpy.exp(1250/TEMP)*0.15' 
		rate_values[430] = 1.36e-13*numpy.exp(1250/TEMP)*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.36e-13*numpy.exp(1250/TEMP)*0.85' 
		rate_values[431] = 1.36e-13*numpy.exp(1250/TEMP)*0.85
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[432] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[433] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*(K298CH3O2*8.0e-12)**0.5*RO2*0.2' 
		rate_values[434] = 2*(K298CH3O2*8.0e-12)**0.5*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*(K298CH3O2*8.0e-12)**0.5*RO2*0.6' 
		rate_values[435] = 2*(K298CH3O2*8.0e-12)**0.5*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*(K298CH3O2*8.0e-12)**0.5*RO2*0.2' 
		rate_values[436] = 2*(K298CH3O2*8.0e-12)**0.5*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[437] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[438] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[439] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.5e-12*numpy.exp(290/TEMP)' 
		rate_values[440] = 7.5e-12*numpy.exp(290/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[441] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0e-12' 
		rate_values[442] = 4.0e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.3*RO2' 
		rate_values[443] = 1.00e-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.7*RO2' 
		rate_values[444] = 1.00e-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[445] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.042' 
		rate_values[446] = KRO2NO*0.042
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.958' 
		rate_values[447] = KRO2NO*0.958
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[448] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.7' 
		rate_values[449] = 9.20e-14*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.3' 
		rate_values[450] = 9.20e-14*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.27e-10' 
		rate_values[451] = 1.27e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[452] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.60e-11' 
		rate_values[453] = 9.60e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[454] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROSEC*O2' 
		rate_values[455] = KROSEC*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.09e-10' 
		rate_values[456] = 1.09e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00e11*numpy.exp(-3160/TEMP)+5.00e-12*O2' 
		rate_values[457] = 7.00e11*numpy.exp(-3160/TEMP)+5.00e-12*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00e-12*O2*3.2*(1-numpy.exp(-550/TEMP))' 
		rate_values[458] = 5.00e-12*O2*3.2*(1-numpy.exp(-550/TEMP))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00e-12*O2*3.2*numpy.exp(-550/TEMP)' 
		rate_values[459] = 5.00e-12*O2*3.2*numpy.exp(-550/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[460] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[461] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[462] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[463] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[464] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[465] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[466] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[467] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[468] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[469] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[470] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[471] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.70e-11' 
		rate_values[472] = 3.70e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]+J[19]' 
		rate_values[473] = J[18]+J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.00e-11' 
		rate_values[474] = 4.00e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[20]*2' 
		rate_values[475] = J[20]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.40' 
		rate_values[476] = KDEC*0.40
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.60' 
		rate_values[477] = KDEC*0.60
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.70e-11' 
		rate_values[478] = 3.70e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[479] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.31e-11' 
		rate_values[480] = 2.31e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[481] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.62e-11' 
		rate_values[482] = 2.62e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[17]' 
		rate_values[483] = J[41]+J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.29e-11' 
		rate_values[484] = 2.29e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[485] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.44e-10' 
		rate_values[486] = 1.44e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[487] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[488] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[489] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.36e-10' 
		rate_values[490] = 1.36e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[491] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.77e-11' 
		rate_values[492] = 5.77e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[493] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[494] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.23e-11' 
		rate_values[495] = 1.23e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[496] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[497] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[498] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[499] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[500] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[501] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.22e-10' 
		rate_values[502] = 1.22e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[503] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[504] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[505] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[506] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.67e-11' 
		rate_values[507] = 3.67e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]+J[15]' 
		rate_values[508] = J[34]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.13e-11' 
		rate_values[509] = 8.13e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[510] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[511] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[512] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.60e-12' 
		rate_values[513] = 3.60e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[514] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.60e-12' 
		rate_values[515] = 2.60e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.47e-12' 
		rate_values[516] = 3.47e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[517] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[518] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[519] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[520] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.00e-13' 
		rate_values[521] = 3.00e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.6e-12' 
		rate_values[522] = 4.6e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.75e-11' 
		rate_values[523] = 2.75e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]' 
		rate_values[524] = J[41]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.25e-11' 
		rate_values[525] = 2.25e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]+J[15]' 
		rate_values[526] = J[55]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[527] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.41e-11' 
		rate_values[528] = 2.41e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[529] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[530] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[531] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[532] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2' 
		rate_values[533] = 6.70e-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.28e-11' 
		rate_values[534] = 6.28e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[35]' 
		rate_values[535] = J[41]+J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.85e-11' 
		rate_values[536] = 2.85e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]+J[35]' 
		rate_values[537] = J[55]+J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[538] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.93e-11' 
		rate_values[539] = 5.93e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[540] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.80' 
		rate_values[541] = KDEC*0.80
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.20' 
		rate_values[542] = KDEC*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.69e-11' 
		rate_values[543] = 2.69e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[544] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.00e-11' 
		rate_values[545] = 3.00e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]' 
		rate_values[546] = J[41]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.52e-11' 
		rate_values[547] = 2.52e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[548] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.16e-12' 
		rate_values[549] = 9.16e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[22]' 
		rate_values[550] = J[41]+J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.70e-12' 
		rate_values[551] = 5.70e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[552] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.56e-12' 
		rate_values[553] = 5.56e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[554] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.36e-11' 
		rate_values[555] = 2.36e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[22]' 
		rate_values[556] = J[41]+J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.20e+10*numpy.exp(-3523/TEMP)' 
		rate_values[557] = 4.20e+10*numpy.exp(-3523/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.05e-11' 
		rate_values[558] = 1.05e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[22]' 
		rate_values[559] = J[41]+J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[560] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.20e-12' 
		rate_values[561] = 7.20e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[562] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.02e-11' 
		rate_values[563] = 1.02e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[22]' 
		rate_values[564] = J[41]+J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.60e-12' 
		rate_values[565] = 6.60e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[566] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.29e-11' 
		rate_values[567] = 1.29e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[22]' 
		rate_values[568] = J[41]+J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[569] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[570] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[571] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.39e-12' 
		rate_values[572] = 8.39e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[573] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[574] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6e-12*numpy.exp(305/TEMP)' 
		rate_values[575] = 1.6e-12*numpy.exp(305/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[576] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[577] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*2.4' 
		rate_values[578] = KNO3AL*2.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.9e-12*numpy.exp(575/TEMP)' 
		rate_values[579] = 1.9e-12*numpy.exp(575/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13' 
		rate_values[580] = 8.00e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.70e-12' 
		rate_values[581] = 3.70e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[582] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-14' 
		rate_values[583] = 3e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[584] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.06e-11' 
		rate_values[585] = 7.06e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[586] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.26e-11' 
		rate_values[587] = 1.26e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[588] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.72e-11' 
		rate_values[589] = 6.72e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[590] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[591] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[592] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[593] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[594] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[595] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.7*RO2' 
		rate_values[596] = 1.00e-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.3*RO2' 
		rate_values[597] = 1.00e-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.18e-12' 
		rate_values[598] = 6.18e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[599] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[600] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[601] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.20e-15' 
		rate_values[602] = 1.20e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-14' 
		rate_values[603] = 1.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-15' 
		rate_values[604] = 1.00e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00e-14' 
		rate_values[605] = 7.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.00e-18*H2O' 
		rate_values[606] = 6.00e-18*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-17*H2O' 
		rate_values[607] = 1.00e-17*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.387' 
		rate_values[608] = KRO2HO2*0.387
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[609] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[610] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.2*RO2' 
		rate_values[611] = 2.00e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.6*RO2' 
		rate_values[612] = 2.00e-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.2*RO2' 
		rate_values[613] = 2.00e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.68e-11' 
		rate_values[614] = 3.68e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[615] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[616] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.78e-11' 
		rate_values[617] = 1.78e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.78e-11' 
		rate_values[618] = 2.78e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4e-12' 
		rate_values[619] = 1.4e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.33e-11' 
		rate_values[620] = 7.33e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[621] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.97e-11' 
		rate_values[622] = 6.97e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[623] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.45e-11' 
		rate_values[624] = 2.45e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]*2' 
		rate_values[625] = J[34]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.25e-15' 
		rate_values[626] = 2.25e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.00e-14' 
		rate_values[627] = 3.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[628] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[629] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[630] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[631] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[632] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[633] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[634] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[635] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[636] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[637] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[638] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[639] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[640] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[641] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[642] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[643] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[644] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[645] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[646] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[647] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[648] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[649] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[650] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[651] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[652] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-11' 
		rate_values[653] = 6.70e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[654] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[655] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.125' 
		rate_values[656] = KRO2NO*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.875' 
		rate_values[657] = KRO2NO*0.875
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[658] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.7*RO2' 
		rate_values[659] = 6.70e-15*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.3*RO2' 
		rate_values[660] = 6.70e-15*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.03e-12' 
		rate_values[661] = 8.03e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[662] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[663] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[664] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.278' 
		rate_values[665] = KRO2NO*0.278
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.722' 
		rate_values[666] = KRO2NO*0.722
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[667] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*0.6*RO2' 
		rate_values[668] = 2.50e-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*0.2*RO2' 
		rate_values[669] = 2.50e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*0.2*RO2' 
		rate_values[670] = 2.50e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[671] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[672] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[673] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[674] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[675] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[676] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[677] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[678] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[679] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.104' 
		rate_values[680] = KRO2NO*0.104
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.896' 
		rate_values[681] = KRO2NO*0.896
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[682] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2*0.7' 
		rate_values[683] = 6.70e-15*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2*0.3' 
		rate_values[684] = 6.70e-15*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KNO3AL*5.5' 
		rate_values[685] = 2*KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.33e-10' 
		rate_values[686] = 1.33e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[687] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[688] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[689] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[690] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2' 
		rate_values[691] = 6.70e-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[692] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.118' 
		rate_values[693] = KRO2NO*0.118
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.882' 
		rate_values[694] = KRO2NO*0.882
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[695] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.7*RO2' 
		rate_values[696] = 6.70e-15*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.3*RO2' 
		rate_values[697] = 6.70e-15*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[698] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[699] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[700] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2' 
		rate_values[701] = 6.70e-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[702] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.92e-11*0.232' 
		rate_values[703] = 8.92e-11*0.232
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.92e-11*0.768' 
		rate_values[704] = 8.92e-11*0.768
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[705] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]' 
		rate_values[706] = J[41]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.58e-11' 
		rate_values[707] = 1.58e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*8.0' 
		rate_values[708] = KNO3AL*8.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.44e-11' 
		rate_values[709] = 3.44e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[710] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.16e-12' 
		rate_values[711] = 1.16e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.19e-11' 
		rate_values[712] = 2.19e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.91e-11' 
		rate_values[713] = 2.91e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[714] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[715] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[716] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[717] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL' 
		rate_values[718] = KNO3AL
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.200' 
		rate_values[719] = 1.00e-11*0.200
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.800' 
		rate_values[720] = 1.00e-11*0.800
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[721] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[722] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[723] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[724] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[725] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.6*RO2' 
		rate_values[726] = 8.80e-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[727] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[728] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[729] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[730] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[731] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[732] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[733] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[734] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[735] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.00e-13' 
		rate_values[736] = 9.00e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[737] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[738] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[739] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[740] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[741] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[742] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[743] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.68e-11' 
		rate_values[744] = 6.68e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[745] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[746] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.23e-10' 
		rate_values[747] = 1.23e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[748] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[749] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.07e-11' 
		rate_values[750] = 6.07e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.18e-11' 
		rate_values[751] = 9.18e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[752] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[753] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[754] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[755] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[756] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[757] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[758] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[759] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[760] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[761] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[762] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[763] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.01e-11' 
		rate_values[764] = 8.01e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]' 
		rate_values[765] = J[41]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.03e-11' 
		rate_values[766] = 7.03e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]+J[15]' 
		rate_values[767] = J[55]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[768] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.66e-11' 
		rate_values[769] = 7.66e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[770] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[771] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[772] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[773] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[774] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-10' 
		rate_values[775] = 2.00e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[35]' 
		rate_values[776] = J[41]+J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.23e-11' 
		rate_values[777] = 2.23e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]+J[35]' 
		rate_values[778] = J[54]+J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROSEC*O2' 
		rate_values[779] = KROSEC*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.26e-10' 
		rate_values[780] = 1.26e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[781] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.04e-11' 
		rate_values[782] = 1.04e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[783] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[784] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.138' 
		rate_values[785] = KRO2NO*0.138
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.862' 
		rate_values[786] = KRO2NO*0.862
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[787] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*RO2*0.2' 
		rate_values[788] = 1.30e-12*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*RO2*0.6' 
		rate_values[789] = 1.30e-12*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*RO2*0.2' 
		rate_values[790] = 1.30e-12*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.29e-12' 
		rate_values[791] = 7.29e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.77e-12' 
		rate_values[792] = 6.77e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[793] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.61e-11' 
		rate_values[794] = 3.61e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]' 
		rate_values[795] = J[41]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.56e-11' 
		rate_values[796] = 2.56e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]+J[15]' 
		rate_values[797] = J[55]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.70e+14*numpy.exp(-6643/TEMP)' 
		rate_values[798] = 2.70e+14*numpy.exp(-6643/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-11' 
		rate_values[799] = 2.86e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[800] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[801] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[802] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[803] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2' 
		rate_values[804] = 2.00e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.29e-11' 
		rate_values[805] = 1.29e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[22]' 
		rate_values[806] = J[41]+J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[807] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.05e-11' 
		rate_values[808] = 2.05e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[35]' 
		rate_values[809] = J[41]+J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.37e-12' 
		rate_values[810] = 5.37e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]+J[35]' 
		rate_values[811] = J[55]+J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[812] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.69e-11' 
		rate_values[813] = 1.69e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[814] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.45e-11' 
		rate_values[815] = 3.45e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]' 
		rate_values[816] = J[41]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[817] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[818] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[819] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[820] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[821] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[822] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[823] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[824] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[825] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[826] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[827] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[828] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2*0.2' 
		rate_values[829] = 2.00e-12*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2*0.6' 
		rate_values[830] = 2.00e-12*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2*0.2' 
		rate_values[831] = 2.00e-12*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[832] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[833] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[834] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[835] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[836] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[837] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.7*RO2' 
		rate_values[838] = 1.00e-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.3*RO2' 
		rate_values[839] = 1.00e-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.66e-11' 
		rate_values[840] = 4.66e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[841] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[842] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.55e-11' 
		rate_values[843] = 2.55e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.68e-12' 
		rate_values[844] = 5.68e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[845] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[846] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[847] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[848] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[849] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[850] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[851] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.70e-11' 
		rate_values[852] = 7.70e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]*2' 
		rate_values[853] = J[34]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[854] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[855] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[856] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[857] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[858] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[859] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.75e-12' 
		rate_values[860] = 4.75e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[35]' 
		rate_values[861] = J[41]+J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[862] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[863] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[864] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2' 
		rate_values[865] = 2.00e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.83e-13' 
		rate_values[866] = 8.83e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[867] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.55e-11' 
		rate_values[868] = 7.55e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[22]+J[15]' 
		rate_values[869] = J[41]+J[22]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.19e-11' 
		rate_values[870] = 7.19e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[871] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[872] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[873] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[874] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.6*RO2' 
		rate_values[875] = 8.80e-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[876] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[877] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.25e-11' 
		rate_values[878] = 3.25e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[879] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.00e+04' 
		rate_values[880] = 4.00e+04
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROSEC*O2' 
		rate_values[881] = KROSEC*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.70e-11' 
		rate_values[882] = 1.70e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[883] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.29e-12' 
		rate_values[884] = 3.29e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[53]' 
		rate_values[885] = J[53]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[886] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*8.5' 
		rate_values[887] = KNO3AL*8.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.63e-11' 
		rate_values[888] = 2.63e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[889] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.89e-12' 
		rate_values[890] = 7.89e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[891] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.104' 
		rate_values[892] = KRO2NO*0.104
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.896' 
		rate_values[893] = KRO2NO*0.896
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[894] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2*0.7' 
		rate_values[895] = 6.70e-15*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2*0.3' 
		rate_values[896] = 6.70e-15*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.33e-11' 
		rate_values[897] = 8.33e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[22]+J[15]' 
		rate_values[898] = J[41]+J[22]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[899] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[900] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[901] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[902] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2' 
		rate_values[903] = 6.70e-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.22e-12' 
		rate_values[904] = 3.22e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[905] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[906] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.098' 
		rate_values[907] = KRO2NO*0.098
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.902' 
		rate_values[908] = KRO2NO*0.902
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[909] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[910] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.6*RO2' 
		rate_values[911] = 8.80e-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[912] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[913] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[914] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[915] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[916] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.39e-11' 
		rate_values[917] = 2.39e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]*2' 
		rate_values[918] = J[22]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.70e-11' 
		rate_values[919] = 2.70e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[22]*2' 
		rate_values[920] = J[41]+J[22]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.29e-11' 
		rate_values[921] = 2.29e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[922] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.23e-11' 
		rate_values[923] = 3.23e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[22]*2' 
		rate_values[924] = J[41]+J[22]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[925] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.55e-11' 
		rate_values[926] = 3.55e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[927] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.54e-11*0.890' 
		rate_values[928] = 2.54e-11*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.54e-11*0.110' 
		rate_values[929] = 2.54e-11*0.110
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]*2' 
		rate_values[930] = J[22]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.73e-12' 
		rate_values[931] = 2.73e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.19e-12' 
		rate_values[932] = 6.19e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[933] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.12e-12' 
		rate_values[934] = 1.12e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[935] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[936] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[937] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[938] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[939] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[940] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[941] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[942] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.01e-11' 
		rate_values[943] = 8.01e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[944] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.66e-11' 
		rate_values[945] = 7.66e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[946] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.01e-11' 
		rate_values[947] = 1.01e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[35]' 
		rate_values[948] = J[41]+J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[949] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[950] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.33e-11' 
		rate_values[951] = 1.33e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[952] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KNO3AL*4.0' 
		rate_values[953] = 2*KNO3AL*4.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.39e-11' 
		rate_values[954] = 3.39e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[955] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[956] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.20e-10' 
		rate_values[957] = 1.20e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]' 
		rate_values[958] = J[41]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[959] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.25e-12' 
		rate_values[960] = 1.25e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]' 
		rate_values[961] = J[55]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[962] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[963] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[964] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.7' 
		rate_values[965] = 9.20e-14*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.3' 
		rate_values[966] = 9.20e-14*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[967] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[968] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[969] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[970] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[971] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[972] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[973] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[974] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[975] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[976] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[977] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*RO2' 
		rate_values[978] = 1.30e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.35e-11' 
		rate_values[979] = 8.35e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]' 
		rate_values[980] = J[41]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.96e-11' 
		rate_values[981] = 4.96e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]+J[15]' 
		rate_values[982] = J[55]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.70e+14*numpy.exp(-6643/TEMP)' 
		rate_values[983] = 2.70e+14*numpy.exp(-6643/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-11' 
		rate_values[984] = 8.00e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[985] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[986] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[987] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[988] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[989] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[990] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[991] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.3*RO2' 
		rate_values[992] = 1.00e-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.7*RO2' 
		rate_values[993] = 1.00e-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.51e-11' 
		rate_values[994] = 1.51e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[22]' 
		rate_values[995] = J[41]+J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[996] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[997] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[998] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[999] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.6*RO2' 
		rate_values[1000] = 2.00e-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.2*RO2' 
		rate_values[1001] = 2.00e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.2*RO2' 
		rate_values[1002] = 2.00e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.69e-11' 
		rate_values[1003] = 8.69e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[35]' 
		rate_values[1004] = J[41]+J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '71.11e-12' 
		rate_values[1005] = 71.11e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[1006] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1007] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.78e-11' 
		rate_values[1008] = 3.78e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[1009] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.49e-11' 
		rate_values[1010] = 7.49e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]' 
		rate_values[1011] = J[41]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1012] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[1013] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[1014] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[1015] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1016] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1017] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1018] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[1019] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[1020] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[1021] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.017' 
		rate_values[1022] = KRO2NO*0.017
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.983' 
		rate_values[1023] = KRO2NO*0.983
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1024] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.2*RO2' 
		rate_values[1025] = 2.00e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.6*RO2' 
		rate_values[1026] = 2.00e-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.2*RO2' 
		rate_values[1027] = 2.00e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.06e-11' 
		rate_values[1028] = 3.06e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[34]' 
		rate_values[1029] = J[41]+J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.74e-11' 
		rate_values[1030] = 2.74e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1031] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[1032] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[1033] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1034] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1035] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1036] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[1037] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[1038] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[1039] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1040] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1041] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1042] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[1043] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.520' 
		rate_values[1044] = KRO2HO2*0.520
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1045] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1046] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2' 
		rate_values[1047] = 2.00e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[1048] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1049] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1050] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[1051] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.09e-11' 
		rate_values[1052] = 1.09e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1053] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1054] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.42e-12' 
		rate_values[1055] = 7.42e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.65e-12' 
		rate_values[1056] = 9.65e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1057] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.57e-12' 
		rate_values[1058] = 6.57e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.96e-12' 
		rate_values[1059] = 2.96e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1060] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.27e-11' 
		rate_values[1061] = 1.27e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1062] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1063] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[1064] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.129' 
		rate_values[1065] = KRO2NO*0.129
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.871' 
		rate_values[1066] = KRO2NO*0.871
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1067] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2*0.6' 
		rate_values[1068] = 2.50e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2*0.2' 
		rate_values[1069] = 2.50e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2*0.2' 
		rate_values[1070] = 2.50e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[1071] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.14e-11' 
		rate_values[1072] = 2.14e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]' 
		rate_values[1073] = J[41]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.49e-11' 
		rate_values[1074] = 2.49e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1075] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.10e-11' 
		rate_values[1076] = 2.10e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[1077] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1078] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1079] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[1080] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[35]' 
		rate_values[1081] = J[41]+J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1082] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.99e-12' 
		rate_values[1083] = 5.99e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1084] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[1085] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.69e-12' 
		rate_values[1086] = 2.69e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[1087] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[1088] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*4.0' 
		rate_values[1089] = KNO3AL*4.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.23e-11' 
		rate_values[1090] = 1.23e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[1091] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[1092] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[1093] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.65e-11' 
		rate_values[1094] = 6.65e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[1095] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KNO3AL*2.4' 
		rate_values[1096] = 2*KNO3AL*2.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.29e-11' 
		rate_values[1097] = 4.29e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.34e-11' 
		rate_values[1098] = 2.34e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[1099] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.65e-11' 
		rate_values[1100] = 2.65e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[22]' 
		rate_values[1101] = J[41]+J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.60e-12' 
		rate_values[1102] = 7.60e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1103] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1104] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.77e-11' 
		rate_values[1105] = 5.77e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[56]*0.91' 
		rate_values[1106] = J[56]*0.91
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.23e-12' 
		rate_values[1107] = 2.23e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1108] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[1109] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*4.0' 
		rate_values[1110] = KNO3AL*4.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.45e-11' 
		rate_values[1111] = 2.45e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[1112] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.88e-11' 
		rate_values[1113] = 1.88e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]' 
		rate_values[1114] = J[41]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.23e-12' 
		rate_values[1115] = 4.23e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1116] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.12e-13' 
		rate_values[1117] = 3.12e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.63e-11' 
		rate_values[1118] = 1.63e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[34]' 
		rate_values[1119] = J[41]+J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.27e-11' 
		rate_values[1120] = 1.27e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1121] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.41e-11' 
		rate_values[1122] = 2.41e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[34]' 
		rate_values[1123] = J[41]+J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1124] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.85e-11' 
		rate_values[1125] = 1.85e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1126] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1127] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1128] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.104' 
		rate_values[1129] = KRO2NO*0.104
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.896' 
		rate_values[1130] = KRO2NO*0.896
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1131] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2*0.7' 
		rate_values[1132] = 6.70e-15*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2*0.3' 
		rate_values[1133] = 6.70e-15*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[1134] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1135] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1136] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2' 
		rate_values[1137] = 6.70e-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.10e-10' 
		rate_values[1138] = 1.10e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]*2' 
		rate_values[1139] = J[41]+J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.33e-11' 
		rate_values[1140] = 4.33e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]+J[15]*2' 
		rate_values[1141] = J[54]+J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.80e-14*numpy.exp(-260/TEMP)*O2' 
		rate_values[1142] = 1.80e-14*numpy.exp(-260/TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.99e-11' 
		rate_values[1143] = 6.99e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[1144] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-10' 
		rate_values[1145] = 1.00e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[22]' 
		rate_values[1146] = J[41]+J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1147] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[1148] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[1149] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1150] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1151] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1152] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[1153] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[1154] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[1155] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1156] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1157] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1158] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[1159] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.86e-11' 
		rate_values[1160] = 1.86e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[34]' 
		rate_values[1161] = J[41]+J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.82e-12' 
		rate_values[1162] = 7.82e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]+J[34]' 
		rate_values[1163] = J[55]+J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1164] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.75e-11' 
		rate_values[1165] = 1.75e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[1166] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.31e-11' 
		rate_values[1167] = 3.31e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1168] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1169] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[1170] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.37e-11' 
		rate_values[1171] = 2.37e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[1172] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[1173] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[1174] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1175] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.34e-12' 
		rate_values[1176] = 7.34e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1177] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.74e-12' 
		rate_values[1178] = 3.74e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.68e-11' 
		rate_values[1179] = 1.68e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1180] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.32e-11' 
		rate_values[1181] = 1.32e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1182] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.69e-11' 
		rate_values[1183] = 6.69e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]+J[15]' 
		rate_values[1184] = J[34]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[1185] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1186] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1187] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.8e-13*RO2' 
		rate_values[1188] = 8.8e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[1189] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1190] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1191] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[1192] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[1193] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[1194] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1195] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1196] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1197] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[1198] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[1199] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1200] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1201] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2' 
		rate_values[1202] = 2.00e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]+J[41]' 
		rate_values[1203] = J[22]+J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.38e-11' 
		rate_values[1204] = 3.38e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1205] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.46e-11' 
		rate_values[1206] = 7.46e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1207] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1208] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.55e-12' 
		rate_values[1209] = 6.55e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[35]' 
		rate_values[1210] = J[41]+J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.92e-12' 
		rate_values[1211] = 2.92e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1212] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.61e-12' 
		rate_values[1213] = 9.61e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[35]' 
		rate_values[1214] = J[41]+J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1215] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.44e-11' 
		rate_values[1216] = 1.44e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]+J[35]' 
		rate_values[1217] = J[34]+J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.32' 
		rate_values[1218] = 0.32
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.1*KDEC' 
		rate_values[1219] = 0.1*KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.022' 
		rate_values[1220] = 0.022
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.07' 
		rate_values[1221] = 0.07
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.25' 
		rate_values[1222] = 0.25
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.07' 
		rate_values[1223] = 0.07
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.25' 
		rate_values[1224] = 0.25
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.5' 
		rate_values[1225] = 1.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.45' 
		rate_values[1226] = 0.45
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.1' 
		rate_values[1227] = 0.1
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.01' 
		rate_values[1228] = 0.01
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.009' 
		rate_values[1229] = KRO2NO*0.009
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.005' 
		rate_values[1230] = KRO2NO*0.005
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.12' 
		rate_values[1231] = KRO2NO*0.12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.002' 
		rate_values[1232] = KRO2NO*0.002
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.006' 
		rate_values[1233] = KRO2NO*0.006
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.01' 
		rate_values[1234] = KRO2NO*0.01
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.04' 
		rate_values[1235] = KRO2NO*0.04
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.004' 
		rate_values[1236] = KRO2NO*0.004
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.008' 
		rate_values[1237] = KRO2NO*0.008
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.02' 
		rate_values[1238] = KRO2NO*0.02
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.02' 
		rate_values[1239] = KRO2NO*0.02
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.897' 
		rate_values[1240] = KRO2NO*0.897
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.978' 
		rate_values[1241] = KRO2NO*0.978
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.93' 
		rate_values[1242] = KRO2NO*0.93
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.898' 
		rate_values[1243] = KRO2NO*0.898
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.894' 
		rate_values[1244] = KRO2NO*0.894
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.59' 
		rate_values[1245] = KRO2NO*0.59
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.96' 
		rate_values[1246] = KRO2NO*0.96
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.998' 
		rate_values[1247] = KRO2NO*0.998
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.996' 
		rate_values[1248] = KRO2NO*0.996
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.98' 
		rate_values[1249] = KRO2NO*0.98
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.98' 
		rate_values[1250] = KRO2NO*0.98
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.90' 
		rate_values[1251] = KDEC*0.90
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.85' 
		rate_values[1252] = KDEC*0.85
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.90' 
		rate_values[1253] = KDEC*0.90
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.90' 
		rate_values[1254] = KDEC*0.90
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.40' 
		rate_values[1255] = KDEC*0.40
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.90' 
		rate_values[1256] = KDEC*0.90
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.90' 
		rate_values[1257] = KDEC*0.90
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.85' 
		rate_values[1258] = KDEC*0.85
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.70' 
		rate_values[1259] = KDEC*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.05' 
		rate_values[1260] = KDEC*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.10' 
		rate_values[1261] = KDEC*0.10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.60' 
		rate_values[1262] = KDEC*0.60
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.35' 
		rate_values[1263] = KDEC*0.35
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.05' 
		rate_values[1264] = KDEC*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.05' 
		rate_values[1265] = KDEC*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.55' 
		rate_values[1266] = KDEC*0.55
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.05' 
		rate_values[1267] = KDEC*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.05' 
		rate_values[1268] = KDEC*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.10' 
		rate_values[1269] = KDEC*0.10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.25' 
		rate_values[1270] = KDEC*0.25
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.50' 
		rate_values[1271] = KDEC*0.50
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.45' 
		rate_values[1272] = KDEC*0.45
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.05' 
		rate_values[1273] = KDEC*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.05' 
		rate_values[1274] = KDEC*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.05' 
		rate_values[1275] = KDEC*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.05' 
		rate_values[1276] = KDEC*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.05' 
		rate_values[1277] = KDEC*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.05' 
		rate_values[1278] = KDEC*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.05' 
		rate_values[1279] = KDEC*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.05' 
		rate_values[1280] = KDEC*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.05' 
		rate_values[1281] = KDEC*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.05' 
		rate_values[1282] = KDEC*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.05' 
		rate_values[1283] = KDEC*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.1' 
		rate_values[1284] = KRO2NO*0.1
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.02' 
		rate_values[1285] = KRO2NO*0.02
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.03' 
		rate_values[1286] = KRO2NO*0.03
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.1' 
		rate_values[1287] = KRO2NO*0.1
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.1' 
		rate_values[1288] = KRO2NO*0.1
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.4' 
		rate_values[1289] = KRO2NO*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.77*0.01' 
		rate_values[1290] = KRO2HO2*0.77*0.01
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.77*0.06' 
		rate_values[1291] = KRO2HO2*0.77*0.06
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.77*0.2' 
		rate_values[1292] = KRO2HO2*0.77*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.77*0.6' 
		rate_values[1293] = KRO2HO2*0.77*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.77*0.04' 
		rate_values[1294] = KRO2HO2*0.77*0.04
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.77*0.55' 
		rate_values[1295] = KRO2HO2*0.77*0.55
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706*0.005' 
		rate_values[1296] = KRO2HO2*0.706*0.005
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706*0.3' 
		rate_values[1297] = KRO2HO2*0.706*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706*0.6' 
		rate_values[1298] = KRO2HO2*0.706*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706*0.1' 
		rate_values[1299] = KRO2HO2*0.706*0.1
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706*0.6' 
		rate_values[1300] = KRO2HO2*0.706*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*2.62e-11*0.6*0.4' 
		rate_values[1301] = RO2*2.62e-11*0.6*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*5.38e-11*0.6*0.4' 
		rate_values[1302] = RO2*5.38e-11*0.6*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*8.14e-11*0.6*0.4' 
		rate_values[1303] = RO2*8.14e-11*0.6*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*1.24e-11*0.6*0.4' 
		rate_values[1304] = RO2*1.24e-11*0.6*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*4e-11*0.6*0.4' 
		rate_values[1305] = RO2*4e-11*0.6*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*6.76e-11*0.6*0.4' 
		rate_values[1306] = RO2*6.76e-11*0.6*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.99*0.77' 
		rate_values[1307] = KRO2HO2*0.99*0.77
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.92*0.77' 
		rate_values[1308] = KRO2HO2*0.92*0.77
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.85*0.77' 
		rate_values[1309] = KRO2HO2*0.85*0.77
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.4*0.77' 
		rate_values[1310] = KRO2HO2*0.4*0.77
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.98*0.77' 
		rate_values[1311] = KRO2HO2*0.98*0.77
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.45*0.77' 
		rate_values[1312] = KRO2HO2*0.45*0.77
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.995*0.77' 
		rate_values[1313] = KRO2HO2*0.995*0.77
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.4*0.77' 
		rate_values[1314] = KRO2HO2*0.4*0.77
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.4*0.77' 
		rate_values[1315] = KRO2HO2*0.4*0.77
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.9*0.77' 
		rate_values[1316] = KRO2HO2*0.9*0.77
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.4*0.77' 
		rate_values[1317] = KRO2HO2*0.4*0.77
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*2.62e-11*0.2*0.4' 
		rate_values[1318] = RO2*2.62e-11*0.2*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*5.38e-11*0.2*0.4' 
		rate_values[1319] = RO2*5.38e-11*0.2*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*8.14e-11*0.2*0.4' 
		rate_values[1320] = RO2*8.14e-11*0.2*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*1.24e-11*0.2*0.4' 
		rate_values[1321] = RO2*1.24e-11*0.2*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*4e-11*0.2*0.4' 
		rate_values[1322] = RO2*4e-11*0.2*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*6.76e-11*0.2*0.4' 
		rate_values[1323] = RO2*6.76e-11*0.2*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*2.97e-11*0.2*0.4' 
		rate_values[1324] = RO2*2.97e-11*0.2*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*4.34e-11*0.2*0.4' 
		rate_values[1325] = RO2*4.34e-11*0.2*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*5.72e-11*0.2*0.4' 
		rate_values[1326] = RO2*5.72e-11*0.2*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*7.1e-11*0.2*0.4' 
		rate_values[1327] = RO2*7.1e-11*0.2*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*8.48e-11*0.2*0.4' 
		rate_values[1328] = RO2*8.48e-11*0.2*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*2.62e-11*0.21*0.4' 
		rate_values[1329] = RO2*2.62e-11*0.21*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*5.38e-11*0.21*0.4' 
		rate_values[1330] = RO2*5.38e-11*0.21*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*8.14e-11*0.21*0.4' 
		rate_values[1331] = RO2*8.14e-11*0.21*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*1.24e-11*0.21*0.4' 
		rate_values[1332] = RO2*1.24e-11*0.21*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*4e-11*0.21*0.4' 
		rate_values[1333] = RO2*4e-11*0.21*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*6.76e-11*0.21*0.4' 
		rate_values[1334] = RO2*6.76e-11*0.21*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*2.97e-11*0.21*0.4' 
		rate_values[1335] = RO2*2.97e-11*0.21*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*4.34e-11*0.21*0.4' 
		rate_values[1336] = RO2*4.34e-11*0.21*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*5.72e-11*0.21*0.4' 
		rate_values[1337] = RO2*5.72e-11*0.21*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*7.1e-11*0.21*0.4' 
		rate_values[1338] = RO2*7.1e-11*0.21*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'RO2*8.48e-11*0.21*0.4' 
		rate_values[1339] = RO2*8.48e-11*0.21*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.17e-11*1.0' 
		rate_values[1340] = 4.17e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.e-10*1.0' 
		rate_values[1341] = 1.e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.69e-11*1.0' 
		rate_values[1342] = 9.69e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.79e-11*1.0' 
		rate_values[1343] = 2.79e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.55e-11*1.0' 
		rate_values[1344] = 5.55e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.31e-11*1.0' 
		rate_values[1345] = 8.31e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.51e-11*1.0' 
		rate_values[1346] = 4.51e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.89e-11*1.0' 
		rate_values[1347] = 5.89e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.27e-11*1.0' 
		rate_values[1348] = 7.27e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.65e-11*1.0' 
		rate_values[1349] = 8.65e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1350] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.41e-11*0.15' 
		rate_values[1351] = 1.41e-11*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.24e-11*1.0' 
		rate_values[1352] = 1.24e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.79e-11*1.0' 
		rate_values[1353] = 2.79e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1354] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.67e-11*1.0' 
		rate_values[1355] = 6.67e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.79e-11*1.0' 
		rate_values[1356] = 2.79e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1357] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1358] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1359] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1360] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1361] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.05e-11*1.0' 
		rate_values[1362] = 3.05e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1363] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1364] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1365] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1366] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1367] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.78e-12*1.0' 
		rate_values[1368] = 9.78e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1369] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.05e-11*1.0' 
		rate_values[1370] = 8.05e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.12e-11*1.0' 
		rate_values[1371] = 5.12e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.24e-11*1.0' 
		rate_values[1372] = 1.24e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.75e-12*1.0' 
		rate_values[1373] = 3.75e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1374] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1375] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1376] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.62e-11*1.0' 
		rate_values[1377] = 2.62e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1378] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.69e-11*1.0' 
		rate_values[1379] = 9.69e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1380] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.55e-11*1.0' 
		rate_values[1381] = 5.55e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.31e-11*1.0' 
		rate_values[1382] = 8.31e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1383] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.27e-11*1.0' 
		rate_values[1384] = 7.27e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.65e-11*1.0' 
		rate_values[1385] = 8.65e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1386] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1387] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1388] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.17e-11*0.5' 
		rate_values[1389] = 4.17e-11*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4e-11*1.0' 
		rate_values[1390] = 4e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.55e-11*1.0' 
		rate_values[1391] = 5.55e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.58e-11*1.0' 
		rate_values[1392] = 1.58e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.43e-11*1.0' 
		rate_values[1393] = 9.43e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.55e-11*1.0' 
		rate_values[1394] = 5.55e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.75e-12*1.0' 
		rate_values[1395] = 3.75e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.75e-11*1.0' 
		rate_values[1396] = 1.75e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1397] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.93e-11*1.0' 
		rate_values[1398] = 1.93e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1399] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.81e-11*1.0' 
		rate_values[1400] = 5.81e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.93e-11*1.0' 
		rate_values[1401] = 1.93e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.2e-12*1.0' 
		rate_values[1402] = 7.2e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.24e-11*1.0' 
		rate_values[1403] = 1.24e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1404] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1405] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.74e-11*1.0' 
		rate_values[1406] = 3.74e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1407] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1408] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.87e-11*1.0' 
		rate_values[1409] = 7.87e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4e-11*1.0' 
		rate_values[1410] = 4e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.13e-11*1.0' 
		rate_values[1411] = 3.13e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1412] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1413] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1414] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.37e-11*1.0' 
		rate_values[1415] = 5.37e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.79e-11*1.0' 
		rate_values[1416] = 2.79e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1417] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.31e-11*1.0' 
		rate_values[1418] = 8.31e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1419] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1420] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1421] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1422] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1423] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1424] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1425] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.e-10*1.0' 
		rate_values[1426] = 1.e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.75e-11*1.0' 
		rate_values[1427] = 6.75e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.31e-11*1.0' 
		rate_values[1428] = 8.31e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.34e-11*1.0' 
		rate_values[1429] = 4.34e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1430] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.31e-11*1.0' 
		rate_values[1431] = 8.31e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.13e-11*1.0' 
		rate_values[1432] = 3.13e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.51e-11*1.0' 
		rate_values[1433] = 4.51e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1e-11*1.0' 
		rate_values[1434] = 2.1e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.69e-11*1.0' 
		rate_values[1435] = 4.69e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.62e-11*1.0' 
		rate_values[1436] = 2.62e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.56e-11*1.0' 
		rate_values[1437] = 8.56e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.69e-11*1.0' 
		rate_values[1438] = 4.69e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.48e-11*1.0' 
		rate_values[1439] = 3.48e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4e-11*1.0' 
		rate_values[1440] = 4e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.92e-12*1.0' 
		rate_values[1441] = 8.92e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1442] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.5e-11*1.0' 
		rate_values[1443] = 6.5e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1444] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1445] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1446] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.75e-11*1.0' 
		rate_values[1447] = 6.75e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.89e-11*1.0' 
		rate_values[1448] = 5.89e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1449] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1450] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.06e-11*1.0' 
		rate_values[1451] = 1.06e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.13e-11*1.0' 
		rate_values[1452] = 8.13e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.55e-11*1.0' 
		rate_values[1453] = 5.55e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.41e-11*1.0' 
		rate_values[1454] = 1.41e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.17e-11*1.0' 
		rate_values[1455] = 4.17e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.e-10*1.0' 
		rate_values[1456] = 1.e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.13e-11*1.0' 
		rate_values[1457] = 3.13e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.51e-11*1.0' 
		rate_values[1458] = 4.51e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.89e-11*1.0' 
		rate_values[1459] = 5.89e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.27e-11*1.0' 
		rate_values[1460] = 7.27e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.65e-11*1.0' 
		rate_values[1461] = 8.65e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1462] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1463] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.41e-11*1.0' 
		rate_values[1464] = 1.41e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1465] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.29e-11*1.0' 
		rate_values[1466] = 5.29e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.41e-11*1.0' 
		rate_values[1467] = 1.41e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1468] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1469] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1470] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1471] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1472] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.67e-11*1.0' 
		rate_values[1473] = 1.67e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1474] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1475] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1476] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1477] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1478] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1479] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1480] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.67e-11*1.0' 
		rate_values[1481] = 6.67e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.74e-11*1.0' 
		rate_values[1482] = 3.74e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1483] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1484] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1485] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.17e-11*1.0' 
		rate_values[1486] = 9.17e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1487] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.24e-11*1.0' 
		rate_values[1488] = 1.24e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1489] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.e-10*1.0' 
		rate_values[1490] = 1.e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.69e-11*1.0' 
		rate_values[1491] = 9.69e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.89e-11*1.0' 
		rate_values[1492] = 5.89e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.27e-11*1.0' 
		rate_values[1493] = 7.27e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.65e-11*1.0' 
		rate_values[1494] = 8.65e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1495] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1496] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.79e-11*0.3' 
		rate_values[1497] = 2.79e-11*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.62e-11*1.0' 
		rate_values[1498] = 2.62e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.17e-11*1.0' 
		rate_values[1499] = 4.17e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.02e-12*1.0' 
		rate_values[1500] = 2.02e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.05e-11*1.0' 
		rate_values[1501] = 8.05e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.17e-11*1.0' 
		rate_values[1502] = 4.17e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1503] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.75e-12*1.0' 
		rate_values[1504] = 3.75e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1505] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.47e-12*1.0' 
		rate_values[1506] = 5.47e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1507] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.43e-11*1.0' 
		rate_values[1508] = 4.43e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.47e-12*1.0' 
		rate_values[1509] = 5.47e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1510] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1511] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1512] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1513] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.36e-11*1.0' 
		rate_values[1514] = 2.36e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1515] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.43e-11*1.0' 
		rate_values[1516] = 9.43e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.5e-11*1.0' 
		rate_values[1517] = 6.5e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.62e-11*1.0' 
		rate_values[1518] = 2.62e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.75e-11*1.0' 
		rate_values[1519] = 1.75e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1520] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1521] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1522] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4e-11*1.0' 
		rate_values[1523] = 4e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.41e-11*1.0' 
		rate_values[1524] = 1.41e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1525] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.65e-11*1.0' 
		rate_values[1526] = 8.65e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1527] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1528] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1529] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1530] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.55e-11*1.0' 
		rate_values[1531] = 5.55e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.37e-11*1.0' 
		rate_values[1532] = 5.37e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.e-10*1.0' 
		rate_values[1533] = 1.e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.96e-11*1.0' 
		rate_values[1534] = 2.96e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1535] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.e-10*1.0' 
		rate_values[1536] = 1.e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.75e-11*1.0' 
		rate_values[1537] = 1.75e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.13e-11*1.0' 
		rate_values[1538] = 3.13e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.2e-12*1.0' 
		rate_values[1539] = 7.2e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.31e-11*1.0' 
		rate_values[1540] = 3.31e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.24e-11*1.0' 
		rate_values[1541] = 1.24e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.19e-11*1.0' 
		rate_values[1542] = 7.19e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.31e-11*1.0' 
		rate_values[1543] = 3.31e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1e-11*1.0' 
		rate_values[1544] = 2.1e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.62e-11*1.0' 
		rate_values[1545] = 2.62e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1546] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1547] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.12e-11*1.0' 
		rate_values[1548] = 5.12e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1549] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1550] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.25e-11*1.0' 
		rate_values[1551] = 9.25e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.37e-11*1.0' 
		rate_values[1552] = 5.37e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.51e-11*1.0' 
		rate_values[1553] = 4.51e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1554] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1555] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1556] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.75e-11*1.0' 
		rate_values[1557] = 6.75e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.17e-11*1.0' 
		rate_values[1558] = 4.17e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.86e-11*1.0' 
		rate_values[1559] = 4.86e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.24e-11*1.0' 
		rate_values[1560] = 6.24e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.62e-11*1.0' 
		rate_values[1561] = 7.62e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9e-11*1.0' 
		rate_values[1562] = 9e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1563] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.75e-11*1.0' 
		rate_values[1564] = 1.75e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.58e-11*1.0' 
		rate_values[1565] = 1.58e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.13e-11*1.0' 
		rate_values[1566] = 3.13e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1567] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.01e-11*1.0' 
		rate_values[1568] = 7.01e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.13e-11*1.0' 
		rate_values[1569] = 3.13e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1570] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1571] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1572] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1573] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1574] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.39e-11*1.0' 
		rate_values[1575] = 3.39e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1576] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1577] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1578] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1579] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1580] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.32e-11*1.0' 
		rate_values[1581] = 1.32e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1582] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.39e-11*1.0' 
		rate_values[1583] = 8.39e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.46e-11*1.0' 
		rate_values[1584] = 5.46e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.58e-11*1.0' 
		rate_values[1585] = 1.58e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.2e-12*1.0' 
		rate_values[1586] = 7.2e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1587] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1588] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1589] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.96e-11*1.0' 
		rate_values[1590] = 2.96e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.75e-12*1.0' 
		rate_values[1591] = 3.75e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.62e-11*1.0' 
		rate_values[1592] = 7.62e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9e-11*1.0' 
		rate_values[1593] = 9e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1594] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1595] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.13e-11*1.0' 
		rate_values[1596] = 3.13e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.96e-11*1.0' 
		rate_values[1597] = 2.96e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.51e-11*1.0' 
		rate_values[1598] = 4.51e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.47e-12*1.0' 
		rate_values[1599] = 5.47e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.39e-11*1.0' 
		rate_values[1600] = 8.39e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.51e-11*1.0' 
		rate_values[1601] = 4.51e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1602] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.2e-12*1.0' 
		rate_values[1603] = 7.2e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1604] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.92e-12*1.0' 
		rate_values[1605] = 8.92e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1606] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.77e-11*1.0' 
		rate_values[1607] = 4.77e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.92e-12*1.0' 
		rate_values[1608] = 8.92e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1609] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.02e-12*1.0' 
		rate_values[1610] = 2.02e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1611] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1612] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.7e-11*1.0' 
		rate_values[1613] = 2.7e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1614] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.77e-11*1.0' 
		rate_values[1615] = 9.77e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.84e-11*1.0' 
		rate_values[1616] = 6.84e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.96e-11*1.0' 
		rate_values[1617] = 2.96e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1e-11*1.0' 
		rate_values[1618] = 2.1e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1619] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1620] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1621] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.34e-11*1.0' 
		rate_values[1622] = 4.34e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.75e-11*1.0' 
		rate_values[1623] = 1.75e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1624] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1625] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1626] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.51e-11*0.3' 
		rate_values[1627] = 4.51e-11*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.34e-11*1.0' 
		rate_values[1628] = 4.34e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.89e-11*1.0' 
		rate_values[1629] = 5.89e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.93e-11*1.0' 
		rate_values[1630] = 1.93e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.77e-11*1.0' 
		rate_values[1631] = 9.77e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.89e-11*1.0' 
		rate_values[1632] = 5.89e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.2e-12*1.0' 
		rate_values[1633] = 7.2e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1e-11*1.0' 
		rate_values[1634] = 2.1e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1635] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.27e-11*1.0' 
		rate_values[1636] = 2.27e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.02e-12*1.0' 
		rate_values[1637] = 2.02e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.15e-11*1.0' 
		rate_values[1638] = 6.15e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.27e-11*1.0' 
		rate_values[1639] = 2.27e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.06e-11*1.0' 
		rate_values[1640] = 1.06e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.58e-11*1.0' 
		rate_values[1641] = 1.58e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1642] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1643] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.08e-11*1.0' 
		rate_values[1644] = 4.08e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1645] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1646] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.22e-11*1.0' 
		rate_values[1647] = 8.22e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.34e-11*1.0' 
		rate_values[1648] = 4.34e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.48e-11*1.0' 
		rate_values[1649] = 3.48e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1650] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1651] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1652] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.72e-11*1.0' 
		rate_values[1653] = 5.72e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.13e-11*1.0' 
		rate_values[1654] = 3.13e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1655] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1656] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.89e-11*1.0' 
		rate_values[1657] = 5.89e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.72e-11*1.0' 
		rate_values[1658] = 5.72e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.27e-11*1.0' 
		rate_values[1659] = 7.27e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.31e-11*1.0' 
		rate_values[1660] = 3.31e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1661] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.27e-11*1.0' 
		rate_values[1662] = 7.27e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1e-11*1.0' 
		rate_values[1663] = 2.1e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.48e-11*1.0' 
		rate_values[1664] = 3.48e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.06e-11*1.0' 
		rate_values[1665] = 1.06e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.65e-11*1.0' 
		rate_values[1666] = 3.65e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.58e-11*1.0' 
		rate_values[1667] = 1.58e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.53e-11*1.0' 
		rate_values[1668] = 7.53e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.65e-11*1.0' 
		rate_values[1669] = 3.65e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.44e-11*1.0' 
		rate_values[1670] = 2.44e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.96e-11*1.0' 
		rate_values[1671] = 2.96e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1672] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1673] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.46e-11*1.0' 
		rate_values[1674] = 5.46e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1675] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1676] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.6e-11*1.0' 
		rate_values[1677] = 9.6e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.72e-11*1.0' 
		rate_values[1678] = 5.72e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.86e-11*1.0' 
		rate_values[1679] = 4.86e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1680] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1681] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1682] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.1e-11*1.0' 
		rate_values[1683] = 7.1e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.51e-11*1.0' 
		rate_values[1684] = 4.51e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1685] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.27e-11*1.0' 
		rate_values[1686] = 7.27e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.1e-11*1.0' 
		rate_values[1687] = 7.1e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.65e-11*1.0' 
		rate_values[1688] = 8.65e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.69e-11*1.0' 
		rate_values[1689] = 4.69e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1690] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.65e-11*1.0' 
		rate_values[1691] = 8.65e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.48e-11*1.0' 
		rate_values[1692] = 3.48e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.86e-11*1.0' 
		rate_values[1693] = 4.86e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.44e-11*1.0' 
		rate_values[1694] = 2.44e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.03e-11*1.0' 
		rate_values[1695] = 5.03e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.96e-11*1.0' 
		rate_values[1696] = 2.96e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.91e-11*1.0' 
		rate_values[1697] = 8.91e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.03e-11*1.0' 
		rate_values[1698] = 5.03e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.82e-11*1.0' 
		rate_values[1699] = 3.82e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.34e-11*1.0' 
		rate_values[1700] = 4.34e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.24e-11*1.0' 
		rate_values[1701] = 1.24e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1702] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.84e-11*1.0' 
		rate_values[1703] = 6.84e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1704] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1705] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1706] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.1e-11*1.0' 
		rate_values[1707] = 7.1e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.24e-11*1.0' 
		rate_values[1708] = 6.24e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1709] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1710] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.41e-11*1.0' 
		rate_values[1711] = 1.41e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.48e-11*1.0' 
		rate_values[1712] = 8.48e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.89e-11*1.0' 
		rate_values[1713] = 5.89e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*0.7' 
		rate_values[1714] = 3e-13*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1715] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1716] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1717] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.91e-11*1.0' 
		rate_values[1718] = 3.91e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1719] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1720] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1721] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1722] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1723] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1724] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.89e-12*1.0' 
		rate_values[1725] = 2.89e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1726] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1727] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1728] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1729] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1730] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1731] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.17e-11*1.0' 
		rate_values[1732] = 9.17e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.29e-11*1.0' 
		rate_values[1733] = 5.29e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.36e-11*1.0' 
		rate_values[1734] = 2.36e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1735] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1736] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1737] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.79e-11*1.0' 
		rate_values[1738] = 7.79e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1739] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1740] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1741] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1742] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1743] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1744] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.74e-11*1.0' 
		rate_values[1745] = 3.74e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1746] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1747] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1748] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1749] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1750] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1751] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.16e-12*1.0' 
		rate_values[1752] = 1.16e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1753] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1754] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1755] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1756] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1757] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1758] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9e-11*1.0' 
		rate_values[1759] = 9e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.12e-11*1.0' 
		rate_values[1760] = 5.12e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.19e-11*1.0' 
		rate_values[1761] = 2.19e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1762] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1763] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1764] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.62e-11*1.0' 
		rate_values[1765] = 7.62e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1766] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1767] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1768] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.41e-11*1.0' 
		rate_values[1769] = 1.41e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1770] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.29e-11*1.0' 
		rate_values[1771] = 5.29e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.41e-11*1.0' 
		rate_values[1772] = 1.41e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1773] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1774] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1775] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1776] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1777] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.67e-11*1.0' 
		rate_values[1778] = 1.67e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1779] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1780] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1781] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1782] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1783] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1784] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1785] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.67e-11*1.0' 
		rate_values[1786] = 6.67e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.74e-11*1.0' 
		rate_values[1787] = 3.74e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1788] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1789] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1790] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.17e-11*1.0' 
		rate_values[1791] = 9.17e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1792] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.24e-11*1.0' 
		rate_values[1793] = 1.24e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1794] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1795] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.32e-11*1.0' 
		rate_values[1796] = 1.32e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1797] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1798] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1799] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1800] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1801] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1802] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1803] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1804] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1805] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1806] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1807] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1808] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1809] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.58e-11*1.0' 
		rate_values[1810] = 6.58e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.7e-11*1.0' 
		rate_values[1811] = 2.7e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1812] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1813] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1814] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.08e-11*1.0' 
		rate_values[1815] = 9.08e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2e-11*1.0' 
		rate_values[1816] = 5.2e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1817] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1818] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1819] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.17e-11*1.0' 
		rate_values[1820] = 9.17e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.29e-11*1.0' 
		rate_values[1821] = 5.29e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.16e-12*1.0' 
		rate_values[1822] = 1.16e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.5e-11*1.0' 
		rate_values[1823] = 1.5e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1824] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.67e-11*1.0' 
		rate_values[1825] = 1.67e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1826] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.55e-11*1.0' 
		rate_values[1827] = 5.55e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.67e-11*1.0' 
		rate_values[1828] = 1.67e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.61e-12*1.0' 
		rate_values[1829] = 4.61e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.78e-12*1.0' 
		rate_values[1830] = 9.78e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1831] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1832] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.48e-11*1.0' 
		rate_values[1833] = 3.48e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1834] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1835] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.62e-11*1.0' 
		rate_values[1836] = 7.62e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.74e-11*1.0' 
		rate_values[1837] = 3.74e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.87e-11*1.0' 
		rate_values[1838] = 2.87e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1839] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1840] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1841] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.12e-11*1.0' 
		rate_values[1842] = 5.12e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.53e-11*1.0' 
		rate_values[1843] = 2.53e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.41e-11*1.0' 
		rate_values[1844] = 1.41e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1845] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1846] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1847] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1848] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1849] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.67e-11*1.0' 
		rate_values[1850] = 1.67e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1851] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1852] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1853] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1854] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1855] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1856] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1857] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.67e-11*1.0' 
		rate_values[1858] = 6.67e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.74e-11*1.0' 
		rate_values[1859] = 3.74e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1860] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1861] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1862] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.17e-11*1.0' 
		rate_values[1863] = 9.17e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1864] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.24e-11*1.0' 
		rate_values[1865] = 1.24e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1866] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1867] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1868] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1869] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1870] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1871] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1872] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1873] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1874] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1875] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1876] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1877] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1878] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.37e-11*1.0' 
		rate_values[1879] = 5.37e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.5e-11*1.0' 
		rate_values[1880] = 1.5e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1881] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1882] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1883] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.87e-11*1.0' 
		rate_values[1884] = 7.87e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4e-11*1.0' 
		rate_values[1885] = 4e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1886] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1887] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1888] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1889] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1890] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1891] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1892] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1893] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1894] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1895] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1896] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1897] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1898] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1899] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.75e-11*1.0' 
		rate_values[1900] = 6.75e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.87e-11*1.0' 
		rate_values[1901] = 2.87e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1902] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1903] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1904] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.25e-11*1.0' 
		rate_values[1905] = 9.25e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.37e-11*1.0' 
		rate_values[1906] = 5.37e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1907] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1908] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1909] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1910] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1911] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1912] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1913] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1914] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1915] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1916] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1917] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1918] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1919] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.34e-11*1.0' 
		rate_values[1920] = 4.34e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.61e-12*1.0' 
		rate_values[1921] = 4.61e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1922] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1923] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1924] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.84e-11*1.0' 
		rate_values[1925] = 6.84e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.96e-11*1.0' 
		rate_values[1926] = 2.96e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1927] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1928] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1929] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1930] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1931] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1932] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1933] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1934] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1935] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1936] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1937] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1938] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.93e-11*1.0' 
		rate_values[1939] = 6.93e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.05e-11*1.0' 
		rate_values[1940] = 3.05e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.16e-12*1.0' 
		rate_values[1941] = 1.16e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1942] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1943] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.43e-11*1.0' 
		rate_values[1944] = 9.43e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.55e-11*1.0' 
		rate_values[1945] = 5.55e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1946] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1947] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1948] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1949] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1950] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1951] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1952] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1953] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1954] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1955] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1956] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.86e-11*1.0' 
		rate_values[1957] = 4.86e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.78e-12*1.0' 
		rate_values[1958] = 9.78e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1959] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1960] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1961] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.36e-11*1.0' 
		rate_values[1962] = 7.36e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.48e-11*1.0' 
		rate_values[1963] = 3.48e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1964] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1965] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1966] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.93e-11*1.0' 
		rate_values[1967] = 1.93e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1968] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1969] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1970] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1971] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1972] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1973] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1974] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.93e-11*1.0' 
		rate_values[1975] = 6.93e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4e-11*1.0' 
		rate_values[1976] = 4e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.16e-12*1.0' 
		rate_values[1977] = 1.16e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1978] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[1979] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.43e-11*1.0' 
		rate_values[1980] = 9.43e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1981] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.5e-11*1.0' 
		rate_values[1982] = 1.5e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1983] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1984] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1985] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1986] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1987] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1988] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1989] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.93e-11*1.0' 
		rate_values[1990] = 6.93e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.05e-11*1.0' 
		rate_values[1991] = 3.05e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.16e-12*1.0' 
		rate_values[1992] = 1.16e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1993] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1994] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.43e-11*1.0' 
		rate_values[1995] = 9.43e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.55e-11*1.0' 
		rate_values[1996] = 5.55e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1997] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1998] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[1999] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2000] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2001] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2002] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2003] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2004] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.72e-11*1.0' 
		rate_values[2005] = 5.72e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.84e-11*1.0' 
		rate_values[2006] = 1.84e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2007] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2008] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2009] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.22e-11*1.0' 
		rate_values[2010] = 8.22e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.34e-11*1.0' 
		rate_values[2011] = 4.34e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2012] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2013] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2014] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2015] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2016] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2017] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2018] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.24e-11*1.0' 
		rate_values[2019] = 6.24e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.36e-11*1.0' 
		rate_values[2020] = 2.36e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2021] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2022] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2023] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.74e-11*1.0' 
		rate_values[2024] = 8.74e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.86e-11*1.0' 
		rate_values[2025] = 4.86e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2026] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2027] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2028] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2029] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2030] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2031] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.13e-11*1.0' 
		rate_values[2032] = 3.13e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2033] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2034] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2035] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2036] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.63e-11*1.0' 
		rate_values[2037] = 5.63e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.75e-11*1.0' 
		rate_values[2038] = 1.75e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2039] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2040] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2041] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2042] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2043] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.93e-11*1.0' 
		rate_values[2044] = 1.93e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2045] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2046] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2047] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2048] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.43e-11*1.0' 
		rate_values[2049] = 4.43e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.47e-12*1.0' 
		rate_values[2050] = 5.47e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2051] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2052] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2053] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2054] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.74e-11*1.0' 
		rate_values[2055] = 8.74e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.86e-11*1.0' 
		rate_values[2056] = 4.86e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.93e-11*1.0' 
		rate_values[2057] = 1.93e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2058] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2059] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2060] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.36e-11*1.0' 
		rate_values[2061] = 7.36e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2062] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2063] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2064] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2065] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2066] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2067] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9e-11*1.0' 
		rate_values[2068] = 9e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.13e-11*1.0' 
		rate_values[2069] = 8.13e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2070] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2071] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.31e-11*1.0' 
		rate_values[2072] = 3.31e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2073] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.79e-11*1.0' 
		rate_values[2074] = 7.79e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2075] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9e-11*1.0' 
		rate_values[2076] = 9e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.12e-11*1.0' 
		rate_values[2077] = 5.12e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.25e-11*1.0' 
		rate_values[2078] = 4.25e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2079] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2080] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2081] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.5e-11*1.0' 
		rate_values[2082] = 6.5e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.91e-11*1.0' 
		rate_values[2083] = 3.91e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.06e-11*1.0' 
		rate_values[2084] = 6.06e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.19e-11*1.0' 
		rate_values[2085] = 2.19e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.32e-11*1.0' 
		rate_values[2086] = 1.32e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2087] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2088] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2089] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.56e-11*1.0' 
		rate_values[2090] = 3.56e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.78e-12*1.0' 
		rate_values[2091] = 9.78e-12*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2092] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2093] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2094] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.62e-11*1.0' 
		rate_values[2095] = 7.62e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2096] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2097] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2098] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2099] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2100] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.75e-11*1.0' 
		rate_values[2101] = 6.75e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2102] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2103] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2104] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2105] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2106] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.81e-11*1.0' 
		rate_values[2107] = 5.81e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2108] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2109] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1e-10*1.0' 
		rate_values[2110] = 1e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.93e-11*1.0' 
		rate_values[2111] = 1.93e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9e-11*1.0' 
		rate_values[2112] = 9e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.41e-11*1.0' 
		rate_values[2113] = 6.41e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2114] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2115] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2116] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.06e-11*1.0' 
		rate_values[2117] = 1.06e-11*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2118] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-13*1.0' 
		rate_values[2119] = 3e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.2' 
		rate_values[2120] = KDEC*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.1' 
		rate_values[2121] = KDEC*0.1
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.45' 
		rate_values[2122] = KDEC*0.45
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.35' 
		rate_values[2123] = KDEC*0.35
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.45' 
		rate_values[2124] = KDEC*0.45
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.004' 
		rate_values[2125] = KDEC*0.004
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.001' 
		rate_values[2126] = KDEC*0.001
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.02' 
		rate_values[2127] = KDEC*0.02
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.45' 
		rate_values[2128] = KDEC*0.45
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.007' 
		rate_values[2129] = KDEC*0.007
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.45' 
		rate_values[2130] = KDEC*0.45
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.8' 
		rate_values[2131] = KDEC*0.8
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.88' 
		rate_values[2132] = KDEC*0.88
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.43' 
		rate_values[2133] = KDEC*0.43
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.1' 
		rate_values[2134] = KDEC*0.1
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.64' 
		rate_values[2135] = KDEC*0.64
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.1' 
		rate_values[2136] = KDEC*0.1
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6' 
		rate_values[2137] = 0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.014' 
		rate_values[2138] = 0.014
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.8' 
		rate_values[2139] = 0.8
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.06' 
		rate_values[2140] = 0.06
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.01' 
		rate_values[2141] = 0.01
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.01' 
		rate_values[2142] = 0.01
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.01' 
		rate_values[2143] = 0.01
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.01' 
		rate_values[2144] = 0.01
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.01' 
		rate_values[2145] = 0.01
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.01' 
		rate_values[2146] = 0.01
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.3' 
		rate_values[2147] = 0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.8' 
		rate_values[2148] = 0.8
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.8' 
		rate_values[2149] = 0.8
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.01' 
		rate_values[2150] = 0.01
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.01' 
		rate_values[2151] = 0.01
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.55*0.90' 
		rate_values[2152] = KDEC*0.55*0.90
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.45*0.90' 
		rate_values[2153] = KDEC*0.45*0.90
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.10' 
		rate_values[2154] = KDEC*0.10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5*0.7' 
		rate_values[2155] = KDEC*0.5*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5*0.7' 
		rate_values[2156] = KDEC*0.5*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.30' 
		rate_values[2157] = KDEC*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.025*1.2E-11*numpy.exp(440/TEMP)' 
		rate_values[2158] = 0.025*1.2E-11*numpy.exp(440/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.01*2.38E-11*numpy.exp(357/TEMP)' 
		rate_values[2159] = 0.01*2.38E-11*numpy.exp(357/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.01*4.28e-11*numpy.exp(401/TEMP)' 
		rate_values[2160] = 0.01*4.28e-11*numpy.exp(401/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.025*1.60e-11*numpy.exp(500./TEMP)' 
		rate_values[2161] = 0.025*1.60e-11*numpy.exp(500./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e17*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[2162] = 1.2e17*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e18*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[2163] = 1.2e18*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e18*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[2164] = 1.2e18*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6e17*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[2165] = 6e17*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2e17*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[2166] = 2e17*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e17*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[2167] = 1.2e17*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e17*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[2168] = 1.2e17*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e16*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[2169] = 1.2e16*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2170] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2171] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2172] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2173] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2174] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2175] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2176] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2177] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2178] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2179] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2180] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*1.0' 
		rate_values[2181] = KRO2NO*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*1.0' 
		rate_values[2182] = KRO2NO*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*1.0' 
		rate_values[2183] = KRO2NO*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.9' 
		rate_values[2184] = KRO2NO*0.9
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.85' 
		rate_values[2185] = KRO2NO*0.85
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.8' 
		rate_values[2186] = KRO2NO*0.8
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.6' 
		rate_values[2187] = KRO2NO*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.5' 
		rate_values[2188] = KRO2NO*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.3' 
		rate_values[2189] = KRO2NO*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2190] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2191] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2192] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2193] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2194] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.1' 
		rate_values[2195] = 0.4*KRO2NO*0.1
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.15' 
		rate_values[2196] = 0.4*KRO2NO*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.2' 
		rate_values[2197] = 0.4*KRO2NO*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.4' 
		rate_values[2198] = 0.4*KRO2NO*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.5' 
		rate_values[2199] = 0.4*KRO2NO*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.7' 
		rate_values[2200] = 0.4*KRO2NO*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*1.0' 
		rate_values[2201] = 0.4*KRO2NO*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*1.0' 
		rate_values[2202] = 0.4*KRO2NO*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2203] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2204] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2205] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.1*7e0/10e0' 
		rate_values[2206] = 0.6*KRO2NO*0.1*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.15*7e0/10e0' 
		rate_values[2207] = 0.6*KRO2NO*0.15*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.2*7e0/10e0' 
		rate_values[2208] = 0.6*KRO2NO*0.2*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.4*7e0/10e0' 
		rate_values[2209] = 0.6*KRO2NO*0.4*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.5*7e0/10e0' 
		rate_values[2210] = 0.6*KRO2NO*0.5*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.7*7e0/10e0' 
		rate_values[2211] = 0.6*KRO2NO*0.7*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*1.0*7e0/10e0' 
		rate_values[2212] = 0.6*KRO2NO*1.0*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*1.0*7e0/10e0' 
		rate_values[2213] = 0.6*KRO2NO*1.0*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2214] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2215] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2216] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.1*3e0/10e0' 
		rate_values[2217] = 0.6*KRO2NO*0.1*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.15*3e0/10e0' 
		rate_values[2218] = 0.6*KRO2NO*0.15*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.2*3e0/10e0' 
		rate_values[2219] = 0.6*KRO2NO*0.2*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.4*3e0/10e0' 
		rate_values[2220] = 0.6*KRO2NO*0.4*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.5*3e0/10e0' 
		rate_values[2221] = 0.6*KRO2NO*0.5*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.7*3e0/10e0' 
		rate_values[2222] = 0.6*KRO2NO*0.7*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*1.0*3e0/10e0' 
		rate_values[2223] = 0.6*KRO2NO*1.0*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*1.0*3e0/10e0' 
		rate_values[2224] = 0.6*KRO2NO*1.0*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2225] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2226] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2227] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2228] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2229] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2230] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2231] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2232] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2233] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2234] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2235] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2236] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2237] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2238] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2239] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2240] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2241] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2242] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2243] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2244] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2245] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2246] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2247] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2248] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2249] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2250] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2251] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2252] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2253] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2254] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2255] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2256] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2257] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2258] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2259] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2260] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2261] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2262] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2263] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2264] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2265] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2266] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2267] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2268] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2269] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2270] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2271] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2272] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2273] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2274] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2275] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2276] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2277] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2278] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2279] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2280] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2281] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2282] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2283] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2284] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2285] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2286] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2287] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2288] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2289] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2290] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2291] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2292] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2293] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2294] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2295] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2296] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2297] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2298] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2299] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2300] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2301] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2302] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2303] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2304] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2305] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2306] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2307] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2308] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2309] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2310] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2311] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2312] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2313] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2314] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2315] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2316] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2317] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2318] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2319] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2320] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2321] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2322] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2323] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2324] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2325] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2326] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2327] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2328] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2329] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2330] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2331] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2332] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2333] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2334] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2335] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2336] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2337] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2338] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2339] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2340] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2341] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2342] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2343] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2344] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2345] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2346] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2347] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2348] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2349] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2350] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2351] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2352] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2353] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2354] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2355] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2356] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2357] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2358] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2359] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2360] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2361] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2362] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2363] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2364] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2365] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2366] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2367] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2368] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2369] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2370] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2371] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2372] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2373] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2374] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2375] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2376] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2377] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2378] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2379] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2380] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2381] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2382] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2383] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2384] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2385] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2386] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2387] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2388] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2389] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2390] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2391] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2392] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2393] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2394] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2395] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2396] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2397] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2398] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2399] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2400] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2401] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2402] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2403] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2404] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2405] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2406] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2407] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2408] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2409] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2410] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2411] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2412] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2413] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2414] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2415] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2416] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2417] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2418] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2419] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2420] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2421] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2422] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2423] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2424] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2425] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2426] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2427] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2428] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2429] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2430] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2431] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2432] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2433] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2434] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2435] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2436] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2437] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2438] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2439] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2440] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2441] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2442] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2443] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2444] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2445] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2446] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2447] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2448] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2449] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2450] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2451] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2452] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2453] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2454] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2455] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2456] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2457] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2458] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2459] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2460] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2461] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2462] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2463] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2464] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2465] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2466] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2467] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2468] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2469] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2470] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2471] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2472] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2473] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2474] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2475] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2476] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2477] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2478] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2479] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2480] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2481] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2482] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2483] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2484] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2485] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2486] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2487] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2488] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2489] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2490] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2491] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2492] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2493] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2494] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2495] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2496] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2497] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2498] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2499] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2500] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2501] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2502] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2503] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2504] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2505] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2506] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2507] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2508] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2509] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2510] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2511] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2512] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2513] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2514] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2515] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2516] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2517] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2518] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2519] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2520] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2521] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2522] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2523] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2524] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2525] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2526] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2527] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2528] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2529] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2530] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2531] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2532] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2533] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2534] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2535] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2536] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2537] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2538] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2539] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2540] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2541] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2542] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2543] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2544] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2545] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2546] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2547] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2548] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2549] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2550] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2551] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2552] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2553] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2554] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2555] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2556] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2557] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2558] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2559] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2560] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2561] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2562] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2563] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2564] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2565] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2566] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2567] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2568] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2569] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2570] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2571] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2572] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2573] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2574] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2575] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2576] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2577] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2578] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2579] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2580] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2581] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2582] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2583] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2584] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2585] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2586] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2587] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2588] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2589] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2590] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2591] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2592] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2593] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2594] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2595] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2596] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2597] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2598] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2599] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2600] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2601] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2602] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2603] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2604] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2605] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2606] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2607] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2608] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2609] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2610] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2611] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2612] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2613] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2614] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2615] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2616] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2617] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2618] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2619] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2620] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2621] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2622] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2623] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2624] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2625] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2626] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2627] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2628] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2629] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2630] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2631] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2632] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2633] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2634] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2635] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2636] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2637] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2638] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2639] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2640] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2641] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2642] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2643] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2644] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2645] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2646] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2647] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2648] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2649] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2650] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2651] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2652] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2653] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2654] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2655] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2656] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2657] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2658] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2659] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2660] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2661] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2662] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2663] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2664] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2665] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2666] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2667] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2668] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2669] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2670] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2671] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2672] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2673] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2674] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2675] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2676] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2677] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2678] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2679] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2680] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2681] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2682] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2683] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2684] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2685] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2686] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2687] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2688] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2689] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2690] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2691] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2692] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2693] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2694] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2695] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2696] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2697] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2698] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2699] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2700] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2701] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2702] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2703] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2704] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2705] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2706] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2707] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2708] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2709] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2710] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2711] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2712] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2713] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2714] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2715] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2716] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2717] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2718] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2719] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2720] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2721] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2722] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2723] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2724] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2725] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2726] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2727] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2728] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2729] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2730] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2731] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2732] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2733] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2734] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2735] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2736] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2737] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2738] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2739] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2740] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2741] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2742] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2743] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2744] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2745] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2746] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2747] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2748] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2749] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2750] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2751] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2752] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2753] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2754] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2755] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2756] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2757] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2758] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2759] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2760] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2761] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2762] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2763] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2764] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2765] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2766] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2767] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2768] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2769] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2770] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2771] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2772] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2773] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2774] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2775] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2776] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2777] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2778] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2779] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2780] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2781] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2782] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2783] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2784] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2785] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2786] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2787] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2788] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2789] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2790] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2791] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2792] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2793] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2794] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2795] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2796] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2797] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2798] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2799] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2800] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2801] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2802] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2803] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2804] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2805] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2806] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2807] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2808] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2809] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2810] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2811] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2812] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2813] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2814] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2815] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2816] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2817] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2818] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2819] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2820] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2821] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2822] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2823] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2824] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2825] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2826] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2827] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2828] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2829] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2830] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2831] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2832] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2833] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2834] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2835] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2836] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2837] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2838] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2839] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2840] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2841] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2842] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2843] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2844] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2845] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2846] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2847] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2848] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2849] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2850] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2851] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2852] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2853] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2854] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2855] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2856] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2857] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2858] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2859] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2860] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2861] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2862] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2863] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2864] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2865] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2866] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2867] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2868] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2869] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2870] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2871] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2872] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2873] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2874] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2875] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2876] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2877] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2878] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2879] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2880] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2881] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2882] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2883] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2884] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2885] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2886] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2887] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2888] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2889] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2890] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2891] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2892] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2893] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2894] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2895] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2896] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2897] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2898] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2899] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2900] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2901] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2902] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2903] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2904] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2905] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2906] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2907] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2908] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2909] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2910] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2911] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2912] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2913] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2914] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2915] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2916] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2917] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2918] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2919] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2920] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2921] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2922] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2923] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2924] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2925] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2926] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2927] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2928] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2929] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2930] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2931] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2932] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2933] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2934] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2935] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2936] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2937] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2938] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2939] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2940] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2941] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2942] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2943] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2944] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2945] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2946] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2947] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2948] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2949] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2950] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2951] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2952] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2953] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2954] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2955] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2956] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2957] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2958] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2959] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2960] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2961] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2962] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2963] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2964] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2965] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2966] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2967] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2968] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2969] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2970] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2971] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2972] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2973] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2974] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2975] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2976] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2977] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2978] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2979] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2980] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2981] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2982] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2983] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2984] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2985] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2986] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2987] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2988] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2989] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2990] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2991] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2992] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2993] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2994] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2995] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2996] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2997] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2998] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[2999] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3000] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3001] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3002] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3003] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3004] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3005] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3006] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3007] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3008] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3009] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3010] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3011] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3012] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3013] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3014] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3015] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3016] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3017] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3018] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3019] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3020] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3021] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3022] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3023] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3024] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3025] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3026] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3027] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3028] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3029] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3030] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3031] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3032] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3033] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3034] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3035] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3036] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3037] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3038] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3039] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3040] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3041] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3042] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3043] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3044] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3045] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3046] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3047] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3048] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3049] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3050] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3051] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3052] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3053] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3054] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3055] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3056] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3057] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3058] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3059] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3060] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3061] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3062] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3063] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3064] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3065] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3066] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3067] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3068] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3069] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3070] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3071] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3072] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3073] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3074] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3075] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3076] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3077] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3078] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3079] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3080] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3081] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3082] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3083] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3084] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3085] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3086] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3087] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3088] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3089] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3090] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3091] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3092] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3093] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3094] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3095] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3096] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3097] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3098] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3099] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3100] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3101] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3102] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3103] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3104] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3105] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3106] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3107] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3108] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3109] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3110] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3111] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3112] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3113] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3114] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3115] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3116] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3117] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3118] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3119] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3120] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3121] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3122] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3123] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3124] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3125] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3126] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3127] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3128] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3129] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3130] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3131] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3132] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3133] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3134] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3135] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3136] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3137] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3138] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3139] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3140] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3141] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3142] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3143] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3144] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3145] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3146] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3147] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3148] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3149] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3150] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3151] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3152] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3153] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3154] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3155] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3156] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3157] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3158] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3159] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3160] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3161] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3162] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3163] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3164] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3165] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3166] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3167] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3168] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3169] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3170] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3171] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3172] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3173] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3174] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3175] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3176] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3177] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3178] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3179] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3180] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3181] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3182] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3183] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3184] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3185] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3186] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3187] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3188] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3189] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3190] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3191] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3192] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3193] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3194] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3195] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3196] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3197] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3198] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3199] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3200] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3201] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3202] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3203] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3204] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3205] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3206] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3207] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3208] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3209] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3210] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3211] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3212] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3213] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3214] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3215] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3216] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3217] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3218] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3219] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3220] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3221] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3222] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3223] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3224] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3225] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3226] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3227] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3228] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3229] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3230] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3231] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3232] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3233] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3234] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3235] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3236] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3237] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3238] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3239] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3240] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3241] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3242] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3243] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3244] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3245] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3246] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3247] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3248] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3249] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3250] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3251] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3252] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3253] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3254] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3255] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3256] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3257] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3258] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3259] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3260] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3261] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3262] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3263] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3264] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3265] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3266] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3267] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3268] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[3269] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-12*0.2*RO2' 
		rate_values[3270] = 1E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-12*0.6*RO2' 
		rate_values[3271] = 1E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-12*0.2*RO2' 
		rate_values[3272] = 1E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[3273] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.6*RO2' 
		rate_values[3274] = 5E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[3275] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[3276] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.6*RO2' 
		rate_values[3277] = 5E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[3278] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[3279] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.6*RO2' 
		rate_values[3280] = 5E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[3281] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7E-12*0.4*RO2' 
		rate_values[3282] = 7E-12*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7E-12*0.4*RO2' 
		rate_values[3283] = 7E-12*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7E-12*0.2*RO2' 
		rate_values[3284] = 7E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8E-12*0.2*RO2' 
		rate_values[3285] = 8E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8E-12*0.4*RO2' 
		rate_values[3286] = 8E-12*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8E-12*0.4*RO2' 
		rate_values[3287] = 8E-12*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9E-12*0.4*RO2' 
		rate_values[3288] = 9E-12*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9E-12*0.2*RO2' 
		rate_values[3289] = 9E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9E-12*0.4*RO2' 
		rate_values[3290] = 9E-12*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.3*RO2' 
		rate_values[3291] = 1E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.2*RO2' 
		rate_values[3292] = 1E-11*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.5*RO2' 
		rate_values[3293] = 1E-11*0.5*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.5*RO2' 
		rate_values[3294] = 1E-11*0.5*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.2*RO2' 
		rate_values[3295] = 1E-11*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.3*RO2' 
		rate_values[3296] = 1E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.4*RO2' 
		rate_values[3297] = 1E-11*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.2*RO2' 
		rate_values[3298] = 1E-11*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.4*RO2' 
		rate_values[3299] = 1E-11*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.5*RO2' 
		rate_values[3300] = 1E-11*0.5*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.0*RO2' 
		rate_values[3301] = 1E-11*0.0*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.5*RO2' 
		rate_values[3302] = 1E-11*0.5*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.1' 
		rate_values[3303] = 9.20e-14*RO2*0.1
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*0.7*RO2*0.9' 
		rate_values[3304] = 9.20e-14*0.7*RO2*0.9
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*0.3*RO2*0.9' 
		rate_values[3305] = 9.20e-14*0.3*RO2*0.9
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2e17*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[3306] = 2e17*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6e16*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[3307] = 6e16*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6e16*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[3308] = 6e16*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e16*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[3309] = 3e16*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[3310] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[3311] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[3312] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[3313] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[3314] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0*KRO2NO' 
		rate_values[3315] = 1.0*KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.9*KRO2NO' 
		rate_values[3316] = 0.9*KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.8*KRO2NO' 
		rate_values[3317] = 0.8*KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.7*KRO2NO' 
		rate_values[3318] = 0.7*KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.5*KRO2NO' 
		rate_values[3319] = 0.5*KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[3320] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.1' 
		rate_values[3321] = 0.4*KRO2NO*0.1
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.2' 
		rate_values[3322] = 0.4*KRO2NO*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.3' 
		rate_values[3323] = 0.4*KRO2NO*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.5' 
		rate_values[3324] = 0.4*KRO2NO*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*1.0' 
		rate_values[3325] = 0.4*KRO2NO*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[3326] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.1*7e0/10e0' 
		rate_values[3327] = 0.6*KRO2NO*0.1*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.2*7e0/10e0' 
		rate_values[3328] = 0.6*KRO2NO*0.2*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.3*7e0/10e0' 
		rate_values[3329] = 0.6*KRO2NO*0.3*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.5*7e0/10e0' 
		rate_values[3330] = 0.6*KRO2NO*0.5*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*1.0*7e0/10e0' 
		rate_values[3331] = 0.6*KRO2NO*1.0*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[3332] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.1*3e0/10e0' 
		rate_values[3333] = 0.6*KRO2NO*0.1*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.2*3e0/10e0' 
		rate_values[3334] = 0.6*KRO2NO*0.2*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.3*3e0/10e0' 
		rate_values[3335] = 0.6*KRO2NO*0.3*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.5*3e0/10e0' 
		rate_values[3336] = 0.6*KRO2NO*0.5*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*1.0*3e0/10e0' 
		rate_values[3337] = 0.6*KRO2NO*1.0*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[3338] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[3339] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[3340] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[3341] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[3342] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[3343] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3344] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3345] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3346] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3347] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3348] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3349] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3350] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3351] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3352] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3353] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3354] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3355] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3356] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3357] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3358] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3359] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3360] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3361] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3362] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3363] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3364] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3365] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3366] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3367] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3368] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3369] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3370] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3371] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3372] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3373] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3374] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3375] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3376] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3377] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3378] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3379] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3380] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3381] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3382] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3383] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3384] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3385] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3386] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3387] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3388] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3389] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3390] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3391] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3392] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3393] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3394] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3395] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3396] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3397] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3398] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3399] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3400] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3401] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3402] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3403] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3404] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3405] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3406] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3407] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3408] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3409] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3410] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3411] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3412] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3413] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3414] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3415] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3416] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3417] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3418] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3419] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3420] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3421] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3422] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3423] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3424] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3425] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3426] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3427] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3428] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3429] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3430] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3431] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3432] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3433] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3434] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3435] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3436] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3437] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3438] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3439] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3440] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3441] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3442] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3443] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3444] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3445] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3446] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3447] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3448] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3449] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3450] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3451] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3452] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3453] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3454] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3455] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3456] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3457] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3458] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3459] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3460] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3461] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3462] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3463] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3464] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3465] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3466] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3467] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3468] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3469] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3470] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3471] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3472] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3473] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3474] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3475] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3476] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3477] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3478] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3479] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3480] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3481] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3482] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3483] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3484] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3485] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3486] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3487] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3488] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3489] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3490] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3491] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3492] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3493] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3494] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3495] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3496] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3497] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3498] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3499] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3500] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3501] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3502] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3503] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3504] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3505] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3506] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3507] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3508] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3509] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3510] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3511] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3512] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3513] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3514] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3515] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3516] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3517] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3518] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3519] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3520] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3521] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3522] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3523] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3524] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3525] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3526] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3527] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3528] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3529] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3530] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[3531] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3532] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3533] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3534] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3535] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3536] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3537] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3538] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3539] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3540] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3541] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3542] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3543] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3544] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3545] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3546] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3547] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3548] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3549] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3550] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3551] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3552] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3553] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3554] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3555] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3556] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3557] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3558] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3559] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3560] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3561] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3562] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3563] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3564] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3565] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3566] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3567] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3568] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3569] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3570] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3571] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3572] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3573] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3574] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3575] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3576] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3577] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3578] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3579] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3580] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3581] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3582] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3583] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3584] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3585] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3586] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3587] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3588] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3589] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3590] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3591] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3592] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3593] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3594] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3595] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3596] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3597] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3598] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3599] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3600] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3601] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3602] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3603] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3604] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3605] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3606] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3607] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3608] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3609] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3610] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3611] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3612] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3613] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3614] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3615] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3616] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3617] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3618] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3619] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3620] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3621] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3622] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3623] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3624] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[3625] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3626] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3627] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3628] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3629] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3630] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3631] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3632] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3633] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3634] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3635] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3636] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3637] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3638] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3639] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3640] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3641] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3642] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3643] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3644] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3645] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3646] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3647] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3648] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3649] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3650] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3651] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3652] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3653] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3654] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3655] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3656] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3657] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3658] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3659] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3660] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3661] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3662] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3663] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3664] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3665] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3666] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3667] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3668] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3669] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3670] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3671] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3672] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3673] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3674] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3675] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3676] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3677] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3678] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3679] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3680] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3681] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3682] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3683] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3684] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3685] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3686] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3687] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3688] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3689] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3690] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3691] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3692] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3693] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3694] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3695] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3696] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3697] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3698] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3699] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3700] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3701] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3702] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3703] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3704] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3705] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3706] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3707] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3708] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3709] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3710] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3711] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3712] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3713] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3714] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3715] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3716] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3717] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3718] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[3719] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3720] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3721] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3722] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3723] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3724] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3725] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3726] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3727] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3728] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3729] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3730] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3731] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3732] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3733] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3734] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3735] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3736] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3737] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3738] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3739] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3740] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3741] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3742] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3743] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3744] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3745] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3746] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3747] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3748] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3749] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3750] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3751] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3752] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3753] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3754] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3755] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3756] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3757] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3758] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3759] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3760] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3761] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3762] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3763] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3764] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3765] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3766] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3767] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3768] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3769] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3770] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3771] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3772] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3773] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3774] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3775] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3776] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3777] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3778] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3779] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3780] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3781] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3782] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3783] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3784] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3785] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3786] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3787] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3788] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3789] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3790] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3791] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3792] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3793] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3794] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3795] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3796] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3797] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3798] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3799] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3800] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3801] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3802] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3803] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3804] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3805] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3806] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3807] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3808] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3809] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3810] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3811] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3812] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3813] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3814] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3815] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3816] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3817] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3818] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3819] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3820] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3821] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3822] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3823] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3824] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3825] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3826] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3827] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3828] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3829] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3830] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3831] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3832] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3833] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3834] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3835] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3836] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3837] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3838] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3839] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3840] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3841] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3842] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3843] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3844] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3845] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3846] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3847] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3848] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3849] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3850] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3851] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3852] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3853] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3854] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3855] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3856] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3857] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3858] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3859] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3860] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3861] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3862] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3863] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3864] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3865] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3866] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3867] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3868] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3869] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3870] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3871] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3872] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3873] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3874] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3875] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3876] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3877] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3878] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3879] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3880] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3881] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3882] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3883] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3884] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3885] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3886] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3887] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3888] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3889] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3890] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3891] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3892] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3893] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3894] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3895] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3896] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3897] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3898] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3899] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3900] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3901] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3902] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3903] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3904] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3905] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3906] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[3907] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[3908] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.6*RO2' 
		rate_values[3909] = 5E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[3910] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[3911] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.6*RO2' 
		rate_values[3912] = 5E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[3913] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8E-12*0.2*RO2' 
		rate_values[3914] = 8E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8E-12*0.6*RO2' 
		rate_values[3915] = 8E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8E-12*0.2*RO2' 
		rate_values[3916] = 8E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.3*RO2' 
		rate_values[3917] = 1E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.4*RO2' 
		rate_values[3918] = 1E-11*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.3*RO2' 
		rate_values[3919] = 1E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.4*RO2' 
		rate_values[3920] = 1E-11*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.2*RO2' 
		rate_values[3921] = 1E-11*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.4*RO2' 
		rate_values[3922] = 1E-11*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.5*RO2' 
		rate_values[3923] = 1E-11*0.5*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.5*RO2' 
		rate_values[3924] = 1E-11*0.5*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3925] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3926] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.7356e-11' 
		rate_values[3927] = 1.7356e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4859e-11' 
		rate_values[3928] = 2.4859e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8872e-11' 
		rate_values[3929] = 2.8872e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8872e-11' 
		rate_values[3930] = 2.8872e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3931] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.1952e-11' 
		rate_values[3932] = 1.1952e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1439e-11' 
		rate_values[3933] = 2.1439e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.7863e-11' 
		rate_values[3934] = 2.7863e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8872e-11' 
		rate_values[3935] = 2.8872e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3936] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.7356e-11' 
		rate_values[3937] = 1.7356e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4859e-11' 
		rate_values[3938] = 2.4859e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3939] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.1952e-11' 
		rate_values[3940] = 1.1952e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1439e-11' 
		rate_values[3941] = 2.1439e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3942] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3943] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3944] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3945] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3946] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3947] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3948] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3949] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3950] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3951] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3952] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3953] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.7356e-11' 
		rate_values[3954] = 1.7356e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.1952e-11' 
		rate_values[3955] = 1.1952e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3956] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3957] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3958] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3959] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3960] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3961] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3962] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3963] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3964] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3965] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3966] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3967] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3968] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3969] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3970] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3971] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0831e-12' 
		rate_values[3972] = 4.0831e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[3973] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[3974] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6457e-11' 
		rate_values[3975] = 2.6457e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.7895e-11' 
		rate_values[3976] = 3.7895e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4012e-11' 
		rate_values[3977] = 4.4012e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4012e-11' 
		rate_values[3978] = 4.4012e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[3979] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8219e-11' 
		rate_values[3980] = 1.8219e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.268e-11' 
		rate_values[3981] = 3.268e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.2475e-11' 
		rate_values[3982] = 4.2475e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4012e-11' 
		rate_values[3983] = 4.4012e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[3984] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6457e-11' 
		rate_values[3985] = 2.6457e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.7895e-11' 
		rate_values[3986] = 3.7895e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[3987] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8219e-11' 
		rate_values[3988] = 1.8219e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.268e-11' 
		rate_values[3989] = 3.268e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[3990] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[3991] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[3992] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[3993] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[3994] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[3995] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[3996] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[3997] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[3998] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[3999] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4000] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4001] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6457e-11' 
		rate_values[4002] = 2.6457e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8219e-11' 
		rate_values[4003] = 1.8219e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4004] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4005] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4006] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4007] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4008] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4009] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4010] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4011] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4012] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4013] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4014] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4015] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4016] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4017] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4018] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4019] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2242e-12' 
		rate_values[4020] = 6.2242e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4021] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4022] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[4023] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8506e-11' 
		rate_values[4024] = 3.8506e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[4025] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[4026] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4027] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[4028] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3207e-11' 
		rate_values[4029] = 3.3207e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.3159e-11' 
		rate_values[4030] = 4.3159e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[4031] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4032] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[4033] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8506e-11' 
		rate_values[4034] = 3.8506e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4035] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[4036] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3207e-11' 
		rate_values[4037] = 3.3207e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4038] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4039] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4040] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4041] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4042] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4043] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4044] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4045] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4046] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4047] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4048] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4049] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[4050] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[4051] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4052] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4053] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4054] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4055] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4056] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4057] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4058] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4059] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4060] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4061] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4062] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4063] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4064] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4065] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4066] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4067] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4068] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4069] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4070] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0092e-11' 
		rate_values[4071] = 1.0092e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4456e-11' 
		rate_values[4072] = 1.4456e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6789e-11' 
		rate_values[4073] = 1.6789e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6789e-11' 
		rate_values[4074] = 1.6789e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4075] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.9499e-12' 
		rate_values[4076] = 6.9499e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2466e-11' 
		rate_values[4077] = 1.2466e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6202e-11' 
		rate_values[4078] = 1.6202e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6789e-11' 
		rate_values[4079] = 1.6789e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4080] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0092e-11' 
		rate_values[4081] = 1.0092e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4456e-11' 
		rate_values[4082] = 1.4456e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4083] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.9499e-12' 
		rate_values[4084] = 6.9499e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2466e-11' 
		rate_values[4085] = 1.2466e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4086] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4087] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4088] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4089] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4090] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4091] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4092] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4093] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4094] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4095] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4096] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4097] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0092e-11' 
		rate_values[4098] = 1.0092e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.9499e-12' 
		rate_values[4099] = 6.9499e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4100] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4101] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4102] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4103] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4104] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4105] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4106] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4107] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4108] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4109] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4110] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4111] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4112] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4113] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4114] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4115] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4116] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4117] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4118] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.2374e-11' 
		rate_values[4119] = 2.2374e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.2047e-11' 
		rate_values[4120] = 3.2047e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.722e-11' 
		rate_values[4121] = 3.722e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.722e-11' 
		rate_values[4122] = 3.722e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4123] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.5408e-11' 
		rate_values[4124] = 1.5408e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.7637e-11' 
		rate_values[4125] = 2.7637e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.592e-11' 
		rate_values[4126] = 3.592e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.722e-11' 
		rate_values[4127] = 3.722e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4128] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.2374e-11' 
		rate_values[4129] = 2.2374e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.2047e-11' 
		rate_values[4130] = 3.2047e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4131] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.5408e-11' 
		rate_values[4132] = 1.5408e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.7637e-11' 
		rate_values[4133] = 2.7637e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4134] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4135] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4136] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4137] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4138] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4139] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4140] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4141] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4142] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4143] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4144] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4145] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.2374e-11' 
		rate_values[4146] = 2.2374e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.5408e-11' 
		rate_values[4147] = 1.5408e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4148] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4149] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4150] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4151] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4152] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4153] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4154] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4155] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4156] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4157] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4158] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4159] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4160] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4161] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4162] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4163] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2637e-12' 
		rate_values[4164] = 5.2637e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4165] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4166] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[4167] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8506e-11' 
		rate_values[4168] = 3.8506e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[4169] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[4170] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4171] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[4172] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3207e-11' 
		rate_values[4173] = 3.3207e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.3159e-11' 
		rate_values[4174] = 4.3159e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[4175] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4176] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[4177] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8506e-11' 
		rate_values[4178] = 3.8506e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4179] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[4180] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3207e-11' 
		rate_values[4181] = 3.3207e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4182] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4183] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4184] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4185] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4186] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4187] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4188] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4189] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4190] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4191] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4192] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4193] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[4194] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[4195] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4196] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4197] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4198] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4199] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4200] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4201] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4202] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4203] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4204] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4205] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4206] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4207] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4208] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4209] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4210] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4211] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4212] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4213] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4214] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4215] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4216] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4217] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4218] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4219] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4220] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4221] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[4222] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4223] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4224] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4225] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4226] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4227] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4228] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4229] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4230] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4231] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4232] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4233] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4234] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4235] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4236] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4237] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4238] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4239] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4240] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4241] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4242] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4243] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4244] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4245] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4246] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4247] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4248] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4249] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4250] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4251] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4252] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4253] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4254] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4255] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4256] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4257] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4258] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4259] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4260] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4261] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4262] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2317e-11' 
		rate_values[4263] = 1.2317e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.7641e-11' 
		rate_values[4264] = 1.7641e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0489e-11' 
		rate_values[4265] = 2.0489e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0489e-11' 
		rate_values[4266] = 2.0489e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4267] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.4816e-12' 
		rate_values[4268] = 8.4816e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.5214e-11' 
		rate_values[4269] = 1.5214e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.9773e-11' 
		rate_values[4270] = 1.9773e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0489e-11' 
		rate_values[4271] = 2.0489e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4272] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2317e-11' 
		rate_values[4273] = 1.2317e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.7641e-11' 
		rate_values[4274] = 1.7641e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4275] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.4816e-12' 
		rate_values[4276] = 8.4816e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.5214e-11' 
		rate_values[4277] = 1.5214e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4278] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4279] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4280] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4281] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4282] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4283] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4284] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4285] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4286] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4287] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4288] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4289] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2317e-11' 
		rate_values[4290] = 1.2317e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.4816e-12' 
		rate_values[4291] = 8.4816e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4292] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4293] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4294] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4295] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4296] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4297] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4298] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4299] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4300] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4301] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4302] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4303] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4304] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4305] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4306] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4307] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8976e-12' 
		rate_values[4308] = 2.8976e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4309] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4310] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8737e-11' 
		rate_values[4311] = 1.8737e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6837e-11' 
		rate_values[4312] = 2.6837e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.1169e-11' 
		rate_values[4313] = 3.1169e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.1169e-11' 
		rate_values[4314] = 3.1169e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4315] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2903e-11' 
		rate_values[4316] = 1.2903e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3144e-11' 
		rate_values[4317] = 2.3144e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.0081e-11' 
		rate_values[4318] = 3.0081e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.1169e-11' 
		rate_values[4319] = 3.1169e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4320] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8737e-11' 
		rate_values[4321] = 1.8737e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6837e-11' 
		rate_values[4322] = 2.6837e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4323] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2903e-11' 
		rate_values[4324] = 1.2903e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3144e-11' 
		rate_values[4325] = 2.3144e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4326] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4327] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4328] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4329] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4330] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4331] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4332] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4333] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4334] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4335] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4336] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4337] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8737e-11' 
		rate_values[4338] = 1.8737e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2903e-11' 
		rate_values[4339] = 1.2903e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4340] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4341] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4342] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4343] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4344] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4345] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4346] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4347] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4348] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4349] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4350] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4351] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4352] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4353] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4354] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4355] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.408e-12' 
		rate_values[4356] = 4.408e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4357] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4358] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3462e-11' 
		rate_values[4359] = 2.3462e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3605e-11' 
		rate_values[4360] = 3.3605e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.9029e-11' 
		rate_values[4361] = 3.9029e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.9029e-11' 
		rate_values[4362] = 3.9029e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4363] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6156e-11' 
		rate_values[4364] = 1.6156e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.898e-11' 
		rate_values[4365] = 2.898e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.7666e-11' 
		rate_values[4366] = 3.7666e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.9029e-11' 
		rate_values[4367] = 3.9029e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4368] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3462e-11' 
		rate_values[4369] = 2.3462e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3605e-11' 
		rate_values[4370] = 3.3605e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4371] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6156e-11' 
		rate_values[4372] = 1.6156e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.898e-11' 
		rate_values[4373] = 2.898e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4374] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4375] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4376] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4377] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4378] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4379] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4380] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4381] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4382] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4383] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4384] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4385] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3462e-11' 
		rate_values[4386] = 2.3462e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6156e-11' 
		rate_values[4387] = 1.6156e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4388] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4389] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4390] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4391] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4392] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4393] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4394] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4395] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4396] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4397] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4398] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4399] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4400] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4401] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4402] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4403] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5195e-12' 
		rate_values[4404] = 5.5195e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4405] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4406] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[4407] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8506e-11' 
		rate_values[4408] = 3.8506e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[4409] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[4410] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4411] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[4412] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3207e-11' 
		rate_values[4413] = 3.3207e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.3159e-11' 
		rate_values[4414] = 4.3159e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[4415] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4416] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[4417] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8506e-11' 
		rate_values[4418] = 3.8506e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4419] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[4420] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3207e-11' 
		rate_values[4421] = 3.3207e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4422] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4423] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4424] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4425] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4426] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4427] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4428] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4429] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4430] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4431] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4432] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4433] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[4434] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[4435] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4436] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4437] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4438] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4439] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4440] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4441] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4442] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4443] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4444] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4445] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4446] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4447] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4448] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4449] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4450] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4451] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[4452] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4453] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4454] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4455] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4456] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4457] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4458] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4459] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4460] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4461] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[4462] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4463] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4464] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4465] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4466] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4467] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4468] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4469] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4470] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4471] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4472] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4473] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4474] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4475] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4476] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4477] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4478] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4479] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4480] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4481] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4482] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4483] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4484] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4485] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4486] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4487] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4488] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4489] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4490] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4491] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4492] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4493] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4494] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4495] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4496] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4497] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4498] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4499] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4500] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4501] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4502] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4503] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4504] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4505] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4506] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4507] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4508] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4509] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[4510] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4511] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4512] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4513] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4514] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4515] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4516] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4517] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4518] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4519] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4520] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4521] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4522] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4523] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4524] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4525] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4526] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4527] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4528] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4529] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4530] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4531] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4532] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4533] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4534] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4535] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4536] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4537] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4538] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4539] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4540] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4541] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4542] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4543] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4544] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4545] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4546] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4547] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4548] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4549] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4550] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0092e-11' 
		rate_values[4551] = 1.0092e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4456e-11' 
		rate_values[4552] = 1.4456e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6789e-11' 
		rate_values[4553] = 1.6789e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6789e-11' 
		rate_values[4554] = 1.6789e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4555] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.9499e-12' 
		rate_values[4556] = 6.9499e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2466e-11' 
		rate_values[4557] = 1.2466e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6202e-11' 
		rate_values[4558] = 1.6202e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6789e-11' 
		rate_values[4559] = 1.6789e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4560] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0092e-11' 
		rate_values[4561] = 1.0092e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4456e-11' 
		rate_values[4562] = 1.4456e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4563] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.9499e-12' 
		rate_values[4564] = 6.9499e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2466e-11' 
		rate_values[4565] = 1.2466e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4566] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4567] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4568] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4569] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4570] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4571] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4572] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4573] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4574] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4575] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4576] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4577] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0092e-11' 
		rate_values[4578] = 1.0092e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.9499e-12' 
		rate_values[4579] = 6.9499e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4580] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4581] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4582] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4583] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4584] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4585] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4586] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4587] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4588] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4589] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4590] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4591] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4592] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4593] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4594] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4595] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4596] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4597] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4598] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4599] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4600] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4601] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4602] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4603] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4604] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4605] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[4606] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4607] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4608] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4609] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4610] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4611] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4612] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4613] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4614] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4615] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4616] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4617] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4618] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4619] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4620] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4621] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4622] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4623] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4624] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4625] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4626] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4627] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4628] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4629] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4630] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4631] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4632] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4633] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4634] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4635] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4636] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4637] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4638] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4639] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4640] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4641] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4642] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4643] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4644] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4645] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4646] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.5741e-11' 
		rate_values[4647] = 2.5741e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.6869e-11' 
		rate_values[4648] = 3.6869e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.282e-11' 
		rate_values[4649] = 4.282e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.282e-11' 
		rate_values[4650] = 4.282e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4651] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.7726e-11' 
		rate_values[4652] = 1.7726e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.1796e-11' 
		rate_values[4653] = 3.1796e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.1325e-11' 
		rate_values[4654] = 4.1325e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.282e-11' 
		rate_values[4655] = 4.282e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4656] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.5741e-11' 
		rate_values[4657] = 2.5741e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.6869e-11' 
		rate_values[4658] = 3.6869e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4659] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.7726e-11' 
		rate_values[4660] = 1.7726e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.1796e-11' 
		rate_values[4661] = 3.1796e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4662] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4663] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4664] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4665] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4666] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4667] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4668] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4669] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4670] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4671] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4672] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4673] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.5741e-11' 
		rate_values[4674] = 2.5741e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.7726e-11' 
		rate_values[4675] = 1.7726e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4676] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4677] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4678] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4679] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4680] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4681] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4682] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4683] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4684] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4685] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4686] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4687] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4688] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4689] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4690] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4691] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0557e-12' 
		rate_values[4692] = 6.0557e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4693] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4694] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0092e-11' 
		rate_values[4695] = 1.0092e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4456e-11' 
		rate_values[4696] = 1.4456e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6789e-11' 
		rate_values[4697] = 1.6789e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6789e-11' 
		rate_values[4698] = 1.6789e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4699] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.9499e-12' 
		rate_values[4700] = 6.9499e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2466e-11' 
		rate_values[4701] = 1.2466e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6202e-11' 
		rate_values[4702] = 1.6202e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6789e-11' 
		rate_values[4703] = 1.6789e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4704] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0092e-11' 
		rate_values[4705] = 1.0092e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4456e-11' 
		rate_values[4706] = 1.4456e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4707] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.9499e-12' 
		rate_values[4708] = 6.9499e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2466e-11' 
		rate_values[4709] = 1.2466e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4710] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4711] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4712] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4713] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4714] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4715] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4716] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4717] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4718] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4719] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4720] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4721] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0092e-11' 
		rate_values[4722] = 1.0092e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.9499e-12' 
		rate_values[4723] = 6.9499e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4724] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4725] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4726] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4727] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4728] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4729] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4730] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4731] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4732] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4733] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4734] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4735] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4736] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4737] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4738] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4739] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3743e-12' 
		rate_values[4740] = 2.3743e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4741] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4742] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4743] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4744] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4745] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4746] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4747] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4748] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4749] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[4750] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4751] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4752] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4753] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4754] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4755] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4756] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4757] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4758] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4759] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4760] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4761] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4762] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4763] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4764] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4765] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4766] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4767] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4768] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4769] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4770] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4771] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4772] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4773] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4774] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4775] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4776] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4777] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4778] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4779] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4780] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4781] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4782] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4783] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4784] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4785] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4786] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4787] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4788] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4789] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4790] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4791] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4792] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4793] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4794] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4795] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4796] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4797] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[4798] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4799] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4800] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4801] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4802] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4803] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4804] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4805] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4806] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4807] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4808] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4809] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4810] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4811] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4812] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4813] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4814] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4815] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4816] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4817] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4818] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4819] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4820] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4821] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4822] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4823] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4824] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4825] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4826] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4827] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4828] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4829] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4830] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4831] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4832] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4833] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4834] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4835] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4836] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4837] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4838] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4839] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4840] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4841] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4842] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4843] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4844] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4845] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[4846] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4847] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4848] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4849] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4850] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4851] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4852] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4853] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4854] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4855] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4856] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4857] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4858] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4859] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4860] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4861] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4862] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4863] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4864] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4865] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4866] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4867] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4868] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4869] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4870] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4871] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4872] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4873] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4874] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4875] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4876] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4877] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4878] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4879] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4880] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4881] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4882] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4883] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4884] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4885] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4886] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4887] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4888] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4889] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4890] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4891] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4892] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4893] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[4894] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4895] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4896] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4897] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4898] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4899] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4900] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4901] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4902] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4903] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4904] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4905] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4906] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4907] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4908] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4909] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4910] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4911] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4912] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4913] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4914] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4915] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4916] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4917] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4918] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4919] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4920] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4921] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4922] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4923] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4924] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4925] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4926] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4927] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4928] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4929] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4930] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4931] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4932] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4933] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4934] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4935] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4936] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4937] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4938] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4939] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4940] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4941] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[4942] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[4943] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4944] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4945] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[4946] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4947] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4948] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[4949] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4950] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4951] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4952] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4953] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4954] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4955] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4956] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4957] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4958] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4959] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4960] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4961] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[4962] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[4963] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4964] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4965] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4966] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4967] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4968] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4969] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4970] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4971] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4972] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4973] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4974] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4975] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4976] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4977] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4978] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4979] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[4980] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[4981] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[4982] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.18e-11' 
		rate_values[4983] = 1.18e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6901e-11' 
		rate_values[4984] = 1.6901e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.9629e-11' 
		rate_values[4985] = 1.9629e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.9629e-11' 
		rate_values[4986] = 1.9629e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[4987] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.1258e-12' 
		rate_values[4988] = 8.1258e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4576e-11' 
		rate_values[4989] = 1.4576e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8944e-11' 
		rate_values[4990] = 1.8944e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.9629e-11' 
		rate_values[4991] = 1.9629e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[4992] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.18e-11' 
		rate_values[4993] = 1.18e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6901e-11' 
		rate_values[4994] = 1.6901e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[4995] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.1258e-12' 
		rate_values[4996] = 8.1258e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4576e-11' 
		rate_values[4997] = 1.4576e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[4998] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[4999] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5000] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5001] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5002] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5003] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5004] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5005] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5006] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5007] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5008] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5009] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.18e-11' 
		rate_values[5010] = 1.18e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.1258e-12' 
		rate_values[5011] = 8.1258e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5012] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5013] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5014] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5015] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5016] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5017] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5018] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5019] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5020] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5021] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5022] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5023] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5024] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5025] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5026] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5027] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.776e-12' 
		rate_values[5028] = 2.776e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5029] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5030] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5031] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5032] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5033] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5034] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5035] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5036] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5037] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[5038] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5039] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5040] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5041] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5042] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5043] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5044] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5045] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5046] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5047] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5048] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5049] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5050] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5051] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5052] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5053] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5054] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5055] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5056] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5057] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5058] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5059] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5060] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5061] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5062] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5063] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5064] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5065] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5066] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5067] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5068] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5069] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5070] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5071] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5072] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5073] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5074] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5075] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5076] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5077] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5078] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5079] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5080] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5081] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5082] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5083] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5084] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5085] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[5086] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5087] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5088] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5089] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5090] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5091] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5092] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5093] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5094] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5095] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5096] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5097] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5098] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5099] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5100] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5101] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5102] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5103] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5104] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5105] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5106] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5107] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5108] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5109] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5110] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5111] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5112] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5113] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5114] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5115] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5116] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5117] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5118] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5119] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5120] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5121] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5122] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5123] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5124] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5125] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5126] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5127] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5128] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5129] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5130] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5131] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5132] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5133] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[5134] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5135] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5136] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5137] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5138] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5139] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5140] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5141] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5142] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5143] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5144] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5145] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5146] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5147] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5148] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5149] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5150] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5151] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5152] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5153] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5154] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5155] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5156] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5157] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5158] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5159] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5160] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5161] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5162] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5163] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5164] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5165] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5166] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5167] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5168] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5169] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5170] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5171] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5172] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5173] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5174] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5175] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5176] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5177] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5178] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5179] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5180] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5181] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[5182] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5183] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5184] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5185] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5186] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5187] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5188] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5189] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5190] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5191] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5192] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5193] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5194] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5195] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5196] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5197] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5198] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5199] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5200] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5201] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5202] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5203] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5204] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5205] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5206] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5207] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5208] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5209] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5210] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5211] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5212] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5213] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5214] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5215] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5216] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5217] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5218] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5219] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5220] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5221] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5222] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5223] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5224] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5225] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5226] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5227] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5228] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5229] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[5230] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5231] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5232] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5233] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5234] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5235] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5236] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5237] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5238] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5239] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5240] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5241] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5242] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5243] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5244] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5245] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5246] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5247] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5248] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5249] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5250] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5251] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5252] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5253] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5254] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5255] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5256] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5257] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5258] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5259] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5260] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5261] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5262] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5263] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5264] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5265] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5266] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5267] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5268] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5269] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5270] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5271] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5272] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5273] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5274] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5275] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5276] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5277] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[5278] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5279] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5280] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5281] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5282] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5283] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5284] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5285] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5286] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5287] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5288] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5289] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5290] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5291] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5292] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5293] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5294] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5295] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5296] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5297] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5298] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5299] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5300] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5301] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5302] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5303] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5304] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5305] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5306] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5307] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5308] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5309] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5310] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5311] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5312] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5313] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5314] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5315] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5316] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5317] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5318] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[5319] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8506e-11' 
		rate_values[5320] = 3.8506e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[5321] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[5322] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5323] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[5324] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3207e-11' 
		rate_values[5325] = 3.3207e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.3159e-11' 
		rate_values[5326] = 4.3159e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[5327] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5328] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[5329] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8506e-11' 
		rate_values[5330] = 3.8506e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5331] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[5332] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3207e-11' 
		rate_values[5333] = 3.3207e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5334] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5335] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5336] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5337] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5338] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5339] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5340] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5341] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5342] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5343] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5344] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5345] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[5346] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[5347] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5348] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5349] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5350] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5351] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5352] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5353] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5354] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5355] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5356] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5357] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5358] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5359] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5360] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5361] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5362] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5363] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5364] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5365] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5366] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[5367] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8506e-11' 
		rate_values[5368] = 3.8506e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[5369] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[5370] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5371] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[5372] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3207e-11' 
		rate_values[5373] = 3.3207e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.3159e-11' 
		rate_values[5374] = 4.3159e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[5375] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5376] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[5377] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8506e-11' 
		rate_values[5378] = 3.8506e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5379] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[5380] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3207e-11' 
		rate_values[5381] = 3.3207e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5382] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5383] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5384] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5385] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5386] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5387] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5388] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5389] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5390] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5391] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5392] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5393] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[5394] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[5395] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5396] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5397] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5398] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5399] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5400] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5401] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5402] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5403] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5404] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5405] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5406] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5407] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5408] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5409] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5410] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5411] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5412] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5413] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5414] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0936e-11' 
		rate_values[5415] = 2.0936e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.9987e-11' 
		rate_values[5416] = 2.9987e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4827e-11' 
		rate_values[5417] = 3.4827e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4827e-11' 
		rate_values[5418] = 3.4827e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5419] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4417e-11' 
		rate_values[5420] = 1.4417e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.586e-11' 
		rate_values[5421] = 2.586e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.361e-11' 
		rate_values[5422] = 3.361e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4827e-11' 
		rate_values[5423] = 3.4827e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5424] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0936e-11' 
		rate_values[5425] = 2.0936e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.9987e-11' 
		rate_values[5426] = 2.9987e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5427] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4417e-11' 
		rate_values[5428] = 1.4417e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.586e-11' 
		rate_values[5429] = 2.586e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5430] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5431] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5432] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5433] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5434] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5435] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5436] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5437] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5438] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5439] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5440] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5441] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0936e-11' 
		rate_values[5442] = 2.0936e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4417e-11' 
		rate_values[5443] = 1.4417e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5444] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5445] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5446] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5447] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5448] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5449] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5450] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5451] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5452] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5453] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5454] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5455] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5456] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5457] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5458] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5459] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9252e-12' 
		rate_values[5460] = 4.9252e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5461] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5462] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5463] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5464] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5465] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5466] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5467] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5468] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5469] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[5470] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5471] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5472] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5473] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5474] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5475] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5476] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5477] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5478] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5479] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5480] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5481] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5482] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5483] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5484] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5485] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5486] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5487] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5488] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5489] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5490] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5491] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5492] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5493] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5494] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5495] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5496] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5497] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5498] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5499] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5500] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5501] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5502] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5503] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5504] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5505] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5506] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5507] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5508] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5509] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5510] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5511] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5512] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5513] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5514] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5515] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5516] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5517] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[5518] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5519] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5520] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5521] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5522] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5523] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5524] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5525] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5526] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5527] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5528] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5529] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5530] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5531] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5532] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5533] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5534] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5535] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5536] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5537] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5538] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5539] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5540] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5541] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5542] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5543] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5544] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5545] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5546] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5547] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5548] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5549] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5550] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5551] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5552] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5553] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5554] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5555] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5556] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5557] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5558] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[5559] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8506e-11' 
		rate_values[5560] = 3.8506e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[5561] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[5562] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5563] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[5564] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3207e-11' 
		rate_values[5565] = 3.3207e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.3159e-11' 
		rate_values[5566] = 4.3159e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[5567] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5568] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[5569] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8506e-11' 
		rate_values[5570] = 3.8506e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5571] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[5572] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3207e-11' 
		rate_values[5573] = 3.3207e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5574] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5575] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5576] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5577] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5578] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5579] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5580] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5581] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5582] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5583] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5584] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5585] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[5586] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[5587] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5588] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5589] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5590] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5591] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5592] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5593] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5594] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5595] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5596] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5597] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5598] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5599] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5600] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5601] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5602] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5603] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5604] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5605] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5606] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[5607] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8506e-11' 
		rate_values[5608] = 3.8506e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[5609] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[5610] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5611] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[5612] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3207e-11' 
		rate_values[5613] = 3.3207e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.3159e-11' 
		rate_values[5614] = 4.3159e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.4721e-11' 
		rate_values[5615] = 4.4721e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5616] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[5617] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8506e-11' 
		rate_values[5618] = 3.8506e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5619] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[5620] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3207e-11' 
		rate_values[5621] = 3.3207e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5622] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5623] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5624] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5625] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5626] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5627] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5628] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5629] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5630] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5631] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5632] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5633] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6884e-11' 
		rate_values[5634] = 2.6884e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8513e-11' 
		rate_values[5635] = 1.8513e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5636] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5637] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5638] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5639] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5640] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5641] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5642] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5643] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5644] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5645] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5646] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5647] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5648] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5649] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5650] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5651] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.3246e-12' 
		rate_values[5652] = 6.3246e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5653] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5654] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5655] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5656] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5657] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5658] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5659] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5660] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5661] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[5662] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5663] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5664] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5665] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5666] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5667] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5668] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5669] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5670] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5671] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5672] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5673] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5674] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5675] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5676] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5677] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5678] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5679] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5680] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5681] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5682] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5683] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5684] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5685] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5686] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5687] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5688] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5689] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5690] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5691] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5692] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5693] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5694] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5695] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5696] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5697] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5698] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5699] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5700] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5701] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5702] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.7712e-12' 
		rate_values[5703] = 8.7712e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2563e-11' 
		rate_values[5704] = 1.2563e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4591e-11' 
		rate_values[5705] = 1.4591e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4591e-11' 
		rate_values[5706] = 1.4591e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5707] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0401e-12' 
		rate_values[5708] = 6.0401e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0834e-11' 
		rate_values[5709] = 1.0834e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4081e-11' 
		rate_values[5710] = 1.4081e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4591e-11' 
		rate_values[5711] = 1.4591e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5712] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.7712e-12' 
		rate_values[5713] = 8.7712e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2563e-11' 
		rate_values[5714] = 1.2563e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5715] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0401e-12' 
		rate_values[5716] = 6.0401e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0834e-11' 
		rate_values[5717] = 1.0834e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5718] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5719] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5720] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5721] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5722] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5723] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5724] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5725] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5726] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5727] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5728] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5729] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.7712e-12' 
		rate_values[5730] = 8.7712e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0401e-12' 
		rate_values[5731] = 6.0401e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5732] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5733] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5734] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5735] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5736] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5737] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5738] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5739] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5740] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5741] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5742] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5743] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5744] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5745] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5746] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5747] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0635e-12' 
		rate_values[5748] = 2.0635e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5749] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5750] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5751] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5752] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5753] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5754] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5755] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5756] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5757] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3639e-12' 
		rate_values[5758] = 2.3639e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4495e-12' 
		rate_values[5759] = 2.4495e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5760] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5761] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.1091e-12' 
		rate_values[5762] = 2.1091e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5763] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5764] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8188e-12' 
		rate_values[5765] = 1.8188e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5766] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5767] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5768] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5769] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5770] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5771] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5772] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5773] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5774] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5775] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5776] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5777] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4725e-12' 
		rate_values[5778] = 1.4725e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.014e-12' 
		rate_values[5779] = 1.014e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5780] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5781] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5782] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5783] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5784] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5785] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5786] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5787] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5788] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5789] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5790] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5791] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5792] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5793] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5794] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5795] = 3.4641e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.4641e-13' 
		rate_values[5796] = 3.4641e-13
	except:
		erf = 1 # flag error
		err_mess = (str('Error: Could not calculate rate coefficient for equation number ' + str(gprn) + ' ' + rc_eq_now + ' (message from rate coeffs.py)'))
	
	# aqueous-phase reactions
	
	# surface (e.g. wall) reactions
	
	return(rate_values, erf, err_mess)
