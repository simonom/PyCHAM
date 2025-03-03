##########################################################################################
#                                                                                        #
#    Copyright (C) 2018-2025 Simon O'Meara : simon.omeara@manchester.ac.uk               #
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
# created at 2025-03-03 17:31:55.089275

import numpy
import photolysisRates

def evaluate_rates(RO2, H2O, RH, TEMP, time, M, N2, O2, Jlen, NO, HO2, NO3, sumt, self):

	# inputs: ------------------------------------------------------------------
	# RO2 - total concentration of alkyl peroxy radicals (# molecules/cm3) 
	# M - third body concentration (# molecules/cm3 (air))
	# N2 - nitrogen concentration (# molecules/cm3 (air))
	# O2 - oxygen concentration (# molecules/cm3 (air))
	# H2O, RH, TEMP: given by the user
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

	# calculate any generic reaction rate 
	# coefficients given by chemical scheme 

	try:
		gprn=0
		# keep count on reaction number 
		gprn += 1 
		KRO2NO=2.7e-12*numpy.exp(360/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KRO2HO2=2.91e-13*numpy.exp(1300/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KAPHO2=5.2e-13*numpy.exp(980/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KAPNO=7.5e-12*numpy.exp(290/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KRO2NO3=2.3e-12 
		# keep count on reaction number 
		gprn += 1 
		KNO3AL=1.4e-12*numpy.exp(-1860/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KDEC=1.00e+06 
		# keep count on reaction number 
		gprn += 1 
		KROPRIM=2.50e-14*numpy.exp(-300/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KROSEC=2.50e-14*numpy.exp(-300/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KCH3O2=1.03e-13*numpy.exp(365/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		K298CH3O2=3.5e-13 
		# keep count on reaction number 
		gprn += 1 
		K14ISOM1=3.00e7*numpy.exp(-5300/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		FCD=0.30 
		# keep count on reaction number 
		gprn += 1 
		KD0=1.10e-05*M*numpy.exp(-10100/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KDI=1.90e17*numpy.exp(-14100/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KRD=KD0/KDI 
		# keep count on reaction number 
		gprn += 1 
		NCD=0.75-1.27*(numpy.log10(FCD)) 
		# keep count on reaction number 
		gprn += 1 
		FD=10**(numpy.log10(FCD)/(1+(numpy.log10(KRD)/NCD)**(2))) 
		# keep count on reaction number 
		gprn += 1 
		KBPAN=(KD0*KDI)*FD/(KD0+KDI) 
		# keep count on reaction number 
		gprn += 1 
		FCC=0.30 
		# keep count on reaction number 
		gprn += 1 
		KC0=3.28e-28*M*(TEMP/300)**(-6.87) 
		# keep count on reaction number 
		gprn += 1 
		KCI=1.125e-11*(TEMP/300)**(-1.105) 
		# keep count on reaction number 
		gprn += 1 
		KRC=KC0/KCI 
		# keep count on reaction number 
		gprn += 1 
		NC=0.75-1.27*(numpy.log10(FCC)) 
		# keep count on reaction number 
		gprn += 1 
		FC=10**(numpy.log10(FCC)/(1+(numpy.log10(KRC)/NC)**(2))) 
		# keep count on reaction number 
		gprn += 1 
		KFPAN=(KC0*KCI)*FC/(KC0+KCI) 
		# keep count on reaction number 
		gprn += 1 
		FC1=0.85 
		# keep count on reaction number 
		gprn += 1 
		K10=1.0e-31*M*(TEMP/300)**(-1.6) 
		# keep count on reaction number 
		gprn += 1 
		K1I=5.0e-11*(TEMP/300)**(-0.3) 
		# keep count on reaction number 
		gprn += 1 
		KR1=K10/K1I 
		# keep count on reaction number 
		gprn += 1 
		NC1=0.75-1.27*(numpy.log10(FC1)) 
		# keep count on reaction number 
		gprn += 1 
		F1=10**(numpy.log10(FC1)/(1+(numpy.log10(KR1)/NC1)**(2))) 
		# keep count on reaction number 
		gprn += 1 
		KMT01=(K10*K1I)*F1/(K10+K1I) 
		# keep count on reaction number 
		gprn += 1 
		FC2=0.6 
		# keep count on reaction number 
		gprn += 1 
		K20=1.3e-31*M*(TEMP/300)**(-1.5) 
		# keep count on reaction number 
		gprn += 1 
		K2I=2.3e-11*(TEMP/300)**(0.24) 
		# keep count on reaction number 
		gprn += 1 
		KR2=K20/K2I 
		# keep count on reaction number 
		gprn += 1 
		NC2=0.75-1.27*(numpy.log10(FC2)) 
		# keep count on reaction number 
		gprn += 1 
		F2=10**(numpy.log10(FC2)/(1+(numpy.log10(KR2)/NC2)**(2))) 
		# keep count on reaction number 
		gprn += 1 
		KMT02=(K20*K2I)*F2/(K20+K2I) 
		# keep count on reaction number 
		gprn += 1 
		FC3=0.35 
		# keep count on reaction number 
		gprn += 1 
		K30=3.6e-30*M*(TEMP/300)**(-4.1) 
		# keep count on reaction number 
		gprn += 1 
		K3I=1.9e-12*(TEMP/300)**(0.2) 
		# keep count on reaction number 
		gprn += 1 
		KR3=K30/K3I 
		# keep count on reaction number 
		gprn += 1 
		NC3=0.75-1.27*(numpy.log10(FC3)) 
		# keep count on reaction number 
		gprn += 1 
		F3=10**(numpy.log10(FC3)/(1+(numpy.log10(KR3)/NC3)**(2))) 
		# keep count on reaction number 
		gprn += 1 
		KMT03=(K30*K3I)*F3/(K30+K3I) 
		# keep count on reaction number 
		gprn += 1 
		FC4=0.35 
		# keep count on reaction number 
		gprn += 1 
		K40=1.3e-3*M*(TEMP/300)**(-3.5)*numpy.exp(-11000/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		K4I=9.7e+14*(TEMP/300)**(0.1)*numpy.exp(-11080/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KR4=K40/K4I 
		# keep count on reaction number 
		gprn += 1 
		NC4=0.75-1.27*(numpy.log10(FC4)) 
		# keep count on reaction number 
		gprn += 1 
		F4=10**(numpy.log10(FC4)/(1+(numpy.log10(KR4)/NC4)**(2))) 
		# keep count on reaction number 
		gprn += 1 
		KMT05=1.44e-13*(1+(M/4.2e+19)) 
		# keep count on reaction number 
		gprn += 1 
		KMT04=(K40*K4I)*F4/(K40+K4I) 
		# keep count on reaction number 
		gprn += 1 
		FC7=0.81 
		# keep count on reaction number 
		gprn += 1 
		K70=7.4e-31*M*(TEMP/300)**(-2.4) 
		# keep count on reaction number 
		gprn += 1 
		K7I=3.3e-11*(TEMP/300)**(-0.3) 
		# keep count on reaction number 
		gprn += 1 
		KR7=K70/K7I 
		# keep count on reaction number 
		gprn += 1 
		NC7=0.75-1.27*(numpy.log10(FC7)) 
		# keep count on reaction number 
		gprn += 1 
		F7=10**(numpy.log10(FC7)/(1+(numpy.log10(KR7)/NC7)**(2))) 
		# keep count on reaction number 
		gprn += 1 
		KMT07=(K70*K7I)*F7/(K70+K7I) 
		# keep count on reaction number 
		gprn += 1 
		KMT06=1+(1.40e-21*numpy.exp(2200/TEMP)*H2O) 
		# keep count on reaction number 
		gprn += 1 
		FC8=0.41 
		# keep count on reaction number 
		gprn += 1 
		K80=3.2e-30*M*(TEMP/300)**(-4.5) 
		# keep count on reaction number 
		gprn += 1 
		K8I=3.0e-11 
		# keep count on reaction number 
		gprn += 1 
		KR8=K80/K8I 
		# keep count on reaction number 
		gprn += 1 
		NC8=0.75-1.27*(numpy.log10(FC8)) 
		# keep count on reaction number 
		gprn += 1 
		F8=10**(numpy.log10(FC8)/(1+(numpy.log10(KR8)/NC8)**(2))) 
		# keep count on reaction number 
		gprn += 1 
		KMT08=(K80*K8I)*F8/(K80+K8I) 
		# keep count on reaction number 
		gprn += 1 
		FC9=0.4 
		# keep count on reaction number 
		gprn += 1 
		K90=1.4e-31*M*(TEMP/300)**(-3.1) 
		# keep count on reaction number 
		gprn += 1 
		K9I=4.0e-12 
		# keep count on reaction number 
		gprn += 1 
		KR9=K90/K9I 
		# keep count on reaction number 
		gprn += 1 
		NC9=0.75-1.27*(numpy.log10(FC9)) 
		# keep count on reaction number 
		gprn += 1 
		F9=10**(numpy.log10(FC9)/(1+(numpy.log10(KR9)/NC9)**(2))) 
		# keep count on reaction number 
		gprn += 1 
		KMT09=(K90*K9I)*F9/(K90+K9I) 
		# keep count on reaction number 
		gprn += 1 
		FC10=0.4 
		# keep count on reaction number 
		gprn += 1 
		K100=4.10e-05*M*numpy.exp(-10650/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		K10I=6.0e+15*numpy.exp(-11170/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KR10=K100/K10I 
		# keep count on reaction number 
		gprn += 1 
		NC10=0.75-1.27*(numpy.log10(FC10)) 
		# keep count on reaction number 
		gprn += 1 
		F10=10**(numpy.log10(FC10)/(1+(numpy.log10(KR10)/NC10)**(2))) 
		# keep count on reaction number 
		gprn += 1 
		KMT10=(K100*K10I)*F10/(K100+K10I) 
		# keep count on reaction number 
		gprn += 1 
		K3=6.50e-34*numpy.exp(1335/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		K4=2.70e-17*numpy.exp(2199/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		K1=2.40e-14*numpy.exp(460/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		K2=(K3*M)/(1+(K3*M/K4)) 
		# keep count on reaction number 
		gprn += 1 
		KMT11=K1+K2 
		# keep count on reaction number 
		gprn += 1 
		FC12=0.53 
		# keep count on reaction number 
		gprn += 1 
		K120=2.5e-31*M*(TEMP/300)**(-2.6) 
		# keep count on reaction number 
		gprn += 1 
		K12I=2.0e-12 
		# keep count on reaction number 
		gprn += 1 
		KR12=K120/K12I 
		# keep count on reaction number 
		gprn += 1 
		NC12=0.75-1.27*(numpy.log10(FC12)) 
		# keep count on reaction number 
		gprn += 1 
		F12=10**(numpy.log10(FC12)/(1.0+(numpy.log10(KR12)/NC12)**(2))) 
		# keep count on reaction number 
		gprn += 1 
		KMT12=(K120*K12I*F12)/(K120+K12I) 
		# keep count on reaction number 
		gprn += 1 
		FC13=0.36 
		# keep count on reaction number 
		gprn += 1 
		K130=2.5e-30*M*(TEMP/300)**(-5.5) 
		# keep count on reaction number 
		gprn += 1 
		K13I=1.8e-11 
		# keep count on reaction number 
		gprn += 1 
		KR13=K130/K13I 
		# keep count on reaction number 
		gprn += 1 
		NC13=0.75-1.27*(numpy.log10(FC13)) 
		# keep count on reaction number 
		gprn += 1 
		F13=10**(numpy.log10(FC13)/(1+(numpy.log10(KR13)/NC13)**(2))) 
		# keep count on reaction number 
		gprn += 1 
		KMT13=(K130*K13I)*F13/(K130+K13I) 
		# keep count on reaction number 
		gprn += 1 
		FC14=0.36 
		# keep count on reaction number 
		gprn += 1 
		K140=9.0e-5*numpy.exp(-9690/TEMP)*M 
		# keep count on reaction number 
		gprn += 1 
		K14I=1.1e+16*numpy.exp(-10560/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KR14=K140/K14I 
		# keep count on reaction number 
		gprn += 1 
		NC14=0.75-1.27*(numpy.log10(FC14)) 
		# keep count on reaction number 
		gprn += 1 
		F14=10**(numpy.log10(FC14)/(1+(numpy.log10(KR14)/NC14)**(2))) 
		# keep count on reaction number 
		gprn += 1 
		KMT14=(K140*K14I)*F14/(K140+K14I) 

	except:
		erf = 1 # flag error
		err_mess = str('Error: generic reaction rates failed to be calculated inside rate_coeffs.py at number ' + str(gprn) + ', please check chemical scheme and associated chemical scheme markers, which are stated in the model variables input file') # error message
		return([], erf, err_mess)
	# estimate and append photolysis rates
	J = photolysisRates.PhotolysisCalculation(TEMP, Jlen, sumt, self)

	if (self.light_stat_now == 0):
		J = [0]*len(J)
	rate_values = numpy.zeros((280))
	
	# if reactions have been found in the chemical scheme
	# gas-phase reactions
	gprn = 0 # keep count on reaction number
	try:
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.6e-34*N2*(TEMP/300)**(-2.6)*O2' 
		rate_values[0] = 5.6e-34*N2*(TEMP/300)**(-2.6)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0e-34*O2*(TEMP/300)**(-2.6)*O2' 
		rate_values[1] = 6.0e-34*O2*(TEMP/300)**(-2.6)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.0e-12*numpy.exp(-2060/TEMP)' 
		rate_values[2] = 8.0e-12*numpy.exp(-2060/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT01' 
		rate_values[3] = KMT01
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5e-12*numpy.exp(188/TEMP)' 
		rate_values[4] = 5.5e-12*numpy.exp(188/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT02' 
		rate_values[5] = KMT02
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.2e-11*numpy.exp(67/TEMP)*O2' 
		rate_values[6] = 3.2e-11*numpy.exp(67/TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0e-11*numpy.exp(130/TEMP)*N2' 
		rate_values[7] = 2.0e-11*numpy.exp(130/TEMP)*N2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4e-12*numpy.exp(-1310/TEMP)' 
		rate_values[8] = 1.4e-12*numpy.exp(-1310/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4e-13*numpy.exp(-2470/TEMP)' 
		rate_values[9] = 1.4e-13*numpy.exp(-2470/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3e-39*numpy.exp(530/TEMP)*O2' 
		rate_values[10] = 3.3e-39*numpy.exp(530/TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8e-11*numpy.exp(110/TEMP)' 
		rate_values[11] = 1.8e-11*numpy.exp(110/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.50e-14*numpy.exp(-1260/TEMP)' 
		rate_values[12] = 4.50e-14*numpy.exp(-1260/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT03' 
		rate_values[13] = KMT03
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.14e-10*H2O' 
		rate_values[14] = 2.14e-10*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.70e-12*numpy.exp(-940/TEMP)' 
		rate_values[15] = 1.70e-12*numpy.exp(-940/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.7e-12*numpy.exp(-2100/TEMP)' 
		rate_values[16] = 7.7e-12*numpy.exp(-2100/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT05' 
		rate_values[17] = KMT05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.9e-12*numpy.exp(-160/TEMP)' 
		rate_values[18] = 2.9e-12*numpy.exp(-160/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.03e-16*(TEMP/300)**(4.57)*numpy.exp(693/TEMP)' 
		rate_values[19] = 2.03e-16*(TEMP/300)**(4.57)*numpy.exp(693/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.8e-11*numpy.exp(250/TEMP)' 
		rate_values[20] = 4.8e-11*numpy.exp(250/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.20e-13*KMT06*numpy.exp(600/TEMP)' 
		rate_values[21] = 2.20e-13*KMT06*numpy.exp(600/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-33*M*KMT06*numpy.exp(980/TEMP)' 
		rate_values[22] = 1.90e-33*M*KMT06*numpy.exp(980/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT07' 
		rate_values[23] = KMT07
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT08' 
		rate_values[24] = KMT08
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0e-11' 
		rate_values[25] = 2.0e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.45e-12*numpy.exp(270/TEMP)' 
		rate_values[26] = 3.45e-12*numpy.exp(270/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT09' 
		rate_values[27] = KMT09
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.2e-13*numpy.exp(690/TEMP)*1.0' 
		rate_values[28] = 3.2e-13*numpy.exp(690/TEMP)*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0e-12' 
		rate_values[29] = 4.0e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.5e-12*numpy.exp(260/TEMP)' 
		rate_values[30] = 2.5e-12*numpy.exp(260/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT11' 
		rate_values[31] = KMT11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0e-32*numpy.exp(-1000/TEMP)*M' 
		rate_values[32] = 4.0e-32*numpy.exp(-1000/TEMP)*M
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT12' 
		rate_values[33] = KMT12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.3e-12*numpy.exp(-330/TEMP)*O2' 
		rate_values[34] = 1.3e-12*numpy.exp(-330/TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.00e-06' 
		rate_values[35] = 6.00e-06
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.00e-04' 
		rate_values[36] = 4.00e-04
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.20e-15*H2O' 
		rate_values[37] = 1.20e-15*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[1]' 
		rate_values[38] = J[1]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[2]' 
		rate_values[39] = J[2]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[3]' 
		rate_values[40] = J[3]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]' 
		rate_values[41] = J[4]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[5]' 
		rate_values[42] = J[5]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[6]' 
		rate_values[43] = J[6]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[7]' 
		rate_values[44] = J[7]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[8]' 
		rate_values[45] = J[8]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT04' 
		rate_values[46] = KMT04
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT10' 
		rate_values[47] = KMT10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.5e-13' 
		rate_values[48] = 4.5e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.99e-11*0.16' 
		rate_values[49] = 6.99e-11*0.16
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.99e-11*0.06' 
		rate_values[50] = 6.99e-11*0.06
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.99e-11*0.119' 
		rate_values[51] = 6.99e-11*0.119
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.99e-11*0.179' 
		rate_values[52] = 6.99e-11*0.179
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.99e-11*0.482' 
		rate_values[53] = 6.99e-11*0.482
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.10e-17' 
		rate_values[54] = 2.10e-17
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[55] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.60e-10' 
		rate_values[56] = 2.60e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[57] = 2.08e-12
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
		rc_eq_now = '4.378e-12*0.8*RO2' 
		rate_values[60] = 4.378e-12*0.8*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.378e-12*0.2*RO2' 
		rate_values[61] = 4.378e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[62] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[63] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.56e-10' 
		rate_values[64] = 2.56e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.3' 
		rate_values[65] = KDEC*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.7' 
		rate_values[66] = KDEC*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[31]' 
		rate_values[67] = J[31]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[32]' 
		rate_values[68] = J[32]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[33]' 
		rate_values[69] = J[33]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.1e-12*numpy.exp(340/TEMP)' 
		rate_values[70] = 3.1e-12*numpy.exp(340/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL' 
		rate_values[71] = KNO3AL
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[72] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[73] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[74] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[75] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[76] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.7*RO2' 
		rate_values[77] = 1.00e-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.3*RO2' 
		rate_values[78] = 1.00e-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.58e-11' 
		rate_values[79] = 1.58e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[80] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[81] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.23e-11' 
		rate_values[82] = 1.23e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.4e-12*numpy.exp(135/TEMP)' 
		rate_values[83] = 5.4e-12*numpy.exp(135/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[84] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[12]' 
		rate_values[85] = J[12]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5e-16' 
		rate_values[86] = 5.5e-16
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[87] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[88] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00e11*numpy.exp(-3160/TEMP)' 
		rate_values[89] = 7.00e11*numpy.exp(-3160/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00e-12*O2' 
		rate_values[90] = 5.00e-12*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00e-12*O2*3.2*numpy.exp(-550/TEMP)' 
		rate_values[91] = 5.00e-12*O2*3.2*numpy.exp(-550/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00e-12*O2*3.2*(1-numpy.exp(-550/TEMP))' 
		rate_values[92] = 5.00e-12*O2*3.2*(1-numpy.exp(-550/TEMP))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.73e-10*0.03' 
		rate_values[93] = 1.73e-10*0.03
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.73e-10*0.618' 
		rate_values[94] = 1.73e-10*0.618
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.73e-10*0.265' 
		rate_values[95] = 1.73e-10*0.265
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.73e-10*0.087' 
		rate_values[96] = 1.73e-10*0.087
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[97] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[98] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[99] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[100] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.25e-13*0.2*RO2' 
		rate_values[101] = 1.25e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.25e-13*0.8*RO2' 
		rate_values[102] = 1.25e-13*0.8*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[103] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.72e-10' 
		rate_values[104] = 2.72e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.3' 
		rate_values[105] = KDEC*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.7' 
		rate_values[106] = KDEC*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.5' 
		rate_values[107] = J[4]*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.5' 
		rate_values[108] = J[4]*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.59e-11' 
		rate_values[109] = 5.59e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[110] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[111] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.5' 
		rate_values[112] = KAPHO2*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.37' 
		rate_values[113] = KAPHO2*0.37
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.13' 
		rate_values[114] = KAPHO2*0.13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.2' 
		rate_values[115] = 1.00e-11*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.8' 
		rate_values[116] = 1.00e-11*RO2*0.8
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[117] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[118] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[119] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[120] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.89e-11' 
		rate_values[121] = 4.89e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.8' 
		rate_values[122] = KDEC*0.8
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.2' 
		rate_values[123] = KDEC*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.73e-9*0.13' 
		rate_values[124] = 1.73e-9*0.13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.73e-9*0.26' 
		rate_values[125] = 1.73e-9*0.26
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.73e-9*0.61' 
		rate_values[126] = 1.73e-9*0.61
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.95e-11' 
		rate_values[127] = 7.95e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[128] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[129] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-14*0.6*RO2' 
		rate_values[130] = 2.00e-14*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-14*0.2*RO2' 
		rate_values[131] = 2.00e-14*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-14*0.2*RO2' 
		rate_values[132] = 2.00e-14*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[133] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[134] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.43e-10' 
		rate_values[135] = 1.43e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.45e-11' 
		rate_values[136] = 7.45e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.32e-11' 
		rate_values[137] = 5.32e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[138] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.62e-11' 
		rate_values[139] = 2.62e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[140] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[141] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[142] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[143] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.01' 
		rate_values[144] = KDEC*0.01
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.99' 
		rate_values[145] = KDEC*0.99
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4e-12' 
		rate_values[146] = 1.4e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[147] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[148] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[149] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[150] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[151] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[152] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[153] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[154] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.97e-11' 
		rate_values[155] = 6.97e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[156] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.33e-11' 
		rate_values[157] = 7.33e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[158] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[159] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.6*RO2' 
		rate_values[160] = 8.80e-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[161] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[162] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[163] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.66e-11' 
		rate_values[164] = 4.66e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.55e-11' 
		rate_values[165] = 2.55e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.68e-12' 
		rate_values[166] = 5.68e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[167] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[168] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[169] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.42' 
		rate_values[170] = KAPHO2*0.42
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[171] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[172] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[173] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.58' 
		rate_values[174] = KAPHO2*0.58
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.74e-11' 
		rate_values[175] = 5.74e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[176] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[177] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.11e-11' 
		rate_values[178] = 1.11e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[179] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.65e-12' 
		rate_values[180] = 1.65e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[181] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[182] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.39e-12*0.8*RO2' 
		rate_values[183] = 4.39e-12*0.8*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.39e-12*0.2*RO2' 
		rate_values[184] = 4.39e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[185] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[186] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[187] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.91e-11' 
		rate_values[188] = 7.91e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.1*0.5' 
		rate_values[189] = J[4]*0.1*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.1*0.5' 
		rate_values[190] = J[4]*0.1*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[191] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.36e-10' 
		rate_values[192] = 1.36e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[193] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[194] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[195] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.520' 
		rate_values[196] = KRO2HO2*0.520
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[197] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[198] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[199] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[200] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[201] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.44e-10' 
		rate_values[202] = 1.44e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[203] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.77e-11' 
		rate_values[204] = 5.77e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]*2' 
		rate_values[205] = J[17]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.32e-11' 
		rate_values[206] = 4.32e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KNO3AL*4.0' 
		rate_values[207] = 2*KNO3AL*4.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[208] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[209] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[210] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[211] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[212] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[213] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[214] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[215] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.29e-11' 
		rate_values[216] = 2.29e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[217] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[218] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.62e-11' 
		rate_values[219] = 2.62e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[220] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.31e-11' 
		rate_values[221] = 2.31e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[222] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.93e-11*0.82' 
		rate_values[223] = 2.93e-11*0.82
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.93e-11*0.18' 
		rate_values[224] = 2.93e-11*0.18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[225] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[226] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[227] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[228] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[229] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[230] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[231] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[232] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[233] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[234] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0e-10' 
		rate_values[235] = 1.0e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.9e-11' 
		rate_values[236] = 9.9e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.2e-18' 
		rate_values[237] = 9.2e-18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[238] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[239] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[240] = 2.86e-13
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
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[243] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[244] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[245] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[246] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.47e-12' 
		rate_values[247] = 3.47e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.60e-12' 
		rate_values[248] = 2.60e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[249] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[250] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[251] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[252] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[253] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[254] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[255] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[256] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[257] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[258] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[259] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[260] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[261] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[262] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[263] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[264] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[265] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[266] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.70e-11' 
		rate_values[267] = 3.70e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.40' 
		rate_values[268] = KDEC*0.40
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.60' 
		rate_values[269] = KDEC*0.60
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'self.JFAC*1.62E7+(1-self.JFAC)*1.40E7' 
		rate_values[270] = self.JFAC*1.62E7+(1-self.JFAC)*1.40E7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.1e17*J[4]*numpy.exp(-5686/TEMP)*3' 
		rate_values[271] = 3.1e17*J[4]*numpy.exp(-5686/TEMP)*3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.3e21*J[4]*numpy.exp(-8945/TEMP)+J[4]*(5.8E8)*(RH**0.36)' 
		rate_values[272] = 7.3e21*J[4]*numpy.exp(-8945/TEMP)+J[4]*(5.8E8)*(RH**0.36)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'self.JFAC*3.86E6' 
		rate_values[273] = self.JFAC*3.86E6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.7e-5' 
		rate_values[274] = 0.7e-5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6e-5' 
		rate_values[275] = 1.6e-5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[8]' 
		rate_values[276] = J[8]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*2.46e9' 
		rate_values[277] = J[4]*2.46e9
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'self.JFAC*5.98e6+(1-self.JFAC)*1.85e6*1.7' 
		rate_values[278] = self.JFAC*5.98e6+(1-self.JFAC)*1.85e6*1.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'self.JFAC*6.4e-5' 
		rate_values[279] = self.JFAC*6.4e-5
	except:
		erf = 1 # flag error
		err_mess = (str('Error: Could not calculate '+ 
		'rate coefficient for equation number ' 
		+ str(gprn) + ' ' + rc_eq_now + 
		' (message from rate coeffs.py)'))
	
	# aqueous-phase reactions
	
	# surface (e.g. wall) reactions
	
	return(rate_values, erf, err_mess)
