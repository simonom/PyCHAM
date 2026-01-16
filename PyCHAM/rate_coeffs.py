##########################################################################################
#                                                                                        #
#    Copyright (C) 2018-2026 Simon O'Meara : simon.omeara@manchester.ac.uk               #
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
# created at 2026-01-16 11:34:04.335530

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
		K298CH3O2=3.5e-13 
		# keep count on reaction number 
		gprn += 1 
		KAPHO2=5.2e-13*numpy.exp(980/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KAPNO=7.5e-12*numpy.exp(290/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KCH3O2=1.03e-13*numpy.exp(365/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KDEC=1.00e+06 
		# keep count on reaction number 
		gprn += 1 
		KMT05=1.44e-13*(1+(M/4.2e+19)) 
		# keep count on reaction number 
		gprn += 1 
		KMT06=1+(1.40e-21*numpy.exp(2200/TEMP)*H2O) 
		# keep count on reaction number 
		gprn += 1 
		KNO3AL=1.44e-12*numpy.exp(-1862/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KRO2HO2=2.91e-13*numpy.exp(1300/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KRO2NO=2.7e-12*numpy.exp(360/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KRO2NO3=2.3e-12 
		# keep count on reaction number 
		gprn += 1 
		KROPRIM=2.50e-14*numpy.exp(-300/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KROSEC=2.50e-14*numpy.exp(-300/TEMP) 
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
		FCPPN=0.36 
		# keep count on reaction number 
		gprn += 1 
		KPPN0=1.7e-03*numpy.exp(-11280/TEMP)*M 
		# keep count on reaction number 
		gprn += 1 
		KPPNI=8.3e+16*numpy.exp(-13940/TEMP) 
		# keep count on reaction number 
		gprn += 1 
		KRPPN=KPPN0/KPPNI 
		# keep count on reaction number 
		gprn += 1 
		NCPPN=0.75-1.27*(numpy.log10(FCPPN)) 
		# keep count on reaction number 
		gprn += 1 
		FPPN=10**(numpy.log10(FCPPN)/(1+(numpy.log10(KRPPN)/NCPPN)**(2))) 
		# keep count on reaction number 
		gprn += 1 
		KBPPN=(KPPN0*KPPNI)*FPPN/(KPPN0+KPPNI) 
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
	rate_values = numpy.zeros((3173))
	
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
		rc_eq_now = '2.3e-12*numpy.exp(360/TEMP)*0.999' 
		rate_values[48] = 2.3e-12*numpy.exp(360/TEMP)*0.999
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3e-12*numpy.exp(360/TEMP)*0.001' 
		rate_values[49] = 2.3e-12*numpy.exp(360/TEMP)*0.001
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.2e-14*numpy.exp(-1080/TEMP)*O2' 
		rate_values[50] = 7.2e-14*numpy.exp(-1080/TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT13' 
		rate_values[51] = KMT13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT14' 
		rate_values[52] = KMT14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-12' 
		rate_values[53] = 1.2e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8e-13*numpy.exp(780/TEMP)*(1-1/(1+498*numpy.exp(-1160/TEMP)))' 
		rate_values[54] = 3.8e-13*numpy.exp(780/TEMP)*(1-1/(1+498*numpy.exp(-1160/TEMP)))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KCH3O2*RO2*7.18*numpy.exp(-885/TEMP)' 
		rate_values[55] = 2*KCH3O2*RO2*7.18*numpy.exp(-885/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KCH3O2*RO2*0.5*(1-7.18*numpy.exp(-885/TEMP))' 
		rate_values[56] = 2*KCH3O2*RO2*0.5*(1-7.18*numpy.exp(-885/TEMP))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KCH3O2*RO2*0.5*(1-7.18*numpy.exp(-885/TEMP))' 
		rate_values[57] = 2*KCH3O2*RO2*0.5*(1-7.18*numpy.exp(-885/TEMP))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0e-13*numpy.exp(-845/TEMP)' 
		rate_values[58] = 4.0e-13*numpy.exp(-845/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[51]' 
		rate_values[59] = J[51]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.3e-12*numpy.exp(190/TEMP)*0.6' 
		rate_values[60] = 5.3e-12*numpy.exp(190/TEMP)*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.3e-12*numpy.exp(190/TEMP)*0.4' 
		rate_values[61] = 5.3e-12*numpy.exp(190/TEMP)*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[62] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.55e-12*numpy.exp(380/TEMP)*0.991' 
		rate_values[63] = 2.55e-12*numpy.exp(380/TEMP)*0.991
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.55e-12*numpy.exp(380/TEMP)*0.009' 
		rate_values[64] = 2.55e-12*numpy.exp(380/TEMP)*0.009
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4e-14*numpy.exp(-325/TEMP)*O2' 
		rate_values[65] = 2.4e-14*numpy.exp(-325/TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[66] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.3e-13*numpy.exp(870/TEMP)' 
		rate_values[67] = 4.3e-13*numpy.exp(870/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*(KCH3O2*6.4e-14)**(0.5)*RO2*0.6' 
		rate_values[68] = 2*(KCH3O2*6.4e-14)**(0.5)*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*(KCH3O2*6.4e-14)**(0.5)*RO2*0.2' 
		rate_values[69] = 2*(KCH3O2*6.4e-14)**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*(KCH3O2*6.4e-14)**(0.5)*RO2*0.2' 
		rate_values[70] = 2*(KCH3O2*6.4e-14)**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.7e-13*numpy.exp(-395/TEMP)' 
		rate_values[71] = 6.7e-13*numpy.exp(-395/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[52]' 
		rate_values[72] = J[52]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[73] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.01e-12' 
		rate_values[74] = 8.01e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[75] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[76] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[77] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.65e-11' 
		rate_values[78] = 6.65e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[79] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[80] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[81] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[82] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[83] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[84] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[85] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.12e-13' 
		rate_values[86] = 3.12e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[87] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[88] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.23e-12' 
		rate_values[89] = 4.23e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[90] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[91] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.29e-11' 
		rate_values[92] = 4.29e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KNO3AL*2.4' 
		rate_values[93] = 2*KNO3AL*2.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[94] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[95] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[96] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[97] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[98] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[99] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.7*RO2' 
		rate_values[100] = 1.00e-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.3*RO2' 
		rate_values[101] = 1.00e-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.10e-11' 
		rate_values[102] = 2.10e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.49e-11' 
		rate_values[103] = 2.49e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[104] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[105] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.14e-11' 
		rate_values[106] = 2.14e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[107] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[108] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.40e-18*TEMP**(2)*numpy.exp(194/TEMP)' 
		rate_values[109] = 1.40e-18*TEMP**(2)*numpy.exp(194/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[110] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[111] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[112] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[113] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.6*RO2' 
		rate_values[114] = 2.00e-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.2*RO2' 
		rate_values[115] = 2.00e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.2*RO2' 
		rate_values[116] = 2.00e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.99e-12' 
		rate_values[117] = 5.99e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[118] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[119] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[120] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.69e-12' 
		rate_values[121] = 2.69e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[122] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.23e-11' 
		rate_values[123] = 1.23e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*4.0' 
		rate_values[124] = KNO3AL*4.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[125] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.39e-11' 
		rate_values[126] = 1.39e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*4.0' 
		rate_values[127] = KNO3AL*4.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.995' 
		rate_values[128] = KRO2NO*0.995
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.005' 
		rate_values[129] = KRO2NO*0.005
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.50e+13*numpy.exp(-5988/TEMP)' 
		rate_values[130] = 9.50e+13*numpy.exp(-5988/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROPRIM*O2' 
		rate_values[131] = KROPRIM*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[132] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.53e-13*numpy.exp(1300/TEMP)' 
		rate_values[133] = 1.53e-13*numpy.exp(1300/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*(KCH3O2*7.8e-14*numpy.exp(1000/TEMP))**(0.5)*RO2*0.6' 
		rate_values[134] = 2*(KCH3O2*7.8e-14*numpy.exp(1000/TEMP))**(0.5)*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*(KCH3O2*7.8e-14*numpy.exp(1000/TEMP))**(0.5)*RO2*0.2' 
		rate_values[135] = 2*(KCH3O2*7.8e-14*numpy.exp(1000/TEMP))**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*(KCH3O2*7.8e-14*numpy.exp(1000/TEMP))**(0.5)*RO2*0.2' 
		rate_values[136] = 2*(KCH3O2*7.8e-14*numpy.exp(1000/TEMP))**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.40e-13' 
		rate_values[137] = 8.40e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[138] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[139] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.38e-11' 
		rate_values[140] = 1.38e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.800' 
		rate_values[141] = 1.00e-11*0.800
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.200' 
		rate_values[142] = 1.00e-11*0.200
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[143] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL' 
		rate_values[144] = KNO3AL
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[145] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[146] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[147] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[148] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[149] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[150] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.7*RO2' 
		rate_values[151] = 1.00e-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.3*RO2' 
		rate_values[152] = 1.00e-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.73e-12' 
		rate_values[153] = 2.73e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[154] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.19e-12' 
		rate_values[155] = 6.19e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6e-12*numpy.exp(305/TEMP)' 
		rate_values[156] = 1.6e-12*numpy.exp(305/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[157] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[158] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.08e-11' 
		rate_values[159] = 3.08e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*4.0' 
		rate_values[160] = KNO3AL*4.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[161] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[162] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[163] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[164] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[165] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.82e-12' 
		rate_values[166] = 3.82e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[167] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[168] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.38e-11' 
		rate_values[169] = 1.38e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[170] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8e-13*numpy.exp(780/TEMP)*(1/(1+498*numpy.exp(-1160/TEMP)))' 
		rate_values[171] = 3.8e-13*numpy.exp(780/TEMP)*(1/(1+498*numpy.exp(-1160/TEMP)))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.7e-12*numpy.exp(345/TEMP)*0.05' 
		rate_values[172] = 4.7e-12*numpy.exp(345/TEMP)*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.36e-13*numpy.exp(1250/TEMP)*0.15' 
		rate_values[173] = 1.36e-13*numpy.exp(1250/TEMP)*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.0e-12*numpy.exp(20/TEMP)*0.05' 
		rate_values[174] = 3.0e-12*numpy.exp(20/TEMP)*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.18' 
		rate_values[175] = KDEC*0.18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.57' 
		rate_values[176] = KDEC*0.57
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.125' 
		rate_values[177] = KDEC*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.125' 
		rate_values[178] = KDEC*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[31]' 
		rate_values[179] = J[31]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[32]' 
		rate_values[180] = J[32]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[33]' 
		rate_values[181] = J[33]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.1e-12*numpy.exp(340/TEMP)' 
		rate_values[182] = 3.1e-12*numpy.exp(340/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL' 
		rate_values[183] = KNO3AL
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[184] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[185] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.9e-12*numpy.exp(575/TEMP)' 
		rate_values[186] = 1.9e-12*numpy.exp(575/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*2.4' 
		rate_values[187] = KNO3AL*2.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[188] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[189] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[190] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[191] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.7*RO2' 
		rate_values[192] = 1.00e-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.3*RO2' 
		rate_values[193] = 1.00e-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.58e-11' 
		rate_values[194] = 1.58e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[195] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[196] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.23e-11' 
		rate_values[197] = 1.23e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.983' 
		rate_values[198] = KRO2NO*0.983
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.017' 
		rate_values[199] = KRO2NO*0.017
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[200] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[201] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[202] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.6*RO2' 
		rate_values[203] = 2.00e-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.2*RO2' 
		rate_values[204] = 2.00e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.2*RO2' 
		rate_values[205] = 2.00e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.23e-12' 
		rate_values[206] = 2.23e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.77e-11' 
		rate_values[207] = 5.77e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[208] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.88e-11' 
		rate_values[209] = 1.88e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[210] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[211] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.45e-11' 
		rate_values[212] = 2.45e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*4.0' 
		rate_values[213] = KNO3AL*4.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[214] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[215] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[216] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[217] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[218] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[219] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.74e-12' 
		rate_values[220] = 3.74e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.34e-12' 
		rate_values[221] = 7.34e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[222] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[223] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.11' 
		rate_values[224] = KDEC*0.11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.89' 
		rate_values[225] = KDEC*0.89
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.0e-14' 
		rate_values[226] = 7.0e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-15' 
		rate_values[227] = 1.2e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0e-14' 
		rate_values[228] = 1.0e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0e-15' 
		rate_values[229] = 1.0e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0e-18*H2O' 
		rate_values[230] = 6.0e-18*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.0e-14' 
		rate_values[231] = 7.0e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-15' 
		rate_values[232] = 1.2e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0e-14' 
		rate_values[233] = 1.0e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0e-15' 
		rate_values[234] = 1.0e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0e-18*H2O' 
		rate_values[235] = 6.0e-18*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0e-17*H2O' 
		rate_values[236] = 1.0e-17*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[237] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.0e-13' 
		rate_values[238] = 8.0e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.11' 
		rate_values[239] = KDEC*0.11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.89' 
		rate_values[240] = KDEC*0.89
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.0e-14' 
		rate_values[241] = 7.0e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-15' 
		rate_values[242] = 1.2e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0e-14' 
		rate_values[243] = 1.0e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0e-15' 
		rate_values[244] = 1.0e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0e-18*H2O' 
		rate_values[245] = 6.0e-18*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0e-17*H2O' 
		rate_values[246] = 1.0e-17*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.4e-12*numpy.exp(135/TEMP)' 
		rate_values[247] = 5.4e-12*numpy.exp(135/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[248] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[12]' 
		rate_values[249] = J[12]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5e-16' 
		rate_values[250] = 5.5e-16
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.7e-12*numpy.exp(345/TEMP)*0.95' 
		rate_values[251] = 4.7e-12*numpy.exp(345/TEMP)*0.95
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4e-12*numpy.exp(-1860/TEMP)' 
		rate_values[252] = 1.4e-12*numpy.exp(-1860/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.5e-12*numpy.exp(290/TEMP)' 
		rate_values[253] = 7.5e-12*numpy.exp(290/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[254] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[255] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0e-12' 
		rate_values[256] = 4.0e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[257] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[258] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.7*RO2' 
		rate_values[259] = 1.00e-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.3*RO2' 
		rate_values[260] = 1.00e-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e-14' 
		rate_values[261] = 3e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[262] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.70e-12' 
		rate_values[263] = 3.70e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]' 
		rate_values[264] = J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9e-12*numpy.exp(405/TEMP)' 
		rate_values[265] = 4.9e-12*numpy.exp(405/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.24e-12*numpy.exp(-1860/TEMP)' 
		rate_values[266] = 3.24e-12*numpy.exp(-1860/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.7e-12*numpy.exp(340/TEMP)' 
		rate_values[267] = 6.7e-12*numpy.exp(340/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[268] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPPN' 
		rate_values[269] = KBPPN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[270] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[271] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[272] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.7*RO2' 
		rate_values[273] = 1.00e-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.3*RO2' 
		rate_values[274] = 1.00e-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.27e-12' 
		rate_values[275] = 1.27e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[276] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.42e-12' 
		rate_values[277] = 4.42e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[14]' 
		rate_values[278] = J[14]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[279] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[280] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[281] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.387' 
		rate_values[282] = KRO2HO2*0.387
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.6*RO2' 
		rate_values[283] = 2.00e-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.2*RO2' 
		rate_values[284] = 2.00e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*0.2*RO2' 
		rate_values[285] = 2.00e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[286] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[287] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[288] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.91e-11' 
		rate_values[289] = 2.91e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[21]' 
		rate_values[290] = J[21]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.8e-12*numpy.exp(-1320/TEMP)+1.7e-14*numpy.exp(423/TEMP)' 
		rate_values[291] = 8.8e-12*numpy.exp(-1320/TEMP)+1.7e-14*numpy.exp(423/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[292] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[293] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[294] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.36e-13*numpy.exp(1250/TEMP)*0.85' 
		rate_values[295] = 1.36e-13*numpy.exp(1250/TEMP)*0.85
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*(K298CH3O2*8.0e-12)**(0.5)*RO2*0.6' 
		rate_values[296] = 2*(K298CH3O2*8.0e-12)**(0.5)*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*(K298CH3O2*8.0e-12)**(0.5)*RO2*0.2' 
		rate_values[297] = 2*(K298CH3O2*8.0e-12)**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*(K298CH3O2*8.0e-12)**(0.5)*RO2*0.2' 
		rate_values[298] = 2*(K298CH3O2*8.0e-12)**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[299] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.39e-12' 
		rate_values[300] = 8.39e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[301] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[302] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[303] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[304] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[305] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[306] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.6*RO2' 
		rate_values[307] = 8.80e-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[308] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[309] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[310] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[311] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.88e-11' 
		rate_values[312] = 1.88e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[313] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[314] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[315] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[316] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[317] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[318] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[319] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[320] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[321] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[322] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[323] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2' 
		rate_values[324] = 2.00e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.19e-11' 
		rate_values[325] = 7.19e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.55e-11' 
		rate_values[326] = 7.55e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[327] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[328] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[329] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.33e-11' 
		rate_values[330] = 8.33e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[331] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[332] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[333] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[334] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[335] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[336] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[337] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[338] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[339] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[340] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[341] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[342] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.520' 
		rate_values[343] = KRO2HO2*0.520
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2' 
		rate_values[344] = 2.00e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.27e-11' 
		rate_values[345] = 1.27e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.63e-11' 
		rate_values[346] = 1.63e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[347] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[348] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.41e-11' 
		rate_values[349] = 2.41e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[350] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[351] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.85e-12*numpy.exp(-345/TEMP)' 
		rate_values[352] = 2.85e-12*numpy.exp(-345/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.0e-12*numpy.exp(20/TEMP)*0.9' 
		rate_values[353] = 3.0e-12*numpy.exp(20/TEMP)*0.9
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.0e-12*numpy.exp(20/TEMP)*0.05' 
		rate_values[354] = 3.0e-12*numpy.exp(20/TEMP)*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.86e-12' 
		rate_values[355] = 5.86e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[356] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[357] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[358] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.39e-11' 
		rate_values[359] = 3.39e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KNO3AL*4.0' 
		rate_values[360] = 2*KNO3AL*4.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.04e-13' 
		rate_values[361] = 2.04e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[362] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[363] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[364] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[365] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[366] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2' 
		rate_values[367] = 2.00e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.59e-12' 
		rate_values[368] = 9.59e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[369] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[370] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.45e-11' 
		rate_values[371] = 1.45e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[372] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[373] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[374] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[375] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[376] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[377] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-11' 
		rate_values[378] = 2.00e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2' 
		rate_values[379] = 2.00e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.82e-11' 
		rate_values[380] = 9.82e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[381] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[382] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[383] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[384] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[385] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.59e-12' 
		rate_values[386] = 3.59e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2' 
		rate_values[387] = 2.00e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[388] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.09e-13' 
		rate_values[389] = 1.09e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[390] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[391] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[392] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[393] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.00e-13*RO2' 
		rate_values[394] = 6.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.19e-11' 
		rate_values[395] = 1.19e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[396] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[397] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.09e-11' 
		rate_values[398] = 7.09e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[399] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[400] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[401] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[402] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[403] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[404] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[405] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[406] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-14' 
		rate_values[407] = 1.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.59e-12' 
		rate_values[408] = 3.59e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[409] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13' 
		rate_values[410] = 8.00e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-12' 
		rate_values[411] = 1.2e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.44e-11' 
		rate_values[412] = 1.44e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[413] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[414] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-11*numpy.exp(440/TEMP)*0.572' 
		rate_values[415] = 1.2e-11*numpy.exp(440/TEMP)*0.572
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-11*numpy.exp(440/TEMP)*0.353' 
		rate_values[416] = 1.2e-11*numpy.exp(440/TEMP)*0.353
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-11*numpy.exp(440/TEMP)*0.075' 
		rate_values[417] = 1.2e-11*numpy.exp(440/TEMP)*0.075
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.770' 
		rate_values[418] = KRO2NO*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.230' 
		rate_values[419] = KRO2NO*0.230
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[420] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[421] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[422] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.7' 
		rate_values[423] = 9.20e-14*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.3' 
		rate_values[424] = 9.20e-14*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.50e-12' 
		rate_values[425] = 5.50e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[426] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.83e-11' 
		rate_values[427] = 1.83e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[428] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2e-12*numpy.exp(600/TEMP)*0.772' 
		rate_values[429] = 5.2e-12*numpy.exp(600/TEMP)*0.772
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2e-12*numpy.exp(600/TEMP)*0.228' 
		rate_values[430] = 5.2e-12*numpy.exp(600/TEMP)*0.228
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0e-14' 
		rate_values[431] = 2.0e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[432] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.843' 
		rate_values[433] = KRO2NO*0.843
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.157' 
		rate_values[434] = KRO2NO*0.157
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.20e+10*numpy.exp(-3523/TEMP)' 
		rate_values[435] = 4.20e+10*numpy.exp(-3523/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[436] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[437] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.882' 
		rate_values[438] = KRO2NO*0.882
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.118' 
		rate_values[439] = KRO2NO*0.118
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[440] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.902' 
		rate_values[441] = KRO2NO*0.902
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.098' 
		rate_values[442] = KRO2NO*0.098
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[443] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.950' 
		rate_values[444] = KRO2NO*0.950
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.050' 
		rate_values[445] = KRO2NO*0.050
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[446] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.875' 
		rate_values[447] = KRO2NO*0.875
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.125' 
		rate_values[448] = KRO2NO*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[449] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[450] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[451] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[452] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[453] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[454] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[455] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[456] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.7*RO2' 
		rate_values[457] = 1.00e-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.3*RO2' 
		rate_values[458] = 1.00e-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[459] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[460] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.6*RO2' 
		rate_values[461] = 1.30e-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.2*RO2' 
		rate_values[462] = 1.30e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.2*RO2' 
		rate_values[463] = 1.30e-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[464] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[465] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.7*RO2' 
		rate_values[466] = 6.70e-15*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.3*RO2' 
		rate_values[467] = 6.70e-15*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[468] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[469] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.7*RO2' 
		rate_values[470] = 6.70e-15*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.3*RO2' 
		rate_values[471] = 6.70e-15*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[472] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[473] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.6*RO2' 
		rate_values[474] = 8.80e-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[475] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[476] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[477] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[478] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.7*RO2' 
		rate_values[479] = 6.70e-15*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.3*RO2' 
		rate_values[480] = 6.70e-15*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[481] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[482] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.7*RO2' 
		rate_values[483] = 6.70e-15*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.3*RO2' 
		rate_values[484] = 6.70e-15*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[485] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[486] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.6*RO2' 
		rate_values[487] = 8.80e-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[488] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[489] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[53]' 
		rate_values[490] = J[53]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[491] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.88e-12' 
		rate_values[492] = 2.88e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]' 
		rate_values[493] = J[55]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[494] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.37e-12' 
		rate_values[495] = 5.37e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[496] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '71.11e-12' 
		rate_values[497] = 71.11e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]' 
		rate_values[498] = J[55]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[499] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.25e-11' 
		rate_values[500] = 2.25e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]' 
		rate_values[501] = J[55]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[502] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.03e-11' 
		rate_values[503] = 7.03e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.66e-12' 
		rate_values[504] = 3.66e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[505] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[506] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.73e-12' 
		rate_values[507] = 9.73e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[508] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[509] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[510] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-11' 
		rate_values[511] = 1.30e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[512] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[513] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.05e-11' 
		rate_values[514] = 1.05e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[515] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[516] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.05e-11' 
		rate_values[517] = 2.05e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[518] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[519] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.69e-11' 
		rate_values[520] = 8.69e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[521] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[522] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.75e-11' 
		rate_values[523] = 2.75e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[524] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[525] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.01e-11' 
		rate_values[526] = 8.01e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[527] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[528] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.20e-10' 
		rate_values[529] = 1.20e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[530] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.65e-12' 
		rate_values[531] = 6.65e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[532] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.67e-12' 
		rate_values[533] = 7.67e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[534] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.20e-12' 
		rate_values[535] = 7.20e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[536] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.69e-11' 
		rate_values[537] = 1.69e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[538] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.78e-11' 
		rate_values[539] = 3.78e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[540] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.41e-11' 
		rate_values[541] = 2.41e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[542] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.66e-11' 
		rate_values[543] = 7.66e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[544] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.92e-11*0.232' 
		rate_values[545] = 8.92e-11*0.232
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.92e-11*0.768' 
		rate_values[546] = 8.92e-11*0.768
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[547] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[548] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[549] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[550] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[551] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[552] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[553] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[554] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[555] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[556] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[557] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[558] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.29e-11' 
		rate_values[559] = 2.29e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[560] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]*2' 
		rate_values[561] = J[22]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.70e-11' 
		rate_values[562] = 2.70e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]*2' 
		rate_values[563] = J[22]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.39e-11' 
		rate_values[564] = 2.39e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[565] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]*2' 
		rate_values[566] = J[22]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.23e-11' 
		rate_values[567] = 3.23e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]*2' 
		rate_values[568] = J[22]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.54e-11*0.890' 
		rate_values[569] = 2.54e-11*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.54e-11*0.110' 
		rate_values[570] = 2.54e-11*0.110
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[571] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.55e-11' 
		rate_values[572] = 3.55e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.60e-12' 
		rate_values[573] = 7.60e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[574] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[575] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.65e-11' 
		rate_values[576] = 2.65e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[577] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.34e-11' 
		rate_values[578] = 2.34e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[579] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[580] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[581] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[582] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[583] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[584] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[585] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[586] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[587] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2*0.6' 
		rate_values[588] = 2.00e-12*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2*0.2' 
		rate_values[589] = 2.00e-12*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2*0.2' 
		rate_values[590] = 2.00e-12*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[591] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.64e-11' 
		rate_values[592] = 2.64e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*8.5' 
		rate_values[593] = KNO3AL*8.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[594] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[595] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[596] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[597] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[598] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[599] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[600] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[601] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.22e-12' 
		rate_values[602] = 3.22e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[603] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.33e-11' 
		rate_values[604] = 1.33e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[605] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[606] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-11' 
		rate_values[607] = 6.70e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[608] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[609] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[610] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[611] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[612] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[613] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[614] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[615] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[616] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[617] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[618] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2' 
		rate_values[619] = 6.70e-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[620] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[621] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2' 
		rate_values[622] = 6.70e-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[623] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[624] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[625] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[626] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[627] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[628] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[629] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[630] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[631] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[632] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2' 
		rate_values[633] = 2.00e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.60e-12' 
		rate_values[634] = 6.60e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[635] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[636] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.02e-11' 
		rate_values[637] = 1.02e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[638] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[639] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.29e-11' 
		rate_values[640] = 1.29e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[641] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[642] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.45e-11' 
		rate_values[643] = 3.45e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[644] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[645] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.49e-11' 
		rate_values[646] = 7.49e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.83e-13' 
		rate_values[647] = 8.83e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[648] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[649] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.75e-12' 
		rate_values[650] = 4.75e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[651] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[652] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.01e-11' 
		rate_values[653] = 1.01e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[654] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.33e-10' 
		rate_values[655] = 1.33e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KNO3AL*5.5' 
		rate_values[656] = 2*KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.770' 
		rate_values[657] = KRO2NO*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.230' 
		rate_values[658] = KRO2NO*0.230
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[659] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[660] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[661] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[662] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[663] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[664] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.64e-12' 
		rate_values[665] = 3.64e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[666] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.28e-11' 
		rate_values[667] = 3.28e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.49e-11' 
		rate_values[668] = 1.49e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.18e-12' 
		rate_values[669] = 8.18e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-12*numpy.exp(490/TEMP)*0.65' 
		rate_values[670] = 1.2e-12*numpy.exp(490/TEMP)*0.65
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e-12*numpy.exp(490/TEMP)*0.35' 
		rate_values[671] = 1.2e-12*numpy.exp(490/TEMP)*0.35
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[672] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[673] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[674] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[675] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.9*RO2' 
		rate_values[676] = 6.70e-15*0.9*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.1*RO2' 
		rate_values[677] = 6.70e-15*0.1*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[678] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.87e-12' 
		rate_values[679] = 6.87e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[680] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.00e+05' 
		rate_values[681] = 4.00e+05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROSEC*O2' 
		rate_values[682] = KROSEC*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[683] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[684] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*0.8*RO2' 
		rate_values[685] = 2.50e-13*0.8*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*0.1*RO2' 
		rate_values[686] = 2.50e-13*0.1*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*0.1*RO2' 
		rate_values[687] = 2.50e-13*0.1*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[688] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[689] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.23e-11' 
		rate_values[690] = 1.23e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]' 
		rate_values[691] = J[55]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.55e-12' 
		rate_values[692] = 5.55e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[693] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[694] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[695] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[696] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[697] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.00e+04' 
		rate_values[698] = 4.00e+04
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROSEC*O2' 
		rate_values[699] = KROSEC*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[700] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[701] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2' 
		rate_values[702] = 6.70e-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[703] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[704] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2' 
		rate_values[705] = 6.70e-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[706] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[707] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[708] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[709] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.94e-12' 
		rate_values[710] = 5.94e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[711] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.03e-12' 
		rate_values[712] = 8.03e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[713] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.25e-11' 
		rate_values[714] = 3.25e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]' 
		rate_values[715] = J[55]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.25e-12' 
		rate_values[716] = 1.25e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[717] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[718] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[719] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[720] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[721] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[722] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[723] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[724] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[725] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[726] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[727] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[728] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.85e-11' 
		rate_values[729] = 1.85e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.32e-11' 
		rate_values[730] = 1.32e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[731] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.68e-11' 
		rate_values[732] = 1.68e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.875' 
		rate_values[733] = KRO2NO*0.875
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.125' 
		rate_values[734] = KRO2NO*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[735] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[736] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[737] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2*0.7' 
		rate_values[738] = 6.70e-15*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2*0.3' 
		rate_values[739] = 6.70e-15*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]' 
		rate_values[740] = J[55]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.87e-11' 
		rate_values[741] = 9.87e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[742] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.03e-10' 
		rate_values[743] = 1.03e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.91e-11' 
		rate_values[744] = 9.91e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.722' 
		rate_values[745] = KRO2NO*0.722
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.278' 
		rate_values[746] = KRO2NO*0.278
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROSEC*O2' 
		rate_values[747] = KROSEC*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[748] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[749] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2*0.6' 
		rate_values[750] = 2.50e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2*0.2' 
		rate_values[751] = 2.50e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2*0.2' 
		rate_values[752] = 2.50e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[753] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.60e-11' 
		rate_values[754] = 9.60e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[755] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.27e-10' 
		rate_values[756] = 1.27e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.09e-10' 
		rate_values[757] = 1.09e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.19e-10' 
		rate_values[758] = 1.19e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.958' 
		rate_values[759] = KRO2NO*0.958
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.042' 
		rate_values[760] = KRO2NO*0.042
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[761] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[762] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[763] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.7' 
		rate_values[764] = 9.20e-14*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.3' 
		rate_values[765] = 9.20e-14*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.26e-11' 
		rate_values[766] = 1.26e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[767] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.06e-11' 
		rate_values[768] = 7.06e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.72e-11' 
		rate_values[769] = 6.72e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.05e-16*numpy.exp(-640/TEMP)*0.6' 
		rate_values[770] = 8.05e-16*numpy.exp(-640/TEMP)*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.05e-16*numpy.exp(-640/TEMP)*0.4' 
		rate_values[771] = 8.05e-16*numpy.exp(-640/TEMP)*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.55' 
		rate_values[772] = KDEC*0.55
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.45' 
		rate_values[773] = KDEC*0.45
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[774] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[775] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.875' 
		rate_values[776] = KRO2NO*0.875
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.125' 
		rate_values[777] = KRO2NO*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[778] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.722' 
		rate_values[779] = KRO2NO*0.722
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.278' 
		rate_values[780] = KRO2NO*0.278
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROSEC*O2' 
		rate_values[781] = KROSEC*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[782] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[783] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*0.7*RO2' 
		rate_values[784] = 9.20e-14*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*0.3*RO2' 
		rate_values[785] = 9.20e-14*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[786] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[787] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.7*RO2' 
		rate_values[788] = 6.70e-15*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*0.3*RO2' 
		rate_values[789] = 6.70e-15*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[790] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[791] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*0.6*RO2' 
		rate_values[792] = 2.50e-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*0.2*RO2' 
		rate_values[793] = 2.50e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*0.2*RO2' 
		rate_values[794] = 2.50e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]' 
		rate_values[795] = J[55]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[796] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.85e-11' 
		rate_values[797] = 2.85e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[798] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[799] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.23e-11' 
		rate_values[800] = 2.23e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[801] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[802] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.01e-11' 
		rate_values[803] = 3.01e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[804] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[805] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.28e-11' 
		rate_values[806] = 6.28e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[807] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[808] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-10' 
		rate_values[809] = 2.00e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[810] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.66e-11' 
		rate_values[811] = 2.66e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[812] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.93e-11' 
		rate_values[813] = 5.93e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[814] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.26e-10' 
		rate_values[815] = 1.26e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[816] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.80' 
		rate_values[817] = KDEC*0.80
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.20' 
		rate_values[818] = KDEC*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[819] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[820] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2*0.90' 
		rate_values[821] = 2.00e-12*RO2*0.90
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2*0.05' 
		rate_values[822] = 2.00e-12*RO2*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2*0.05' 
		rate_values[823] = 2.00e-12*RO2*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[824] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[825] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.20e+10*numpy.exp(-3523/TEMP)' 
		rate_values[826] = 4.20e+10*numpy.exp(-3523/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[827] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[828] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[829] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[830] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[831] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[832] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[833] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[834] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.37e-11' 
		rate_values[835] = 2.37e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[836] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[837] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[838] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[839] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[840] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[841] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[842] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[843] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[844] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.92e-12' 
		rate_values[845] = 2.92e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[846] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[847] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.55e-12' 
		rate_values[848] = 6.55e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[849] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[850] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2' 
		rate_values[851] = 2.00e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[852] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[853] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.61e-12' 
		rate_values[854] = 9.61e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[855] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[856] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[857] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[858] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[859] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[860] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[861] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[862] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[863] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*RO2' 
		rate_values[864] = 1.30e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[865] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[866] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2' 
		rate_values[867] = 6.70e-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[868] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[869] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2' 
		rate_values[870] = 6.70e-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[871] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[872] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[873] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[874] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[875] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[876] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.47e-11' 
		rate_values[877] = 5.47e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[878] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[879] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.45e-11' 
		rate_values[880] = 4.45e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[881] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[882] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.47e-11' 
		rate_values[883] = 5.47e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.56e-12' 
		rate_values[884] = 5.56e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[885] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.70e-12' 
		rate_values[886] = 5.70e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[887] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[888] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.16e-12' 
		rate_values[889] = 9.16e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[890] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[891] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.36e-11' 
		rate_values[892] = 2.36e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[893] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[894] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.29e-11' 
		rate_values[895] = 1.29e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[896] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[897] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.51e-11' 
		rate_values[898] = 1.51e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[899] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[900] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-10' 
		rate_values[901] = 1.00e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.50' 
		rate_values[902] = KDEC*0.50
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.50' 
		rate_values[903] = KDEC*0.50
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00e-14' 
		rate_values[904] = 7.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.20e-15' 
		rate_values[905] = 1.20e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-14' 
		rate_values[906] = 1.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-15' 
		rate_values[907] = 1.00e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.40e-17*H2O' 
		rate_values[908] = 1.40e-17*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-18*H2O' 
		rate_values[909] = 2.00e-18*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[910] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[911] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.3e-11*0.352' 
		rate_values[912] = 1.3e-11*0.352
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.3e-11*0.118' 
		rate_values[913] = 1.3e-11*0.118
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.3e-11*0.53' 
		rate_values[914] = 1.3e-11*0.53
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.918' 
		rate_values[915] = KRO2NO*0.918
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.082' 
		rate_values[916] = KRO2NO*0.082
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.50' 
		rate_values[917] = KDEC*0.50
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.50' 
		rate_values[918] = KDEC*0.50
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[919] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[920] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[921] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[922] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[923] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[924] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.30e-11' 
		rate_values[925] = 7.30e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[926] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.77e-11' 
		rate_values[927] = 9.77e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.21e-10' 
		rate_values[928] = 1.21e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[929] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.16e-11' 
		rate_values[930] = 8.16e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.20e-11*0.17' 
		rate_values[931] = 5.20e-11*0.17
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.20e-11*0.83' 
		rate_values[932] = 5.20e-11*0.83
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.14*0.4' 
		rate_values[933] = J[4]*0.14*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.14*0.6' 
		rate_values[934] = J[4]*0.14*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[935] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[936] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[937] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.60*RO2' 
		rate_values[938] = 8.80e-13*0.60*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.20*RO2' 
		rate_values[939] = 8.80e-13*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.20*RO2' 
		rate_values[940] = 8.80e-13*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[941] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[942] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[943] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[944] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[945] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[946] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[947] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.70*RO2' 
		rate_values[948] = 1.00e-11*0.70*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.30*RO2' 
		rate_values[949] = 1.00e-11*0.30*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[950] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[951] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[952] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.22e-10' 
		rate_values[953] = 1.22e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[954] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[955] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.67e-11' 
		rate_values[956] = 3.67e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]*2' 
		rate_values[957] = J[34]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.45e-11' 
		rate_values[958] = 2.45e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[959] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[960] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[961] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[962] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[963] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[964] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.97e-11' 
		rate_values[965] = 6.97e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[966] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.33e-11' 
		rate_values[967] = 7.33e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[968] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.13e-11' 
		rate_values[969] = 8.13e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[970] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.36e-10' 
		rate_values[971] = 1.36e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-18' 
		rate_values[972] = 2.00e-18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4e-12' 
		rate_values[973] = 1.4e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[974] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[975] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[976] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[977] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.6*RO2' 
		rate_values[978] = 8.80e-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[979] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[980] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[981] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.66e-11' 
		rate_values[982] = 4.66e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.55e-11' 
		rate_values[983] = 2.55e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.68e-12' 
		rate_values[984] = 5.68e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-18' 
		rate_values[985] = 2.00e-18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KNO3AL*2.0' 
		rate_values[986] = 2*KNO3AL*2.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.1*0.5' 
		rate_values[987] = J[4]*0.1*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.1*0.5' 
		rate_values[988] = J[4]*0.1*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KNO3AL*2.75' 
		rate_values[989] = 2*KNO3AL*2.75
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.08e-11*0.69' 
		rate_values[990] = 6.08e-11*0.69
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.08e-11*0.31' 
		rate_values[991] = 6.08e-11*0.31
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.895' 
		rate_values[992] = KRO2NO*0.895
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.105' 
		rate_values[993] = KRO2NO*0.105
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[994] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[995] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[996] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.77' 
		rate_values[997] = KRO2HO2*0.77
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.6*RO2' 
		rate_values[998] = 8.80e-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[999] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[1000] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.38e-11' 
		rate_values[1001] = 4.38e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1002] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[1003] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.31e-10' 
		rate_values[1004] = 1.31e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[1005] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.23e-11' 
		rate_values[1006] = 8.23e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[1007] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[1008] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1009] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1010] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1011] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.520' 
		rate_values[1012] = KRO2HO2*0.520
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[1013] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[1014] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[1015] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1016] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[1017] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.44e-10' 
		rate_values[1018] = 1.44e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[1019] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.77e-11' 
		rate_values[1020] = 5.77e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1021] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1022] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1023] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1024] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1025] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1026] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[1027] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[1028] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[1029] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[1030] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1031] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[1032] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[1033] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[1034] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[1035] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.05e-11' 
		rate_values[1036] = 4.05e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1037] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[1038] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[1039] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.37e-11' 
		rate_values[1040] = 4.37e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[1041] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[1042] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.06e-11' 
		rate_values[1043] = 4.06e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1044] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[1045] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[1046] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.52e-11' 
		rate_values[1047] = 7.52e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[1048] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[1049] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.75e-11' 
		rate_values[1050] = 7.75e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[1051] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[1052] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[1053] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.90e-11' 
		rate_values[1054] = 4.90e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]*2' 
		rate_values[1055] = J[17]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.32e-11' 
		rate_values[1056] = 4.32e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KNO3AL*4.0' 
		rate_values[1057] = 2*KNO3AL*4.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1058] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1059] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1060] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1061] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[1062] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[1063] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[1064] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[1065] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.29e-11' 
		rate_values[1066] = 2.29e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1067] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[1068] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.62e-11' 
		rate_values[1069] = 2.62e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[1070] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.31e-11' 
		rate_values[1071] = 2.31e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.6e-12' 
		rate_values[1072] = 4.6e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1073] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1074] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1075] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[1076] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[1077] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[1078] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[1079] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1080] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.23e-10' 
		rate_values[1081] = 1.23e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.18e-11' 
		rate_values[1082] = 9.18e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.07e-11' 
		rate_values[1083] = 6.07e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.45e-11' 
		rate_values[1084] = 4.45e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1085] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1086] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1087] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[1088] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[1089] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[1090] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[1091] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1092] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.68e-11' 
		rate_values[1093] = 3.68e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.78e-11' 
		rate_values[1094] = 2.78e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.78e-11' 
		rate_values[1095] = 1.78e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[1096] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.44e-11' 
		rate_values[1097] = 3.44e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*8.0' 
		rate_values[1098] = KNO3AL*8.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.7e-13*numpy.exp(1220/TEMP)*0.06' 
		rate_values[1099] = 4.7e-13*numpy.exp(1220/TEMP)*0.06
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.7e-13*numpy.exp(1220/TEMP)*0.14' 
		rate_values[1100] = 4.7e-13*numpy.exp(1220/TEMP)*0.14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.7e-13*numpy.exp(1220/TEMP)*0.8' 
		rate_values[1101] = 4.7e-13*numpy.exp(1220/TEMP)*0.8
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8e-12*0.742' 
		rate_values[1102] = 3.8e-12*0.742
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8e-12*0.258' 
		rate_values[1103] = 3.8e-12*0.258
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.29' 
		rate_values[1104] = KDEC*0.29
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[1105] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1106] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.00e-14' 
		rate_values[1107] = 9.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[1108] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[1109] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1110] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1111] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[1112] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[1113] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1114] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.00e-13' 
		rate_values[1115] = 9.00e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0e-10' 
		rate_values[1116] = 1.0e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.9e-11' 
		rate_values[1117] = 9.9e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.2e-18' 
		rate_values[1118] = 9.2e-18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1119] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[1120] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[1121] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1122] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1123] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[1124] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[1125] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1126] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1127] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.47e-12' 
		rate_values[1128] = 3.47e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.60e-12' 
		rate_values[1129] = 2.60e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1130] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1131] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1132] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[1133] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[1134] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1135] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1136] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1137] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1138] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1139] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[1140] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[1141] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1142] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1143] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1144] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1145] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.00e-14' 
		rate_values[1146] = 3.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.25e-15' 
		rate_values[1147] = 2.25e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1148] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1149] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1150] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[1151] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[1152] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1153] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1154] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1155] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1156] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1157] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[1158] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[1159] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1160] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1161] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1162] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.71' 
		rate_values[1163] = KDEC*0.71
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1164] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[1165] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.7' 
		rate_values[1166] = 8.00e-13*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.3' 
		rate_values[1167] = 8.00e-13*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1168] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.16e-10' 
		rate_values[1169] = 1.16e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.13e-10' 
		rate_values[1170] = 1.13e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1171] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1172] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1173] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[1174] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.7' 
		rate_values[1175] = 8.00e-13*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.3' 
		rate_values[1176] = 8.00e-13*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1177] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.07e-10' 
		rate_values[1178] = 1.07e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.04e-10' 
		rate_values[1179] = 1.04e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[1180] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.895' 
		rate_values[1181] = KRO2NO*0.895
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.105' 
		rate_values[1182] = KRO2NO*0.105
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1183] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.5e-13*numpy.exp(1310/TEMP)' 
		rate_values[1184] = 1.5e-13*numpy.exp(1310/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*(KCH3O2*2.4e-14*numpy.exp(1620/TEMP))**(0.5)*RO2*0.6' 
		rate_values[1185] = 2*(KCH3O2*2.4e-14*numpy.exp(1620/TEMP))**(0.5)*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*(KCH3O2*2.4e-14*numpy.exp(1620/TEMP))**(0.5)*RO2*0.2' 
		rate_values[1186] = 2*(KCH3O2*2.4e-14*numpy.exp(1620/TEMP))**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*(KCH3O2*2.4e-14*numpy.exp(1620/TEMP))**(0.5)*RO2*0.2' 
		rate_values[1187] = 2*(KCH3O2*2.4e-14*numpy.exp(1620/TEMP))**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROPRIM*O2' 
		rate_values[1188] = KROPRIM*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1189] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.05e-11' 
		rate_values[1190] = 2.05e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.7e-11' 
		rate_values[1191] = 2.7e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[53]' 
		rate_values[1192] = J[53]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.03e-12' 
		rate_values[1193] = 6.03e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[1194] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.99e-11' 
		rate_values[1195] = 7.99e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.2*0.4' 
		rate_values[1196] = J[4]*0.2*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.2*0.6' 
		rate_values[1197] = J[4]*0.2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.35' 
		rate_values[1198] = KDEC*0.35
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2e-11*0.52' 
		rate_values[1199] = 6.2e-11*0.52
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2e-11*0.48' 
		rate_values[1200] = 6.2e-11*0.48
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1201] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1202] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1203] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1204] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[1205] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[1206] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.70*RO2' 
		rate_values[1207] = 1.00e-11*0.70*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.30*RO2' 
		rate_values[1208] = 1.00e-11*0.30*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1209] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1210] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[1211] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.40e-12*0.60*RO2' 
		rate_values[1212] = 1.40e-12*0.60*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.40e-12*0.20*RO2' 
		rate_values[1213] = 1.40e-12*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.40e-12*0.20*RO2' 
		rate_values[1214] = 1.40e-12*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1215] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1216] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[1217] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[1218] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1219] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.21e-10' 
		rate_values[1220] = 1.21e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[1221] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[1222] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.48e-11' 
		rate_values[1223] = 7.48e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[1224] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[1225] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.26e-11' 
		rate_values[1226] = 6.26e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.59e-12' 
		rate_values[1227] = 3.59e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[20]*2' 
		rate_values[1228] = J[20]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.44e-11' 
		rate_values[1229] = 5.44e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[1230] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[1231] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.43e-11' 
		rate_values[1232] = 5.43e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-18*0.5' 
		rate_values[1233] = 2.00e-18*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-18*0.5' 
		rate_values[1234] = 2.00e-18*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.18' 
		rate_values[1235] = KDEC*0.18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.57' 
		rate_values[1236] = KDEC*0.57
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.125' 
		rate_values[1237] = KDEC*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.125' 
		rate_values[1238] = KDEC*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*2.75' 
		rate_values[1239] = KNO3AL*2.75
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[1240] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]*2' 
		rate_values[1241] = J[35]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.72e-11' 
		rate_values[1242] = 1.72e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[1243] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[1244] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.17' 
		rate_values[1245] = KDEC*0.17
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.83' 
		rate_values[1246] = KDEC*0.83
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.35' 
		rate_values[1247] = KDEC*0.35
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.65' 
		rate_values[1248] = KDEC*0.65
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.65' 
		rate_values[1249] = KDEC*0.65
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.50e-12' 
		rate_values[1250] = 1.50e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1251] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1252] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1253] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[1254] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.70' 
		rate_values[1255] = 9.20e-14*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.30' 
		rate_values[1256] = 9.20e-14*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1257] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.69e-11' 
		rate_values[1258] = 1.69e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.34e-11' 
		rate_values[1259] = 1.34e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.2*0.3' 
		rate_values[1260] = J[4]*0.2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.2*0.4' 
		rate_values[1261] = J[4]*0.2*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.2*0.3' 
		rate_values[1262] = J[4]*0.2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.41e-11*0.23' 
		rate_values[1263] = 4.41e-11*0.23
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.41e-11*0.385' 
		rate_values[1264] = 4.41e-11*0.385
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.41e-11*0.385' 
		rate_values[1265] = 4.41e-11*0.385
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*4.25' 
		rate_values[1266] = KNO3AL*4.25
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*4.25' 
		rate_values[1267] = KNO3AL*4.25
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1268] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1269] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[1270] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.70' 
		rate_values[1271] = 9.20e-14*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.30' 
		rate_values[1272] = 9.20e-14*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1273] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1274] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]*2' 
		rate_values[1275] = J[17]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.58e-11' 
		rate_values[1276] = 4.58e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]*2' 
		rate_values[1277] = J[17]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.70e-11' 
		rate_values[1278] = 5.70e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1279] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1280] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1281] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1282] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[1283] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[1284] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.70*RO2' 
		rate_values[1285] = 1.00e-11*0.70*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.30*RO2' 
		rate_values[1286] = 1.00e-11*0.30*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.37e-11' 
		rate_values[1287] = 4.37e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[1288] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[1289] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.38e-11' 
		rate_values[1290] = 4.38e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1291] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1292] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1293] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1294] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[1295] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[1296] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.37e-11' 
		rate_values[1297] = 4.37e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00e-18*0.50' 
		rate_values[1298] = 5.00e-18*0.50
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00e-18*0.50' 
		rate_values[1299] = 5.00e-18*0.50
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.90e-11' 
		rate_values[1300] = 6.90e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1301] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1302] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1303] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[1304] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.7' 
		rate_values[1305] = 9.20e-14*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.3' 
		rate_values[1306] = 9.20e-14*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1307] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.53e-11' 
		rate_values[1308] = 2.53e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.19e-11' 
		rate_values[1309] = 2.19e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[1310] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.13e-13*1.0' 
		rate_values[1311] = 3.13e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[1312] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[1313] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1314] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1315] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[1316] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[1317] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1318] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.53e-12' 
		rate_values[1319] = 1.53e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.10e-14' 
		rate_values[1320] = 5.10e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.83e-15' 
		rate_values[1321] = 7.83e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1322] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1323] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1324] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[1325] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[1326] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1327] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1328] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1329] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1330] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1331] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[1332] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[1333] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1334] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1335] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[1336] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.8e-12' 
		rate_values[1337] = 2.8e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[1338] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[1339] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.70e-11' 
		rate_values[1340] = 3.70e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[20]*2' 
		rate_values[1341] = J[20]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.00e-11' 
		rate_values[1342] = 4.00e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1343] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[1344] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.69e-12' 
		rate_values[1345] = 3.69e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1346] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1347] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[1348] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[1349] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1350] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1351] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.02e-13' 
		rate_values[1352] = 1.02e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.12e-12' 
		rate_values[1353] = 1.12e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1354] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1355] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1356] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1357] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[1358] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[1359] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.66e-11' 
		rate_values[1360] = 7.66e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.01e-11' 
		rate_values[1361] = 8.01e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1362] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-11' 
		rate_values[1363] = 9.20e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.31e-11*0.50' 
		rate_values[1364] = 2.31e-11*0.50
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.31e-11*0.29' 
		rate_values[1365] = 2.31e-11*0.29
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.31e-11*0.17' 
		rate_values[1366] = 2.31e-11*0.17
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.31e-11*0.04' 
		rate_values[1367] = 2.31e-11*0.04
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.60e-16' 
		rate_values[1368] = 2.60e-16
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.895' 
		rate_values[1369] = KRO2NO*0.895
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.105' 
		rate_values[1370] = KRO2NO*0.105
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1371] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1372] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.60*RO2' 
		rate_values[1373] = 1.30e-12*0.60*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.20*RO2' 
		rate_values[1374] = 1.30e-12*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.20*RO2' 
		rate_values[1375] = 1.30e-12*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROPRIM*O2' 
		rate_values[1376] = KROPRIM*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1377] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.05e-11' 
		rate_values[1378] = 2.05e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.93e-12' 
		rate_values[1379] = 8.93e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[53]' 
		rate_values[1380] = J[53]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.03e-12' 
		rate_values[1381] = 6.03e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[1382] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[1383] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL' 
		rate_values[1384] = KNO3AL
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.29e-11' 
		rate_values[1385] = 1.29e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1386] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1387] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1388] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1389] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.35' 
		rate_values[1390] = KDEC*0.35
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.65' 
		rate_values[1391] = KDEC*0.65
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1392] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN*0.67' 
		rate_values[1393] = KBPAN*0.67
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[1394] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[1395] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.70*RO2' 
		rate_values[1396] = 1.00e-11*0.70*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.30*RO2' 
		rate_values[1397] = 1.00e-11*0.30*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[1398] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[1399] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.66e-12' 
		rate_values[1400] = 4.66e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1401] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.10e-12' 
		rate_values[1402] = 1.10e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1403] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1404] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.862' 
		rate_values[1405] = KRO2NO*0.862
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.138' 
		rate_values[1406] = KRO2NO*0.138
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.26' 
		rate_values[1407] = KDEC*0.26
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.11' 
		rate_values[1408] = KDEC*0.11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.26' 
		rate_values[1409] = KDEC*0.26
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.11' 
		rate_values[1410] = KDEC*0.11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.26' 
		rate_values[1411] = KDEC*0.26
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1412] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1413] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[1414] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[1415] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[1416] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[1417] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.15e-11' 
		rate_values[1418] = 7.15e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1419] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.17e-10' 
		rate_values[1420] = 1.17e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-10' 
		rate_values[1421] = 2.50e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]*0.5' 
		rate_values[1422] = J[22]*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]*0.5' 
		rate_values[1423] = J[22]*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.93e-11*0.5' 
		rate_values[1424] = 7.93e-11*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.93e-11*0.5' 
		rate_values[1425] = 7.93e-11*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.2*0.6' 
		rate_values[1426] = J[4]*0.2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.2*0.4' 
		rate_values[1427] = J[4]*0.2*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00e-11*0.5' 
		rate_values[1428] = 7.00e-11*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00e-11*0.5' 
		rate_values[1429] = 7.00e-11*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*4.25' 
		rate_values[1430] = KNO3AL*4.25
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1431] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1432] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1433] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1434] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[1435] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.70*RO2' 
		rate_values[1436] = 1.00e-11*0.70*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.30*RO2' 
		rate_values[1437] = 1.00e-11*0.30*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.54e-11' 
		rate_values[1438] = 3.54e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[20]*2' 
		rate_values[1439] = J[20]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.90e-11' 
		rate_values[1440] = 3.90e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1441] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1442] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[1443] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*0.70*RO2' 
		rate_values[1444] = 9.20e-14*0.70*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*0.30*RO2' 
		rate_values[1445] = 9.20e-14*0.30*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1446] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1447] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[1448] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[1449] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.46e-11' 
		rate_values[1450] = 4.46e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0e-18*0.5' 
		rate_values[1451] = 2.0e-18*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0e-18*0.5' 
		rate_values[1452] = 2.0e-18*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.11' 
		rate_values[1453] = KDEC*0.11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.89' 
		rate_values[1454] = KDEC*0.89
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[1455] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[1456] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.13e-11' 
		rate_values[1457] = 4.13e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]*2' 
		rate_values[1458] = J[22]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.06e-12' 
		rate_values[1459] = 4.06e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[1460] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.66e-11' 
		rate_values[1461] = 6.66e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1462] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1463] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1464] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1465] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[1466] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[1467] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.99e-11' 
		rate_values[1468] = 6.99e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1469] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.64e-11' 
		rate_values[1470] = 6.64e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1471] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1472] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1473] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1474] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[1475] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[1476] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.66e-13' 
		rate_values[1477] = 5.66e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[1478] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1479] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[1480] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.16e-12' 
		rate_values[1481] = 4.16e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.028*0.5' 
		rate_values[1482] = J[4]*0.028*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.028*0.5' 
		rate_values[1483] = J[4]*0.028*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.02e-11*0.74' 
		rate_values[1484] = 8.02e-11*0.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.02e-11*0.26' 
		rate_values[1485] = 8.02e-11*0.26
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*4.25' 
		rate_values[1486] = KNO3AL*4.25
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.827' 
		rate_values[1487] = KRO2NO*0.827
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.173' 
		rate_values[1488] = KRO2NO*0.173
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[1489] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[1490] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1491] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1492] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[1493] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[1494] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[1495] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.08e-11' 
		rate_values[1496] = 3.08e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1497] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[1498] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[1499] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.41e-11' 
		rate_values[1500] = 7.41e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[1501] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[1502] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.20e-11' 
		rate_values[1503] = 6.20e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[1504] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[1505] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.08e-11' 
		rate_values[1506] = 4.08e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1507] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1508] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1509] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1510] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[1511] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[1512] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.70' 
		rate_values[1513] = 1.00e-11*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.30' 
		rate_values[1514] = 1.00e-11*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1515] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1516] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1517] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[1518] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[1519] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[1520] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[1521] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.99e-11' 
		rate_values[1522] = 5.99e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1523] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[1524] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[1525] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.31e-11' 
		rate_values[1526] = 6.31e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[1527] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[1528] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.00e-11' 
		rate_values[1529] = 6.00e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1530] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[1531] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[1532] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.36e-11' 
		rate_values[1533] = 7.36e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[1534] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[1535] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.76e-11' 
		rate_values[1536] = 9.76e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[1537] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[1538] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[1539] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.22e-11' 
		rate_values[1540] = 5.22e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]*2' 
		rate_values[1541] = J[17]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.28e-11' 
		rate_values[1542] = 4.28e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KNO3AL*8.5' 
		rate_values[1543] = 2*KNO3AL*8.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1544] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1545] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1546] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1547] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[1548] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[1549] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.70' 
		rate_values[1550] = 1.00e-11*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.30' 
		rate_values[1551] = 1.00e-11*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.25e-11' 
		rate_values[1552] = 2.25e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1553] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[1554] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.57e-11' 
		rate_values[1555] = 2.57e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[1556] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.27e-11' 
		rate_values[1557] = 2.27e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1558] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1559] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1560] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[1561] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.7' 
		rate_values[1562] = 9.20e-14*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.3' 
		rate_values[1563] = 9.20e-14*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1564] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]*2' 
		rate_values[1565] = J[17]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.35e-10' 
		rate_values[1566] = 1.35e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.35e-11' 
		rate_values[1567] = 2.35e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1568] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1569] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1570] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1571] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.7' 
		rate_values[1572] = 9.20e-14*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.3' 
		rate_values[1573] = 9.20e-14*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1574] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.48e-11' 
		rate_values[1575] = 9.48e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.13e-11' 
		rate_values[1576] = 9.13e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1577] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1578] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1579] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1580] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[1581] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[1582] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.42e-11' 
		rate_values[1583] = 2.42e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.75e-11' 
		rate_values[1584] = 2.75e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1585] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.15e-11' 
		rate_values[1586] = 2.15e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[1587] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[1588] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.10e-11*0.07' 
		rate_values[1589] = 9.10e-11*0.07
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.10e-11*0.51' 
		rate_values[1590] = 9.10e-11*0.51
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.10e-11*0.42' 
		rate_values[1591] = 9.10e-11*0.42
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.50e-11*0.39' 
		rate_values[1592] = 4.50e-11*0.39
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.50e-11*0.51' 
		rate_values[1593] = 4.50e-11*0.51
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.50e-11*0.10' 
		rate_values[1594] = 4.50e-11*0.10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[1595] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1596] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1597] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1598] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[1599] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1600] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.65e-11' 
		rate_values[1601] = 4.65e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.15' 
		rate_values[1602] = KDEC*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.05e-10' 
		rate_values[1603] = 2.05e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.01e-10' 
		rate_values[1604] = 2.01e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.81e-17' 
		rate_values[1605] = 2.81e-17
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[1606] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[1607] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.83e-12' 
		rate_values[1608] = 6.83e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.03e-12' 
		rate_values[1609] = 5.03e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1610] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1611] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1612] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1613] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[1614] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1615] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1616] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1617] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1618] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1619] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1620] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[1621] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1622] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1623] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1624] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1625] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1626] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[1627] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1628] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.05e-10' 
		rate_values[1629] = 2.05e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1630] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.13e-13' 
		rate_values[1631] = 3.13e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[1632] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[1633] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1634] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1635] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1636] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[1637] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1638] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.53e-12' 
		rate_values[1639] = 1.53e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.10e-14' 
		rate_values[1640] = 5.10e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.83e-15' 
		rate_values[1641] = 7.83e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1642] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1643] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1644] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1645] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[1646] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1647] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1648] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1649] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1650] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1651] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1652] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[1653] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1654] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1655] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1656] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.85' 
		rate_values[1657] = KDEC*0.85
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1658] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1659] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.70' 
		rate_values[1660] = 8.00e-13*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.30' 
		rate_values[1661] = 8.00e-13*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1662] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.13e-10' 
		rate_values[1663] = 1.13e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.10e-10' 
		rate_values[1664] = 1.10e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1665] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1666] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1667] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1668] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.70' 
		rate_values[1669] = 8.00e-13*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.30' 
		rate_values[1670] = 8.00e-13*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1671] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.07e-10' 
		rate_values[1672] = 1.07e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.03e-10' 
		rate_values[1673] = 1.03e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[1674] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.53e-12' 
		rate_values[1675] = 1.53e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[1676] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.54e-11' 
		rate_values[1677] = 3.54e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.67e-11' 
		rate_values[1678] = 2.67e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1679] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1680] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1681] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[1682] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.70' 
		rate_values[1683] = 9.20e-14*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.30' 
		rate_values[1684] = 9.20e-14*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1685] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.09e-11' 
		rate_values[1686] = 3.09e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.74e-11' 
		rate_values[1687] = 2.74e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[1688] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.68e-11' 
		rate_values[1689] = 2.68e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*8.5' 
		rate_values[1690] = KNO3AL*8.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1691] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1692] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1693] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1694] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1695] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[1696] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[1697] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[1698] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[1699] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1700] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[1701] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[1702] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[1703] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[1704] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.46e-12' 
		rate_values[1705] = 6.46e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1706] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.69e-12' 
		rate_values[1707] = 9.69e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.58e-12' 
		rate_values[1708] = 6.58e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1709] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.01e-11' 
		rate_values[1710] = 6.01e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1711] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.54e-11' 
		rate_values[1712] = 2.54e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1713] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-18' 
		rate_values[1714] = 2.00e-18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.40' 
		rate_values[1715] = KDEC*0.40
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.60' 
		rate_values[1716] = KDEC*0.60
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.00e-13' 
		rate_values[1717] = 3.00e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1718] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[1719] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[1720] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.16e-12' 
		rate_values[1721] = 1.16e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1722] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[1723] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[1724] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1725] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.18e-12' 
		rate_values[1726] = 6.18e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.20e-19' 
		rate_values[1727] = 2.20e-19
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[1728] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[1729] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00e-14' 
		rate_values[1730] = 7.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.20e-15' 
		rate_values[1731] = 1.20e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-14' 
		rate_values[1732] = 1.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-15' 
		rate_values[1733] = 1.00e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.00e-18*H2O' 
		rate_values[1734] = 6.00e-18*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-17*H2O' 
		rate_values[1735] = 1.00e-17*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.19e-11' 
		rate_values[1736] = 2.19e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.00e-13' 
		rate_values[1737] = 3.00e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1738] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1739] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1740] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[1741] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[1742] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1743] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.68e-11' 
		rate_values[1744] = 6.68e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]*2' 
		rate_values[1745] = J[34]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.70e-11' 
		rate_values[1746] = 7.70e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1747] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1748] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1749] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1750] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[1751] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[1752] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1753] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[1754] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.06e-11' 
		rate_values[1755] = 3.06e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-12' 
		rate_values[1756] = 1.00e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1757] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1758] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1759] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[1760] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2' 
		rate_values[1761] = 9.20e-14*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1762] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.44e-12' 
		rate_values[1763] = 4.44e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-19' 
		rate_values[1764] = 8.00e-19
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[1765] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[1766] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00e-14' 
		rate_values[1767] = 7.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.20e-15' 
		rate_values[1768] = 1.20e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-14' 
		rate_values[1769] = 1.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-15' 
		rate_values[1770] = 1.00e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.00e-18*H2O' 
		rate_values[1771] = 6.00e-18*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-17*H2O' 
		rate_values[1772] = 1.00e-17*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.76e-13' 
		rate_values[1773] = 8.76e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1774] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1775] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1776] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1777] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[1778] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[1779] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1780] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[1781] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.60e-11' 
		rate_values[1782] = 3.60e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-12' 
		rate_values[1783] = 1.00e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1784] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1785] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1786] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1787] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2' 
		rate_values[1788] = 9.20e-14*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1789] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.54e-11' 
		rate_values[1790] = 7.54e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[1791] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[1792] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.04e-11' 
		rate_values[1793] = 7.04e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-12' 
		rate_values[1794] = 1.00e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1795] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1796] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1797] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[1798] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2' 
		rate_values[1799] = 9.20e-14*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1800] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.80e-12' 
		rate_values[1801] = 5.80e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-19' 
		rate_values[1802] = 8.00e-19
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.90e-13' 
		rate_values[1803] = 3.90e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.43e-11*0.625' 
		rate_values[1804] = 1.43e-11*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.43e-11*0.155' 
		rate_values[1805] = 1.43e-11*0.155
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.43e-11*0.12' 
		rate_values[1806] = 1.43e-11*0.12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.43e-11*0.10' 
		rate_values[1807] = 1.43e-11*0.10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00e-16' 
		rate_values[1808] = 5.00e-16
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.895' 
		rate_values[1809] = KRO2NO*0.895
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.105' 
		rate_values[1810] = KRO2NO*0.105
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1811] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[1812] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.60*RO2' 
		rate_values[1813] = 1.30e-12*0.60*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.20*RO2' 
		rate_values[1814] = 1.30e-12*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.20*RO2' 
		rate_values[1815] = 1.30e-12*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROPRIM*O2' 
		rate_values[1816] = KROPRIM*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1817] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.05e-11' 
		rate_values[1818] = 2.05e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.93e-12' 
		rate_values[1819] = 8.93e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[53]' 
		rate_values[1820] = J[53]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.03e-12' 
		rate_values[1821] = 6.03e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[1822] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[1823] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL' 
		rate_values[1824] = KNO3AL
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.29e-11' 
		rate_values[1825] = 1.29e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1826] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1827] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1828] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1829] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.90e-13' 
		rate_values[1830] = 3.90e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1831] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN*0.67' 
		rate_values[1832] = KBPAN*0.67
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[1833] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[1834] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.70*RO2' 
		rate_values[1835] = 1.00e-11*0.70*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.30*RO2' 
		rate_values[1836] = 1.00e-11*0.30*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[1837] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[1838] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1839] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.66e-12' 
		rate_values[1840] = 4.66e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.10e-12' 
		rate_values[1841] = 1.10e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1842] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1843] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.6e-12' 
		rate_values[1844] = 3.6e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1845] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1846] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1847] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[1848] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2' 
		rate_values[1849] = 9.20e-14*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1850] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.52e-10' 
		rate_values[1851] = 1.52e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.862' 
		rate_values[1852] = KRO2NO*0.862
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.138' 
		rate_values[1853] = KRO2NO*0.138
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.18' 
		rate_values[1854] = KDEC*0.18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.64' 
		rate_values[1855] = KDEC*0.64
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.18' 
		rate_values[1856] = KDEC*0.18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1857] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1858] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.70' 
		rate_values[1859] = 9.20e-14*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.30' 
		rate_values[1860] = 9.20e-14*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]' 
		rate_values[1861] = J[55]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.16e-11' 
		rate_values[1862] = 7.16e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1863] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.65e-11' 
		rate_values[1864] = 7.65e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.12e-10' 
		rate_values[1865] = 1.12e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.0057' 
		rate_values[1866] = J[4]*0.0057
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.90e-11' 
		rate_values[1867] = 4.90e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1868] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1869] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[1870] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.60*RO2' 
		rate_values[1871] = 8.80e-13*0.60*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.20*RO2' 
		rate_values[1872] = 8.80e-13*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.20*RO2' 
		rate_values[1873] = 8.80e-13*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.3' 
		rate_values[1874] = KDEC*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.7' 
		rate_values[1875] = KDEC*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1876] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]*2' 
		rate_values[1877] = J[22]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[1878] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.75e-11' 
		rate_values[1879] = 9.75e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[1880] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[1881] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.04e-12' 
		rate_values[1882] = 9.04e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]*2' 
		rate_values[1883] = J[22]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.00e-11' 
		rate_values[1884] = 4.00e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00e-18' 
		rate_values[1885] = 5.00e-18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.1*0.5' 
		rate_values[1886] = J[4]*0.1*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.1*0.5' 
		rate_values[1887] = J[4]*0.1*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.88e-11*0.74' 
		rate_values[1888] = 7.88e-11*0.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.88e-11*0.26' 
		rate_values[1889] = 7.88e-11*0.26
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*4.25' 
		rate_values[1890] = KNO3AL*4.25
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00e-18' 
		rate_values[1891] = 5.00e-18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.827' 
		rate_values[1892] = KRO2NO*0.827
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.173' 
		rate_values[1893] = KRO2NO*0.173
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[1894] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[1895] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1896] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1897] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-12*RO2*0.6' 
		rate_values[1898] = 8.80e-12*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-12*RO2*0.2' 
		rate_values[1899] = 8.80e-12*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-12*RO2*0.2' 
		rate_values[1900] = 8.80e-12*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1901] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[1902] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[1903] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1904] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[1905] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-12*RO2' 
		rate_values[1906] = 8.80e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1907] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[1908] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[1909] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.83e-11' 
		rate_values[1910] = 7.83e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1911] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[1912] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[1913] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.26e-11' 
		rate_values[1914] = 7.26e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[1915] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[1916] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.03e-11' 
		rate_values[1917] = 6.03e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[1918] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[1919] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.58e-11' 
		rate_values[1920] = 3.58e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1921] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1922] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1923] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1924] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[1925] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[1926] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.70' 
		rate_values[1927] = 1.00e-11*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.30' 
		rate_values[1928] = 1.00e-11*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1929] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1930] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1931] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[1932] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.70' 
		rate_values[1933] = 9.20e-14*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.30' 
		rate_values[1934] = 9.20e-14*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.85e-11' 
		rate_values[1935] = 5.85e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1936] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[1937] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[1938] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.18e-11' 
		rate_values[1939] = 6.18e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[1940] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[1941] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.87e-11' 
		rate_values[1942] = 5.87e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1943] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[1944] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[1945] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[1946] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.25e-11' 
		rate_values[1947] = 8.25e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[1948] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[1949] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[1950] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.14e-10' 
		rate_values[1951] = 1.14e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]*2' 
		rate_values[1952] = J[17]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.18e-11' 
		rate_values[1953] = 4.18e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.35e-11' 
		rate_values[1954] = 2.35e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1955] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1956] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1957] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1958] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.70' 
		rate_values[1959] = 9.20e-14*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.30' 
		rate_values[1960] = 9.20e-14*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1961] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.48e-11' 
		rate_values[1962] = 9.48e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.13e-11' 
		rate_values[1963] = 9.13e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1964] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1965] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1966] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1967] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[1968] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0e-11*RO2' 
		rate_values[1969] = 1.0e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.42e-11' 
		rate_values[1970] = 2.42e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1971] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.75e-11' 
		rate_values[1972] = 2.75e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-12' 
		rate_values[1973] = 1.00e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[1974] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[1975] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[1976] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[1977] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2' 
		rate_values[1978] = 9.20e-14*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1979] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.54e-11' 
		rate_values[1980] = 7.54e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[1981] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[1982] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.35e-11' 
		rate_values[1983] = 8.35e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[1984] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[1985] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[1986] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[1987] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[1988] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[1989] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[1990] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[1991] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.60e-11' 
		rate_values[1992] = 3.60e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.28e-11' 
		rate_values[1993] = 3.28e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.70e-11' 
		rate_values[1994] = 3.70e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.06e-12' 
		rate_values[1995] = 1.06e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.74e-11' 
		rate_values[1996] = 2.74e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.28e-11' 
		rate_values[1997] = 3.28e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.06e-12' 
		rate_values[1998] = 1.06e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-11*0.07' 
		rate_values[1999] = 8.00e-11*0.07
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-11*0.2' 
		rate_values[2000] = 8.00e-11*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-11*0.73' 
		rate_values[2001] = 8.00e-11*0.73
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.48e-11*0.39' 
		rate_values[2002] = 3.48e-11*0.39
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.48e-11*0.51' 
		rate_values[2003] = 3.48e-11*0.51
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.48e-11*0.10' 
		rate_values[2004] = 3.48e-11*0.10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[2005] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2006] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2007] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2008] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[2009] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2010] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.65e-11' 
		rate_values[2011] = 4.65e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.32' 
		rate_values[2012] = KDEC*0.32
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.05e-10' 
		rate_values[2013] = 2.05e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.01e-10' 
		rate_values[2014] = 2.01e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.81e-17' 
		rate_values[2015] = 2.81e-17
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[2016] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[2017] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.83e-12' 
		rate_values[2018] = 6.83e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.03e-12' 
		rate_values[2019] = 5.03e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2020] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2021] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2022] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2023] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[2024] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2025] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[2026] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2027] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2028] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2029] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2030] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[2031] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2032] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[2033] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2034] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2035] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2036] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[2037] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2038] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.05e-10' 
		rate_values[2039] = 2.05e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2040] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.13e-13' 
		rate_values[2041] = 3.13e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[2042] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[2043] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2044] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2045] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2046] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[2047] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2048] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.53e-12' 
		rate_values[2049] = 1.53e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.10e-14' 
		rate_values[2050] = 5.10e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.83e-15' 
		rate_values[2051] = 7.83e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2052] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2053] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2054] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2055] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[2056] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2057] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[2058] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2059] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2060] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2061] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2062] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[2063] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2064] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[2065] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[2066] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2067] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2068] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.68' 
		rate_values[2069] = KDEC*0.68
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[2070] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[2071] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.66e-11' 
		rate_values[2072] = 4.66e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2073] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2074] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.70' 
		rate_values[2075] = 8.00e-13*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.30' 
		rate_values[2076] = 8.00e-13*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2077] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.51e-10' 
		rate_values[2078] = 1.51e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.48e-10' 
		rate_values[2079] = 1.48e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2080] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2081] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2082] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2083] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.70' 
		rate_values[2084] = 8.00e-13*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.30' 
		rate_values[2085] = 8.00e-13*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2086] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.43e-10' 
		rate_values[2087] = 1.43e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.40e-10' 
		rate_values[2088] = 1.40e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[2089] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.53e-12' 
		rate_values[2090] = 1.53e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.42e-11' 
		rate_values[2091] = 2.42e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2092] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2093] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2094] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[2095] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.70' 
		rate_values[2096] = 9.20e-14*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.30' 
		rate_values[2097] = 9.20e-14*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2098] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.78e-11' 
		rate_values[2099] = 2.78e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.44e-11' 
		rate_values[2100] = 2.44e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[2101] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.15e-11' 
		rate_values[2102] = 2.15e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[2103] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[2104] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[2105] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2106] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[2107] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[2108] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[2109] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[2110] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.36e-13' 
		rate_values[2111] = 7.36e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2112] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.34e-12' 
		rate_values[2113] = 4.34e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.76e-13' 
		rate_values[2114] = 8.76e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2115] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.27e-11' 
		rate_values[2116] = 1.27e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*4.0' 
		rate_values[2117] = KNO3AL*4.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[2118] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-12' 
		rate_values[2119] = 1.00e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2120] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2121] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2122] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[2123] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2' 
		rate_values[2124] = 9.20e-14*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2125] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.16e-12' 
		rate_values[2126] = 5.16e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-19' 
		rate_values[2127] = 8.00e-19
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.36e-11*0.55' 
		rate_values[2128] = 1.36e-11*0.55
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.36e-11*0.24' 
		rate_values[2129] = 1.36e-11*0.24
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.36e-11*0.16' 
		rate_values[2130] = 1.36e-11*0.16
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.36e-11*0.05' 
		rate_values[2131] = 1.36e-11*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.10e-16' 
		rate_values[2132] = 4.10e-16
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.895' 
		rate_values[2133] = KRO2NO*0.895
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.105' 
		rate_values[2134] = KRO2NO*0.105
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2135] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2136] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.60*RO2' 
		rate_values[2137] = 1.30e-12*0.60*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.20*RO2' 
		rate_values[2138] = 1.30e-12*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.20*RO2' 
		rate_values[2139] = 1.30e-12*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROPRIM*O2' 
		rate_values[2140] = KROPRIM*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2141] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.75e-11' 
		rate_values[2142] = 2.75e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.60e-11' 
		rate_values[2143] = 1.60e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[53]' 
		rate_values[2144] = J[53]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.31e-11' 
		rate_values[2145] = 1.31e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[2146] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[2147] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL' 
		rate_values[2148] = KNO3AL
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.29e-11' 
		rate_values[2149] = 1.29e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[2150] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2151] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2152] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2153] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[2154] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN*0.67' 
		rate_values[2155] = KBPAN*0.67
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[2156] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[2157] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.70*RO2' 
		rate_values[2158] = 1.00e-11*0.70*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.30*RO2' 
		rate_values[2159] = 1.00e-11*0.30*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[2160] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[2161] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.66e-12' 
		rate_values[2162] = 4.66e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2163] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.10e-12' 
		rate_values[2164] = 1.10e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.65e-11' 
		rate_values[2165] = 4.65e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2166] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.862' 
		rate_values[2167] = KRO2NO*0.862
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.138' 
		rate_values[2168] = KRO2NO*0.138
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.245' 
		rate_values[2169] = KDEC*0.245
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.175' 
		rate_values[2170] = KDEC*0.175
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.08' 
		rate_values[2171] = KDEC*0.08
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.175' 
		rate_values[2172] = KDEC*0.175
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.08' 
		rate_values[2173] = KDEC*0.08
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.245' 
		rate_values[2174] = KDEC*0.245
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2175] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2176] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.70' 
		rate_values[2177] = 9.20e-14*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.30' 
		rate_values[2178] = 9.20e-14*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]' 
		rate_values[2179] = J[55]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.28e-11' 
		rate_values[2180] = 7.28e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2181] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.06e-11' 
		rate_values[2182] = 8.06e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2183] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.62e-11' 
		rate_values[2184] = 7.62e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.11e-10' 
		rate_values[2185] = 1.11e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.06e-12' 
		rate_values[2186] = 1.06e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.05*0.5' 
		rate_values[2187] = J[4]*0.05*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.05*0.5' 
		rate_values[2188] = J[4]*0.05*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.99e-11*0.368' 
		rate_values[2189] = 5.99e-11*0.368
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.99e-11*0.632' 
		rate_values[2190] = 5.99e-11*0.632
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KNO3AL*4.25' 
		rate_values[2191] = 2*KNO3AL*4.25
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-18' 
		rate_values[2192] = 2.00e-18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.827' 
		rate_values[2193] = KRO2NO*0.827
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.173' 
		rate_values[2194] = KRO2NO*0.173
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[2195] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[2196] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2197] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2198] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-12*RO2*0.60' 
		rate_values[2199] = 8.80e-12*RO2*0.60
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-12*RO2*0.20' 
		rate_values[2200] = 8.80e-12*RO2*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-12*RO2*0.20' 
		rate_values[2201] = 8.80e-12*RO2*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2202] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[2203] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[2204] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2205] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[2206] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2' 
		rate_values[2207] = 9.20e-14*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2208] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2209] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2210] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.97e-11' 
		rate_values[2211] = 6.97e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.92e-11' 
		rate_values[2212] = 2.92e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.27e-11' 
		rate_values[2213] = 3.27e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2214] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]*2' 
		rate_values[2215] = J[17]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.33e-11' 
		rate_values[2216] = 9.33e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]*2' 
		rate_values[2217] = J[17]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.11e-11' 
		rate_values[2218] = 8.11e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2219] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[2220] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.84e-11' 
		rate_values[2221] = 4.84e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[2222] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[2223] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[2224] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2225] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[2226] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[2227] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.70' 
		rate_values[2228] = 1.00e-11*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.30' 
		rate_values[2229] = 1.00e-11*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2230] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2231] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2232] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[2233] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.70' 
		rate_values[2234] = 9.20e-14*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.30' 
		rate_values[2235] = 9.20e-14*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.96e-11' 
		rate_values[2236] = 3.96e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2237] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[2238] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[2239] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.28e-11' 
		rate_values[2240] = 4.28e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[2241] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[2242] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.98e-11' 
		rate_values[2243] = 3.98e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2244] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2245] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[2246] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[2247] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.28e-11' 
		rate_values[2248] = 4.28e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2249] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[2250] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[2251] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.24e-11' 
		rate_values[2252] = 5.24e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2253] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2254] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.63e-11' 
		rate_values[2255] = 6.63e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]*2' 
		rate_values[2256] = J[17]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.24e-11' 
		rate_values[2257] = 4.24e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2*KNO3AL*8.5' 
		rate_values[2258] = 2*KNO3AL*8.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[2259] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[2260] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[2261] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2262] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[2263] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[2264] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.70' 
		rate_values[2265] = 1.00e-11*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.30' 
		rate_values[2266] = 1.00e-11*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.21e-11' 
		rate_values[2267] = 2.21e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2268] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2269] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.53e-11' 
		rate_values[2270] = 2.53e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2271] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.22e-11' 
		rate_values[2272] = 2.22e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.2*0.4' 
		rate_values[2273] = J[4]*0.2*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.2*0.6' 
		rate_values[2274] = J[4]*0.2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.10e-11*0.72' 
		rate_values[2275] = 7.10e-11*0.72
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.10e-11*0.28' 
		rate_values[2276] = 7.10e-11*0.28
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2277] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2278] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2279] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[2280] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*0.7*RO2' 
		rate_values[2281] = 9.20e-14*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*0.3*RO2' 
		rate_values[2282] = 9.20e-14*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2283] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2284] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2285] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.12e-11' 
		rate_values[2286] = 3.12e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2287] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2288] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.77e-11' 
		rate_values[2289] = 2.77e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[2290] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2291] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[2292] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[2293] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[2294] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[2295] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.70*RO2' 
		rate_values[2296] = 1.00e-11*0.70*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.30*RO2' 
		rate_values[2297] = 1.00e-11*0.30*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[20]*2' 
		rate_values[2298] = J[20]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.44e-11' 
		rate_values[2299] = 5.44e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.08e-11' 
		rate_values[2300] = 5.08e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[2301] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.08e-11' 
		rate_values[2302] = 5.08e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*4.25' 
		rate_values[2303] = KNO3AL*4.25
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00e-18*0.5' 
		rate_values[2304] = 5.00e-18*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00e-18*0.5' 
		rate_values[2305] = 5.00e-18*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.22' 
		rate_values[2306] = KDEC*0.22
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.27' 
		rate_values[2307] = KDEC*0.27
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00e-14' 
		rate_values[2308] = 7.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.20e-15' 
		rate_values[2309] = 1.20e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-14' 
		rate_values[2310] = 1.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-15' 
		rate_values[2311] = 1.00e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.00e-18*H2O' 
		rate_values[2312] = 6.00e-18*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.65' 
		rate_values[2313] = KDEC*0.65
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.35' 
		rate_values[2314] = KDEC*0.35
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.35e-10' 
		rate_values[2315] = 2.35e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2316] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2317] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2318] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2319] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[2320] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[2321] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[2322] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2323] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.61e-10' 
		rate_values[2324] = 1.61e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-10' 
		rate_values[2325] = 1.30e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.85e-11' 
		rate_values[2326] = 9.85e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[2327] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[2328] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[2329] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2330] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[2331] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[2332] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.51e-11' 
		rate_values[2333] = 9.51e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2334] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.83e-11' 
		rate_values[2335] = 9.83e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2336] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[2337] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[2338] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-12' 
		rate_values[2339] = 1.00e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2340] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2341] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2342] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2343] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[2344] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2345] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.05e-10' 
		rate_values[2346] = 1.05e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]*2' 
		rate_values[2347] = J[34]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.15e-10' 
		rate_values[2348] = 1.15e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[2349] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[2350] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[2351] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2352] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[2353] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[2354] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.82e-11' 
		rate_values[2355] = 3.82e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2356] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[2357] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.15e-11' 
		rate_values[2358] = 4.15e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-11*0.07' 
		rate_values[2359] = 8.00e-11*0.07
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-11*0.20' 
		rate_values[2360] = 8.00e-11*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-11*0.73' 
		rate_values[2361] = 8.00e-11*0.73
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.20e-11*0.39' 
		rate_values[2362] = 3.20e-11*0.39
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.20e-11*0.51' 
		rate_values[2363] = 3.20e-11*0.51
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.20e-11*0.10' 
		rate_values[2364] = 3.20e-11*0.10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[2365] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.48' 
		rate_values[2366] = KDEC*0.48
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.05e-10' 
		rate_values[2367] = 2.05e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.01e-10' 
		rate_values[2368] = 2.01e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.81e-17' 
		rate_values[2369] = 2.81e-17
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[2370] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[2371] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.83e-12' 
		rate_values[2372] = 6.83e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.03e-12' 
		rate_values[2373] = 5.03e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2374] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2375] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2376] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2377] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[2378] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2379] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[2380] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2381] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2382] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2383] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2384] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[2385] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2386] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[2387] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2388] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2389] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2390] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[2391] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2392] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.05e-10' 
		rate_values[2393] = 2.05e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2394] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.13e-13' 
		rate_values[2395] = 3.13e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[2396] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[2397] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2398] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2399] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2400] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[2401] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2402] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.53e-12' 
		rate_values[2403] = 1.53e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.10e-14' 
		rate_values[2404] = 5.10e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.83e-15' 
		rate_values[2405] = 7.83e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2406] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2407] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2408] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2409] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[2410] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2411] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[2412] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2413] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2414] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2415] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2416] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[2417] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2418] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[2419] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2420] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.52' 
		rate_values[2421] = KDEC*0.52
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2422] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2423] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.70' 
		rate_values[2424] = 8.00e-13*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.30' 
		rate_values[2425] = 8.00e-13*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2426] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.51e-10' 
		rate_values[2427] = 1.51e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.48e-10' 
		rate_values[2428] = 1.48e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2429] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2430] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2431] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2432] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.70' 
		rate_values[2433] = 8.00e-13*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.30' 
		rate_values[2434] = 8.00e-13*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2435] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.43e-10' 
		rate_values[2436] = 1.43e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.40e-10' 
		rate_values[2437] = 1.40e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[2438] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.53e-12' 
		rate_values[2439] = 1.53e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.01e-11' 
		rate_values[2440] = 3.01e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2441] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2442] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2443] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[2444] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.70' 
		rate_values[2445] = 9.20e-14*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.30' 
		rate_values[2446] = 9.20e-14*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2447] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.12e-12' 
		rate_values[2448] = 7.12e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.76e-12' 
		rate_values[2449] = 3.76e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2450] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.90e-12' 
		rate_values[2451] = 2.90e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2452] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2453] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2454] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[2455] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.6' 
		rate_values[2456] = 8.80e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[2457] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.2' 
		rate_values[2458] = 8.80e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2459] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.42e-11' 
		rate_values[2460] = 4.42e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[2461] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.91e-11' 
		rate_values[2462] = 1.91e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2463] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.65e-13' 
		rate_values[2464] = 4.65e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2465] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2466] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2467] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-12*RO2' 
		rate_values[2468] = 2.00e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2469] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2470] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.72e-11' 
		rate_values[2471] = 1.72e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-12' 
		rate_values[2472] = 1.00e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2473] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2474] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2475] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[2476] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2' 
		rate_values[2477] = 9.20e-14*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2478] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.95e-12' 
		rate_values[2479] = 4.95e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.40e-18' 
		rate_values[2480] = 1.40e-18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.66e-11' 
		rate_values[2481] = 2.66e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2482] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2483] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2484] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[2485] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.70' 
		rate_values[2486] = 9.20e-14*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.30' 
		rate_values[2487] = 9.20e-14*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2488] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.37e-11' 
		rate_values[2489] = 1.37e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.02e-11' 
		rate_values[2490] = 1.02e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2491] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[2492] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.85e-11' 
		rate_values[2493] = 1.85e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-12' 
		rate_values[2494] = 1.00e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2495] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2496] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2497] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[2498] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2' 
		rate_values[2499] = 9.20e-14*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2500] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.93e-11' 
		rate_values[2501] = 9.93e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-17' 
		rate_values[2502] = 8.00e-17
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2503] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2504] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2505] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[2506] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2507] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.65e-11' 
		rate_values[2508] = 4.65e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2509] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2510] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[2511] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[2512] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[2513] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[2514] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2515] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2516] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]*2' 
		rate_values[2517] = J[17]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[2518] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[2519] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00e-12*0.07' 
		rate_values[2520] = 7.00e-12*0.07
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00e-12*0.65' 
		rate_values[2521] = 7.00e-12*0.65
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00e-12*0.10' 
		rate_values[2522] = 7.00e-12*0.10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00e-12*0.18' 
		rate_values[2523] = 7.00e-12*0.18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.20e-16' 
		rate_values[2524] = 1.20e-16
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.900' 
		rate_values[2525] = KRO2NO*0.900
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.100' 
		rate_values[2526] = KRO2NO*0.100
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2527] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2528] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.60*RO2' 
		rate_values[2529] = 1.30e-12*0.60*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.20*RO2' 
		rate_values[2530] = 1.30e-12*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*0.20*RO2' 
		rate_values[2531] = 1.30e-12*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROPRIM*O2' 
		rate_values[2532] = KROPRIM*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2533] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.17e-11' 
		rate_values[2534] = 2.17e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.01e-11' 
		rate_values[2535] = 1.01e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[53]' 
		rate_values[2536] = J[53]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.18e-12' 
		rate_values[2537] = 7.18e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2538] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL' 
		rate_values[2539] = KNO3AL
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-11' 
		rate_values[2540] = 2.08e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[2541] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[2542] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN*0.67' 
		rate_values[2543] = KBPAN*0.67
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2544] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[2545] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[2546] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.70*RO2' 
		rate_values[2547] = 1.00e-11*0.70*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.30*RO2' 
		rate_values[2548] = 1.00e-11*0.30*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2549] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.66e-12' 
		rate_values[2550] = 4.66e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.10e-12' 
		rate_values[2551] = 1.10e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.06e-12' 
		rate_values[2552] = 1.06e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.90' 
		rate_values[2553] = KRO2NO*0.90
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.10' 
		rate_values[2554] = KRO2NO*0.10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.2' 
		rate_values[2555] = KDEC*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.3' 
		rate_values[2556] = KDEC*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.3' 
		rate_values[2557] = KDEC*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.2' 
		rate_values[2558] = KDEC*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2559] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2560] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.60' 
		rate_values[2561] = 8.80e-13*RO2*0.60
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.20' 
		rate_values[2562] = 8.80e-13*RO2*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.20' 
		rate_values[2563] = 8.80e-13*RO2*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[2564] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.28e-11' 
		rate_values[2565] = 7.28e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2566] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.75e-11' 
		rate_values[2567] = 9.75e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.21e-10' 
		rate_values[2568] = 1.21e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2569] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.11e-11' 
		rate_values[2570] = 8.11e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.2*0.4' 
		rate_values[2571] = J[4]*0.2*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.2*0.6' 
		rate_values[2572] = J[4]*0.2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*2.75' 
		rate_values[2573] = KNO3AL*2.75
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.58e-11*0.5' 
		rate_values[2574] = 5.58e-11*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.58e-11*0.5' 
		rate_values[2575] = 5.58e-11*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[2576] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[2577] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[2578] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2579] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[2580] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[2581] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.70*RO2' 
		rate_values[2582] = 1.00e-11*0.70*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*0.30*RO2' 
		rate_values[2583] = 1.00e-11*0.30*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2584] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2585] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[2586] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.60*RO2' 
		rate_values[2587] = 8.80e-13*0.60*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.20*RO2' 
		rate_values[2588] = 8.80e-13*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.20*RO2' 
		rate_values[2589] = 8.80e-13*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2590] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.73e-11' 
		rate_values[2591] = 3.73e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[20]*2' 
		rate_values[2592] = J[20]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.09e-11' 
		rate_values[2593] = 4.09e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[2594] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.73e-11' 
		rate_values[2595] = 3.73e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2596] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2597] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[2598] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.16e-10' 
		rate_values[2599] = 1.16e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[2600] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2601] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.21e-11' 
		rate_values[2602] = 6.21e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[2603] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.71e-11' 
		rate_values[2604] = 7.71e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]*2' 
		rate_values[2605] = J[35]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[2606] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.14e-11' 
		rate_values[2607] = 3.14e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[2608] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[2609] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[2610] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2611] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[2612] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[2613] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2614] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]*2' 
		rate_values[2615] = J[35]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.09e-12' 
		rate_values[2616] = 6.09e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.49e-12' 
		rate_values[2617] = 2.49e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-18*0.5' 
		rate_values[2618] = 2.00e-18*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00e-18*0.5' 
		rate_values[2619] = 2.00e-18*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.18' 
		rate_values[2620] = KDEC*0.18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.57' 
		rate_values[2621] = KDEC*0.57
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.125' 
		rate_values[2622] = KDEC*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.125' 
		rate_values[2623] = KDEC*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2624] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2625] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.35' 
		rate_values[2626] = KDEC*0.35
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.65' 
		rate_values[2627] = KDEC*0.65
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.29e-11' 
		rate_values[2628] = 2.29e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2629] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2630] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2631] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[2632] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.60' 
		rate_values[2633] = 8.80e-13*RO2*0.60
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.20' 
		rate_values[2634] = 8.80e-13*RO2*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.20' 
		rate_values[2635] = 8.80e-13*RO2*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2636] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.23e-11' 
		rate_values[2637] = 4.23e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.33e-11' 
		rate_values[2638] = 3.33e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.24e-11' 
		rate_values[2639] = 3.24e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[2640] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.35e-11' 
		rate_values[2641] = 4.35e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*8.5' 
		rate_values[2642] = KNO3AL*8.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2643] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2644] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[2645] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.6*RO2' 
		rate_values[2646] = 8.80e-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[2647] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*0.2*RO2' 
		rate_values[2648] = 8.80e-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2649] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[2650] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.95e-11' 
		rate_values[2651] = 3.95e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.09e-11' 
		rate_values[2652] = 7.09e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[2653] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[2654] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[2655] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2656] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[2657] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[2658] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.79e-12' 
		rate_values[2659] = 2.79e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.68e-12' 
		rate_values[2660] = 8.68e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2661] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2662] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-19' 
		rate_values[2663] = 8.00e-19
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-12' 
		rate_values[2664] = 1.00e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2665] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[2666] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[2667] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.73e-12' 
		rate_values[2668] = 6.73e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2669] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[2670] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[2671] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2672] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.18e-11' 
		rate_values[2673] = 8.18e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.35e-11' 
		rate_values[2674] = 2.35e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2675] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2676] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2677] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2678] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.60' 
		rate_values[2679] = 8.80e-13*RO2*0.60
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.20' 
		rate_values[2680] = 8.80e-13*RO2*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.20' 
		rate_values[2681] = 8.80e-13*RO2*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2682] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.46e-10' 
		rate_values[2683] = 1.46e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.15e-10' 
		rate_values[2684] = 1.15e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.32e-11' 
		rate_values[2685] = 8.32e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[2686] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[2687] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[2688] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2689] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[2690] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[2691] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.14e-11' 
		rate_values[2692] = 9.14e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2693] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.46e-11' 
		rate_values[2694] = 9.46e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2695] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2696] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-12' 
		rate_values[2697] = 1.00e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2698] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2699] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2700] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2701] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[2702] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2703] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.74e-11' 
		rate_values[2704] = 8.74e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]*2' 
		rate_values[2705] = J[34]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.76e-11' 
		rate_values[2706] = 9.76e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[2707] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[2708] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[2709] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2710] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[2711] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2' 
		rate_values[2712] = 1.00e-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2713] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[2714] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.78e-11' 
		rate_values[2715] = 3.78e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.45e-11' 
		rate_values[2716] = 3.45e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.028*0.5' 
		rate_values[2717] = J[4]*0.028*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.028*0.5' 
		rate_values[2718] = J[4]*0.028*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*2.75' 
		rate_values[2719] = KNO3AL*2.75
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.10e-11*0.76' 
		rate_values[2720] = 8.10e-11*0.76
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.10e-11*0.26' 
		rate_values[2721] = 8.10e-11*0.26
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.827' 
		rate_values[2722] = KRO2NO*0.827
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.173' 
		rate_values[2723] = KRO2NO*0.173
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[2724] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[2725] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2726] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2727] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.60' 
		rate_values[2728] = 8.80e-13*RO2*0.60
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.20' 
		rate_values[2729] = 8.80e-13*RO2*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.20' 
		rate_values[2730] = 8.80e-13*RO2*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.17e-11' 
		rate_values[2731] = 3.17e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2732] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2733] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2734] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.51e-11' 
		rate_values[2735] = 7.51e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2736] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2737] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.29e-11' 
		rate_values[2738] = 6.29e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2739] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[35]' 
		rate_values[2740] = J[35]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.25e-11' 
		rate_values[2741] = 4.25e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[2742] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[2743] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[2744] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2745] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[2746] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[2747] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.70' 
		rate_values[2748] = 1.00e-11*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.30' 
		rate_values[2749] = 1.00e-11*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.08e-11' 
		rate_values[2750] = 6.08e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2751] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[2752] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[2753] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.40e-11' 
		rate_values[2754] = 6.40e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[2755] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[2756] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.09e-11' 
		rate_values[2757] = 6.09e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2758] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2759] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2760] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[2761] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.60' 
		rate_values[2762] = 8.80e-13*RO2*0.60
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.20' 
		rate_values[2763] = 8.80e-13*RO2*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2*0.20' 
		rate_values[2764] = 8.80e-13*RO2*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2765] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2766] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[2767] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[2768] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.53e-11' 
		rate_values[2769] = 9.53e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2770] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[2771] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[2772] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.19e-11' 
		rate_values[2773] = 1.19e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[2774] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.78e-11' 
		rate_values[2775] = 6.78e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2776] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2777] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[2778] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[2779] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2780] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2781] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[2782] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.73e-11' 
		rate_values[2783] = 5.73e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.10e-11*0.073' 
		rate_values[2784] = 5.10e-11*0.073
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.10e-11*0.200' 
		rate_values[2785] = 5.10e-11*0.200
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.10e-11*0.727' 
		rate_values[2786] = 5.10e-11*0.727
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.25e-11*0.391' 
		rate_values[2787] = 1.25e-11*0.391
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.25e-11*0.506' 
		rate_values[2788] = 1.25e-11*0.506
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.25e-11*0.103' 
		rate_values[2789] = 1.25e-11*0.103
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[2790] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2791] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2792] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2793] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[2794] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2795] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.65e-11' 
		rate_values[2796] = 4.65e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.32' 
		rate_values[2797] = KDEC*0.32
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.05e-10*1.0' 
		rate_values[2798] = 2.05e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.01e-10*1.0' 
		rate_values[2799] = 2.01e-10*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.81e-17' 
		rate_values[2800] = 2.81e-17
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[2801] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[2802] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2803] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2804] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2805] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[2806] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2807] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.05e-10' 
		rate_values[2808] = 2.05e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*1.0' 
		rate_values[2809] = KDEC*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]' 
		rate_values[2810] = J[24]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.89e-11' 
		rate_values[2811] = 4.89e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.13e-13*1.0' 
		rate_values[2812] = 3.13e-13*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[2813] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-13' 
		rate_values[2814] = 2.86e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2815] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2816] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2817] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[2818] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2819] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.53e-12' 
		rate_values[2820] = 1.53e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.83e-12' 
		rate_values[2821] = 6.83e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.03e-12' 
		rate_values[2822] = 5.03e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2823] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2824] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2825] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2826] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[2827] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2828] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[2829] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2830] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2831] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2832] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2833] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[2834] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2835] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[2836] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.10e-14' 
		rate_values[2837] = 5.10e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.83e-15' 
		rate_values[2838] = 7.83e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2839] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2840] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2841] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2842] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[2843] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2844] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[2845] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2846] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2847] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2848] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2849] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2' 
		rate_values[2850] = 8.00e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2851] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[2852] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90e-12*numpy.exp(190/TEMP)' 
		rate_values[2853] = 1.90e-12*numpy.exp(190/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2854] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.50e-12' 
		rate_values[2855] = 1.50e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2856] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2857] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2858] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[2859] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.70' 
		rate_values[2860] = 9.20e-14*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.30' 
		rate_values[2861] = 9.20e-14*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2862] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.69e-11' 
		rate_values[2863] = 1.69e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.34e-11' 
		rate_values[2864] = 1.34e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2865] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.68' 
		rate_values[2866] = KDEC*0.68
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2867] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2868] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.70' 
		rate_values[2869] = 8.00e-13*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.30' 
		rate_values[2870] = 8.00e-13*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2871] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.15e-10' 
		rate_values[2872] = 1.15e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.11e-10' 
		rate_values[2873] = 1.11e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2874] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2875] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2876] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2877] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.70' 
		rate_values[2878] = 8.00e-13*RO2*0.70
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-13*RO2*0.30' 
		rate_values[2879] = 8.00e-13*RO2*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2880] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[2881] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.07e-10' 
		rate_values[2882] = 1.07e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[2883] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.04e-10' 
		rate_values[2884] = 1.04e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08e-12' 
		rate_values[2885] = 2.08e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.53e-12' 
		rate_values[2886] = 1.53e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00e-18' 
		rate_values[2887] = 5.00e-18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[2888] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.53e-12' 
		rate_values[2889] = 2.53e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[2890] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[2891] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00e-14' 
		rate_values[2892] = 7.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.20e-15' 
		rate_values[2893] = 1.20e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-14' 
		rate_values[2894] = 1.00e-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-15' 
		rate_values[2895] = 1.00e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.00e-18*H2O' 
		rate_values[2896] = 6.00e-18*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-17*H2O' 
		rate_values[2897] = 1.00e-17*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.255' 
		rate_values[2898] = KDEC*0.255
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.255' 
		rate_values[2899] = KDEC*0.255
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[2900] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.896' 
		rate_values[2901] = KRO2NO*0.896
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.104' 
		rate_values[2902] = KRO2NO*0.104
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.70e+14*numpy.exp(-6643/TEMP)' 
		rate_values[2903] = 2.70e+14*numpy.exp(-6643/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.896' 
		rate_values[2904] = KRO2NO*0.896
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.104' 
		rate_values[2905] = KRO2NO*0.104
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.70e+14*numpy.exp(-6643/TEMP)' 
		rate_values[2906] = 2.70e+14*numpy.exp(-6643/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.871' 
		rate_values[2907] = KRO2NO*0.871
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.129' 
		rate_values[2908] = KRO2NO*0.129
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.80e-14*numpy.exp(-260/TEMP)*O2' 
		rate_values[2909] = 1.80e-14*numpy.exp(-260/TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2910] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[2911] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[2912] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[2913] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[2914] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[2915] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[2916] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2917] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2918] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2*0.7' 
		rate_values[2919] = 6.70e-15*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2*0.3' 
		rate_values[2920] = 6.70e-15*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2921] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[2922] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2*0.7' 
		rate_values[2923] = 6.70e-15*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2*0.3' 
		rate_values[2924] = 6.70e-15*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2925] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[2926] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2*0.6' 
		rate_values[2927] = 2.50e-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2*0.2' 
		rate_values[2928] = 2.50e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2*0.2' 
		rate_values[2929] = 2.50e-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.52e-11' 
		rate_values[2930] = 2.52e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2931] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[2932] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.00e-11' 
		rate_values[2933] = 3.00e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[2934] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.69e-11' 
		rate_values[2935] = 2.69e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]' 
		rate_values[2936] = J[55]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[2937] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.56e-11' 
		rate_values[2938] = 2.56e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2939] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[2940] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.61e-11' 
		rate_values[2941] = 3.61e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[2942] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86e-11' 
		rate_values[2943] = 2.86e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]' 
		rate_values[2944] = J[55]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[2945] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.96e-11' 
		rate_values[2946] = 4.96e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2947] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[2948] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.35e-11' 
		rate_values[2949] = 8.35e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[2950] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00e-11' 
		rate_values[2951] = 8.00e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[2952] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[2953] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.33e-11' 
		rate_values[2954] = 4.33e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2955] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[2956] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.10e-10' 
		rate_values[2957] = 1.10e-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2' 
		rate_values[2958] = J[15]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.99e-11' 
		rate_values[2959] = 6.99e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.20' 
		rate_values[2960] = KDEC*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.80' 
		rate_values[2961] = KDEC*0.80
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[2962] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.862' 
		rate_values[2963] = KRO2NO*0.862
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2964] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[2965] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2966] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.896' 
		rate_values[2967] = KRO2NO*0.896
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.104' 
		rate_values[2968] = KRO2NO*0.104
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2969] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[2970] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[2971] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[2972] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[2973] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[2974] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[2975] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[2976] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2977] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2978] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*RO2*0.6' 
		rate_values[2979] = 1.30e-12*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*RO2*0.2' 
		rate_values[2980] = 1.30e-12*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*RO2*0.2' 
		rate_values[2981] = 1.30e-12*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2982] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2983] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.7' 
		rate_values[2984] = 9.20e-14*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.3' 
		rate_values[2985] = 9.20e-14*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[2986] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[2987] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2*0.7' 
		rate_values[2988] = 6.70e-15*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2*0.3' 
		rate_values[2989] = 6.70e-15*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[55]' 
		rate_values[2990] = J[55]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[2991] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.82e-12' 
		rate_values[2992] = 7.82e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2993] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.04e-11' 
		rate_values[2994] = 1.04e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2995] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.70e-11' 
		rate_values[2996] = 1.70e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2997] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.09e-11' 
		rate_values[2998] = 1.09e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[2999] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[3000] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.86e-11' 
		rate_values[3001] = 1.86e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.29e-12' 
		rate_values[3002] = 7.29e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.89e-12' 
		rate_values[3003] = 7.89e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.42e-12' 
		rate_values[3004] = 7.42e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.75e-11' 
		rate_values[3005] = 1.75e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[3006] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.77e-12' 
		rate_values[3007] = 6.77e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[3008] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[3009] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.69e-11' 
		rate_values[3010] = 6.69e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[3011] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.63e-11' 
		rate_values[3012] = 2.63e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*8.5' 
		rate_values[3013] = KNO3AL*8.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[3014] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[3015] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[3016] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[3017] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[3018] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[3019] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[3020] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[3021] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[3022] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[3023] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[3024] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[3025] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[3026] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[3027] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[3028] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[3029] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30e-12*RO2' 
		rate_values[3030] = 1.30e-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[3031] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[3032] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70e-15*RO2' 
		rate_values[3033] = 6.70e-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[3034] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[3035] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80e-13*RO2' 
		rate_values[3036] = 8.80e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.96e-12' 
		rate_values[3037] = 2.96e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[3038] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.65e-12' 
		rate_values[3039] = 9.65e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.57e-12' 
		rate_values[3040] = 6.57e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[3041] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.27e-11' 
		rate_values[3042] = 1.27e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[3043] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.31e-11' 
		rate_values[3044] = 3.31e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[3045] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.46e-11' 
		rate_values[3046] = 7.46e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.138' 
		rate_values[3047] = KRO2NO*0.138
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[53]' 
		rate_values[3048] = J[53]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.29e-12' 
		rate_values[3049] = 3.29e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.9e-12*numpy.exp(225/TEMP)' 
		rate_values[3050] = 5.9e-12*numpy.exp(225/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.40e-15' 
		rate_values[3051] = 2.40e-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[3052] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[3053] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[3054] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[3055] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[3056] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[3057] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.00e-13' 
		rate_values[3058] = 9.00e-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[3059] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN*0.67' 
		rate_values[3060] = KBPAN*0.67
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[3061] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[3062] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.7' 
		rate_values[3063] = 1.00e-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00e-11*RO2*0.3' 
		rate_values[3064] = 1.00e-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[3065] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50e-13*RO2' 
		rate_values[3066] = 2.50e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[3067] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.66e-12' 
		rate_values[3068] = 4.66e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.10e-12' 
		rate_values[3069] = 1.10e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.60e-12' 
		rate_values[3070] = 3.60e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.06e-12' 
		rate_values[3071] = 1.06e-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3072] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3073] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3074] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3075] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3076] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3077] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3078] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3079] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3080] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3081] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3082] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3083] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3084] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3085] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3086] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3087] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3088] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3089] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3090] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3091] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3092] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3093] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3094] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3095] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3096] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3097] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3098] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3099] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3100] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3101] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3102] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3103] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3104] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3105] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3106] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3107] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3108] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3109] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3110] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3111] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3112] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3113] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3114] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3115] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3116] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3117] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3118] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3119] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3120] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3121] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3122] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3123] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3124] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3125] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3126] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3127] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3128] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3129] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3130] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[3131] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[3132] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[3133] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[3134] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.8e-13*RO2' 
		rate_values[3135] = 8.8e-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[3136] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[3137] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[3138] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.38e-11' 
		rate_values[3139] = 3.38e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[20]*2' 
		rate_values[3140] = J[20]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.72e-11' 
		rate_values[3141] = 4.72e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[3142] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[3143] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[3144] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2' 
		rate_values[3145] = 9.20e-14*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[3146] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[3147] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[3148] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[3149] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.42e-11' 
		rate_values[3150] = 3.42e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.40e-17*0.27' 
		rate_values[3151] = 2.40e-17*0.27
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.40e-17*0.73' 
		rate_values[3152] = 2.40e-17*0.73
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[3153] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[20]*2' 
		rate_values[3154] = J[20]*2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.72e-11' 
		rate_values[3155] = 4.72e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[3156] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[3157] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[3158] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2' 
		rate_values[3159] = 9.20e-14*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[3160] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[3161] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[3162] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[3163] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.42e-11' 
		rate_values[3164] = 3.42e-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.40e-17*0.73' 
		rate_values[3165] = 2.40e-17*0.73
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.40e-17*0.27' 
		rate_values[3166] = 2.40e-17*0.27
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[3167] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00e11*numpy.exp(-3160/TEMP)' 
		rate_values[3168] = 7.00e11*numpy.exp(-3160/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00e-12*O2' 
		rate_values[3169] = 5.00e-12*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00e-12*O2*3.2*numpy.exp(-550/TEMP)' 
		rate_values[3170] = 5.00e-12*O2*3.2*numpy.exp(-550/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00e-12*O2*3.2*(1-numpy.exp(-550/TEMP))' 
		rate_values[3171] = 5.00e-12*O2*3.2*(1-numpy.exp(-550/TEMP))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[56]*0.91' 
		rate_values[3172] = J[56]*0.91
	except:
		erf = 1 # flag error
		err_mess = (str('Error: Could not calculate '+ 
		'rate coefficient for equation number ' 
		+ str(gprn) + ' ' + rc_eq_now + 
		' (message from rate coeffs.py)'))
	
	# aqueous-phase reactions
	
	# surface (e.g. wall) reactions
	
	return(rate_values, erf, err_mess)
