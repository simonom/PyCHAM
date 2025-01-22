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
# created at 2025-01-22 12:02:23.647768

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
		K298CH3O2=3.5E-13 
		gprn += 1 # keep count on reaction number
		KAPHO2=5.2E-13*numpy.exp(980./TEMP) 
		gprn += 1 # keep count on reaction number
		KAPNO=7.5E-12*numpy.exp(290./TEMP) 
		gprn += 1 # keep count on reaction number
		KCH3O2=1.03E-13*numpy.exp(365./TEMP) 
		gprn += 1 # keep count on reaction number
		KDEC=1.00E+06 
		gprn += 1 # keep count on reaction number
		KMT05=1.44E-13*(1.+(M/4.2E+19)) 
		gprn += 1 # keep count on reaction number
		KMT06=1.+(1.40E-21*numpy.exp(2200./TEMP)*H2O) 
		gprn += 1 # keep count on reaction number
		KNO3AL=1.44E-12*numpy.exp(-1862./TEMP) 
		gprn += 1 # keep count on reaction number
		KRO2HO2=2.91E-13*numpy.exp(1300./TEMP) 
		gprn += 1 # keep count on reaction number
		KRO2NO=2.7E-12*numpy.exp(360./TEMP) 
		gprn += 1 # keep count on reaction number
		KRO2NO3=2.3E-12 
		gprn += 1 # keep count on reaction number
		KROSEC=2.50E-14*numpy.exp(-300./TEMP) 
		gprn += 1 # keep count on reaction number
		FCD=0.30 
		gprn += 1 # keep count on reaction number
		KD0=1.10E-05*M*numpy.exp(-10100./TEMP) 
		gprn += 1 # keep count on reaction number
		KDI=1.90E17*numpy.exp(-14100./TEMP) 
		gprn += 1 # keep count on reaction number
		KRD=KD0/KDI 
		gprn += 1 # keep count on reaction number
		NCD=0.75-1.27*(numpy.log10(FCD)) 
		gprn += 1 # keep count on reaction number
		FD=10.**(numpy.log10(FCD)/(1.+(numpy.log10(KRD)/NCD)**(2.))) 
		gprn += 1 # keep count on reaction number
		KBPAN=(KD0*KDI)*FD/(KD0+KDI) 
		gprn += 1 # keep count on reaction number
		FCC=0.30 
		gprn += 1 # keep count on reaction number
		KC0=3.28E-28*M*(TEMP/300.)**(-6.87) 
		gprn += 1 # keep count on reaction number
		KCI=1.125E-11*(TEMP/300.)**(-1.105) 
		gprn += 1 # keep count on reaction number
		KRC=KC0/KCI 
		gprn += 1 # keep count on reaction number
		NC=0.75-1.27*(numpy.log10(FCC)) 
		gprn += 1 # keep count on reaction number
		FC=10.**(numpy.log10(FCC)/(1.+(numpy.log10(KRC)/NC)**(2.))) 
		gprn += 1 # keep count on reaction number
		KFPAN=(KC0*KCI)*FC/(KC0+KCI) 
		gprn += 1 # keep count on reaction number
		FC1=0.85 
		gprn += 1 # keep count on reaction number
		K10=1.0E-31*M*(TEMP/300.)**(-1.6) 
		gprn += 1 # keep count on reaction number
		K1I=5.0E-11*(TEMP/300.)**(-0.3) 
		gprn += 1 # keep count on reaction number
		KR1=K10/K1I 
		gprn += 1 # keep count on reaction number
		NC1=0.75-1.27*(numpy.log10(FC1)) 
		gprn += 1 # keep count on reaction number
		F1=10.**(numpy.log10(FC1)/(1.+(numpy.log10(KR1)/NC1)**(2.))) 
		gprn += 1 # keep count on reaction number
		KMT01=(K10*K1I)*F1/(K10+K1I) 
		gprn += 1 # keep count on reaction number
		FC2=0.6 
		gprn += 1 # keep count on reaction number
		K20=1.3E-31*M*(TEMP/300.)**(-1.5) 
		gprn += 1 # keep count on reaction number
		K2I=2.3E-11*(TEMP/300.)**(0.24) 
		gprn += 1 # keep count on reaction number
		KR2=K20/K2I 
		gprn += 1 # keep count on reaction number
		NC2=0.75-1.27*(numpy.log10(FC2)) 
		gprn += 1 # keep count on reaction number
		F2=10.**(numpy.log10(FC2)/(1.+(numpy.log10(KR2)/NC2)**(2.))) 
		gprn += 1 # keep count on reaction number
		KMT02=(K20*K2I)*F2/(K20+K2I) 
		gprn += 1 # keep count on reaction number
		FC3=0.35 
		gprn += 1 # keep count on reaction number
		K30=3.6E-30*M*(TEMP/300.)**(-4.1) 
		gprn += 1 # keep count on reaction number
		K3I=1.9E-12*(TEMP/300.)**(0.2) 
		gprn += 1 # keep count on reaction number
		KR3=K30/K3I 
		gprn += 1 # keep count on reaction number
		NC3=0.75-1.27*(numpy.log10(FC3)) 
		gprn += 1 # keep count on reaction number
		F3=10.**(numpy.log10(FC3)/(1.+(numpy.log10(KR3)/NC3)**(2.))) 
		gprn += 1 # keep count on reaction number
		KMT03=(K30*K3I)*F3/(K30+K3I) 
		gprn += 1 # keep count on reaction number
		FC4=0.35 
		gprn += 1 # keep count on reaction number
		K40=1.3E-3*M*(TEMP/300.)**(-3.5)*numpy.exp(-11000./TEMP) 
		gprn += 1 # keep count on reaction number
		K4I=9.7E+14*(TEMP/300.)**(0.1)*numpy.exp(-11080./TEMP) 
		gprn += 1 # keep count on reaction number
		KR4=K40/K4I 
		gprn += 1 # keep count on reaction number
		NC4=0.75-1.27*(numpy.log10(FC4)) 
		gprn += 1 # keep count on reaction number
		F4=10.**(numpy.log10(FC4)/(1.+(numpy.log10(KR4)/NC4)**(2.))) 
		gprn += 1 # keep count on reaction number
		KMT04=(K40*K4I)*F4/(K40+K4I) 
		gprn += 1 # keep count on reaction number
		FC7=0.81 
		gprn += 1 # keep count on reaction number
		K70=7.4E-31*M*(TEMP/300.)**(-2.4) 
		gprn += 1 # keep count on reaction number
		K7I=3.3E-11*(TEMP/300.)**(-0.3) 
		gprn += 1 # keep count on reaction number
		KR7=K70/K7I 
		gprn += 1 # keep count on reaction number
		NC7=0.75-1.27*(numpy.log10(FC7)) 
		gprn += 1 # keep count on reaction number
		F7=10.**(numpy.log10(FC7)/(1.+(numpy.log10(KR7)/NC7)**(2.))) 
		gprn += 1 # keep count on reaction number
		KMT07=(K70*K7I)*F7/(K70+K7I) 
		gprn += 1 # keep count on reaction number
		FC8=0.41 
		gprn += 1 # keep count on reaction number
		K80=3.2E-30*M*(TEMP/300.)**(-4.5) 
		gprn += 1 # keep count on reaction number
		K8I=3.0E-11 
		gprn += 1 # keep count on reaction number
		KR8=K80/K8I 
		gprn += 1 # keep count on reaction number
		NC8=0.75-1.27*(numpy.log10(FC8)) 
		gprn += 1 # keep count on reaction number
		F8=10.**(numpy.log10(FC8)/(1.+(numpy.log10(KR8)/NC8)**(2.))) 
		gprn += 1 # keep count on reaction number
		KMT08=(K80*K8I)*F8/(K80+K8I) 
		gprn += 1 # keep count on reaction number
		FC9=0.4 
		gprn += 1 # keep count on reaction number
		K90=1.4E-31*M*(TEMP/300.)**(-3.1) 
		gprn += 1 # keep count on reaction number
		K9I=4.0E-12 
		gprn += 1 # keep count on reaction number
		KR9=K90/K9I 
		gprn += 1 # keep count on reaction number
		NC9=0.75-1.27*(numpy.log10(FC9)) 
		gprn += 1 # keep count on reaction number
		F9=10.**(numpy.log10(FC9)/(1.+(numpy.log10(KR9)/NC9)**(2.))) 
		gprn += 1 # keep count on reaction number
		KMT09=(K90*K9I)*F9/(K90+K9I) 
		gprn += 1 # keep count on reaction number
		FC10=0.4 
		gprn += 1 # keep count on reaction number
		K100=4.10E-05*M*numpy.exp(-10650./TEMP) 
		gprn += 1 # keep count on reaction number
		K10I=6.0E+15*numpy.exp(-11170./TEMP) 
		gprn += 1 # keep count on reaction number
		KR10=K100/K10I 
		gprn += 1 # keep count on reaction number
		NC10=0.75-1.27*(numpy.log10(FC10)) 
		gprn += 1 # keep count on reaction number
		F10=10.**(numpy.log10(FC10)/(1.+(numpy.log10(KR10)/NC10)**(2.))) 
		gprn += 1 # keep count on reaction number
		KMT10=(K100*K10I)*F10/(K100+K10I) 
		gprn += 1 # keep count on reaction number
		K3=6.50E-34*numpy.exp(1335./TEMP) 
		gprn += 1 # keep count on reaction number
		K4=2.70E-17*numpy.exp(2199./TEMP) 
		gprn += 1 # keep count on reaction number
		K1=2.40E-14*numpy.exp(460./TEMP) 
		gprn += 1 # keep count on reaction number
		K2=(K3*M)/(1.+(K3*M/K4)) 
		gprn += 1 # keep count on reaction number
		KMT11=K1+K2 
		gprn += 1 # keep count on reaction number
		FC12=0.53 
		gprn += 1 # keep count on reaction number
		K120=2.5E-31*M*(TEMP/300.)**(-2.6) 
		gprn += 1 # keep count on reaction number
		K12I=2.0E-12 
		gprn += 1 # keep count on reaction number
		KR12=K120/K12I 
		gprn += 1 # keep count on reaction number
		NC12=0.75-1.27*(numpy.log10(FC12)) 
		gprn += 1 # keep count on reaction number
		F12=10.**(numpy.log10(FC12)/(1.0+(numpy.log10(KR12)/NC12)**(2.))) 
		gprn += 1 # keep count on reaction number
		KMT12=(K120*K12I*F12)/(K120+K12I) 
		gprn += 1 # keep count on reaction number
		FC13=0.36 
		gprn += 1 # keep count on reaction number
		K130=2.5E-30*M*(TEMP/300.)**(-5.5) 
		gprn += 1 # keep count on reaction number
		K13I=1.8E-11 
		gprn += 1 # keep count on reaction number
		KR13=K130/K13I 
		gprn += 1 # keep count on reaction number
		NC13=0.75-1.27*(numpy.log10(FC13)) 
		gprn += 1 # keep count on reaction number
		F13=10.**(numpy.log10(FC13)/(1.+(numpy.log10(KR13)/NC13)**(2.))) 
		gprn += 1 # keep count on reaction number
		KMT13=(K130*K13I)*F13/(K130+K13I) 
		gprn += 1 # keep count on reaction number
		FC14=0.36 
		gprn += 1 # keep count on reaction number
		K140=9.0E-5*numpy.exp(-9690./TEMP)*M 
		gprn += 1 # keep count on reaction number
		K14I=1.1E+16*numpy.exp(-10560./TEMP) 
		gprn += 1 # keep count on reaction number
		KR14=K140/K14I 
		gprn += 1 # keep count on reaction number
		NC14=0.75-1.27*(numpy.log10(FC14)) 
		gprn += 1 # keep count on reaction number
		F14=10.**(numpy.log10(FC14)/(1.+(numpy.log10(KR14)/NC14)**(2.))) 
		gprn += 1 # keep count on reaction number
		KMT14=(K140*K14I)*F14/(K140+K14I) 
		gprn += 1 # keep count on reaction number
		#{1.}APINOOA=C107O2+OH:KDEC*0.55*0.90; 
		gprn += 1 # keep count on reaction number
		#{2.}APINOOA=C109O2+OH:KDEC*0.45*0.90; 
		gprn += 1 # keep count on reaction number
		#RO2+NO-}R=O+NO2+HO2 
		gprn += 1 # keep count on reaction number
		#RO2(MCMpool)+RO2(PRAM)-}RC=O+ROH,orRO2(MCMpool)+RO2(PRAM)-}RO+RO+O2,or 
		gprn += 1 # keep count on reaction number
		#RO2(MCMpool)+RO2(PRAM)-}ROH+RC=O 
		gprn += 1 # keep count on reaction number
		#{1153.}C107O2=C107O:	9.20e-14*0.7*RO2*0.9; 
		gprn += 1 # keep count on reaction number
		#{1154.}C107O2=C107OH:	9.20e-14*0.3*RO2*0.9; 
		gprn += 1 # keep count on reaction number
		#RO2+NO-}R=O+NO2+HO2 
		gprn += 1 # keep count on reaction number
		#RO2(MCMpool)+RO2(PRAM)-}RC=O+ROH,orRO2(MCMpool)+RO2(PRAM)-}RO+RO+O2,or 
		gprn += 1 # keep count on reaction number
		#RO2(MCMpool)+RO2(PRAM)-}ROH+RC=O 

	except:
		erf = 1 # flag error
		err_mess = str('Error: generic reaction rates failed to be calculated inside rate_coeffs.py at number ' + str(gprn) + ', please check chemical scheme and associated chemical scheme markers, which are stated in the model variables input file') # error message
		return([], erf, err_mess)
	# estimate and append photolysis rates
	J = photolysisRates.PhotolysisCalculation(TEMP, Jlen, sumt, self)

	if (self.light_stat_now == 0):
		J = [0]*len(J)
	rate_values = numpy.zeros((2667))
	
	# if reactions have been found in the chemical scheme
	# gas-phase reactions
	gprn = 0 # keep count on reaction number
	try:
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.6E-34*N2*(TEMP/300.)**(-2.6)*O2+6.0E-34*O2*(TEMP/300.)**(-2.6)*O2' 
		rate_values[0] = 5.6E-34*N2*(TEMP/300.)**(-2.6)*O2+6.0E-34*O2*(TEMP/300.)**(-2.6)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.0E-12*numpy.exp(-2060./TEMP)' 
		rate_values[1] = 8.0E-12*numpy.exp(-2060./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT01' 
		rate_values[2] = KMT01
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5E-12*numpy.exp(188./TEMP)' 
		rate_values[3] = 5.5E-12*numpy.exp(188./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT02' 
		rate_values[4] = KMT02
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.2E-11*numpy.exp(67./TEMP)*O2+2.0E-11*numpy.exp(130./TEMP)*N2' 
		rate_values[5] = 3.2E-11*numpy.exp(67./TEMP)*O2+2.0E-11*numpy.exp(130./TEMP)*N2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4E-12*numpy.exp(-1310./TEMP)' 
		rate_values[6] = 1.4E-12*numpy.exp(-1310./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4E-13*numpy.exp(-2470./TEMP)' 
		rate_values[7] = 1.4E-13*numpy.exp(-2470./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.3E-39*numpy.exp(530./TEMP)*O2' 
		rate_values[8] = 3.3E-39*numpy.exp(530./TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.8E-11*numpy.exp(110./TEMP)' 
		rate_values[9] = 1.8E-11*numpy.exp(110./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.50E-14*numpy.exp(-1260./TEMP)' 
		rate_values[10] = 4.50E-14*numpy.exp(-1260./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT03' 
		rate_values[11] = KMT03
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.14E-10*H2O' 
		rate_values[12] = 2.14E-10*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.70E-12*numpy.exp(-940./TEMP)' 
		rate_values[13] = 1.70E-12*numpy.exp(-940./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.7E-12*numpy.exp(-2100./TEMP)' 
		rate_values[14] = 7.7E-12*numpy.exp(-2100./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT05' 
		rate_values[15] = KMT05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.9E-12*numpy.exp(-160./TEMP)' 
		rate_values[16] = 2.9E-12*numpy.exp(-160./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.03E-16*(TEMP/300.)**(4.57)*numpy.exp(693./TEMP)' 
		rate_values[17] = 2.03E-16*(TEMP/300.)**(4.57)*numpy.exp(693./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.8E-11*numpy.exp(250./TEMP)' 
		rate_values[18] = 4.8E-11*numpy.exp(250./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.20E-13*KMT06*numpy.exp(600./TEMP)+1.90E-33*M*KMT06*numpy.exp(980./TEMP)' 
		rate_values[19] = 2.20E-13*KMT06*numpy.exp(600./TEMP)+1.90E-33*M*KMT06*numpy.exp(980./TEMP)
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
		rc_eq_now = '2.0E-11' 
		rate_values[22] = 2.0E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.45E-12*numpy.exp(270./TEMP)' 
		rate_values[23] = 3.45E-12*numpy.exp(270./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT09' 
		rate_values[24] = KMT09
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.2E-13*numpy.exp(690./TEMP)*1.0' 
		rate_values[25] = 3.2E-13*numpy.exp(690./TEMP)*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0E-12' 
		rate_values[26] = 4.0E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.5E-12*numpy.exp(260./TEMP)' 
		rate_values[27] = 2.5E-12*numpy.exp(260./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT11' 
		rate_values[28] = KMT11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0E-32*numpy.exp(-1000./TEMP)*M' 
		rate_values[29] = 4.0E-32*numpy.exp(-1000./TEMP)*M
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT12' 
		rate_values[30] = KMT12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.3E-12*numpy.exp(-330./TEMP)*O2' 
		rate_values[31] = 1.3E-12*numpy.exp(-330./TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.00E-06' 
		rate_values[32] = 6.00E-06
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.00E-04' 
		rate_values[33] = 4.00E-04
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.20E-15*H2O' 
		rate_values[34] = 1.20E-15*H2O
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
		rc_eq_now = '1.85E-12*numpy.exp(-1690./TEMP)' 
		rate_values[45] = 1.85E-12*numpy.exp(-1690./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3E-12*numpy.exp(360./TEMP)*0.999' 
		rate_values[46] = 2.3E-12*numpy.exp(360./TEMP)*0.999
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3E-12*numpy.exp(360./TEMP)*0.001' 
		rate_values[47] = 2.3E-12*numpy.exp(360./TEMP)*0.001
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.2E-14*numpy.exp(-1080./TEMP)*O2' 
		rate_values[48] = 7.2E-14*numpy.exp(-1080./TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT13' 
		rate_values[49] = KMT13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT14' 
		rate_values[50] = KMT14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2E-12' 
		rate_values[51] = 1.2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8E-13*numpy.exp(780./TEMP)*(1.-1./(1.+498.*numpy.exp(-1160./TEMP)))' 
		rate_values[52] = 3.8E-13*numpy.exp(780./TEMP)*(1.-1./(1.+498.*numpy.exp(-1160./TEMP)))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*KCH3O2*RO2*7.18*numpy.exp(-885./TEMP)' 
		rate_values[53] = 2.*KCH3O2*RO2*7.18*numpy.exp(-885./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*KCH3O2*RO2*0.5*(1.-7.18*numpy.exp(-885./TEMP))' 
		rate_values[54] = 2.*KCH3O2*RO2*0.5*(1.-7.18*numpy.exp(-885./TEMP))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*KCH3O2*RO2*0.5*(1.-7.18*numpy.exp(-885./TEMP))' 
		rate_values[55] = 2.*KCH3O2*RO2*0.5*(1.-7.18*numpy.exp(-885./TEMP))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0E-13*numpy.exp(-845./TEMP)' 
		rate_values[56] = 4.0E-13*numpy.exp(-845./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[20]' 
		rate_values[57] = J[20]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.3E-12*numpy.exp(190./TEMP)*0.6' 
		rate_values[58] = 5.3E-12*numpy.exp(190./TEMP)*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.3E-12*numpy.exp(190./TEMP)*0.4' 
		rate_values[59] = 5.3E-12*numpy.exp(190./TEMP)*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[60] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[61] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[62] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.65E-11' 
		rate_values[63] = 6.65E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[64] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[65] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[66] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[67] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[68] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[69] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2' 
		rate_values[70] = 1.00E-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.12E-13' 
		rate_values[71] = 3.12E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[11]' 
		rate_values[72] = J[19]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.23E-12' 
		rate_values[73] = 4.23E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]*2.' 
		rate_values[74] = J[11]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.29E-11' 
		rate_values[75] = 4.29E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*KNO3AL*2.4' 
		rate_values[76] = 2.*KNO3AL*2.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[77] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[78] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[79] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[80] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[81] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[82] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.7*RO2' 
		rate_values[83] = 1.00E-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.3*RO2' 
		rate_values[84] = 1.00E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.10E-11' 
		rate_values[85] = 2.10E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.49E-11' 
		rate_values[86] = 2.49E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[11]' 
		rate_values[87] = J[19]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.14E-11' 
		rate_values[88] = 2.14E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[89] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[90] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[91] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[92] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[93] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*0.6*RO2' 
		rate_values[94] = 2.00E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*0.2*RO2' 
		rate_values[95] = 2.00E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*0.2*RO2' 
		rate_values[96] = 2.00E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.99E-12' 
		rate_values[97] = 5.99E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[98] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[18]' 
		rate_values[99] = J[19]+J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.69E-12' 
		rate_values[100] = 2.69E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[101] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.23E-11' 
		rate_values[102] = 1.23E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*4.0' 
		rate_values[103] = KNO3AL*4.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.6E-12*numpy.exp(-1240./TEMP)' 
		rate_values[104] = 6.6E-12*numpy.exp(-1240./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.800' 
		rate_values[105] = 1.00E-11*0.800
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.200' 
		rate_values[106] = 1.00E-11*0.200
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[107] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL' 
		rate_values[108] = KNO3AL
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[109] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[110] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[111] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[112] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[113] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[114] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.7*RO2' 
		rate_values[115] = 1.00E-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.3*RO2' 
		rate_values[116] = 1.00E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.73E-12' 
		rate_values[117] = 2.73E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[118] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.19E-12' 
		rate_values[119] = 6.19E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6E-12*numpy.exp(305./TEMP)' 
		rate_values[120] = 1.6E-12*numpy.exp(305./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]' 
		rate_values[121] = J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8E-13*numpy.exp(780./TEMP)*(1./(1.+498.*numpy.exp(-1160./TEMP)))' 
		rate_values[122] = 3.8E-13*numpy.exp(780./TEMP)*(1./(1.+498.*numpy.exp(-1160./TEMP)))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.36E-13*numpy.exp(1250./TEMP)*0.15' 
		rate_values[123] = 1.36E-13*numpy.exp(1250./TEMP)*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[14]' 
		rate_values[124] = J[14]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[125] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[16]' 
		rate_values[126] = J[16]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.1E-12*numpy.exp(340./TEMP)' 
		rate_values[127] = 3.1E-12*numpy.exp(340./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL' 
		rate_values[128] = KNO3AL
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[129] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[130] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.9E-12*numpy.exp(575./TEMP)' 
		rate_values[131] = 1.9E-12*numpy.exp(575./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*2.4' 
		rate_values[132] = KNO3AL*2.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[133] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[134] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[135] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[136] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.7*RO2' 
		rate_values[137] = 1.00E-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.3*RO2' 
		rate_values[138] = 1.00E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.58E-11' 
		rate_values[139] = 1.58E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[11]' 
		rate_values[140] = J[19]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.23E-11' 
		rate_values[141] = 1.23E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.983' 
		rate_values[142] = KRO2NO*0.983
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.017' 
		rate_values[143] = KRO2NO*0.017
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[144] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[145] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[146] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*0.6*RO2' 
		rate_values[147] = 2.00E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*0.2*RO2' 
		rate_values[148] = 2.00E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*0.2*RO2' 
		rate_values[149] = 2.00E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.23E-12' 
		rate_values[150] = 2.23E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.77E-11' 
		rate_values[151] = 5.77E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[152] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.88E-11' 
		rate_values[153] = 1.88E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]' 
		rate_values[154] = J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[155] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.45E-11' 
		rate_values[156] = 2.45E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*4.0' 
		rate_values[157] = KNO3AL*4.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[158] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[159] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[160] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[161] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[162] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2' 
		rate_values[163] = 1.00E-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.74E-12' 
		rate_values[164] = 3.74E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.34E-12' 
		rate_values[165] = 7.34E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[166] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]' 
		rate_values[167] = J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.4E-12*numpy.exp(135./TEMP)' 
		rate_values[168] = 5.4E-12*numpy.exp(135./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[9]' 
		rate_values[169] = J[9]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[10]' 
		rate_values[170] = J[10]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5E-16' 
		rate_values[171] = 5.5E-16
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.5E-12*numpy.exp(290./TEMP)' 
		rate_values[172] = 7.5E-12*numpy.exp(290./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[173] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[174] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.0E-12' 
		rate_values[175] = 4.0E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[176] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[177] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.7*RO2' 
		rate_values[178] = 1.00E-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.3*RO2' 
		rate_values[179] = 1.00E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.E-14' 
		rate_values[180] = 3.E-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[181] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.70E-12' 
		rate_values[182] = 3.70E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[183] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[184] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[185] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.387' 
		rate_values[186] = KRO2HO2*0.387
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*0.6*RO2' 
		rate_values[187] = 2.00E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*0.2*RO2' 
		rate_values[188] = 2.00E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*0.2*RO2' 
		rate_values[189] = 2.00E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[190] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[191] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[192] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.91E-11' 
		rate_values[193] = 2.91E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[12]' 
		rate_values[194] = J[12]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.8E-12*numpy.exp(-1320./TEMP)+1.7E-14*numpy.exp(423./TEMP)' 
		rate_values[195] = 8.8E-12*numpy.exp(-1320./TEMP)+1.7E-14*numpy.exp(423./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[196] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[197] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[198] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.36E-13*numpy.exp(1250./TEMP)*0.85' 
		rate_values[199] = 1.36E-13*numpy.exp(1250./TEMP)*0.85
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(K298CH3O2*8.0E-12)**(0.5)*RO2*0.6' 
		rate_values[200] = 2.*(K298CH3O2*8.0E-12)**(0.5)*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(K298CH3O2*8.0E-12)**(0.5)*RO2*0.2' 
		rate_values[201] = 2.*(K298CH3O2*8.0E-12)**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(K298CH3O2*8.0E-12)**(0.5)*RO2*0.2' 
		rate_values[202] = 2.*(K298CH3O2*8.0E-12)**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[203] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.39E-12' 
		rate_values[204] = 8.39E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[205] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]' 
		rate_values[206] = J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[207] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[208] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[209] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[210] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[211] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[212] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[213] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2' 
		rate_values[214] = 1.00E-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[215] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[216] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*RO2' 
		rate_values[217] = 2.00E-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.19E-11' 
		rate_values[218] = 7.19E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.55E-11' 
		rate_values[219] = 7.55E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[13]+J[11]' 
		rate_values[220] = J[19]+J[13]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.33E-11' 
		rate_values[221] = 8.33E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[13]+J[11]' 
		rate_values[222] = J[19]+J[13]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[223] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[224] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[225] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[226] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[227] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[228] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[229] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2' 
		rate_values[230] = 1.00E-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[231] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.520' 
		rate_values[232] = KRO2HO2*0.520
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*RO2' 
		rate_values[233] = 2.00E-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.27E-11' 
		rate_values[234] = 1.27E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.63E-11' 
		rate_values[235] = 1.63E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[17]' 
		rate_values[236] = J[19]+J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.41E-11' 
		rate_values[237] = 2.41E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[17]' 
		rate_values[238] = J[19]+J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.85E-12*numpy.exp(-345./TEMP)' 
		rate_values[239] = 2.85E-12*numpy.exp(-345./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[240] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[241] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.39E-11' 
		rate_values[242] = 3.39E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*KNO3AL*4.0' 
		rate_values[243] = 2.*KNO3AL*4.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00E-13' 
		rate_values[244] = 8.00E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.44E-11' 
		rate_values[245] = 1.44E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]+J[18]' 
		rate_values[246] = J[17]+J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2E-11*numpy.exp(440./TEMP)*0.5577' 
		rate_values[247] = 1.2E-11*numpy.exp(440./TEMP)*0.5577
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2E-11*numpy.exp(440./TEMP)*0.344' 
		rate_values[248] = 1.2E-11*numpy.exp(440./TEMP)*0.344
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2E-11*numpy.exp(440./TEMP)*0.073125' 
		rate_values[249] = 1.2E-11*numpy.exp(440./TEMP)*0.073125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.770' 
		rate_values[250] = KRO2NO*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.230' 
		rate_values[251] = KRO2NO*0.230
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[252] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[253] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[254] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20E-14*RO2*0.7' 
		rate_values[255] = 9.20E-14*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20E-14*RO2*0.3' 
		rate_values[256] = 9.20E-14*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.50E-12' 
		rate_values[257] = 5.50E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[258] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.83E-11' 
		rate_values[259] = 1.83E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[260] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2E-12*numpy.exp(600./TEMP)*0.772' 
		rate_values[261] = 5.2E-12*numpy.exp(600./TEMP)*0.772
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.2E-12*numpy.exp(600./TEMP)*0.228' 
		rate_values[262] = 5.2E-12*numpy.exp(600./TEMP)*0.228
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.0E-14' 
		rate_values[263] = 2.0E-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[264] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.843' 
		rate_values[265] = KRO2NO*0.843
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.157' 
		rate_values[266] = KRO2NO*0.157
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.20E+10*numpy.exp(-3523./TEMP)' 
		rate_values[267] = 4.20E+10*numpy.exp(-3523./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[268] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[269] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.882' 
		rate_values[270] = KRO2NO*0.882
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.118' 
		rate_values[271] = KRO2NO*0.118
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[272] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.902' 
		rate_values[273] = KRO2NO*0.902
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.098' 
		rate_values[274] = KRO2NO*0.098
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[275] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.950' 
		rate_values[276] = KRO2NO*0.950
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.050' 
		rate_values[277] = KRO2NO*0.050
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[278] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.875' 
		rate_values[279] = KRO2NO*0.875
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.125' 
		rate_values[280] = KRO2NO*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[281] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[282] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[283] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[284] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[285] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[286] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[287] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[288] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.7*RO2' 
		rate_values[289] = 1.00E-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.3*RO2' 
		rate_values[290] = 1.00E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[291] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[292] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30E-12*0.6*RO2' 
		rate_values[293] = 1.30E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30E-12*0.2*RO2' 
		rate_values[294] = 1.30E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30E-12*0.2*RO2' 
		rate_values[295] = 1.30E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[296] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[297] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*0.7*RO2' 
		rate_values[298] = 6.70E-15*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*0.3*RO2' 
		rate_values[299] = 6.70E-15*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[300] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[301] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*0.7*RO2' 
		rate_values[302] = 6.70E-15*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*0.3*RO2' 
		rate_values[303] = 6.70E-15*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[304] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[305] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.6*RO2' 
		rate_values[306] = 8.80E-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.2*RO2' 
		rate_values[307] = 8.80E-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.2*RO2' 
		rate_values[308] = 8.80E-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[309] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[310] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*0.7*RO2' 
		rate_values[311] = 6.70E-15*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*0.3*RO2' 
		rate_values[312] = 6.70E-15*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[313] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[314] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*0.7*RO2' 
		rate_values[315] = 6.70E-15*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*0.3*RO2' 
		rate_values[316] = 6.70E-15*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[317] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[318] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.6*RO2' 
		rate_values[319] = 8.80E-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.2*RO2' 
		rate_values[320] = 8.80E-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.2*RO2' 
		rate_values[321] = 8.80E-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[21]+J[13]' 
		rate_values[322] = J[21]+J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.88E-12' 
		rate_values[323] = 2.88E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]+J[18]' 
		rate_values[324] = J[23]+J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.37E-12' 
		rate_values[325] = 5.37E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[326] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '71.11E-12' 
		rate_values[327] = 71.11E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]+J[11]' 
		rate_values[328] = J[23]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.25E-11' 
		rate_values[329] = 2.25E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]+J[11]' 
		rate_values[330] = J[23]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.03E-11' 
		rate_values[331] = 7.03E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.66E-12' 
		rate_values[332] = 3.66E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[13]' 
		rate_values[333] = J[19]+J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.73E-12' 
		rate_values[334] = 9.73E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[13]' 
		rate_values[335] = J[19]+J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[336] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30E-11' 
		rate_values[337] = 1.30E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[13]' 
		rate_values[338] = J[19]+J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.05E-11' 
		rate_values[339] = 1.05E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[18]' 
		rate_values[340] = J[19]+J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.05E-11' 
		rate_values[341] = 2.05E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[18]' 
		rate_values[342] = J[19]+J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.69E-11' 
		rate_values[343] = 8.69E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[11]' 
		rate_values[344] = J[19]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.75E-11' 
		rate_values[345] = 2.75E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[11]' 
		rate_values[346] = J[19]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.01E-11' 
		rate_values[347] = 8.01E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[11]' 
		rate_values[348] = J[19]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.20E-10' 
		rate_values[349] = 1.20E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]' 
		rate_values[350] = J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.65E-12' 
		rate_values[351] = 6.65E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]' 
		rate_values[352] = J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.67E-12' 
		rate_values[353] = 7.67E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]' 
		rate_values[354] = J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.20E-12' 
		rate_values[355] = 7.20E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[356] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.69E-11' 
		rate_values[357] = 1.69E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[358] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.78E-11' 
		rate_values[359] = 3.78E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]' 
		rate_values[360] = J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.41E-11' 
		rate_values[361] = 2.41E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[362] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.66E-11' 
		rate_values[363] = 7.66E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[364] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.92E-11*0.232' 
		rate_values[365] = 8.92E-11*0.232
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.92E-11*0.768' 
		rate_values[366] = 8.92E-11*0.768
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[367] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[368] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[369] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[370] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[371] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[372] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[373] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[374] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[375] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[376] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2*0.7' 
		rate_values[377] = 1.00E-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2*0.3' 
		rate_values[378] = 1.00E-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.29E-11' 
		rate_values[379] = 2.29E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[13]*2.' 
		rate_values[380] = J[19]+J[13]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.70E-11' 
		rate_values[381] = 2.70E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]*2.' 
		rate_values[382] = J[13]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.39E-11' 
		rate_values[383] = 2.39E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[13]*2.' 
		rate_values[384] = J[19]+J[13]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.23E-11' 
		rate_values[385] = 3.23E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]*2.' 
		rate_values[386] = J[13]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.54E-11*0.890' 
		rate_values[387] = 2.54E-11*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.54E-11*0.110' 
		rate_values[388] = 2.54E-11*0.110
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[389] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.55E-11' 
		rate_values[390] = 3.55E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.60E-12' 
		rate_values[391] = 7.60E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[13]' 
		rate_values[392] = J[19]+J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.65E-11' 
		rate_values[393] = 2.65E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]' 
		rate_values[394] = J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.34E-11' 
		rate_values[395] = 2.34E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[396] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[397] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[398] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[399] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[400] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2*0.7' 
		rate_values[401] = 1.00E-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2*0.3' 
		rate_values[402] = 1.00E-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[403] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[404] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*RO2*0.6' 
		rate_values[405] = 2.00E-12*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*RO2*0.2' 
		rate_values[406] = 2.00E-12*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*RO2*0.2' 
		rate_values[407] = 2.00E-12*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[408] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.64E-11' 
		rate_values[409] = 2.64E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*8.5' 
		rate_values[410] = KNO3AL*8.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[411] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[412] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[413] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[414] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[415] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[416] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[417] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]' 
		rate_values[418] = J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.22E-12' 
		rate_values[419] = 3.22E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[420] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.33E-11' 
		rate_values[421] = 1.33E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[422] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[423] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-11' 
		rate_values[424] = 6.70E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[425] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[426] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[427] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[428] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[429] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[430] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[431] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[432] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2' 
		rate_values[433] = 1.00E-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[434] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[435] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*RO2' 
		rate_values[436] = 6.70E-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[437] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[438] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*RO2' 
		rate_values[439] = 6.70E-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[440] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[441] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2' 
		rate_values[442] = 8.80E-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[443] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[444] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[445] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[446] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2' 
		rate_values[447] = 1.00E-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[448] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[449] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*RO2' 
		rate_values[450] = 2.00E-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.60E-12' 
		rate_values[451] = 6.60E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[13]' 
		rate_values[452] = J[19]+J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.02E-11' 
		rate_values[453] = 1.02E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[13]' 
		rate_values[454] = J[19]+J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.29E-11' 
		rate_values[455] = 1.29E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[11]' 
		rate_values[456] = J[19]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.45E-11' 
		rate_values[457] = 3.45E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[11]' 
		rate_values[458] = J[19]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.49E-11' 
		rate_values[459] = 7.49E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.83E-13' 
		rate_values[460] = 8.83E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[18]' 
		rate_values[461] = J[19]+J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.75E-12' 
		rate_values[462] = 4.75E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[18]' 
		rate_values[463] = J[19]+J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.01E-11' 
		rate_values[464] = 1.01E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]*2.' 
		rate_values[465] = J[11]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.33E-10' 
		rate_values[466] = 1.33E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*KNO3AL*5.5' 
		rate_values[467] = 2.*KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.770' 
		rate_values[468] = KRO2NO*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.230' 
		rate_values[469] = KRO2NO*0.230
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[470] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[471] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[472] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.6' 
		rate_values[473] = 8.80E-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.2' 
		rate_values[474] = 8.80E-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.2' 
		rate_values[475] = 8.80E-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.64E-12' 
		rate_values[476] = 3.64E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[477] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.28E-11' 
		rate_values[478] = 3.28E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.49E-11' 
		rate_values[479] = 1.49E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.18E-12' 
		rate_values[480] = 8.18E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2E-12*numpy.exp(490./TEMP)*0.65' 
		rate_values[481] = 1.2E-12*numpy.exp(490./TEMP)*0.65
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2E-12*numpy.exp(490./TEMP)*0.35' 
		rate_values[482] = 1.2E-12*numpy.exp(490./TEMP)*0.35
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[483] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[484] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[485] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[486] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*0.9*RO2' 
		rate_values[487] = 6.70E-15*0.9*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*0.1*RO2' 
		rate_values[488] = 6.70E-15*0.1*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[489] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.87E-12' 
		rate_values[490] = 6.87E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[491] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.00E+05' 
		rate_values[492] = 4.00E+05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROSEC*O2' 
		rate_values[493] = KROSEC*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[494] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[495] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50E-13*0.8*RO2' 
		rate_values[496] = 2.50E-13*0.8*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50E-13*0.1*RO2' 
		rate_values[497] = 2.50E-13*0.1*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50E-13*0.1*RO2' 
		rate_values[498] = 2.50E-13*0.1*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[499] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[500] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.23E-11' 
		rate_values[501] = 1.23E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[502] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.55E-12' 
		rate_values[503] = 5.55E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[504] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[505] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[506] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[507] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[508] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.00E+04' 
		rate_values[509] = 4.00E+04
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROSEC*O2' 
		rate_values[510] = KROSEC*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[511] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[512] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*RO2' 
		rate_values[513] = 6.70E-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[514] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[515] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*RO2' 
		rate_values[516] = 6.70E-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[517] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[518] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50E-13*RO2' 
		rate_values[519] = 2.50E-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[520] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.94E-12' 
		rate_values[521] = 5.94E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[522] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.03E-12' 
		rate_values[523] = 8.03E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[524] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.25E-11' 
		rate_values[525] = 3.25E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[526] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.25E-12' 
		rate_values[527] = 1.25E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[528] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[529] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[530] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[531] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[532] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2' 
		rate_values[533] = 8.80E-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[534] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[535] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[536] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[537] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2' 
		rate_values[538] = 1.00E-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[539] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.85E-11' 
		rate_values[540] = 1.85E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.32E-11' 
		rate_values[541] = 1.32E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[542] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.68E-11' 
		rate_values[543] = 1.68E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.875' 
		rate_values[544] = KRO2NO*0.875
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.125' 
		rate_values[545] = KRO2NO*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[546] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[547] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[548] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*RO2*0.7' 
		rate_values[549] = 6.70E-15*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*RO2*0.3' 
		rate_values[550] = 6.70E-15*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]' 
		rate_values[551] = J[23]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.87E-11' 
		rate_values[552] = 9.87E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[553] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.03E-10' 
		rate_values[554] = 1.03E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.91E-11' 
		rate_values[555] = 9.91E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.722' 
		rate_values[556] = KRO2NO*0.722
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.278' 
		rate_values[557] = KRO2NO*0.278
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROSEC*O2' 
		rate_values[558] = KROSEC*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[559] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[560] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50E-13*RO2*0.6' 
		rate_values[561] = 2.50E-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50E-13*RO2*0.2' 
		rate_values[562] = 2.50E-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50E-13*RO2*0.2' 
		rate_values[563] = 2.50E-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[564] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.60E-11' 
		rate_values[565] = 9.60E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[566] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.27E-10' 
		rate_values[567] = 1.27E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.09E-10' 
		rate_values[568] = 1.09E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.19E-10' 
		rate_values[569] = 1.19E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.958' 
		rate_values[570] = KRO2NO*0.958
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.042' 
		rate_values[571] = KRO2NO*0.042
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[572] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[573] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[574] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20E-14*RO2*0.7' 
		rate_values[575] = 9.20E-14*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20E-14*RO2*0.3' 
		rate_values[576] = 9.20E-14*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.26E-11' 
		rate_values[577] = 1.26E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[578] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.06E-11' 
		rate_values[579] = 7.06E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.72E-11' 
		rate_values[580] = 6.72E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.05E-16*numpy.exp(-640./TEMP)*0.6' 
		rate_values[581] = 8.05E-16*numpy.exp(-640./TEMP)*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.05E-16*numpy.exp(-640./TEMP)*0.4' 
		rate_values[582] = 8.05E-16*numpy.exp(-640./TEMP)*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.55*0.90' 
		rate_values[583] = KDEC*0.55*0.90
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.45*0.90' 
		rate_values[584] = KDEC*0.45*0.90
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[585] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[586] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.875' 
		rate_values[587] = KRO2NO*0.875
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.125' 
		rate_values[588] = KRO2NO*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[589] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.722' 
		rate_values[590] = KRO2NO*0.722
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.278' 
		rate_values[591] = KRO2NO*0.278
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROSEC*O2' 
		rate_values[592] = KROSEC*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[593] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[594] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20E-14*0.7*RO2*0.9' 
		rate_values[595] = 9.20E-14*0.7*RO2*0.9
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20E-14*0.3*RO2*0.9' 
		rate_values[596] = 9.20E-14*0.3*RO2*0.9
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[597] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[598] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*0.7*RO2' 
		rate_values[599] = 6.70E-15*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*0.3*RO2' 
		rate_values[600] = 6.70E-15*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[601] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[602] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50E-13*0.6*RO2' 
		rate_values[603] = 2.50E-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50E-13*0.2*RO2' 
		rate_values[604] = 2.50E-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50E-13*0.2*RO2' 
		rate_values[605] = 2.50E-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]+J[18]' 
		rate_values[606] = J[23]+J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.85E-11' 
		rate_values[607] = 2.85E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]+J[18]' 
		rate_values[608] = J[22]+J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.23E-11' 
		rate_values[609] = 2.23E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[11]' 
		rate_values[610] = J[19]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.01E-11' 
		rate_values[611] = 3.01E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[18]' 
		rate_values[612] = J[19]+J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.28E-11' 
		rate_values[613] = 6.28E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[18]' 
		rate_values[614] = J[19]+J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-10' 
		rate_values[615] = 2.00E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[616] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.66E-11' 
		rate_values[617] = 2.66E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[618] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.93E-11' 
		rate_values[619] = 5.93E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[620] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.26E-10' 
		rate_values[621] = 1.26E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[622] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.80' 
		rate_values[623] = KDEC*0.80
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.20' 
		rate_values[624] = KDEC*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[625] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[626] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*RO2*0.90' 
		rate_values[627] = 2.00E-12*RO2*0.90
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*RO2*0.05' 
		rate_values[628] = 2.00E-12*RO2*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*RO2*0.05' 
		rate_values[629] = 2.00E-12*RO2*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[630] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[631] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.20E+10*numpy.exp(-3523./TEMP)' 
		rate_values[632] = 4.20E+10*numpy.exp(-3523./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[633] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[634] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[635] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[636] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[637] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[638] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[639] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[640] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.37E-11' 
		rate_values[641] = 2.37E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*5.5' 
		rate_values[642] = KNO3AL*5.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[643] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[644] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[645] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[646] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[647] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[648] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[649] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2' 
		rate_values[650] = 1.00E-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.92E-12' 
		rate_values[651] = 2.92E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[18]' 
		rate_values[652] = J[19]+J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.55E-12' 
		rate_values[653] = 6.55E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[654] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[655] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*RO2' 
		rate_values[656] = 2.00E-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[18]' 
		rate_values[657] = J[19]+J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.61E-12' 
		rate_values[658] = 9.61E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[659] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[660] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[661] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[662] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[663] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2*0.7' 
		rate_values[664] = 1.00E-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2*0.3' 
		rate_values[665] = 1.00E-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[666] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[667] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30E-12*RO2' 
		rate_values[668] = 1.30E-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[669] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[670] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*RO2' 
		rate_values[671] = 6.70E-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[672] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.890' 
		rate_values[673] = KRO2HO2*0.890
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*RO2' 
		rate_values[674] = 6.70E-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[675] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[676] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2' 
		rate_values[677] = 8.80E-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[11]' 
		rate_values[678] = J[19]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]' 
		rate_values[679] = J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.47E-11' 
		rate_values[680] = 5.47E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]' 
		rate_values[681] = J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[682] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.45E-11' 
		rate_values[683] = 4.45E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]+J[11]' 
		rate_values[684] = J[17]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.47E-11' 
		rate_values[685] = 5.47E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.56E-12' 
		rate_values[686] = 5.56E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]' 
		rate_values[687] = J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.70E-12' 
		rate_values[688] = 5.70E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[13]' 
		rate_values[689] = J[19]+J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.16E-12' 
		rate_values[690] = 9.16E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[13]' 
		rate_values[691] = J[19]+J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.36E-11' 
		rate_values[692] = 2.36E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[13]' 
		rate_values[693] = J[19]+J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.29E-11' 
		rate_values[694] = 1.29E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[13]' 
		rate_values[695] = J[19]+J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.51E-11' 
		rate_values[696] = 1.51E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[13]' 
		rate_values[697] = J[19]+J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-10' 
		rate_values[698] = 1.00E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.50' 
		rate_values[699] = KDEC*0.50
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.50' 
		rate_values[700] = KDEC*0.50
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00E-14' 
		rate_values[701] = 7.00E-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.20E-15' 
		rate_values[702] = 1.20E-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-14' 
		rate_values[703] = 1.00E-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-15' 
		rate_values[704] = 1.00E-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.40E-17*H2O' 
		rate_values[705] = 1.40E-17*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-18*H2O' 
		rate_values[706] = 2.00E-18*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[707] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]' 
		rate_values[708] = J[18]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]*2.' 
		rate_values[709] = J[11]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.77E-11' 
		rate_values[710] = 5.77E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[711] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.12E-12' 
		rate_values[712] = 1.12E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[713] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.896' 
		rate_values[714] = KRO2NO*0.896
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.104' 
		rate_values[715] = KRO2NO*0.104
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.70E+14*numpy.exp(-6643./TEMP)' 
		rate_values[716] = 2.70E+14*numpy.exp(-6643./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.896' 
		rate_values[717] = KRO2NO*0.896
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.104' 
		rate_values[718] = KRO2NO*0.104
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.70E+14*numpy.exp(-6643./TEMP)' 
		rate_values[719] = 2.70E+14*numpy.exp(-6643./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.871' 
		rate_values[720] = KRO2NO*0.871
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.129' 
		rate_values[721] = KRO2NO*0.129
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.80E-14*numpy.exp(-260./TEMP)*O2' 
		rate_values[722] = 1.80E-14*numpy.exp(-260./TEMP)*O2
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
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[726] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[727] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2*0.7' 
		rate_values[728] = 1.00E-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2*0.3' 
		rate_values[729] = 1.00E-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[730] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[731] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*RO2*0.7' 
		rate_values[732] = 6.70E-15*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*RO2*0.3' 
		rate_values[733] = 6.70E-15*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[734] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.914' 
		rate_values[735] = KRO2HO2*0.914
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*RO2*0.7' 
		rate_values[736] = 6.70E-15*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*RO2*0.3' 
		rate_values[737] = 6.70E-15*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[738] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[739] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50E-13*RO2*0.6' 
		rate_values[740] = 2.50E-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50E-13*RO2*0.2' 
		rate_values[741] = 2.50E-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50E-13*RO2*0.2' 
		rate_values[742] = 2.50E-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.52E-11' 
		rate_values[743] = 2.52E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[11]' 
		rate_values[744] = J[19]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.00E-11' 
		rate_values[745] = 3.00E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[746] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.69E-11' 
		rate_values[747] = 2.69E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]+J[11]' 
		rate_values[748] = J[23]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.56E-11' 
		rate_values[749] = 2.56E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[11]' 
		rate_values[750] = J[19]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.61E-11' 
		rate_values[751] = 3.61E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[752] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86E-11' 
		rate_values[753] = 2.86E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]+J[11]' 
		rate_values[754] = J[23]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.96E-11' 
		rate_values[755] = 4.96E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[11]' 
		rate_values[756] = J[19]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.35E-11' 
		rate_values[757] = 8.35E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[758] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00E-11' 
		rate_values[759] = 8.00E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]+J[11]*2.' 
		rate_values[760] = J[22]+J[11]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.33E-11' 
		rate_values[761] = 4.33E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[11]*2.' 
		rate_values[762] = J[19]+J[11]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.10E-10' 
		rate_values[763] = 1.10E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]*2.' 
		rate_values[764] = J[11]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.99E-11' 
		rate_values[765] = 6.99E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.20' 
		rate_values[766] = KDEC*0.20
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.80' 
		rate_values[767] = KDEC*0.80
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[768] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.862' 
		rate_values[769] = KRO2NO*0.862
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[770] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[771] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[772] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.896' 
		rate_values[773] = KRO2NO*0.896
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.104' 
		rate_values[774] = KRO2NO*0.104
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[775] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[776] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[777] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[778] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2*0.7' 
		rate_values[779] = 1.00E-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2*0.3' 
		rate_values[780] = 1.00E-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[781] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[782] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[783] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[784] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30E-12*RO2*0.6' 
		rate_values[785] = 1.30E-12*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30E-12*RO2*0.2' 
		rate_values[786] = 1.30E-12*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30E-12*RO2*0.2' 
		rate_values[787] = 1.30E-12*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[788] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[789] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20E-14*RO2*0.7' 
		rate_values[790] = 9.20E-14*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20E-14*RO2*0.3' 
		rate_values[791] = 9.20E-14*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[792] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.859' 
		rate_values[793] = KRO2HO2*0.859
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*RO2*0.7' 
		rate_values[794] = 6.70E-15*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*RO2*0.3' 
		rate_values[795] = 6.70E-15*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[23]+J[17]' 
		rate_values[796] = J[23]+J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.82E-12' 
		rate_values[797] = 7.82E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[798] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.04E-11' 
		rate_values[799] = 1.04E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[800] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.70E-11' 
		rate_values[801] = 1.70E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[802] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.09E-11' 
		rate_values[803] = 1.09E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[17]' 
		rate_values[804] = J[19]+J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.86E-11' 
		rate_values[805] = 1.86E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.29E-12' 
		rate_values[806] = 7.29E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.89E-12' 
		rate_values[807] = 7.89E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.42E-12' 
		rate_values[808] = 7.42E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.75E-11' 
		rate_values[809] = 1.75E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[810] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.77E-12' 
		rate_values[811] = 6.77E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]+J[11]' 
		rate_values[812] = J[17]+J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.69E-11' 
		rate_values[813] = 6.69E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[814] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.63E-11' 
		rate_values[815] = 2.63E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*8.5' 
		rate_values[816] = KNO3AL*8.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[817] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[818] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[819] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[820] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[821] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[822] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[823] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[824] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[825] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[826] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[827] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[828] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2*0.7' 
		rate_values[829] = 1.00E-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2*0.3' 
		rate_values[830] = 1.00E-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[831] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[832] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.30E-12*RO2' 
		rate_values[833] = 1.30E-12*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[834] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.820' 
		rate_values[835] = KRO2HO2*0.820
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.70E-15*RO2' 
		rate_values[836] = 6.70E-15*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[837] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[838] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2' 
		rate_values[839] = 8.80E-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.96E-12' 
		rate_values[840] = 2.96E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[841] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.65E-12' 
		rate_values[842] = 9.65E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.57E-12' 
		rate_values[843] = 6.57E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[844] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.27E-11' 
		rate_values[845] = 1.27E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[846] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.31E-11' 
		rate_values[847] = 3.31E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]' 
		rate_values[848] = J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.46E-11' 
		rate_values[849] = 7.46E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.138' 
		rate_values[850] = KRO2NO*0.138
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[21]' 
		rate_values[851] = J[21]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.29E-12' 
		rate_values[852] = 3.29E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[853] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[854] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[855] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[856] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[857] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[858] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[859] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[860] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[861] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[862] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[863] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[864] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[865] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[866] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[867] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[868] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[869] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[870] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[871] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[872] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[873] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[874] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.8E-13*RO2' 
		rate_values[875] = 8.8E-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[876] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[19]+J[13]' 
		rate_values[877] = J[19]+J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.38E-11' 
		rate_values[878] = 3.38E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00E11*numpy.exp(-3160./TEMP)+5.00E-12*O2' 
		rate_values[879] = 7.00E11*numpy.exp(-3160./TEMP)+5.00E-12*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00E-12*O2*3.2*numpy.exp(-550./TEMP)' 
		rate_values[880] = 5.00E-12*O2*3.2*numpy.exp(-550./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00E-12*O2*3.2*(1.-numpy.exp(-550./TEMP))' 
		rate_values[881] = 5.00E-12*O2*3.2*(1.-numpy.exp(-550./TEMP))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[24]*0.91' 
		rate_values[882] = J[24]*0.91
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.' 
		rate_values[883] = 5.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.10' 
		rate_values[884] = KDEC*0.10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5*0.7' 
		rate_values[885] = KDEC*0.5*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5*0.7' 
		rate_values[886] = KDEC*0.5*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.30' 
		rate_values[887] = KDEC*0.30
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.025*1.2E-11*numpy.exp(440/TEMP)' 
		rate_values[888] = 0.025*1.2E-11*numpy.exp(440/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.01*2.38E-11*numpy.exp(357/TEMP)' 
		rate_values[889] = 0.01*2.38E-11*numpy.exp(357/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.01*4.28e-11*numpy.exp(401/TEMP)' 
		rate_values[890] = 0.01*4.28e-11*numpy.exp(401/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.025*1.60e-11*numpy.exp(500./TEMP)' 
		rate_values[891] = 0.025*1.60e-11*numpy.exp(500./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e17*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[892] = 1.2e17*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e18*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[893] = 1.2e18*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e18*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[894] = 1.2e18*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6e17*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[895] = 6e17*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2e17*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[896] = 2e17*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e17*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[897] = 1.2e17*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e17*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[898] = 1.2e17*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2e16*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[899] = 1.2e16*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[900] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[901] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[902] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[903] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[904] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[905] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[906] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[907] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[908] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[909] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[910] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*1.0' 
		rate_values[911] = KRO2NO*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*1.0' 
		rate_values[912] = KRO2NO*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*1.0' 
		rate_values[913] = KRO2NO*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.9' 
		rate_values[914] = KRO2NO*0.9
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.85' 
		rate_values[915] = KRO2NO*0.85
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.8' 
		rate_values[916] = KRO2NO*0.8
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.6' 
		rate_values[917] = KRO2NO*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.5' 
		rate_values[918] = KRO2NO*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.3' 
		rate_values[919] = KRO2NO*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[920] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[921] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[922] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[923] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[924] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.1' 
		rate_values[925] = 0.4*KRO2NO*0.1
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.15' 
		rate_values[926] = 0.4*KRO2NO*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.2' 
		rate_values[927] = 0.4*KRO2NO*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.4' 
		rate_values[928] = 0.4*KRO2NO*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.5' 
		rate_values[929] = 0.4*KRO2NO*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.7' 
		rate_values[930] = 0.4*KRO2NO*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*1.0' 
		rate_values[931] = 0.4*KRO2NO*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*1.0' 
		rate_values[932] = 0.4*KRO2NO*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[933] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[934] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[935] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.1*7e0/10e0' 
		rate_values[936] = 0.6*KRO2NO*0.1*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.15*7e0/10e0' 
		rate_values[937] = 0.6*KRO2NO*0.15*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.2*7e0/10e0' 
		rate_values[938] = 0.6*KRO2NO*0.2*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.4*7e0/10e0' 
		rate_values[939] = 0.6*KRO2NO*0.4*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.5*7e0/10e0' 
		rate_values[940] = 0.6*KRO2NO*0.5*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.7*7e0/10e0' 
		rate_values[941] = 0.6*KRO2NO*0.7*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*1.0*7e0/10e0' 
		rate_values[942] = 0.6*KRO2NO*1.0*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*1.0*7e0/10e0' 
		rate_values[943] = 0.6*KRO2NO*1.0*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[944] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[945] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[946] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.1*3e0/10e0' 
		rate_values[947] = 0.6*KRO2NO*0.1*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.15*3e0/10e0' 
		rate_values[948] = 0.6*KRO2NO*0.15*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.2*3e0/10e0' 
		rate_values[949] = 0.6*KRO2NO*0.2*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.4*3e0/10e0' 
		rate_values[950] = 0.6*KRO2NO*0.4*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.5*3e0/10e0' 
		rate_values[951] = 0.6*KRO2NO*0.5*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.7*3e0/10e0' 
		rate_values[952] = 0.6*KRO2NO*0.7*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*1.0*3e0/10e0' 
		rate_values[953] = 0.6*KRO2NO*1.0*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*1.0*3e0/10e0' 
		rate_values[954] = 0.6*KRO2NO*1.0*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[955] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[956] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[957] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[958] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[959] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[960] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[961] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[962] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[963] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[964] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[965] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[966] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[967] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[968] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[969] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[970] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[971] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[972] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[973] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[974] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[975] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[976] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[977] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[978] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[979] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[980] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[981] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[982] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[983] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[984] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[985] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[986] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[987] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[988] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[989] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[990] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[991] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[992] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[993] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[994] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[995] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[996] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[997] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[998] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[999] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1000] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1001] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1002] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1003] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1004] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1005] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1006] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1007] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1008] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1009] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1010] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1011] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1012] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1013] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1014] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1015] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1016] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1017] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1018] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1019] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1020] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1021] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1022] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1023] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1024] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1025] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1026] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1027] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1028] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1029] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1030] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1031] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1032] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1033] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1034] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1035] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1036] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1037] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1038] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1039] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1040] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1041] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1042] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1043] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1044] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1045] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1046] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1047] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1048] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1049] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1050] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1051] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1052] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1053] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1054] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1055] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1056] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1057] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1058] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1059] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1060] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1061] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1062] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1063] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1064] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1065] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1066] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1067] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1068] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1069] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1070] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1071] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1072] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1073] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1074] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1075] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1076] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1077] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1078] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1079] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1080] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1081] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1082] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1083] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1084] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1085] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1086] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1087] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1088] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1089] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1090] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1091] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1092] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1093] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1094] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1095] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1096] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1097] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1098] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1099] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1100] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1101] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1102] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1103] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1104] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1105] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1106] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1107] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1108] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1109] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1110] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1111] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1112] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1113] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1114] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1115] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1116] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1117] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1118] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1119] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1120] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1121] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1122] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1123] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1124] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1125] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1126] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1127] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1128] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1129] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1130] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1131] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1132] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1133] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1134] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1135] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1136] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1137] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1138] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1139] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1140] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1141] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1142] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1143] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1144] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1145] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1146] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1147] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1148] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1149] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1150] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1151] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1152] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1153] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1154] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1155] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1156] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1157] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1158] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1159] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1160] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1161] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1162] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1163] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1164] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1165] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1166] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1167] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1168] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1169] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1170] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1171] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1172] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1173] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1174] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1175] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1176] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1177] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1178] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1179] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1180] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1181] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1182] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1183] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1184] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1185] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1186] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1187] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1188] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1189] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1190] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1191] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1192] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1193] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1194] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1195] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1196] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1197] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1198] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1199] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1200] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1201] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1202] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1203] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1204] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1205] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1206] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1207] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1208] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1209] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1210] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1211] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1212] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1213] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1214] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1215] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1216] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1217] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1218] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1219] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1220] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1221] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1222] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1223] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1224] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1225] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1226] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1227] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1228] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1229] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1230] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1231] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1232] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1233] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1234] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1235] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1236] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1237] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1238] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1239] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1240] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1241] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1242] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1243] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1244] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1245] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1246] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[1247] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1248] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1249] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1250] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1251] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1252] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1253] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1254] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1255] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1256] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1257] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1258] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1259] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1260] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1261] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1262] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1263] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1264] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1265] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1266] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1267] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1268] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1269] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1270] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1271] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1272] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1273] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1274] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1275] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1276] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1277] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1278] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1279] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1280] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1281] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1282] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1283] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1284] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1285] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1286] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1287] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1288] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1289] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1290] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1291] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1292] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1293] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1294] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1295] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1296] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1297] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1298] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1299] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1300] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1301] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1302] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1303] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1304] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1305] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1306] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1307] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1308] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1309] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1310] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1311] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1312] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1313] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1314] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1315] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1316] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1317] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1318] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1319] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1320] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1321] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1322] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1323] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1324] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1325] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1326] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1327] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1328] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1329] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1330] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1331] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1332] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1333] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1334] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1335] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1336] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1337] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1338] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1339] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1340] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[1341] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1342] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1343] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1344] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1345] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1346] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1347] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1348] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1349] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1350] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1351] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1352] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1353] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1354] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1355] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1356] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1357] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1358] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1359] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1360] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1361] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1362] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1363] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1364] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1365] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1366] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1367] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1368] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1369] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1370] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1371] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1372] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1373] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1374] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1375] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1376] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1377] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1378] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1379] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1380] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1381] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1382] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1383] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1384] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1385] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1386] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1387] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1388] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1389] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1390] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1391] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1392] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1393] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1394] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1395] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1396] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1397] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1398] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1399] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1400] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1401] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1402] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1403] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1404] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1405] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1406] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1407] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1408] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1409] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1410] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1411] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1412] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1413] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1414] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1415] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1416] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1417] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1418] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1419] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1420] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1421] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1422] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1423] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1424] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1425] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1426] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1427] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1428] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1429] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1430] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1431] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1432] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1433] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1434] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[1435] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1436] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1437] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1438] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1439] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1440] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1441] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1442] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1443] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1444] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1445] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1446] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1447] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1448] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1449] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1450] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1451] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1452] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1453] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1454] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1455] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1456] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1457] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1458] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1459] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1460] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1461] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1462] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1463] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1464] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1465] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1466] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1467] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1468] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1469] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1470] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1471] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1472] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1473] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1474] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1475] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1476] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1477] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1478] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1479] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1480] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1481] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1482] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1483] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1484] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1485] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1486] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1487] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1488] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1489] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1490] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1491] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1492] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1493] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1494] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1495] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1496] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1497] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1498] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1499] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1500] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1501] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1502] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1503] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1504] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1505] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1506] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1507] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1508] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1509] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1510] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1511] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1512] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1513] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1514] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1515] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1516] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1517] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1518] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1519] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1520] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1521] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1522] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1523] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1524] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1525] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1526] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1527] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1528] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1529] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1530] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1531] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1532] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1533] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1534] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1535] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1536] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1537] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1538] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1539] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1540] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1541] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1542] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1543] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1544] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1545] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1546] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1547] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1548] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1549] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1550] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1551] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1552] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1553] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1554] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1555] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1556] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1557] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1558] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1559] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1560] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1561] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1562] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1563] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1564] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1565] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1566] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1567] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1568] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1569] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1570] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1571] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1572] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1573] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1574] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1575] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1576] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1577] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1578] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1579] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1580] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1581] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1582] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1583] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1584] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1585] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1586] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1587] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1588] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1589] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1590] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1591] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1592] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1593] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1594] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1595] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1596] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1597] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1598] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1599] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1600] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1601] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1602] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1603] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1604] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1605] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1606] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1607] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1608] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1609] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1610] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1611] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1612] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1613] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1614] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1615] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1616] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1617] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1618] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1619] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1620] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1621] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1622] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[1623] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1624] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1625] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1626] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1627] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1628] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1629] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1630] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1631] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1632] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1633] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1634] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1635] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1636] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1637] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1638] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1639] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1640] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1641] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1642] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1643] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1644] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1645] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1646] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1647] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1648] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1649] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1650] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1651] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1652] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1653] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1654] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1655] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1656] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1657] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1658] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1659] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1660] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1661] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1662] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1663] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1664] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1665] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1666] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1667] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1668] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1669] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1670] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1671] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1672] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1673] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1674] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1675] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1676] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1677] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1678] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1679] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1680] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1681] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1682] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1683] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1684] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1685] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1686] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1687] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1688] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1689] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1690] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1691] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1692] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1693] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1694] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1695] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1696] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1697] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1698] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1699] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1700] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1701] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1702] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1703] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1704] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1705] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1706] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1707] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1708] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1709] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1710] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1711] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1712] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1713] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1714] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1715] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1716] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1717] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1718] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1719] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1720] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1721] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1722] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1723] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1724] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1725] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1726] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1727] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1728] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1729] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1730] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1731] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1732] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1733] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1734] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1735] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1736] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1737] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1738] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1739] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1740] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1741] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1742] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1743] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1744] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1745] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1746] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1747] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1748] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1749] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1750] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1751] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1752] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1753] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1754] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1755] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1756] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1757] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1758] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1759] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1760] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1761] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1762] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1763] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1764] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1765] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1766] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1767] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1768] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1769] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1770] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1771] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1772] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1773] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1774] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1775] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1776] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1777] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1778] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1779] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1780] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1781] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1782] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1783] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1784] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1785] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1786] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1787] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1788] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1789] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1790] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1791] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1792] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1793] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1794] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1795] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1796] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1797] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1798] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1799] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1800] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1801] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1802] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1803] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1804] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1805] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1806] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1807] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1808] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1809] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1810] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1811] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1812] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1813] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1814] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1815] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1816] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1817] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1818] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1819] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1820] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1821] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1822] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1823] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1824] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1825] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1826] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1827] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1828] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1829] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1830] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1831] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1832] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1833] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1834] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1835] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1836] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1837] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1838] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1839] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1840] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1841] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1842] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1843] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1844] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1845] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1846] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1847] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1848] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1849] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1850] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1851] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1852] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1853] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1854] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1855] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1856] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1857] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1858] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1859] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1860] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1861] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1862] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1863] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1864] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1865] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1866] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1867] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1868] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1869] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1870] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1871] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1872] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1873] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1874] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1875] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1876] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1877] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1878] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1879] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1880] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1881] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1882] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1883] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1884] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1885] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1886] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1887] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1888] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1889] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1890] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1891] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1892] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1893] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1894] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1895] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1896] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1897] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1898] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1899] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1900] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1901] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1902] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1903] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1904] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1905] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1906] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1907] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1908] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1909] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1910] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1911] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1912] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1913] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1914] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1915] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1916] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1917] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1918] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1919] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1920] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1921] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1922] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1923] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1924] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1925] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1926] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1927] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1928] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1929] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1930] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1931] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1932] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1933] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1934] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1935] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1936] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1937] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1938] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1939] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1940] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1941] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1942] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1943] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1944] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1945] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1946] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1947] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1948] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1949] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1950] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1951] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1952] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1953] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1954] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1955] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1956] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1957] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1958] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1959] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1960] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1961] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1962] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1963] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1964] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1965] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1966] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1967] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1968] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1969] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1970] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1971] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1972] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1973] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1974] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1975] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1976] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1977] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1978] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1979] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1980] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1981] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1982] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1983] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1984] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1985] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1986] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1987] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1988] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1989] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1990] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1991] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1992] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1993] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1994] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1995] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1996] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1997] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1998] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12' 
		rate_values[1999] = 5E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-12*0.2*RO2' 
		rate_values[2000] = 1E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-12*0.6*RO2' 
		rate_values[2001] = 1E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-12*0.2*RO2' 
		rate_values[2002] = 1E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[2003] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.6*RO2' 
		rate_values[2004] = 5E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[2005] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[2006] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.6*RO2' 
		rate_values[2007] = 5E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[2008] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[2009] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.6*RO2' 
		rate_values[2010] = 5E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[2011] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7E-12*0.4*RO2' 
		rate_values[2012] = 7E-12*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7E-12*0.4*RO2' 
		rate_values[2013] = 7E-12*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7E-12*0.2*RO2' 
		rate_values[2014] = 7E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8E-12*0.2*RO2' 
		rate_values[2015] = 8E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8E-12*0.4*RO2' 
		rate_values[2016] = 8E-12*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8E-12*0.4*RO2' 
		rate_values[2017] = 8E-12*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9E-12*0.4*RO2' 
		rate_values[2018] = 9E-12*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9E-12*0.2*RO2' 
		rate_values[2019] = 9E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9E-12*0.4*RO2' 
		rate_values[2020] = 9E-12*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.3*RO2' 
		rate_values[2021] = 1E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.2*RO2' 
		rate_values[2022] = 1E-11*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.5*RO2' 
		rate_values[2023] = 1E-11*0.5*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.5*RO2' 
		rate_values[2024] = 1E-11*0.5*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.2*RO2' 
		rate_values[2025] = 1E-11*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.3*RO2' 
		rate_values[2026] = 1E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.4*RO2' 
		rate_values[2027] = 1E-11*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.2*RO2' 
		rate_values[2028] = 1E-11*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.4*RO2' 
		rate_values[2029] = 1E-11*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.5*RO2' 
		rate_values[2030] = 1E-11*0.5*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.0*RO2' 
		rate_values[2031] = 1E-11*0.0*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.5*RO2' 
		rate_values[2032] = 1E-11*0.5*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20e-14*RO2*0.1' 
		rate_values[2033] = 9.20e-14*RO2*0.1
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2e17*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[2034] = 2e17*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6e16*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[2035] = 6e16*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6e16*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[2036] = 6e16*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3e16*numpy.exp(-1.2077e4/TEMP)' 
		rate_values[2037] = 3e16*numpy.exp(-1.2077e4/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2038] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2039] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2040] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2041] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[2042] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0*KRO2NO' 
		rate_values[2043] = 1.0*KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.9*KRO2NO' 
		rate_values[2044] = 0.9*KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.8*KRO2NO' 
		rate_values[2045] = 0.8*KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.7*KRO2NO' 
		rate_values[2046] = 0.7*KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.5*KRO2NO' 
		rate_values[2047] = 0.5*KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2048] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.1' 
		rate_values[2049] = 0.4*KRO2NO*0.1
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.2' 
		rate_values[2050] = 0.4*KRO2NO*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.3' 
		rate_values[2051] = 0.4*KRO2NO*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*0.5' 
		rate_values[2052] = 0.4*KRO2NO*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.4*KRO2NO*1.0' 
		rate_values[2053] = 0.4*KRO2NO*1.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2054] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.1*7e0/10e0' 
		rate_values[2055] = 0.6*KRO2NO*0.1*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.2*7e0/10e0' 
		rate_values[2056] = 0.6*KRO2NO*0.2*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.3*7e0/10e0' 
		rate_values[2057] = 0.6*KRO2NO*0.3*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.5*7e0/10e0' 
		rate_values[2058] = 0.6*KRO2NO*0.5*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*1.0*7e0/10e0' 
		rate_values[2059] = 0.6*KRO2NO*1.0*7e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2060] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.1*3e0/10e0' 
		rate_values[2061] = 0.6*KRO2NO*0.1*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.2*3e0/10e0' 
		rate_values[2062] = 0.6*KRO2NO*0.2*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.3*3e0/10e0' 
		rate_values[2063] = 0.6*KRO2NO*0.3*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*0.5*3e0/10e0' 
		rate_values[2064] = 0.6*KRO2NO*0.5*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.6*KRO2NO*1.0*3e0/10e0' 
		rate_values[2065] = 0.6*KRO2NO*1.0*3e0/10e0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2066] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2067] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2068] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2069] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2070] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2' 
		rate_values[2071] = KRO2HO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2072] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2073] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2074] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2075] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2076] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2077] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2078] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2079] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2080] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2081] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2082] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2083] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2084] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2085] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2086] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2087] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2088] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2089] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2090] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2091] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2092] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2093] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2094] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2095] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2096] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2097] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2098] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2099] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2100] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2101] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2102] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2103] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2104] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2105] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2106] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2107] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2108] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2109] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2110] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2111] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2112] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2113] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2114] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2115] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2116] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2117] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2118] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2119] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2120] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2121] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2122] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2123] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2124] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2125] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2126] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2127] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2128] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2129] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2130] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2131] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2132] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2133] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2134] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2135] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2136] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2137] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2138] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2139] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2140] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2141] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2142] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2143] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2144] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2145] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2146] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2147] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2148] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2149] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2150] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2151] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2152] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2153] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2154] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2155] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2156] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2157] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2158] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2159] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2160] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2161] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2162] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2163] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2164] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2165] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2166] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2167] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2168] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2169] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2170] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2171] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2172] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2173] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2174] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2175] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2176] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2177] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2178] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2179] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2180] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2181] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2182] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2183] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2184] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2185] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2186] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2187] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2188] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2189] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2190] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2191] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2192] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2193] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2194] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2195] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2196] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2197] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2198] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2199] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2200] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2201] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2202] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2203] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2204] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2205] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2206] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2207] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2208] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2209] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2210] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2211] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2212] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2213] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2214] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2215] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2216] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2217] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2218] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2219] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2220] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2221] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2222] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2223] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2224] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2225] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2226] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2227] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2228] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2229] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2230] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2231] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2232] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2233] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2234] = 1E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-13' 
		rate_values[2235] = 1E-13
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
		rc_eq_now = '5E-13' 
		rate_values[2260] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2261] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2262] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2263] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2264] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2265] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2266] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2267] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2268] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2269] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2270] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2271] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2272] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2273] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2274] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2275] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2276] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2277] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2278] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2279] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2280] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2281] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2282] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2283] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2284] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2285] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2286] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2287] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2288] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2289] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2290] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2291] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2292] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2293] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2294] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2295] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2296] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2297] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2298] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2299] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2300] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2301] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2302] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2303] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2304] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2305] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2306] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2307] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2308] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2309] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2310] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2311] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2312] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2313] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2314] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2315] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2316] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2317] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2318] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2319] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2320] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2321] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2322] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2323] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2324] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2325] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2326] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2327] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2328] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2329] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2330] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2331] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2332] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2333] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2334] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2335] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2336] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2337] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2338] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2339] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2340] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2341] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2342] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2343] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2344] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2345] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2346] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2347] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2348] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2349] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2350] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2351] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2352] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-13' 
		rate_values[2353] = 5E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2354] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2355] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2356] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2357] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2358] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2359] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2360] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2361] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2362] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2363] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2364] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2365] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2366] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2367] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2368] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2369] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2370] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2371] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2372] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2373] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2374] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2375] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2376] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2377] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2378] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2379] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2380] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2381] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2382] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2383] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2384] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2385] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2386] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2387] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2388] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2389] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2390] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2391] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2392] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2393] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2394] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2395] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2396] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2397] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2398] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2399] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2400] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2401] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2402] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2403] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2404] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2405] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2406] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2407] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2408] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2409] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2410] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2411] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2412] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2413] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2414] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2415] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2416] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2417] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2418] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2419] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2420] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2421] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2422] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2423] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2424] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2425] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2426] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2427] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2428] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2429] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2430] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2431] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2432] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2433] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2434] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2435] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2436] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2437] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2438] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2439] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2440] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2441] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2442] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2443] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2444] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2445] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2446] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2E-12' 
		rate_values[2447] = 2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2448] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2449] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2450] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2451] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2452] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2453] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2454] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2455] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2456] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2457] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2458] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2459] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2460] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2461] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2462] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2463] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2464] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2465] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2466] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2467] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2468] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2469] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2470] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2471] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2472] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2473] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2474] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2475] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2476] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2477] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2478] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2479] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2480] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2481] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2482] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2483] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2484] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2485] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2486] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2487] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2488] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2489] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2490] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2491] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2492] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2493] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2494] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2495] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2496] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2497] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2498] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2499] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2500] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2501] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2502] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2503] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2504] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2505] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2506] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2507] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2508] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2509] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2510] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2511] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2512] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2513] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2514] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2515] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2516] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2517] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2518] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2519] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2520] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2521] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2522] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2523] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2524] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2525] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2526] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2527] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2528] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2529] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2530] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2531] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2532] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2533] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2534] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2535] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2536] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2537] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2538] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2539] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2540] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2541] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2542] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2543] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2544] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2545] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2546] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2547] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2548] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2549] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2550] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2551] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2552] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2553] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2554] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2555] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2556] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2557] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2558] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2559] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2560] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2561] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2562] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2563] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2564] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2565] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2566] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2567] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2568] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2569] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2570] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2571] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2572] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2573] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2574] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2575] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2576] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2577] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2578] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2579] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2580] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2581] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2582] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2583] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2584] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2585] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2586] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2587] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2588] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2589] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2590] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2591] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2592] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2593] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2594] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2595] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2596] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2597] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2598] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2599] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2600] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2601] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2602] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2603] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2604] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2605] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2606] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2607] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2608] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2609] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2610] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2611] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2612] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2613] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2614] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2615] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2616] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2617] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2618] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2619] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2620] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2621] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2622] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2623] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2624] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2625] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2626] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2627] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2628] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2629] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2630] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2631] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2632] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2633] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2634] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3E-12' 
		rate_values[2635] = 3E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[2636] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.6*RO2' 
		rate_values[2637] = 5E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[2638] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[2639] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.6*RO2' 
		rate_values[2640] = 5E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5E-12*0.2*RO2' 
		rate_values[2641] = 5E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8E-12*0.2*RO2' 
		rate_values[2642] = 8E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8E-12*0.6*RO2' 
		rate_values[2643] = 8E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8E-12*0.2*RO2' 
		rate_values[2644] = 8E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.3*RO2' 
		rate_values[2645] = 1E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.4*RO2' 
		rate_values[2646] = 1E-11*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.3*RO2' 
		rate_values[2647] = 1E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.4*RO2' 
		rate_values[2648] = 1E-11*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.2*RO2' 
		rate_values[2649] = 1E-11*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.4*RO2' 
		rate_values[2650] = 1E-11*0.4*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.5*RO2' 
		rate_values[2651] = 1E-11*0.5*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1E-11*0.5*RO2' 
		rate_values[2652] = 1E-11*0.5*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2653] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2654] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2655] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2656] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2657] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2658] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2659] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2660] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2661] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2662] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2663] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2664] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2665] = 0.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '0.0' 
		rate_values[2666] = 0.0
	except:
		erf = 1 # flag error
		err_mess = (str('Error: Could not calculate rate coefficient for equation number ' + str(gprn) + ' ' + rc_eq_now + ' (message from rate coeffs.py)'))
	
	# aqueous-phase reactions
	
	# surface (e.g. wall) reactions
	
	return(rate_values, erf, err_mess)
