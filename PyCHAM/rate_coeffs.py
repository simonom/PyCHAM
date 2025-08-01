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
# created at 2025-08-01 14:55:45.036456

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
		K14ISOM1=3.00E7*numpy.exp(-5300./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		K298CH3O2=3.5E-13
 
		# keep count on reaction number 
		gprn += 1 
		KAPHO2=5.2E-13*numpy.exp(980./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		KAPNO=7.5E-12*numpy.exp(290./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		KCH3O2=1.03E-13*numpy.exp(365./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		KDEC=1.00E+06
 
		# keep count on reaction number 
		gprn += 1 
		KMT05=1.44E-13*(1.+(M/4.2E+19))
 
		# keep count on reaction number 
		gprn += 1 
		KMT06=1.+(1.40E-21*numpy.exp(2200./TEMP)*H2O)
 
		# keep count on reaction number 
		gprn += 1 
		KMT18=9.5E-39*O2*numpy.exp(5270./TEMP)/(1.+7.5E-29*O2*numpy.exp(5610./TEMP))
 
		# keep count on reaction number 
		gprn += 1 
		KNO3AL=1.44E-12*numpy.exp(-1862./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		KRO2HO2=2.91E-13*numpy.exp(1300./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		KRO2NO=2.7E-12*numpy.exp(360./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		KRO2NO3=2.3E-12
 
		# keep count on reaction number 
		gprn += 1 
		KROPRIM=2.50E-14*numpy.exp(-300./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		KROSEC=2.50E-14*numpy.exp(-300./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		FCD=0.30
 
		# keep count on reaction number 
		gprn += 1 
		KD0=1.10E-05*M*numpy.exp(-10100./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		KDI=1.90E17*numpy.exp(-14100./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		KRD=KD0/KDI
 
		# keep count on reaction number 
		gprn += 1 
		NCD=0.75-1.27*(numpy.log10(FCD))
 
		# keep count on reaction number 
		gprn += 1 
		FD=10.**(numpy.log10(FCD)/(1.+(numpy.log10(KRD)/NCD)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KBPAN=(KD0*KDI)*FD/(KD0+KDI)
 
		# keep count on reaction number 
		gprn += 1 
		FCPPN=0.36
 
		# keep count on reaction number 
		gprn += 1 
		KPPN0=1.7E-03*numpy.exp(-11280./TEMP)*M
 
		# keep count on reaction number 
		gprn += 1 
		KPPNI=8.3E+16*numpy.exp(-13940./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		KRPPN=KPPN0/KPPNI
 
		# keep count on reaction number 
		gprn += 1 
		NCPPN=0.75-1.27*(numpy.log10(FCPPN))
 
		# keep count on reaction number 
		gprn += 1 
		FPPN=10.**(numpy.log10(FCPPN)/(1.+(numpy.log10(KRPPN)/NCPPN)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KBPPN=(KPPN0*KPPNI)*FPPN/(KPPN0+KPPNI)
 
		# keep count on reaction number 
		gprn += 1 
		FCC=0.30
 
		# keep count on reaction number 
		gprn += 1 
		KC0=3.28E-28*M*(TEMP/300.)**(-6.87)
 
		# keep count on reaction number 
		gprn += 1 
		KCI=1.125E-11*(TEMP/300.)**(-1.105)
 
		# keep count on reaction number 
		gprn += 1 
		KRC=KC0/KCI
 
		# keep count on reaction number 
		gprn += 1 
		NC=0.75-1.27*(numpy.log10(FCC))
 
		# keep count on reaction number 
		gprn += 1 
		FC=10.**(numpy.log10(FCC)/(1.+(numpy.log10(KRC)/NC)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KFPAN=(KC0*KCI)*FC/(KC0+KCI)
 
		# keep count on reaction number 
		gprn += 1 
		FC1=0.85
 
		# keep count on reaction number 
		gprn += 1 
		K10=1.0E-31*M*(TEMP/300.)**(-1.6)
 
		# keep count on reaction number 
		gprn += 1 
		K1I=5.0E-11*(TEMP/300.)**(-0.3)
 
		# keep count on reaction number 
		gprn += 1 
		KR1=K10/K1I
 
		# keep count on reaction number 
		gprn += 1 
		NC1=0.75-1.27*(numpy.log10(FC1))
 
		# keep count on reaction number 
		gprn += 1 
		F1=10.**(numpy.log10(FC1)/(1.+(numpy.log10(KR1)/NC1)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KMT01=(K10*K1I)*F1/(K10+K1I)
 
		# keep count on reaction number 
		gprn += 1 
		FC2=0.6
 
		# keep count on reaction number 
		gprn += 1 
		K20=1.3E-31*M*(TEMP/300.)**(-1.5)
 
		# keep count on reaction number 
		gprn += 1 
		K2I=2.3E-11*(TEMP/300.)**(0.24)
 
		# keep count on reaction number 
		gprn += 1 
		KR2=K20/K2I
 
		# keep count on reaction number 
		gprn += 1 
		NC2=0.75-1.27*(numpy.log10(FC2))
 
		# keep count on reaction number 
		gprn += 1 
		F2=10.**(numpy.log10(FC2)/(1.+(numpy.log10(KR2)/NC2)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KMT02=(K20*K2I)*F2/(K20+K2I)
 
		# keep count on reaction number 
		gprn += 1 
		FC3=0.35
 
		# keep count on reaction number 
		gprn += 1 
		K30=3.6E-30*M*(TEMP/300.)**(-4.1)
 
		# keep count on reaction number 
		gprn += 1 
		K3I=1.9E-12*(TEMP/300.)**(0.2)
 
		# keep count on reaction number 
		gprn += 1 
		KR3=K30/K3I
 
		# keep count on reaction number 
		gprn += 1 
		NC3=0.75-1.27*(numpy.log10(FC3))
 
		# keep count on reaction number 
		gprn += 1 
		F3=10.**(numpy.log10(FC3)/(1.+(numpy.log10(KR3)/NC3)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KMT03=(K30*K3I)*F3/(K30+K3I)
 
		# keep count on reaction number 
		gprn += 1 
		FC4=0.35
 
		# keep count on reaction number 
		gprn += 1 
		K40=1.3E-3*M*(TEMP/300.)**(-3.5)*numpy.exp(-11000./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		K4I=9.7E+14*(TEMP/300.)**(0.1)*numpy.exp(-11080./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		KR4=K40/K4I
 
		# keep count on reaction number 
		gprn += 1 
		NC4=0.75-1.27*(numpy.log10(FC4))
 
		# keep count on reaction number 
		gprn += 1 
		F4=10.**(numpy.log10(FC4)/(1.+(numpy.log10(KR4)/NC4)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KMT04=(K40*K4I)*F4/(K40+K4I)
 
		# keep count on reaction number 
		gprn += 1 
		FC7=0.81
 
		# keep count on reaction number 
		gprn += 1 
		K70=7.4E-31*M*(TEMP/300.)**(-2.4)
 
		# keep count on reaction number 
		gprn += 1 
		K7I=3.3E-11*(TEMP/300.)**(-0.3)
 
		# keep count on reaction number 
		gprn += 1 
		KR7=K70/K7I
 
		# keep count on reaction number 
		gprn += 1 
		NC7=0.75-1.27*(numpy.log10(FC7))
 
		# keep count on reaction number 
		gprn += 1 
		F7=10.**(numpy.log10(FC7)/(1.+(numpy.log10(KR7)/NC7)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KMT07=(K70*K7I)*F7/(K70+K7I)
 
		# keep count on reaction number 
		gprn += 1 
		FC8=0.41
 
		# keep count on reaction number 
		gprn += 1 
		K80=3.2E-30*M*(TEMP/300.)**(-4.5)
 
		# keep count on reaction number 
		gprn += 1 
		K8I=3.0E-11
 
		# keep count on reaction number 
		gprn += 1 
		KR8=K80/K8I
 
		# keep count on reaction number 
		gprn += 1 
		NC8=0.75-1.27*(numpy.log10(FC8))
 
		# keep count on reaction number 
		gprn += 1 
		F8=10.**(numpy.log10(FC8)/(1.+(numpy.log10(KR8)/NC8)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KMT08=(K80*K8I)*F8/(K80+K8I)
 
		# keep count on reaction number 
		gprn += 1 
		FC9=0.4
 
		# keep count on reaction number 
		gprn += 1 
		K90=1.4E-31*M*(TEMP/300.)**(-3.1)
 
		# keep count on reaction number 
		gprn += 1 
		K9I=4.0E-12
 
		# keep count on reaction number 
		gprn += 1 
		KR9=K90/K9I
 
		# keep count on reaction number 
		gprn += 1 
		NC9=0.75-1.27*(numpy.log10(FC9))
 
		# keep count on reaction number 
		gprn += 1 
		F9=10.**(numpy.log10(FC9)/(1.+(numpy.log10(KR9)/NC9)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KMT09=(K90*K9I)*F9/(K90+K9I)
 
		# keep count on reaction number 
		gprn += 1 
		FC10=0.4
 
		# keep count on reaction number 
		gprn += 1 
		K100=4.10E-05*M*numpy.exp(-10650./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		K10I=6.0E+15*numpy.exp(-11170./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		KR10=K100/K10I
 
		# keep count on reaction number 
		gprn += 1 
		NC10=0.75-1.27*(numpy.log10(FC10))
 
		# keep count on reaction number 
		gprn += 1 
		F10=10.**(numpy.log10(FC10)/(1.+(numpy.log10(KR10)/NC10)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KMT10=(K100*K10I)*F10/(K100+K10I)
 
		# keep count on reaction number 
		gprn += 1 
		K3=6.50E-34*numpy.exp(1335./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		K4=2.70E-17*numpy.exp(2199./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		K1=2.40E-14*numpy.exp(460./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		K2=(K3*M)/(1.+(K3*M/K4))
 
		# keep count on reaction number 
		gprn += 1 
		KMT11=K1+K2
 
		# keep count on reaction number 
		gprn += 1 
		FC12=0.53
 
		# keep count on reaction number 
		gprn += 1 
		K120=2.5E-31*M*(TEMP/300.)**(-2.6)
 
		# keep count on reaction number 
		gprn += 1 
		K12I=2.0E-12
 
		# keep count on reaction number 
		gprn += 1 
		KR12=K120/K12I
 
		# keep count on reaction number 
		gprn += 1 
		NC12=0.75-1.27*(numpy.log10(FC12))
 
		# keep count on reaction number 
		gprn += 1 
		F12=10.**(numpy.log10(FC12)/(1.0+(numpy.log10(KR12)/NC12)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KMT12=(K120*K12I*F12)/(K120+K12I)
 
		# keep count on reaction number 
		gprn += 1 
		FC13=0.36
 
		# keep count on reaction number 
		gprn += 1 
		K130=2.5E-30*M*(TEMP/300.)**(-5.5)
 
		# keep count on reaction number 
		gprn += 1 
		K13I=1.8E-11
 
		# keep count on reaction number 
		gprn += 1 
		KR13=K130/K13I
 
		# keep count on reaction number 
		gprn += 1 
		NC13=0.75-1.27*(numpy.log10(FC13))
 
		# keep count on reaction number 
		gprn += 1 
		F13=10.**(numpy.log10(FC13)/(1.+(numpy.log10(KR13)/NC13)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KMT13=(K130*K13I)*F13/(K130+K13I)
 
		# keep count on reaction number 
		gprn += 1 
		FC14=0.36
 
		# keep count on reaction number 
		gprn += 1 
		K140=9.0E-5*numpy.exp(-9690./TEMP)*M
 
		# keep count on reaction number 
		gprn += 1 
		K14I=1.1E+16*numpy.exp(-10560./TEMP)
 
		# keep count on reaction number 
		gprn += 1 
		KR14=K140/K14I
 
		# keep count on reaction number 
		gprn += 1 
		NC14=0.75-1.27*(numpy.log10(FC14))
 
		# keep count on reaction number 
		gprn += 1 
		F14=10.**(numpy.log10(FC14)/(1.+(numpy.log10(KR14)/NC14)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KMT14=(K140*K14I)*F14/(K140+K14I)
 
		# keep count on reaction number 
		gprn += 1 
		FC15=0.48
 
		# keep count on reaction number 
		gprn += 1 
		K150=8.6E-29*M*(TEMP/300.)**(-3.1)
 
		# keep count on reaction number 
		gprn += 1 
		K15I=9.0E-12*(TEMP/300.)**(-0.85)
 
		# keep count on reaction number 
		gprn += 1 
		KR15=K150/K15I
 
		# keep count on reaction number 
		gprn += 1 
		NC15=0.75-1.27*(numpy.log10(FC15))
 
		# keep count on reaction number 
		gprn += 1 
		F15=10.**(numpy.log10(FC15)/(1.+(numpy.log10(KR15)/NC15)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KMT15=(K150*K15I)*F15/(K150+K15I)
 
		# keep count on reaction number 
		gprn += 1 
		FC16=0.5
 
		# keep count on reaction number 
		gprn += 1 
		K160=8.E-27*M*(TEMP/300.)**(-3.5)
 
		# keep count on reaction number 
		gprn += 1 
		K16I=3.0E-11*(TEMP/300.)**(-1.)
 
		# keep count on reaction number 
		gprn += 1 
		KR16=K160/K16I
 
		# keep count on reaction number 
		gprn += 1 
		NC16=0.75-1.27*(numpy.log10(FC16))
 
		# keep count on reaction number 
		gprn += 1 
		F16=10.**(numpy.log10(FC16)/(1.+(numpy.log10(KR16)/NC16)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KMT16=(K160*K16I)*F16/(K160+K16I)
 
		# keep count on reaction number 
		gprn += 1 
		FC17=0.17*numpy.exp(-51./TEMP)+numpy.exp(-TEMP/204.)
 
		# keep count on reaction number 
		gprn += 1 
		K170=5.0E-30*M*(TEMP/300.)**(-1.5)
 
		# keep count on reaction number 
		gprn += 1 
		K17I=1.0E-12
 
		# keep count on reaction number 
		gprn += 1 
		KR17=K170/K17I
 
		# keep count on reaction number 
		gprn += 1 
		NC17=0.75-1.27*(numpy.log10(FC17))
 
		# keep count on reaction number 
		gprn += 1 
		F17=10.**(numpy.log10(FC17)/(1.0+(numpy.log10(KR17)/NC17)**(2.)))
 
		# keep count on reaction number 
		gprn += 1 
		KMT17=(K170*K17I*F17)/(K170+K17I)
 

	except:
		erf = 1 # flag error
		err_mess = str('Error: generic reaction rates failed to be calculated inside rate_coeffs.py at number ' + str(gprn) + ', please check chemical scheme and associated chemical scheme markers, which are stated in the model variables input file') # error message
		return([], erf, err_mess)
	# estimate and append photolysis rates
	J = photolysisRates.PhotolysisCalculation(TEMP, Jlen, sumt, self)

	if (self.light_stat_now == 0):
		J = [0]*len(J)
	rate_values = numpy.zeros((641))
	
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
		rc_eq_now = 'J[51]' 
		rate_values[57] = J[51]
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
		rc_eq_now = 'J[41]' 
		rate_values[60] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.55E-12*numpy.exp(380./TEMP)*0.991' 
		rate_values[61] = 2.55E-12*numpy.exp(380./TEMP)*0.991
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.55E-12*numpy.exp(380./TEMP)*0.009' 
		rate_values[62] = 2.55E-12*numpy.exp(380./TEMP)*0.009
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.4E-14*numpy.exp(-325./TEMP)*O2' 
		rate_values[63] = 2.4E-14*numpy.exp(-325./TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[64] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.3E-13*numpy.exp(870./TEMP)' 
		rate_values[65] = 4.3E-13*numpy.exp(870./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(KCH3O2*6.4E-14)**(0.5)*RO2*0.6' 
		rate_values[66] = 2.*(KCH3O2*6.4E-14)**(0.5)*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(KCH3O2*6.4E-14)**(0.5)*RO2*0.2' 
		rate_values[67] = 2.*(KCH3O2*6.4E-14)**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(KCH3O2*6.4E-14)**(0.5)*RO2*0.2' 
		rate_values[68] = 2.*(KCH3O2*6.4E-14)**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.7E-13*numpy.exp(-395./TEMP)' 
		rate_values[69] = 6.7E-13*numpy.exp(-395./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[52]' 
		rate_values[70] = J[52]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[71] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.01E-12' 
		rate_values[72] = 8.01E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[73] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.6E-12*numpy.exp(-585./TEMP)*0.264' 
		rate_values[74] = 7.6E-12*numpy.exp(-585./TEMP)*0.264
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.6E-12*numpy.exp(-585./TEMP)*0.736' 
		rate_values[75] = 7.6E-12*numpy.exp(-585./TEMP)*0.736
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.9E-12*numpy.exp(350./TEMP)*0.980' 
		rate_values[76] = 2.9E-12*numpy.exp(350./TEMP)*0.980
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.9E-12*numpy.exp(350./TEMP)*0.020' 
		rate_values[77] = 2.9E-12*numpy.exp(350./TEMP)*0.020
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6E-14*numpy.exp(-255./TEMP)*O2' 
		rate_values[78] = 2.6E-14*numpy.exp(-255./TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[79] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.520' 
		rate_values[80] = KRO2HO2*0.520
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(K298CH3O2*3.E-13)**(0.5)*RO2*0.6' 
		rate_values[81] = 2.*(K298CH3O2*3.E-13)**(0.5)*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(K298CH3O2*3.E-13)**(0.5)*RO2*0.2' 
		rate_values[82] = 2.*(K298CH3O2*3.E-13)**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(K298CH3O2*3.E-13)**(0.5)*RO2*0.2' 
		rate_values[83] = 2.*(K298CH3O2*3.E-13)**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.8E-13' 
		rate_values[84] = 5.8E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[53]' 
		rate_values[85] = J[53]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[86] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.10E-11' 
		rate_values[87] = 1.10E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[88] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.7E-12*numpy.exp(360./TEMP)*0.958' 
		rate_values[89] = 2.7E-12*numpy.exp(360./TEMP)*0.958
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.7E-12*numpy.exp(360./TEMP)*0.042' 
		rate_values[90] = 2.7E-12*numpy.exp(360./TEMP)*0.042
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.5E-14*numpy.exp(-230./TEMP)*O2' 
		rate_values[91] = 1.5E-14*numpy.exp(-230./TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[92] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.520' 
		rate_values[93] = KRO2HO2*0.520
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(KCH3O2*1.6E-12*numpy.exp(-2200./TEMP))**(0.5)*RO2*0.6' 
		rate_values[94] = 2.*(KCH3O2*1.6E-12*numpy.exp(-2200./TEMP))**(0.5)*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(KCH3O2*1.6E-12*numpy.exp(-2200./TEMP))**(0.5)*RO2*0.2' 
		rate_values[95] = 2.*(KCH3O2*1.6E-12*numpy.exp(-2200./TEMP))**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(KCH3O2*1.6E-12*numpy.exp(-2200./TEMP))**(0.5)*RO2*0.2' 
		rate_values[96] = 2.*(KCH3O2*1.6E-12*numpy.exp(-2200./TEMP))**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2E-13*numpy.exp(-230./TEMP)' 
		rate_values[97] = 6.2E-13*numpy.exp(-230./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[98] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[99] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.66E-11' 
		rate_values[100] = 1.66E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[101] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.981' 
		rate_values[102] = KRO2NO*0.981
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.019' 
		rate_values[103] = KRO2NO*0.019
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROPRIM*O2' 
		rate_values[104] = KROPRIM*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[105] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.520' 
		rate_values[106] = KRO2HO2*0.520
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.00E-13*0.6*RO2' 
		rate_values[107] = 6.00E-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.00E-13*0.2*RO2' 
		rate_values[108] = 6.00E-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.00E-13*0.2*RO2' 
		rate_values[109] = 6.00E-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.44E-12' 
		rate_values[110] = 4.44E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[53]' 
		rate_values[111] = J[53]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[112] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.52E-11' 
		rate_values[113] = 1.52E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[114] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.73E-12' 
		rate_values[115] = 9.73E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[116] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.06E-11' 
		rate_values[117] = 3.06E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*2.4' 
		rate_values[118] = KNO3AL*2.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[119] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[120] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[121] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[122] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[123] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[124] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.7*RO2' 
		rate_values[125] = 1.00E-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.3*RO2' 
		rate_values[126] = 1.00E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.51E-12' 
		rate_values[127] = 4.51E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.73E-11' 
		rate_values[128] = 1.73E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[129] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.39E-11' 
		rate_values[130] = 1.39E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.6E-12*numpy.exp(-1240./TEMP)' 
		rate_values[131] = 6.6E-12*numpy.exp(-1240./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4E-10*0.59*numpy.exp(-90./TEMP)' 
		rate_values[132] = 1.4E-10*0.59*numpy.exp(-90./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4E-10*0.43*numpy.exp(75./TEMP)' 
		rate_values[133] = 1.4E-10*0.43*numpy.exp(75./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.995' 
		rate_values[134] = KRO2NO*0.995
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.005' 
		rate_values[135] = KRO2NO*0.005
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.50E+13*numpy.exp(-5988./TEMP)' 
		rate_values[136] = 9.50E+13*numpy.exp(-5988./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KROPRIM*O2' 
		rate_values[137] = KROPRIM*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[138] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.53E-13*numpy.exp(1300./TEMP)' 
		rate_values[139] = 1.53E-13*numpy.exp(1300./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(KCH3O2*7.8E-14*numpy.exp(1000./TEMP))**(0.5)*RO2*0.6' 
		rate_values[140] = 2.*(KCH3O2*7.8E-14*numpy.exp(1000./TEMP))**(0.5)*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(KCH3O2*7.8E-14*numpy.exp(1000./TEMP))**(0.5)*RO2*0.2' 
		rate_values[141] = 2.*(KCH3O2*7.8E-14*numpy.exp(1000./TEMP))**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(KCH3O2*7.8E-14*numpy.exp(1000./TEMP))**(0.5)*RO2*0.2' 
		rate_values[142] = 2.*(KCH3O2*7.8E-14*numpy.exp(1000./TEMP))**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.40E-13' 
		rate_values[143] = 8.40E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[144] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[145] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.38E-11' 
		rate_values[146] = 1.38E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.800' 
		rate_values[147] = 1.00E-11*0.800
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.200' 
		rate_values[148] = 1.00E-11*0.200
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[149] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL' 
		rate_values[150] = KNO3AL
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[151] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[152] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[153] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[154] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[155] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[156] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.7*RO2' 
		rate_values[157] = 1.00E-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.3*RO2' 
		rate_values[158] = 1.00E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.73E-12' 
		rate_values[159] = 2.73E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[160] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.19E-12' 
		rate_values[161] = 6.19E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.977' 
		rate_values[162] = KRO2NO*0.977
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.023' 
		rate_values[163] = KRO2NO*0.023
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E+14*numpy.exp(-6410./TEMP)' 
		rate_values[164] = 2.00E+14*numpy.exp(-6410./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[165] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.520' 
		rate_values[166] = KRO2HO2*0.520
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.6*RO2' 
		rate_values[167] = 8.80E-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.2*RO2' 
		rate_values[168] = 8.80E-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.2*RO2' 
		rate_values[169] = 8.80E-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.16E-13' 
		rate_values[170] = 9.16E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[171] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[172] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.44E-11' 
		rate_values[173] = 2.44E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.6E-12*numpy.exp(305./TEMP)' 
		rate_values[174] = 1.6E-12*numpy.exp(305./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[175] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8E-13*numpy.exp(780./TEMP)*(1./(1.+498.*numpy.exp(-1160./TEMP)))' 
		rate_values[176] = 3.8E-13*numpy.exp(780./TEMP)*(1./(1.+498.*numpy.exp(-1160./TEMP)))
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.7E-12*numpy.exp(345./TEMP)*0.05' 
		rate_values[177] = 4.7E-12*numpy.exp(345./TEMP)*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.36E-13*numpy.exp(1250./TEMP)*0.15' 
		rate_values[178] = 1.36E-13*numpy.exp(1250./TEMP)*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.0E-12*numpy.exp(20./TEMP)*0.05' 
		rate_values[179] = 3.0E-12*numpy.exp(20./TEMP)*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.18' 
		rate_values[180] = KDEC*0.18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.57' 
		rate_values[181] = KDEC*0.57
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.125' 
		rate_values[182] = KDEC*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.125' 
		rate_values[183] = KDEC*0.125
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[31]' 
		rate_values[184] = J[31]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[32]' 
		rate_values[185] = J[32]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[33]' 
		rate_values[186] = J[33]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.1E-12*numpy.exp(340./TEMP)' 
		rate_values[187] = 3.1E-12*numpy.exp(340./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL' 
		rate_values[188] = KNO3AL
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[189] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[190] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.9E-12*numpy.exp(575./TEMP)' 
		rate_values[191] = 1.9E-12*numpy.exp(575./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*2.4' 
		rate_values[192] = KNO3AL*2.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[193] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[194] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[195] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[196] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.7*RO2' 
		rate_values[197] = 1.00E-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.3*RO2' 
		rate_values[198] = 1.00E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.58E-11' 
		rate_values[199] = 1.58E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]' 
		rate_values[200] = J[41]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.23E-11' 
		rate_values[201] = 1.23E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.0E-14' 
		rate_values[202] = 7.0E-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2E-15' 
		rate_values[203] = 1.2E-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0E-14' 
		rate_values[204] = 1.0E-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0E-15' 
		rate_values[205] = 1.0E-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.0E-18*H2O' 
		rate_values[206] = 6.0E-18*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0E-17*H2O' 
		rate_values[207] = 1.0E-17*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.4E-12*numpy.exp(135./TEMP)' 
		rate_values[208] = 5.4E-12*numpy.exp(135./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[11]' 
		rate_values[209] = J[11]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[12]' 
		rate_values[210] = J[12]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.5E-16' 
		rate_values[211] = 5.5E-16
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.7E-12*numpy.exp(345./TEMP)*0.95' 
		rate_values[212] = 4.7E-12*numpy.exp(345./TEMP)*0.95
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4E-12*numpy.exp(-1860./TEMP)' 
		rate_values[213] = 1.4E-12*numpy.exp(-1860./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.5E-12*numpy.exp(290./TEMP)' 
		rate_values[214] = 7.5E-12*numpy.exp(290./TEMP)
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
		rc_eq_now = '4.0E-12' 
		rate_values[217] = 4.0E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[218] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[219] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.7*RO2' 
		rate_values[220] = 1.00E-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.3*RO2' 
		rate_values[221] = 1.00E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.E-14' 
		rate_values[222] = 3.E-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[223] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.70E-12' 
		rate_values[224] = 3.70E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[13]' 
		rate_values[225] = J[13]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.9E-12*numpy.exp(405./TEMP)' 
		rate_values[226] = 4.9E-12*numpy.exp(405./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.24E-12*numpy.exp(-1860./TEMP)' 
		rate_values[227] = 3.24E-12*numpy.exp(-1860./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.7E-12*numpy.exp(340./TEMP)' 
		rate_values[228] = 6.7E-12*numpy.exp(340./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[229] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPPN' 
		rate_values[230] = KBPPN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[231] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[232] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[233] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.7*RO2' 
		rate_values[234] = 1.00E-11*0.7*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.3*RO2' 
		rate_values[235] = 1.00E-11*0.3*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.27E-12' 
		rate_values[236] = 1.27E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[237] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.42E-12' 
		rate_values[238] = 4.42E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[14]' 
		rate_values[239] = J[14]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[240] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[241] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[242] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.387' 
		rate_values[243] = KRO2HO2*0.387
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*0.6*RO2' 
		rate_values[244] = 2.00E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*0.2*RO2' 
		rate_values[245] = 2.00E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*0.2*RO2' 
		rate_values[246] = 2.00E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[247] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]' 
		rate_values[248] = J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[249] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.91E-11' 
		rate_values[250] = 2.91E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[21]' 
		rate_values[251] = J[21]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.8E-12*numpy.exp(-1320./TEMP)+1.7E-14*numpy.exp(423./TEMP)' 
		rate_values[252] = 8.8E-12*numpy.exp(-1320./TEMP)+1.7E-14*numpy.exp(423./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[253] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[254] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[255] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.36E-13*numpy.exp(1250./TEMP)*0.85' 
		rate_values[256] = 1.36E-13*numpy.exp(1250./TEMP)*0.85
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(K298CH3O2*8.0E-12)**(0.5)*RO2*0.6' 
		rate_values[257] = 2.*(K298CH3O2*8.0E-12)**(0.5)*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(K298CH3O2*8.0E-12)**(0.5)*RO2*0.2' 
		rate_values[258] = 2.*(K298CH3O2*8.0E-12)**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*(K298CH3O2*8.0E-12)**(0.5)*RO2*0.2' 
		rate_values[259] = 2.*(K298CH3O2*8.0E-12)**(0.5)*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[260] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.39E-12' 
		rate_values[261] = 8.39E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[262] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[263] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.85E-12*numpy.exp(-345./TEMP)' 
		rate_values[264] = 2.85E-12*numpy.exp(-345./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.0E-12*numpy.exp(20./TEMP)*0.9' 
		rate_values[265] = 3.0E-12*numpy.exp(20./TEMP)*0.9
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.0E-12*numpy.exp(20./TEMP)*0.05' 
		rate_values[266] = 3.0E-12*numpy.exp(20./TEMP)*0.05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.6E-12*numpy.exp(70./TEMP)*0.494' 
		rate_values[267] = 4.6E-12*numpy.exp(70./TEMP)*0.494
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.6E-12*numpy.exp(70./TEMP)*0.443' 
		rate_values[268] = 4.6E-12*numpy.exp(70./TEMP)*0.443
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6E-12*numpy.exp(200./TEMP)*0.861' 
		rate_values[269] = 2.6E-12*numpy.exp(200./TEMP)*0.861
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.6E-12*numpy.exp(200./TEMP)*0.139' 
		rate_values[270] = 2.6E-12*numpy.exp(200./TEMP)*0.139
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.991' 
		rate_values[271] = KRO2NO*0.991
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.009' 
		rate_values[272] = KRO2NO*0.009
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E+14*numpy.exp(-5505./TEMP)' 
		rate_values[273] = 2.00E+14*numpy.exp(-5505./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[274] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.520' 
		rate_values[275] = KRO2HO2*0.520
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*0.6*RO2' 
		rate_values[276] = 2.00E-12*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*0.2*RO2' 
		rate_values[277] = 2.00E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-12*0.2*RO2' 
		rate_values[278] = 2.00E-12*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.71E-12' 
		rate_values[279] = 1.71E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[280] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.83E-11' 
		rate_values[281] = 1.83E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[282] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.7E-11' 
		rate_values[283] = 1.7E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[284] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*2.4' 
		rate_values[285] = KNO3AL*2.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[286] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[287] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[288] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[289] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[290] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2' 
		rate_values[291] = 1.00E-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.34E-12' 
		rate_values[292] = 2.34E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.34E-12' 
		rate_values[293] = 9.34E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[294] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.45E-11' 
		rate_values[295] = 1.45E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.10E-11*0.387' 
		rate_values[296] = 2.10E-11*0.387
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.10E-11*0.613' 
		rate_values[297] = 2.10E-11*0.613
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00E-13' 
		rate_values[298] = 8.00E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.2E-12' 
		rate_values[299] = 1.2E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.6E-12*numpy.exp(70./TEMP)*0.063' 
		rate_values[300] = 4.6E-12*numpy.exp(70./TEMP)*0.063
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3E-12*numpy.exp(-190./TEMP)*0.352' 
		rate_values[301] = 2.3E-12*numpy.exp(-190./TEMP)*0.352
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3E-12*numpy.exp(-190./TEMP)*0.118' 
		rate_values[302] = 2.3E-12*numpy.exp(-190./TEMP)*0.118
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.3E-12*numpy.exp(-190./TEMP)*0.53' 
		rate_values[303] = 2.3E-12*numpy.exp(-190./TEMP)*0.53
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.918' 
		rate_values[304] = KRO2NO*0.918
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.082' 
		rate_values[305] = KRO2NO*0.082
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.50' 
		rate_values[306] = KDEC*0.50
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.50' 
		rate_values[307] = KDEC*0.50
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[308] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[309] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.6' 
		rate_values[310] = 8.80E-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.2' 
		rate_values[311] = 8.80E-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.2' 
		rate_values[312] = 8.80E-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[313] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.30E-11' 
		rate_values[314] = 7.30E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[315] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.77E-11' 
		rate_values[316] = 9.77E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.21E-10' 
		rate_values[317] = 1.21E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[22]' 
		rate_values[318] = J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.16E-11' 
		rate_values[319] = 8.16E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.20E-11*0.17' 
		rate_values[320] = 5.20E-11*0.17
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.20E-11*0.83' 
		rate_values[321] = 5.20E-11*0.83
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.14*0.4' 
		rate_values[322] = J[4]*0.14*0.4
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.14*0.6' 
		rate_values[323] = J[4]*0.14*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[324] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[325] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[326] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.60*RO2' 
		rate_values[327] = 8.80E-13*0.60*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.20*RO2' 
		rate_values[328] = 8.80E-13*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.20*RO2' 
		rate_values[329] = 8.80E-13*0.20*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[330] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[331] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[332] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[333] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[334] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[335] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[336] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.70*RO2' 
		rate_values[337] = 1.00E-11*0.70*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*0.30*RO2' 
		rate_values[338] = 1.00E-11*0.30*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[339] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2.' 
		rate_values[340] = J[15]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[341] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.22E-10' 
		rate_values[342] = 1.22E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]+J[15]' 
		rate_values[343] = J[34]+J[15]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.67E-11' 
		rate_values[344] = 3.67E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]*2.' 
		rate_values[345] = J[34]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.45E-11' 
		rate_values[346] = 2.45E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[347] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[348] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[349] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[350] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[351] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2' 
		rate_values[352] = 1.00E-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.97E-11' 
		rate_values[353] = 6.97E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[354] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.33E-11' 
		rate_values[355] = 7.33E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2.' 
		rate_values[356] = J[15]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.13E-11' 
		rate_values[357] = 8.13E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2.' 
		rate_values[358] = J[15]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.36E-10' 
		rate_values[359] = 1.36E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-18' 
		rate_values[360] = 2.00E-18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.4E-12' 
		rate_values[361] = 1.4E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[362] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[363] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[364] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[365] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.6*RO2' 
		rate_values[366] = 8.80E-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.2*RO2' 
		rate_values[367] = 8.80E-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.2*RO2' 
		rate_values[368] = 8.80E-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[369] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.66E-11' 
		rate_values[370] = 4.66E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.55E-11' 
		rate_values[371] = 2.55E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.68E-12' 
		rate_values[372] = 5.68E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.00E-18' 
		rate_values[373] = 2.00E-18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*KNO3AL*2.0' 
		rate_values[374] = 2.*KNO3AL*2.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.1*0.5' 
		rate_values[375] = J[4]*0.1*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[4]*0.1*0.5' 
		rate_values[376] = J[4]*0.1*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*KNO3AL*2.75' 
		rate_values[377] = 2.*KNO3AL*2.75
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.08E-11*0.69' 
		rate_values[378] = 6.08E-11*0.69
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.08E-11*0.31' 
		rate_values[379] = 6.08E-11*0.31
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.895' 
		rate_values[380] = KRO2NO*0.895
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO*0.105' 
		rate_values[381] = KRO2NO*0.105
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[382] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[383] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[384] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.77' 
		rate_values[385] = KRO2HO2*0.77
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.6*RO2' 
		rate_values[386] = 8.80E-13*0.6*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.2*RO2' 
		rate_values[387] = 8.80E-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*0.2*RO2' 
		rate_values[388] = 8.80E-13*0.2*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.38E-11' 
		rate_values[389] = 4.38E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[15]*2.' 
		rate_values[390] = J[41]+J[15]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.31E-10' 
		rate_values[391] = 1.31E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2.' 
		rate_values[392] = J[15]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.23E-11' 
		rate_values[393] = 8.23E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2.+J[22]' 
		rate_values[394] = J[15]*2.+J[22]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[395] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[396] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[397] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.520' 
		rate_values[398] = KRO2HO2*0.520
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.6' 
		rate_values[399] = 8.80E-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.2' 
		rate_values[400] = 8.80E-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.2' 
		rate_values[401] = 8.80E-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[402] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2.' 
		rate_values[403] = J[15]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.44E-10' 
		rate_values[404] = 1.44E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[15]*2.' 
		rate_values[405] = J[15]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.77E-11' 
		rate_values[406] = 5.77E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[407] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[408] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[409] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[410] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[411] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[412] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[413] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[414] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2*0.7' 
		rate_values[415] = 1.00E-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2*0.3' 
		rate_values[416] = 1.00E-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[417] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[418] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.6' 
		rate_values[419] = 8.80E-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.2' 
		rate_values[420] = 8.80E-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.2' 
		rate_values[421] = 8.80E-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.05E-11' 
		rate_values[422] = 4.05E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[18]+J[19]' 
		rate_values[423] = J[41]+J[18]+J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.37E-11' 
		rate_values[424] = 4.37E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]+J[19]' 
		rate_values[425] = J[18]+J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.06E-11' 
		rate_values[426] = 4.06E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[427] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]+J[19]' 
		rate_values[428] = J[18]+J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.52E-11' 
		rate_values[429] = 7.52E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]+J[19]' 
		rate_values[430] = J[18]+J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.75E-11' 
		rate_values[431] = 7.75E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]+J[18]+J[19]' 
		rate_values[432] = J[34]+J[18]+J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.90E-11' 
		rate_values[433] = 4.90E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]*2.' 
		rate_values[434] = J[17]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.32E-11' 
		rate_values[435] = 4.32E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.*KNO3AL*4.0' 
		rate_values[436] = 2.*KNO3AL*4.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[437] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[438] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[439] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[440] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.41' 
		rate_values[441] = KAPHO2*0.41
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.15' 
		rate_values[442] = KAPHO2*0.15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2*0.7' 
		rate_values[443] = 1.00E-11*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2*0.3' 
		rate_values[444] = 1.00E-11*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.29E-11' 
		rate_values[445] = 2.29E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[17]' 
		rate_values[446] = J[41]+J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.62E-11' 
		rate_values[447] = 2.62E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[448] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.31E-11' 
		rate_values[449] = 2.31E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.6E-12' 
		rate_values[450] = 4.6E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[451] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[452] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[453] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[454] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.6' 
		rate_values[455] = 8.80E-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.2' 
		rate_values[456] = 8.80E-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.2' 
		rate_values[457] = 8.80E-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[458] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.23E-10' 
		rate_values[459] = 1.23E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.18E-11' 
		rate_values[460] = 9.18E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.07E-11' 
		rate_values[461] = 6.07E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.45E-11' 
		rate_values[462] = 4.45E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[463] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[464] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[465] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.706' 
		rate_values[466] = KRO2HO2*0.706
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.6' 
		rate_values[467] = 8.80E-13*RO2*0.6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.2' 
		rate_values[468] = 8.80E-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2*0.2' 
		rate_values[469] = 8.80E-13*RO2*0.2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[470] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.68E-11' 
		rate_values[471] = 3.68E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.78E-11' 
		rate_values[472] = 2.78E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.78E-11' 
		rate_values[473] = 1.78E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[474] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.44E-11' 
		rate_values[475] = 3.44E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KNO3AL*8.0' 
		rate_values[476] = KNO3AL*8.0
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.7E-13*numpy.exp(1220./TEMP)*0.06' 
		rate_values[477] = 4.7E-13*numpy.exp(1220./TEMP)*0.06
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.7E-13*numpy.exp(1220./TEMP)*0.14' 
		rate_values[478] = 4.7E-13*numpy.exp(1220./TEMP)*0.14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.7E-13*numpy.exp(1220./TEMP)*0.8' 
		rate_values[479] = 4.7E-13*numpy.exp(1220./TEMP)*0.8
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8E-12*0.742' 
		rate_values[480] = 3.8E-12*0.742
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.8E-12*0.258' 
		rate_values[481] = 3.8E-12*0.258
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.29' 
		rate_values[482] = KDEC*0.29
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86E-13' 
		rate_values[483] = 2.86E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[484] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.00E-14' 
		rate_values[485] = 9.00E-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08E-12' 
		rate_values[486] = 2.08E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86E-13' 
		rate_values[487] = 2.86E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[488] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[489] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[490] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50E-13*RO2' 
		rate_values[491] = 2.50E-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[492] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.00E-13' 
		rate_values[493] = 9.00E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.0E-10' 
		rate_values[494] = 1.0E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.9E-11' 
		rate_values[495] = 9.9E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.2E-18' 
		rate_values[496] = 9.2E-18
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[497] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08E-12' 
		rate_values[498] = 2.08E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.86E-13' 
		rate_values[499] = 2.86E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[500] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[501] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[502] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2' 
		rate_values[503] = 8.80E-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[504] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[505] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.47E-12' 
		rate_values[506] = 3.47E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.60E-12' 
		rate_values[507] = 2.60E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[508] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[509] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[510] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[511] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00E-13*RO2' 
		rate_values[512] = 8.00E-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[513] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[514] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[515] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[516] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[517] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[518] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00E-13*RO2' 
		rate_values[519] = 8.00E-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[520] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[521] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[522] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[523] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.00E-14' 
		rate_values[524] = 3.00E-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.25E-15' 
		rate_values[525] = 2.25E-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[526] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[527] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[528] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[529] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00E-13*RO2' 
		rate_values[530] = 8.00E-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[531] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[532] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[533] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[534] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[535] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[536] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00E-13*RO2' 
		rate_values[537] = 8.00E-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[538] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.90E-12*numpy.exp(190./TEMP)' 
		rate_values[539] = 1.90E-12*numpy.exp(190./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[540] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.71' 
		rate_values[541] = KDEC*0.71
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[542] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[543] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00E-13*RO2*0.7' 
		rate_values[544] = 8.00E-13*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00E-13*RO2*0.3' 
		rate_values[545] = 8.00E-13*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[546] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.16E-10' 
		rate_values[547] = 1.16E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.13E-10' 
		rate_values[548] = 1.13E-10
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
		rc_eq_now = 'KRO2NO3' 
		rate_values[551] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[552] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00E-13*RO2*0.7' 
		rate_values[553] = 8.00E-13*RO2*0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.00E-13*RO2*0.3' 
		rate_values[554] = 8.00E-13*RO2*0.3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[555] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.07E-10' 
		rate_values[556] = 1.07E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.04E-10' 
		rate_values[557] = 1.04E-10
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.08E-12' 
		rate_values[558] = 2.08E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]' 
		rate_values[559] = J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[18]+J[19]' 
		rate_values[560] = J[18]+J[19]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.70E-11' 
		rate_values[561] = 3.70E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[20]*2.' 
		rate_values[562] = J[20]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.00E-11' 
		rate_values[563] = 4.00E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.12E-12' 
		rate_values[564] = 1.12E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[565] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[566] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[567] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[568] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[569] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2' 
		rate_values[570] = 1.00E-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.66E-11' 
		rate_values[571] = 7.66E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.01E-11' 
		rate_values[572] = 8.01E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[573] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.20E-11' 
		rate_values[574] = 9.20E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.40' 
		rate_values[575] = KDEC*0.40
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.60' 
		rate_values[576] = KDEC*0.60
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.00E-13' 
		rate_values[577] = 3.00E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[578] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[579] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[580] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.16E-12' 
		rate_values[581] = 1.16E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[582] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.625' 
		rate_values[583] = KRO2HO2*0.625
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2' 
		rate_values[584] = 8.80E-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[585] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.18E-12' 
		rate_values[586] = 6.18E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.20E-19' 
		rate_values[587] = 2.20E-19
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[588] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC*0.5' 
		rate_values[589] = KDEC*0.5
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00E-14' 
		rate_values[590] = 7.00E-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.20E-15' 
		rate_values[591] = 1.20E-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-14' 
		rate_values[592] = 1.00E-14
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-15' 
		rate_values[593] = 1.00E-15
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.00E-18*H2O' 
		rate_values[594] = 6.00E-18*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-17*H2O' 
		rate_values[595] = 1.00E-17*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.19E-11' 
		rate_values[596] = 2.19E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.00E-13' 
		rate_values[597] = 3.00E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[598] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KDEC' 
		rate_values[599] = KDEC
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[600] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[601] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '8.80E-13*RO2' 
		rate_values[602] = 8.80E-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]' 
		rate_values[603] = J[41]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.68E-11' 
		rate_values[604] = 6.68E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[34]*2.' 
		rate_values[605] = J[34]*2.
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.70E-11' 
		rate_values[606] = 7.70E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPNO' 
		rate_values[607] = KAPNO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KFPAN' 
		rate_values[608] = KFPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KBPAN' 
		rate_values[609] = KBPAN
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3*1.74' 
		rate_values[610] = KRO2NO3*1.74
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.56' 
		rate_values[611] = KAPHO2*0.56
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.00E-11*RO2' 
		rate_values[612] = 1.00E-11*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[41]+J[34]' 
		rate_values[613] = J[41]+J[34]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.06E-11' 
		rate_values[614] = 3.06E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.70E-11' 
		rate_values[615] = 3.70E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.74E-11' 
		rate_values[616] = 2.74E-11
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[17]' 
		rate_values[617] = J[17]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[618] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'J[54]' 
		rate_values[619] = J[54]
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO' 
		rate_values[620] = KRO2NO
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2NO3' 
		rate_values[621] = KRO2NO3
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '9.00E-13' 
		rate_values[622] = 9.00E-13
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KRO2HO2*0.770' 
		rate_values[623] = KRO2HO2*0.770
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.50E-13*RO2' 
		rate_values[624] = 2.50E-13*RO2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.60E-12' 
		rate_values[625] = 3.60E-12
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[626] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[627] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[628] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[629] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[630] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[631] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[632] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[633] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[634] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[635] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[636] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KAPHO2*0.44' 
		rate_values[637] = KAPHO2*0.44
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '7.00E11*numpy.exp(-3160./TEMP)+5.00E-12*O2' 
		rate_values[638] = 7.00E11*numpy.exp(-3160./TEMP)+5.00E-12*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00E-12*O2*3.2*numpy.exp(-550./TEMP)' 
		rate_values[639] = 5.00E-12*O2*3.2*numpy.exp(-550./TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.00E-12*O2*3.2*(1.-numpy.exp(-550./TEMP))' 
		rate_values[640] = 5.00E-12*O2*3.2*(1.-numpy.exp(-550./TEMP))
	except:
		erf = 1 # flag error
		err_mess = (str('Error: Could not calculate '+ 
		'rate coefficient for equation number ' 
		+ str(gprn) + ' ' + rc_eq_now + 
		' (message from rate coeffs.py)'))
	
	# aqueous-phase reactions
	
	# surface (e.g. wall) reactions
	
	return(rate_values, erf, err_mess)
