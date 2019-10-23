# -*- coding: utf-8 -*-
# function to calculate reaction rate coefficients

import numpy
import PhotolysisRates
import ipdb

def mcm_constants(TEMP, H2O, M, N2, O2, time, lat, lon):
    
	# write a name list of MCM constants; Automatically generated
	MCMConstNameList = [
		'KRO2NO','KRO2HO2','KAPHO2','KAPNO','KRO2NO3','KNO3AL','KDEC','KROPRIM','KROSEC',
		'KCH3O2','K298CH3O2','KD0','KDI','KRD','FCD','NCD','FD','KBPAN','KC0','KCI',
		'KRC','FCC','NC','FC','KFPAN','K10','K1I','KR1','FC1','NC1',
		'F1','KMT01','K20','K2I','KR2','FC2','NC2','F2','KMT02','K30',
		'K3I','KR3','FC3','NC3','F3','KMT03','K40','K4I','KR4','FC4',
		'NC4','F4','KMT04','KMT05','KMT06','K70','K7I','KR7','FC7','NC7',
		'F7','KMT07','K80','K8I','KR8','FC8','NC8','F8','KMT08','K90',
		'K9I','KR9','FC9','NC9','F9','KMT09','K100','K10I','KR10','FC10',
		'NC10','F10','KMT10','K1','K3','K4','K2','KMT11','K120','K12I',
		'KR12','FC12','NC12','F12','KMT12','K130','K13I','KR13','FC13','NC13',
		'F13','KMT13','K140','K14I','KR14','FC14','NC14','F14','KMT14','K150',
		'K15I','KR15','FC15','NC15','F15','KMT15','K160','K16I','KR16','FC16',
		'NC16','F16','KMT16','K170','K17I','KR17','FC17','NC17','F17','KMT17',
		'KMT18','KPPN0','KPPNI','KRPPN','FCPPN','NCPPN','FPPN','KBPPN',
	]
    
	# Rate coefficients from the MCM_KPP file
	KRO2NO = 2.7e-12*numpy.exp(360/TEMP)
	KRO2HO2 = 2.91e-13*numpy.exp(1300/TEMP)
	KAPHO2 = 5.2e-13*numpy.exp(980/TEMP)
	KAPNO = 7.5e-12*numpy.exp(290/TEMP)
	KRO2NO3 = 2.3e-12
	KNO3AL = 1.4e-12*numpy.exp(-1860/TEMP)
	KDEC = 1.00e+06
	KROPRIM = 2.50e-14*numpy.exp(-300/TEMP)
	KROSEC = 2.50e-14*numpy.exp(-300/TEMP)
	KCH3O2 = 1.03e-13*numpy.exp(365/TEMP)
	K298CH3O2 = 3.5e-13
	KD0 = 4.90e-3*numpy.exp(-12100/TEMP)*M
	KDI = 5.4e+16*numpy.exp(-13830/TEMP)
	KRD = KD0/KDI
	FCD = 0.30
	NCD = 0.75-1.27*(numpy.log10(FCD))
	FD = 10**(numpy.log10(FCD)/(1+(numpy.log10(KRD)/NCD)**2))
	KBPAN = (KD0*KDI)*FD/(KD0+KDI)
	KC0 = 2.7e-28*M*(TEMP/300)**-7.1
	KCI = 1.2e-11*(TEMP/300)**-0.9
	KRC = KC0/KCI
	FCC = 0.30
	NC = 0.75-1.27*(numpy.log10(FCC))
	FC = 10**(numpy.log10(FCC)/(1+(numpy.log10(KRC)/NC)**2))
	KFPAN = (KC0*KCI)*FC/(KC0+KCI)
	K10 = 1.0e-31*M*(TEMP/300)**-1.6
	K1I = 3.00e-11*(TEMP/300)**0.3
	KR1 = K10/K1I
	FC1 = 0.85
	NC1 = 0.75-1.27*(numpy.log10(FC1))
	F1 = 10**(numpy.log10(FC1)/(1+(numpy.log10(KR1)/NC1)**2))
	KMT01 = (K10*K1I)*F1/(K10+K1I)
	K20 = 1.3e-31*M*(TEMP/300)**-1.5
	K2I = 2.3e-11*(TEMP/300)**0.24
	KR2 = K20/K2I
	FC2 = 0.6
	NC2 = 0.75-1.27*(numpy.log10(FC2))
	F2 = 10**(numpy.log10(FC2)/(1+(numpy.log10(KR2)/NC2)**2))
	KMT02 = (K20*K2I)*F2/(K20+K2I)
	K30 = 3.6e-30*M*(TEMP/300)**-4.1
	K3I = 1.9e-12*(TEMP/300)**0.2
	KR3 = K30/K3I
	FC3 = 0.35
	NC3 = 0.75-1.27*(numpy.log10(FC3))
	F3 = 10**(numpy.log10(FC3)/(1+(numpy.log10(KR3)/NC3)**2))
	KMT03 = (K30*K3I)*F3/(K30+K3I)
	K40 = 1.3e-3*M*(TEMP/300)**-3.5*numpy.exp(-11000/TEMP)
	K4I = 9.7e+14*(TEMP/300)**0.1*numpy.exp(-11080/TEMP)
	KR4 = K40/K4I
	FC4 = 0.35
	NC4 = 0.75-1.27*(numpy.log10(FC4))
	F4 = 10**(numpy.log10(FC4)/(1+(numpy.log10(KR4)/NC4)**2))
	KMT04 = (K40*K4I)*F4/(K40+K4I)
	KMT05 = 1.44e-13*(1+(M/4.2e+19))
	KMT06 = 1 + (1.40e-21*numpy.exp(2200/TEMP)*H2O)
	K70 = 7.4e-31*M*(TEMP/300)**-2.4
	K7I = 3.3e-11*(TEMP/300)**-0.3
	KR7 = K70/K7I
	FC7 = numpy.exp(-TEMP/1420)
	NC7 = 0.75-1.27*(numpy.log10(FC7))
	F7 = 10**(numpy.log10(FC7)/(1+(numpy.log10(KR7)/NC7)**2))
	KMT07 = (K70*K7I)*F7/(K70+K7I)
	K80 = 3.3e-30*M*(TEMP/300)**-3.0
	K8I = 4.1e-11
	KR8 = K80/K8I
	FC8 = 0.4
	NC8 = 0.75-1.27*(numpy.log10(FC8))
	F8 = 10**(numpy.log10(FC8)/(1+(numpy.log10(KR8)/NC8)**2))
	KMT08 = (K80*K8I)*F8/(K80+K8I)
	K90 = 1.8e-31*M*(TEMP/300)**-3.2
	K9I = 4.7e-12
	KR9 = K90/K9I
	FC9 = 0.6
	NC9 = 0.75-1.27*(numpy.log10(FC9))
	F9 = 10**(numpy.log10(FC9)/(1+(numpy.log10(KR9)/NC9)**2))
	KMT09 = (K90*K9I)*F9/(K90+K9I)
	K100 = 4.10e-05*M*numpy.exp(-10650/TEMP)
	K10I = 4.8e+15*numpy.exp(-11170/TEMP)
	KR10 = K100/K10I
	FC10 = 0.6
	NC10 = 0.75-1.27*(numpy.log10(FC10))
	F10 = 10**(numpy.log10(FC10)/(1+(numpy.log10(KR10)/NC10)**2))
	KMT10 = (K100*K10I)*F10/(K100+K10I)
	K1 = 2.40e-14*numpy.exp(460/TEMP)
	K3 = 6.50e-34*numpy.exp(1335/TEMP)
	K4 = 2.70e-17*numpy.exp(2199/TEMP)
	K2 = (K3*M)/(1+(K3*M/K4))
	KMT11 = K1 + K2
	K120 = 4.5e-31*M*(TEMP/300)**-3.9
	K12I = 1.3e-12*(TEMP/300)**-0.7
	KR12 = K120/K12I
	FC12 = 0.525
	NC12 = 0.75-1.27*(numpy.log10(FC12))
	F12 = 10**(numpy.log10(FC12)/(1.0+(numpy.log10(KR12)/NC12)**2))
	KMT12 = (K120*K12I*F12)/(K120+K12I)
	K130 = 2.5e-30*M*(TEMP/300)**-5.5
	K13I = 1.8e-11
	KR13 = K130/K13I
	FC13 = 0.36
	NC13 = 0.75-1.27*(numpy.log10(FC13))
	F13 = 10**(numpy.log10(FC13)/(1+(numpy.log10(KR13)/NC13)**2))
	KMT13 = (K130*K13I)*F13/(K130+K13I)
	K140 = 9.0e-5*numpy.exp(-9690/TEMP)*M
	K14I = 1.1e+16*numpy.exp(-10560/TEMP)
	KR14 = K140/K14I
	FC14 = 0.4
	NC14 = 0.75-1.27*(numpy.log10(FC14))
	F14 = 10**(numpy.log10(FC14)/(1+(numpy.log10(KR14)/NC14)**2))
	KMT14 = (K140*K14I)*F14/(K140+K14I)
	K150 = 8.6e-29*M*(TEMP/300)**-3.1
	K15I = 9.0e-12*(TEMP/300)**-0.85
	KR15 = K150/K15I
	FC15 = 0.48
	NC15 = 0.75-1.27*(numpy.log10(FC15))
	F15 = 10**(numpy.log10(FC15)/(1+(numpy.log10(KR15)/NC15)**2))
	KMT15 = (K150*K15I)*F15/(K150+K15I)
	K160 = 8e-27*M*(TEMP/300)**-3.5
	K16I = 3.0e-11*(TEMP/300)**-1
	KR16 = K160/K16I
	FC16 = 0.5
	NC16 = 0.75-1.27*(numpy.log10(FC16))
	F16 = 10**(numpy.log10(FC16)/(1+(numpy.log10(KR16)/NC16)**2))
	KMT16 = (K160*K16I)*F16/(K160+K16I)
	K170 = 5.0e-30*M*(TEMP/300)**-1.5
	K17I = 1.0e-12
	KR17 = K170/K17I
	FC17 = 0.17*numpy.exp(-51/TEMP)+numpy.exp(-TEMP/204)
	NC17 = 0.75-1.27*(numpy.log10(FC17))
	F17 = 10**(numpy.log10(FC17)/(1.0+(numpy.log10(KR17)/NC17)**2))
	KMT17 = (K170*K17I*F17)/(K170+K17I)
	KMT18 = 9.5e-39*O2*numpy.exp(5270/TEMP)/(1+7.5e-29*O2*numpy.exp(5610/TEMP))
	KPPN0 = 1.7e-03*numpy.exp(-11280/TEMP)*M
	KPPNI = 8.3e+16*numpy.exp(-13940/TEMP)
	KRPPN = KPPN0/KPPNI
	FCPPN = 0.36
	NCPPN = 0.75-1.27*(numpy.log10(FCPPN))
	FPPN = 10**(numpy.log10(FCPPN)/(1+(numpy.log10(KRPPN)/NCPPN)**2))
	KBPPN = (KPPN0*KPPNI)*FCPPN/(KPPN0+KPPNI)
	
	# empty list of required constants
	mcm_constants = []
	
	for i in MCMConstNameList:
		mcm_constants.append(locals()[i])
       
	# estimate and append photolysis rates
	j = PhotolysisRates.PhotolysisCalculation(time, lat, lon)

	mcm_constants.append(j)
	
	# append photolysis rate label to MCMConstNameList
	MCMConstNameList.append('J')
    
    
	return mcm_constants, MCMConstNameList