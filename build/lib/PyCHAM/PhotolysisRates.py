'''module for calculating photolysis rates for MCM (Master Chemical Mechanism) equations'''

# can use ambient sunlight expression or artificial chamber lights for calculating a
# photolysis rate that is fixed when the lights are on

import scipy
import numpy
import os
import numpy as np

def PhotolysisCalculation(time, lat, lon):
	
	(theta, secx, cosx) = zenith(time, lat, lon)

	J = [0.0 for i in range(62)] # modified, so that we dont have to check None type
    
    #J          L           M          N
#	J[1]=6.073E-05*cosx**(1.743)*numpy.exp(-1.0*0.474*secx)
#	J[2]=4.775E-04*cosx**(0.298)*numpy.exp(-1.0*0.080*secx)
#	J[3]=1.041E-05*cosx**(0.723)*numpy.exp(-1.0*0.279*secx)
#	J[4]=1.165E-02*cosx**(0.244)*numpy.exp(-1.0*0.267*secx)
#	J[5]=2.485E-02*cosx**(0.168)*numpy.exp(-1.0*0.108*secx)
#	J[6]=1.747E-01*cosx**(0.155)*numpy.exp(-1.0*0.125*secx)
#	J[7]=2.644E-03*cosx**(0.261)*numpy.exp(-1.0*0.288*secx)
#	J[8]=9.312E-07*cosx**(1.230)*numpy.exp(-1.0*0.307*secx)
#	J[11]=4.642E-05*cosx**(0.762)*numpy.exp(-1.0*0.353*secx)
#	J[12]=6.853E-05*cosx**(0.477)*numpy.exp(-1.0*0.323*secx)
#	J[13]=7.344E-06*cosx**(1.202)*numpy.exp(-1.0*0.417*secx)
#	J[14]=2.879E-05*cosx**(1.067)*numpy.exp(-1.0*0.358*secx)
#	J[15]=2.792E-05*cosx**(0.805)*numpy.exp(-1.0*0.338*secx)
#	J[16]=1.675E-05*cosx**(0.805)*numpy.exp(-1.0*0.338*secx)
#	J[17]=7.914E-05*cosx**(0.764)*numpy.exp(-1.0*0.364*secx)
#	J[18]=1.140E-05*cosx**(0.396)*numpy.exp(-1.0*0.298*secx)
#	J[19]=1.140E-05*cosx**(0.396)*numpy.exp(-1.0*0.298*secx)
#	J[21]=7.992E-07*cosx**(1.578)*numpy.exp(-1.0*0.271*secx)
#	J[22]=5.804E-06*cosx**(1.092)*numpy.exp(-1.0*0.377*secx)
#	J[23]=1.836E-05*cosx**(0.395)*numpy.exp(-1.0*0.296*secx)
#	J[24]=1.836E-05*cosx**(0.395)*numpy.exp(-1.0*0.296*secx)
#	J[31]=6.845E-05*cosx**(0.130)*numpy.exp(-1.0*0.201*secx)
#	J[32]=1.032E-05*cosx**(0.130)*numpy.exp(-1.0*0.201*secx)
#	J[33]=3.802E-05*cosx**(0.644)*numpy.exp(-1.0*0.312*secx)
#	J[34]=1.537E-04*cosx**(0.170)*numpy.exp(-1.0*0.208*secx)
#	J[35]=3.326E-04*cosx**(0.148)*numpy.exp(-1.0*0.215*secx)
#	J[41]=7.649E-06*cosx**(0.682)*numpy.exp(-1.0*0.279*secx)
#	J[51]=1.588E-06*cosx**(1.154)*numpy.exp(-1.0*0.318*secx)
#	J[52]=1.907E-06*cosx**(1.244)*numpy.exp(-1.0*0.335*secx)
#	J[53]=2.485E-06*cosx**(1.196)*numpy.exp(-1.0*0.328*secx)
#	J[54]=4.095E-06*cosx**(1.111)*numpy.exp(-1.0*0.316*secx)
#	J[55]=1.135E-05*cosx**(0.974)*numpy.exp(-1.0*0.309*secx)
#	J[56]=7.549E-06*cosx**(1.015)*numpy.exp(-1.0*0.324*secx)
#	J[57]=3.363E-06*cosx**(1.296)*numpy.exp(-1.0*0.322*secx)
#	J[61]=7.537E-04*cosx**(0.499)*numpy.exp(-1.0*0.266*secx)
	
	# from MAC spectral analysis and Mainz database (xsproc.py)
	J[1] = 2.3706768705670786e-05
	J[2] = 0.0
	J[4] = 0.0011216654096221574
	J[5] = 0.0024080663555999995
	J[6] = 0.02122656504799999
	J[7] = 0.0002830012299999995
	J[12] = 6.155370761783531e-06	
	
	return J

def zenith(time, lat, lon):
    
    # time is a float in s
    # lat is latitude and lon is longitude in degrees

    # solar declination angle
    dec = 23.79 # in [\deg]
    pi = 4.0*numpy.arctan(1.0) # 1pi=180deg=3.14. For future calculation
    radian = 1.80e+02/pi # 1 unit rad; use this to convert [\deg] to [\rad]
    dec = dec/radian # in [\rad]

    # latitude and longitude conversion from degrees to radian
    lat = lat/radian # in [\rad]
    lon = lon/radian # in [\rad]

    # local hour angle (lha): representing the time of the day
    lha = (1.0+time/4.32e+04)*pi + lon

    # need to check this
    theta = (numpy.arccos(numpy.cos(lha)*numpy.cos(dec)*numpy.cos(lat)+
    		numpy.sin(dec)*numpy.sin(lat)))
    sinld = numpy.sin(lat)*numpy.sin(dec)
    cosld = numpy.cos(lat)*numpy.cos(dec)
    cosx = (numpy.cos(lha)*cosld)+sinld
    cosx = float(numpy.cos(theta))
    secx = 1.0E+0/(cosx+1.0E-30)

    
    return (theta, secx, cosx)