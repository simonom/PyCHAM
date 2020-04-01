'''module to calculate the integration needed for coagulation, called on by coag'''

import numpy as np
import scipy.constants as si
import scipy.integrate as integ
import pdb

def W_k_int(sbr_m, sbr_mt, sbr_sum_m, T, A_H, sbn):
	
	# -----------------------------------------------------------
	# inputs:
	
	# sbr_m - size bin radii spead across shells (m)
	# sbr_mt - transpose of sbr_m (m)
	# sbr_sum_m - sum of sbr_m and sbr_mt (m)
	# T - temperature (K)
	# A_H - Hamaker constant (dimensionless)
	# sbn - number of size bins
	# -----------------------------------------------------------
	
	# define the integrand in 15.43
	def integrand(x, a, b, c, m, A_H, rj, T):
	
		# differential of Ep with r=x/b substituted in prior to 
		# differentiation

		# van der Waals interaction potential (15.45)
		Ep0 = (-A_H/6.0)*((a)/((x/rj)**m-b)+
			(a)/((x/rj)**m-c)+
			np.log(((x/rj)**m-b)/((x/rj)**m-c)))
	
		# Ep with r=x/rj substituted in
		# first differential of Ep, all terms	
		Ep1 =  (-A_H/6.0)*(((-m*a*(x/(rj**m)))/(((x/rj)**m-b)**m)-
			(m*a*(x/(rj**m)))/(((x/rj)**m-c)**m))+
			((((m*x)/(rj**m))*((x**m)/(rj**m)-c)-
			((m*x)/(rj**m))*((x**m)/(rj**m)-b))/
			(((x**m)/(rj**m)-c)*((x**m)/(rj**m)-b))))	
	
		# second differential of third term
		f = (m*x*b)/(rj**m)-(m*x*c)/(rj**m)
		fd = (m*b)/(rj**m)-(m*c)/(rj**m)
		g = (x/rj)**4.0-((x/rj)**m)*b-((x/rj)**m)*c+c*b
		gd = (4.0*x**3.0)/(rj**4.0)-(m*x*b)/(rj**m)-(m*x*c)/(rj**m)
		
		# second differential of all terms
		Ep2 =  ((-A_H/6.0)*
			((((6.0*a*x**4.0)/(rj**6.0)-(4.0*x**m*a*b)/
			(rj**4.0)-(m*a*b**m)/(rj**m))/
			(((x/rj)**m-b)**4.0))+ 
			(((6.0*a*x**4.0)/(rj**6.0)-(4.0*x**m*a*c)/(rj**4.0)-
			(m*a*c**m)/(rj**m))/
			(((x/rj)**m-c)**4.0))+
			(fd*g-gd*f)/(g**m))) 
		# integration part of eq. 15.43
		return (((Ep1+(x/rj)*Ep2)*
			np.exp((-1.0/(si.k*T))*(((x/(rj*m))*Ep1+Ep0))))*
			((x/rj)**m))

	# define the integrand in 15.44 
	def integrand2(x, a, b, d, rj, e, m, c, T):
	
		# van der Waals interaction potential (15.45)
		Ep0 = (-A_H/6.0)*((m*a)/((x/rj)**m-b)+
			(m*a)/((x/rj)**m-c)+
			np.log(((x/rj)**m-b)/((x/rj)**m-c)))

		return (1.0+((2.6*a)/(b))*(((a)/(e*((x/rj)-d)))**0.5)+
		(a)/(e*((x/rj)-d)))*np.exp(Ep0/(si.k*T))*(1.0/((x/rj)**m))

	
	# numerator for first two terms in 15.45
	a = 2.0*(sbr_m*sbr_mt) 
	# square of sum of size bin radii (m2)
	b = sbr_sum_m**2.0
	# square of difference in size bin radii (m2)
	c = (sbr_m-sbr_mt)**2.0
	# difference in size bin radii (m2)
	d = sbr_m-sbr_mt
	# integration limits for 15.43
	ilu = 1.0+sbr_m/sbr_mt
	ill = 0.0	
	# constant for integration in eq. 15.43
	m = 2.0
	# empty matrices for the integration in 15.43 and the eq. 15.43 
	W_k_int = np.zeros((sbn, sbn))	
	W_k = np.zeros((sbn, sbn))

	# empty matrices for the integration in 15.44 and the eq. 15.44
	W_c_int = np.zeros((sbn, sbn))
	W_c = np.zeros((sbn, sbn))
		
	# loop through size bins
	for i in range(0, sbn):
		for j in range(i, sbn):	
			
			W_k_int[i,j] = integ.quad(integrand, ill, ilu[i,j], 
			args=(a[i,j], b[i,j], c[i,j], m, A_H, sbr_m[0,j], T))[0]
			W_k[i,j] = (-1.0/(m*b[i,j]*si.k*T))*W_k_int[i,j]

			W_c_int[i,j] = integ.quad(integrand2, ill, ilu[i,j], 
			args=(a[i,j]/2.0, b[i,j], d[i,j], sbr_m[0,j], 
			sbr_sum_m[i,j], m, c[i,j], T))[0]
			W_c[i,j] = 1.0/(sbr_sum_m[i,j]*W_c_int[i,j])

	return W_k, W_c
