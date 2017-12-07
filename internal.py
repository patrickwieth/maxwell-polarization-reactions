'''
	  k_aab
A + A <===> B, +dP
	  k_baa

	  k_abc
A + B <===> C, +dP
	  k_cab

	  
	 C ===> A + A + A, -2 dP
	  k_caaa

r_1 = k_aab * n_a^2 - k_baa * n_b
r_2 = k_abc * n_a * n_b - k_cab n_c
r_3 = k_caaa * n_c - k_aaac * n_a * n_a * n_a

dna = -2 r_1 - r_2 + 3 r_3
dnb = r_1 - r_2
dnc = r_2 - r_3

dPa = 1/tau_a * epsilon_ra/epsilon_rges * (epsilon_0 * n_a * epsilon_ra E - P_a) + dp(2 r_3 - r_1 - 0.5 * r_2)
dPb = 1/tau_b * epsilon_rb/epsilon_rges * (epsilon_0 * n_b * epsilon_rb E - P_b) + dp(r_1 - 0.5 * r_2)
dPc = 1/tau_c * epsilon_rc/epsilon_rges * (epsilon_0 * n_c * epsilon_rc E - P_c) + dp(r_2 - 2 * r_3)

P = P_a + P_b + P_c

k_aab = k_aab0 + alpha_aab P^2 + beta_aab P_a^2
k_baa = k_baa0 + alpha_baa P^2
k_abc = k_abc0 + alpha_abc P^2 + beta_abc P_a P_b
k_cab = k_cab0 + alpha_cab P^2
'''

import numpy as np
from functools import partial

class material_constants:
	def __init__(self, k_aab0, k_baa0, k_abc0, k_cab0, k_caaa0, k_aaac0, alpha_aab, alpha_baa, alpha_abc, alpha_cab, beta_aab, beta_abc, dP, epsilon_ra, epsilon_rb, epsilon_rc, tau_a, tau_b, tau_c, D_na, D_nb, D_nc, D_Pa, D_Pb, D_Pc):
		self.k_aab0 =  k_aab0
		self.k_baa0 = k_baa0
		self.k_abc0 = k_abc0
		self.k_cab0 = k_cab0
		self.k_caaa0 = k_caaa0
		self.k_aaac0 = k_aaac0
		self.alpha_aab = alpha_aab
		self.alpha_baa = alpha_baa
		self.alpha_abc = alpha_abc
		self.alpha_cab = alpha_cab
		self.beta_aab = beta_aab
		self.beta_abc = beta_abc
		self.dP = dP
		self.epsilon_ra = epsilon_ra
		self.epsilon_rb = epsilon_rb
		self.epsilon_rc = epsilon_rc
		self.tau_a = tau_a
		self.tau_b = tau_b
		self.tau_c = tau_c
		self.D_na = D_na
		self.D_Pa = D_Pa
		self.D_nb = D_nb
		self.D_Pb = D_Pb
		self.D_nc = D_nc
		self.D_Pc = D_Pc
		
		if self.tau_a == 0:
			print("please don't set relaxation times to 0, if you want instantaneous processes set tau = dt - setting tau_a = 1")
			self.tau_a = 1
		if self.tau_b == 0:
			print("please don't set relaxation times to 0, if you want instantaneous processes set tau = dt - setting tau_b = 1")
			self.tau_b = 1
		if self.tau_c == 0:
			print("please don't set relaxation times to 0, if you want instantaneous processes set tau = dt - setting tau_c = 1")
			self.tau_c = 1
		

class simulation_parameters:
	def __init__(self, epsilon_0, c, dx, dt, wave_length):
		self.C2 = (c * dt / dx)**2 # courant number
		#print('Courant number for this run is: ', self.C2)
		if self.C2 > 1:	print("Courant number is greater than 1, numerical shittyness unpreventable!")

		self.wave_length = wave_length
		self.epsilon_0 = epsilon_0
		self.dt = dt
		self.dx = dx
		self.c = c


class cell (object):
	def __init__(self, constants, parameter, initial_condition):
		self.E = 0.0 #np.zeros(3)

		self.P = 0.0 #np.zeros(3)
		#self.prev_P = 0.0 #np.zeros(3)
		#self.prev2_P = 0.0 #np.zeros(3)

		self.P_a = initial_condition.P_a #np.zeros(3)
		self.P_b = initial_condition.P_b #np.zeros(3)
		self.P_c = initial_condition.P_c #np.zeros(3)

		self.n_a = initial_condition.n_a
		self.n_b = initial_condition.n_b
		self.n_c = initial_condition.n_c

		self.constants = constants
		self.parameter = parameter

		self.internal_update = partial(self.internal_update_in_general, 
			parameter.epsilon_0, 
			constants.k_aab0, constants.k_baa0, constants.k_abc0, constants.k_cab0, constants.k_caaa0, constants.k_aaac0,
			constants.alpha_aab, constants.alpha_baa, constants.alpha_abc, constants.alpha_cab, 
			constants.beta_aab, constants.beta_abc, constants.dP,
			constants.epsilon_ra, constants.epsilon_rb, constants.epsilon_rc, 
			constants.tau_a, constants.tau_b, constants.tau_c, 
			parameter.dt)

		self.neighbors = []


	def internal_update_in_general(self, epsilon_0, k_aab0, k_baa0, k_abc0, k_cab0, k_caaa0, k_aaac0, alpha_aab, alpha_baa, alpha_abc, alpha_cab, beta_aab, beta_abc, dP, epsilon_ra, epsilon_rb, epsilon_rc, tau_a, tau_b, tau_c, dt):
		

		epsilon_rges = self.n_a * epsilon_ra + self.n_b * epsilon_rb + self.n_c * epsilon_rc
		#self.E = q - self.P / epsilon_0 * epsilon_rges

		# calculate updates
		
		P_squared = np.dot(self.P, self.P)

		k_aab = k_aab0 + alpha_aab * P_squared + beta_aab * np.dot(self.P_a, self.P_a)
		k_baa = k_baa0 + alpha_baa * P_squared
		k_abc = k_abc0 + alpha_abc * P_squared + beta_abc * np.dot(self.P_a, self.P_b)
		k_cab = k_cab0 + alpha_cab * P_squared
		k_caaa = k_caaa0
		k_aaac = k_aaac0

		if k_aab < 0:
			k_aab = 0
			print("k_aab < 0")
		if k_baa < 0:
			k_baa = 0
			print("k_baa < 0")
		if k_abc < 0:
			k_abc = 0
			print("k_abc < 0")
		if k_cab < 0:
			k_cab = 0
			print("k_cab < 0")
		
		r_1 = k_aab * self.n_a**2 - k_baa * self.n_b + k_caaa * self.n_c
		r_2 = k_abc * self.n_a * self.n_b - k_cab * self.n_c
		r_3 = k_caaa * self.n_c - k_aaac * self.n_a * self.n_a * self.n_a

		dn_a = -2*r_1 - r_2 + 3*r_3
		dn_b = r_1 - r_2
		dn_c = r_2 - r_3

		'''	normaler turing
		k = -0.005
		l = 0.1
		s = 0.5
		tau = 0.1
		r_1 = l*self.n_a - self.n_a**3 - s*self.n_c + k + alpha_aab * P_squared
		r_2 = (self.n_a - self.n_c) / tau 				+ alpha_abc * P_squared
		r_3 = 0

		dn_a = r_1
		dn_b = - r_1 - r_2
		dn_c = r_2
		'''
		
		dP_a = (self.n_a*epsilon_ra/epsilon_rges * epsilon_0 * self.E - self.P_a)/tau_a + dP * (2*r_1 - r_1 - 0.5 * r_2)
		dP_b = (self.n_b*epsilon_rb/epsilon_rges * epsilon_0 * self.E - self.P_b)/tau_b + dP * (r_1 - 0.5 * r_2)
		dP_c = (self.n_c*epsilon_rc/epsilon_rges * epsilon_0 * self.E - self.P_c)/tau_c + dP * (r_2 - 2 * r_3)


		# apply updates
		self.n_a += dn_a * dt
		self.n_b += dn_b * dt
		self.n_c += dn_c * dt

		self.P_a += dP_a * dt
		self.P_b += dP_b * dt
		self.P_c += dP_c * dt

		self.P = self.P_a + self.P_b + self.P_c