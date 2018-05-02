'''
	       k_aab
     A + A <===> B, +dP
	       k_baa

	       k_abc
     A + B <===> C, +dP
	       k_cab

	      k_aaac
A + A + A <====> C, 2 dP
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
  
k_aab = k_aab0 + alpha_aab P^2
k_baa = k_baa0 + alpha_baa P^2
k_abc = k_abc0 + alpha_abc P^2
k_cab = k_cab0 + alpha_cab P^2
'''

import numpy as np
from functools import partial


class cell (object):
	def __init__(self, constants, parameter, initial_condition):
		self.E = 0.0
		self.P = 0.0

		self.P_a = initial_condition.P_a
		self.P_b = initial_condition.P_b
		self.P_c = initial_condition.P_c

		self.n_a = initial_condition.n_a
		self.n_b = initial_condition.n_b
		self.n_c = initial_condition.n_c

		self.constants = constants
		self.parameter = parameter

		self.internal_update = partial(self.internal_update_in_general, 
			parameter.epsilon_0, 
			constants.k_aab0, constants.k_baa0, constants.k_abc0, constants.k_cab0, constants.k_caaa0, constants.k_aaac0,
			constants.alpha_aab, constants.alpha_baa, constants.alpha_abc, constants.alpha_cab, constants.alpha_aaac, constants.alpha_caaa, constants.dP,
			constants.epsilon_ra, constants.epsilon_rb, constants.epsilon_rc, 
			constants.tau_a, constants.tau_b, constants.tau_c, 
			parameter.dt)

		self.neighbors = []


	def internal_update_in_general(self, epsilon_0, k_aab0, k_baa0, k_abc0, k_cab0, k_caaa0, k_aaac0, alpha_aab, alpha_baa, alpha_abc, alpha_cab, alpha_aaac, alpha_caaa, dP, epsilon_ra, epsilon_rb, epsilon_rc, tau_a, tau_b, tau_c, dt, q):
		
		epsilon_rges = self.n_a * epsilon_ra + self.n_b * epsilon_rb + self.n_c * epsilon_rc
		self.E = (q - self.P * epsilon_rges)/ epsilon_0

		# calculate updates
		
		P_squared = np.dot(self.P, self.P)

		k_aab = k_aab0 + alpha_aab * P_squared
		k_baa = k_baa0 + alpha_baa * P_squared
		k_abc = k_abc0 + alpha_abc * P_squared
		k_cab = k_cab0 + alpha_cab * P_squared
		k_caaa = k_caaa0 + alpha_caaa * P_squared
		k_aaac = k_aaac0 + alpha_aaac * P_squared

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
		
		r_1 = k_aab * self.n_a**2 - k_baa * self.n_b
		r_2 = k_abc * self.n_a * self.n_b - k_cab * self.n_c
		r_3 = k_aaac * self.n_a**3 - k_caaa * self.n_c

		dn_a = -2*r_1 - r_2 - 3*r_3
		dn_b = r_1 - r_2
		dn_c = r_2 + r_3

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