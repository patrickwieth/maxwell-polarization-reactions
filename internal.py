'''
	  k_aab
A + A <===> B, dP_aab
	  k_baa

	  k_abc
A + B <===> C, dP_abc
	  k_cab

r_1 = k_aab * n_a^2 - k_baa * n_b
r_2 = k_abc * n_a * n_b - k_cab n_c

dna = -2 r_1 - r_2
dnb = r_1 - r_2
dnc = r_2

dPa = 1/tau_a (epsilon_0 * n_a * epsilon_ra E - P_a)
dPb = 1/tau_b (epsilon_0 * n_b * epsilon_rb E - P_b) + r_1 * dP_aab
dPc = 1/tau_c (epsilon_0 * n_c * epsilon_rc E - P_c) + r_2 * dP_abc

P = P_a + P_b + P_c

k_aab = k_aab0 + alpha_aab P^2 + beta_aab P_a^2
k_baa = k_baa0 + alpha_baa P^2
k_abc = k_abc0 + alpha_abc P^2 + beta_abc P_a P_b
k_cab = k_cab0 + alpha_cab P^2
'''

import numpy as np
from functools import partial

class material_constants:
	def __init__(self, k_aab0, k_baa0, k_abc0, k_cab0, alpha_aab, alpha_baa, alpha_abc, alpha_cab, beta_aab, beta_abc, epsilon_ra, epsilon_rb, epsilon_rc, tau_a, tau_b, tau_c, D_a, D_b, D_c):
		self.k_aab0 =  k_aab0
		self.k_baa0 = k_baa0
		self.k_abc0 = k_abc0
		self.k_cab0 = k_cab0
		self.alpha_aab = alpha_aab
		self.alpha_baa = alpha_baa
		self.alpha_abc = alpha_abc
		self.alpha_cab = alpha_cab
		self.beta_aab = beta_aab
		self.beta_abc = beta_abc
		self.epsilon_ra = epsilon_ra
		self.epsilon_rb = epsilon_rb
		self.epsilon_rc = epsilon_rc
		self.tau_a = tau_a
		self.tau_b = tau_b
		self.tau_c = tau_c
		self.D_a = D_a
		self.D_b = D_b
		self.D_c = D_c
		
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
		print('Courant number for this run is: ', self.C2)
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
		self.prev_P = 0.0 #np.zeros(3)
		self.prev2_P = 0.0 #np.zeros(3)

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
			constants.k_aab0, constants.k_baa0, constants.k_abc0, constants.k_cab0, 
			constants.alpha_aab, constants.alpha_baa, constants.alpha_abc, constants.alpha_cab, 
			constants.beta_aab, constants.beta_abc, 
			constants.epsilon_ra, constants.epsilon_rb, constants.epsilon_rc, 
			constants.tau_a, constants.tau_b, constants.tau_c, 
			parameter.dt)

		self.neighbors = []


	def internal_update_in_general(self, epsilon_0, k_aab0, k_baa0, k_abc0, k_cab0, alpha_aab, alpha_baa, alpha_abc, alpha_cab, beta_aab, beta_abc, epsilon_ra, epsilon_rb, epsilon_rc, tau_a, tau_b, tau_c, dt):
		# write history
		self.prev2_P = self.prev_P
		self.prev_P = self.P

		epsilon_rges = self.n_a * epsilon_ra + self.n_b * epsilon_rb + self.n_c * epsilon_rc

		# calculate updates
		
		P_squared = np.dot(self.P, self.P)

		k_aab = k_aab0 + alpha_aab * P_squared + beta_aab * np.dot(self.P_a, self.P_a)
		k_baa = k_baa0 + alpha_baa * P_squared
		k_abc = k_abc0 + alpha_abc * P_squared + beta_abc * np.dot(self.P_a, self.P_b)
		k_cab = k_cab0 + alpha_cab * P_squared

		if k_aab < 0:
			print("k_aab < 0")
		if k_baa < 0:
			print("k_baa < 0")
		if k_abc < 0:
			print("k_abc < 0")
		if k_cab < 0:
			print("k_cab < 0")

		r_1 = k_aab * self.n_a**2 - k_baa * self.n_b
		r_2 = k_abc * self.n_a * self.n_b - k_cab * self.n_c

		dn_a = -2*r_1 - r_2
		dn_b = r_1 - r_2
		dn_c = r_2
		
		dP_a = self.n_a*epsilon_ra/epsilon_rges * (epsilon_0 * self.E - self.P)/tau_a
		dP_b = self.n_b*epsilon_rb/epsilon_rges * (epsilon_0 * self.E - self.P)/tau_b
		dP_c = self.n_c*epsilon_rc/epsilon_rges * (epsilon_0 * self.E - self.P)/tau_c

		# apply updates
		self.n_a += dn_a * dt
		self.n_b += dn_b * dt
		self.n_c += dn_c * dt

		self.P_a += dP_a * dt
		self.P_b += dP_b * dt
		self.P_c += dP_c * dt

		self.P = self.P_a + self.P_b + self.P_c

		#print(self.n_a, self.n_b, self.n_c)




'''
arbitrary_material = material_constants(
					k_aab0 = 1, k_baa0 = 1, k_abc0 = 1, k_cab0 = 1, 
					alpha_aab = 0, alpha_baa = 0, alpha_abc = 0, alpha_cab = 0, beta_aab = 0, beta_abc = 0,
					epsilon_ra = 1, epsilon_rb = 1, epsilon_rc = 1,
					tau_a = 1, tau_b = 1, tau_c = 1)

arbitrary_simulation = simulation_parameters(epsilon_0 = 1, c = 1, dx = 1, dt = 0.1)

test = cell(arbitrary_material, arbitrary_simulation)

test.internal_update(np.empty(3), 0.1)
print(test.n_a, test.n_b)

test.internal_update(np.empty(3), 0.1)
print(test.n_a, test.n_b)

test.internal_update(np.empty(3), 0.1)
print(test.n_a, test.n_b)

test.internal_update(np.empty(3), 0.1)
print(test.n_a, test.n_b)
'''