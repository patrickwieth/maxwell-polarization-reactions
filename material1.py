import numpy as np
from functools import partial

class constants:
	def __init__(self, k_aab0, k_baa0, k_abc0, k_cab0, k_caaa0, k_aaac0, alpha_aab, alpha_baa, alpha_abc, alpha_cab, alpha_aaac, alpha_caaa, epsilon_ra, epsilon_rb, epsilon_rc, tau_a, tau_b, tau_c, D_na, D_nb, D_nc, D_Pa, D_Pb, D_Pc):
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
		self.alpha_aaac = alpha_aaac
		self.alpha_caaa = alpha_caaa
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


simple = constants(
					k_aab0 = 0, k_baa0 = 0, k_abc0 = 0, k_cab0 = 0, k_caaa0 = 0, k_aaac0 = 0,
					alpha_aab = 0.0, alpha_baa = 0.0, alpha_abc = 0.0, alpha_cab = -0.0, alpha_aaac = 0, alpha_caaa = 0,
					epsilon_ra = 1.0, epsilon_rb = 2.0, epsilon_rc = 3.0,
					tau_a = 1.0, tau_b = 1.0, tau_c = 1.0,
					D_na = 0.1, D_nb = 0.1, D_nc = 0.1, D_Pa = 0.1, D_Pb = 0.1, D_Pc = 0.1)

twopeaks = constants(
					k_aab0 = 0, k_baa0 = 0, k_abc0 = 0, k_cab0 = 0, k_caaa0 = 0, k_aaac0 = 0,
					alpha_aab = 0.0, alpha_baa = 0.0, alpha_abc = 0.0, alpha_cab = -0.0, alpha_aaac = 0, alpha_caaa = 0,
					epsilon_ra = 1.0, epsilon_rb = 30.0, epsilon_rc = 0.0,
					tau_a = 1.0, tau_b = 100.0, tau_c = 1.0,
					D_na = 0.1, D_nb = 0.1, D_nc = 0.1, D_Pa = 0.1, D_Pb = 0.1, D_Pc = 0.1)

monoalcohol = constants(
					k_aab0 = 0.004, k_baa0 = 0.5, k_abc0 = 0.005, k_cab0 = 0.005, k_caaa0 = 0.004, k_aaac0 = 0.0001,
					alpha_aab = 0.3, alpha_baa = -0.3, alpha_abc = 0.3, alpha_cab = -0.3, alpha_aaac = 0.3, alpha_caaa = -0.3,
					epsilon_ra = 1.0, epsilon_rb = 3.0, epsilon_rc = 9.0,
					tau_a = 1, tau_b = 3.0, tau_c = 5.0,
					D_na = 0.00028, D_nb = 0.05, D_nc = 0.05, D_Pa = 0.00028, D_Pb = 0.05, D_Pc = 0.05)

arbitrary = constants(
					k_aab0 = 0.004, k_baa0 = 0.5, k_abc0 = 0.005, k_cab0 = 0.005, k_caaa0 = 0.004, k_aaac0 = 0.0001,
					alpha_aab = 0.3, alpha_baa = -0.3, alpha_abc = 0.3, alpha_cab = -0.3, alpha_aaac = 0.3, alpha_caaa = -0.3,
					epsilon_ra = 1.0, epsilon_rb = 10.0, epsilon_rc = 100.0,
					tau_a = 1.0, tau_b = 2.0, tau_c = 3.0,
					D_na = 0.00028, D_nb = 0.05, D_nc = 0.05, D_Pa = 0, D_Pb = 0, D_Pc = 0)


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
			constants.alpha_aab, constants.alpha_baa, constants.alpha_abc, constants.alpha_cab, constants.alpha_aaac, constants.alpha_caaa,
			constants.epsilon_ra, constants.epsilon_rb, constants.epsilon_rc,
			constants.tau_a, constants.tau_b, constants.tau_c,
			parameter.dt)

		self.neighbors = []


	def internal_update_in_general(self, epsilon_0, k_aab0, k_baa0, k_abc0, k_cab0, k_caaa0, k_aaac0, alpha_aab, alpha_baa, alpha_abc, alpha_cab, alpha_aaac, alpha_caaa, epsilon_ra, epsilon_rb, epsilon_rc, tau_a, tau_b, tau_c, dt, q):

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
			#print("k_aab < 0")
		if k_baa < 0:
			k_baa = 0
			#print("k_baa < 0")
		if k_abc < 0:
			k_abc = 0
			#print("k_abc < 0")
		if k_cab < 0:
			k_cab = 0
			#print("k_cab < 0")
		if k_caaa < 0:
			k_caaa = 0
			#print("k_caaa < 0")
		if k_aaac < 0:
			k_aaac = 0
			#print("k_aaac < 0")


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

		dP_a = (self.n_a*epsilon_ra/epsilon_rges * epsilon_0 * self.E - self.P_a)/tau_a
		dP_b = (self.n_b*epsilon_rb/epsilon_rges * epsilon_0 * self.E - self.P_b)/tau_b
		dP_c = (self.n_c*epsilon_rc/epsilon_rges * epsilon_0 * self.E - self.P_c)/tau_c


		# apply updates
		self.n_a += dn_a * dt
		self.n_b += dn_b * dt
		self.n_c += dn_c * dt

		self.P_a += dP_a * dt
		self.P_b += dP_b * dt
		self.P_c += dP_c * dt

		self.P = self.P_a + self.P_b + self.P_c



def load_function(loaded):
	def insert_loaded(cell, idx):
		cell.P_a = loaded[idx, 0]
		cell.P_b = loaded[idx, 1]
		cell.P_c = loaded[idx, 2]
		cell.n_a = loaded[idx, 3]
		cell.n_b = loaded[idx, 4]
		cell.n_c = loaded[idx, 5]

	return insert_loaded

def get_data(grid):
	return np.array([[
			grid.cells[x].P_a,
			grid.cells[x].P_b,
			grid.cells[x].P_c,
			grid.cells[x].n_a,
			grid.cells[x].n_b,
			grid.cells[x].n_c]
			for x in range(grid.size)])
