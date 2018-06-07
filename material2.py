import numpy as np
from functools import partial

class constants:
	def __init__(self, k_aab0, k_baa0, alpha, beta, gamma, delta, epsilon_ra, epsilon_rb, tau_a, tau_b, D_na, D_nb, D_Pa, D_Pb):
		self.k_aab0 =  k_aab0
		self.k_baa0 = k_baa0
		self.alpha = alpha
		self.beta = beta
		self.gamma = gamma
		self.delta = delta
		self.epsilon_ra = epsilon_ra
		self.epsilon_rb = epsilon_rb
		self.tau_a = tau_a
		self.tau_b = tau_b
		self.D_na = D_na
		self.D_Pa = D_Pa
		self.D_nb = D_nb
		self.D_Pb = D_Pb
		
		if self.tau_a == 0:
			print("please don't set relaxation times to 0, if you want instantaneous processes set tau = dt ; setting tau_a = 1")
			self.tau_a = 1
		if self.tau_b == 0:
			print("please don't set relaxation times to 0, if you want instantaneous processes set tau = dt ; setting tau_b = 1")
			self.tau_b = 1


simple = constants(
					k_aab0 = 0, k_baa0 = 0, 
					alpha = 0, beta = 0, gamma = 0, delta = 0, 
					epsilon_ra = 1.0, epsilon_rb = 2.0, 
					tau_a = 1.0, tau_b = 1.0, 
					D_na = 0.1, D_nb = 0.1, D_Pa = 0.1, D_Pb = 0.1)

monoalcohol = constants(
					k_aab0 = 0.01, k_baa0 = 0.2,
					alpha = 0.01, beta = 0.01, gamma = 0.01, delta = 0.01, 
					epsilon_ra = 1.0, epsilon_rb = 3,
					tau_a = 1.0, tau_b = 2.0,
					D_na = 0.00028, D_nb = 0.05, D_Pa = 0.00028, D_Pb = 0.1)

arbitrary = constants(
					k_aab0 = 0, k_baa0 = 0, 
					alpha = 0, beta = 0, gamma = 0, delta = 0, 
					epsilon_ra = 1.0, epsilon_rb = 3.0, 
					tau_a = 1.0, tau_b = 1.0, 
					D_na = 0.1, D_nb = 0.1, D_Pa = 0.1, D_Pb = 0.1)


class cell (object):
	def __init__(self, constants, parameter, initial_condition):
		self.E = 0.0
		self.P = 0.0

		self.P_a = initial_condition.P_a
		self.P_b = initial_condition.P_b

		self.n_a = initial_condition.n_a
		self.n_b = initial_condition.n_b

		self.constants = constants
		self.parameter = parameter

		self.internal_update = partial(self.internal_update_in_general, 
			parameter.epsilon_0, 
			constants.k_aab0, constants.k_baa0,
			constants.alpha, constants.beta, constants.gamma, constants.delta,
			constants.epsilon_ra, constants.epsilon_rb,
			constants.tau_a, constants.tau_b,
			parameter.dt)

		self.neighbors = []


	def internal_update_in_general(self, epsilon_0, k_aab0, k_baa0, alpha, beta, gamma, delta, epsilon_ra, epsilon_rb, tau_a, tau_b, dt, q):
		
		epsilon_rges = self.n_a * epsilon_ra + self.n_b * epsilon_rb
		self.E = (q - self.P * epsilon_rges)/ epsilon_0

		# calculate updates
		P_squared = np.dot(self.P, self.P)

		k_aab = k_aab0 + alpha*P_squared + beta*self.n_a**2
		k_baa = k_baa0 - alpha*P_squared + beta*self.n_b**2

		if k_aab < 0:
			k_aab = 0
			print("k_aab < 0")
		if k_baa < 0:
			k_baa = 0
			print("k_baa < 0")		

		r_1 = k_aab * self.n_a**2 - k_baa * self.n_b

		dn_a = -2*r_1
		dn_b = r_1
		
		dP_a = (self.n_a*epsilon_ra/epsilon_rges*epsilon_0*self.E - self.P_a)/tau_a + gamma*self.n_a*self.E + delta*self.n_a*self.P_a 
		dP_b = (self.n_b*epsilon_rb/epsilon_rges*epsilon_0*self.E - self.P_b)/tau_b + gamma*self.n_b*self.E + delta*self.n_b*self.P_b


		# apply updates
		self.n_a += dn_a * dt
		self.n_b += dn_b * dt

		self.P_a += dP_a * dt
		self.P_b += dP_b * dt

		self.P = self.P_a + self.P_b



def load_function(loaded):
	def insert_loaded(cell, idx):
		cell.P_a = loaded[idx, 0]
		cell.P_b = loaded[idx, 1]
		cell.n_a = loaded[idx, 2]
		cell.n_b = loaded[idx, 3]

	return insert_loaded

def get_data(grid):
	return np.array([[
			grid.cells[x].P_a,
			grid.cells[x].P_b,
			grid.cells[x].n_a,
			grid.cells[x].n_b] 
			for x in range(grid.size)]) 

