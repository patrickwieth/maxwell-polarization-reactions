import numpy as np
from functools import partial

class constants:
	def __init__(self, k_aab0, k_baa0, k_abc0, k_cab0, k_caaa0, k_aaac0, alpha_aab, alpha_baa, alpha_abc, alpha_cab, alpha_aaac, alpha_caaa, dP, epsilon_ra, epsilon_rb, epsilon_rc, tau_a, tau_b, tau_c, I_a, I_b, I_c):
		self.k_aab0 = k_aab0
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
		self.dP = dP
		self.epsilon_ra = epsilon_ra
		self.epsilon_rb = epsilon_rb
		self.epsilon_rc = epsilon_rc
		self.tau_a = tau_a
		self.tau_b = tau_b
		self.tau_c = tau_c
		self.I_a = I_a
		self.I_b = I_b
		self.I_c = I_c
		
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
					epsilon_ra = 1.0, epsilon_rb = 2.0, epsilon_rc = 3.0, dP = 0.0,
					tau_a = 1.0, tau_b = 1.0, tau_c = 1.0,
					I_a = 10, I_b = 10, I_c = 10)

monoalcohol = constants(
					k_aab0 = 0.001, k_baa0 = 0.002, k_abc0 = 0.005, k_cab0 = 0.005, k_caaa0 = 0.004, k_aaac0 = 0.0001,
					alpha_aab = 0.01, alpha_baa = -0.01, alpha_abc = 0.01, alpha_cab = -0.01, alpha_aaac = 0.1, alpha_caaa = -0.01, 
					epsilon_ra = 1.0, epsilon_rb = 3, epsilon_rc = 5, dP = 0.0,
					tau_a = 1.0, tau_b = 2.0, tau_c = 3.0,
					I_a = 10000, I_b = 10000, I_c = 10000)

arbitrary = constants(
					k_aab0 = 0.001, k_baa0 = 0.09, k_abc0 = 0.001, k_cab0 = 0.09, k_aaac0 = 0.001, k_caaa0 = 0.09,
					alpha_aab = 1, alpha_baa = -1, alpha_abc = 1, alpha_cab = -1, alpha_aaac = 1, alpha_caaa = -1, 
					epsilon_ra = 1.0, epsilon_rb = 10.0, epsilon_rc = 100.0, dP = 0.0,
					tau_a = 25.0, tau_b = 5.0, tau_c = 1.0,
					I_a = 1, I_b = 1, I_c = 1)


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
		self.rot_a = initial_condition.rot_a
		self.rot_b = initial_condition.rot_b
		self.rot_c = initial_condition.rot_c

		self.dipolar_field = np.zeros(3)

		self.constants = constants
		self.parameter = parameter

		self.internal_update = partial(self.internal_update_in_general, 
			parameter.epsilon_0, 
			constants.k_aab0, constants.k_baa0, constants.k_abc0, constants.k_cab0, constants.k_caaa0, constants.k_aaac0,
			constants.alpha_aab, constants.alpha_baa, constants.alpha_abc, constants.alpha_cab, constants.alpha_aaac, constants.alpha_caaa, constants.dP,
			constants.epsilon_ra, constants.epsilon_rb, constants.epsilon_rc, 
			constants.tau_a, constants.tau_b, constants.tau_c, 
			constants.I_a, constants.I_b, constants.I_c,
			parameter.dt)

		self.neighbors = []


	def internal_update_in_general(self, epsilon_0, k_aab0, k_baa0, k_abc0, k_cab0, k_caaa0, k_aaac0, alpha_aab, alpha_baa, alpha_abc, alpha_cab, alpha_aaac, alpha_caaa, dP, epsilon_ra, epsilon_rb, epsilon_rc, tau_a, tau_b, tau_c, I_a, I_b, I_c, dt, q):
		# permittivity of whole cell
		epsilon_rges = self.n_a * epsilon_ra + self.n_b * epsilon_rb + self.n_c * epsilon_rc

		# electrical field at cell
		self.E = q #(q - self.P * epsilon_rges)/ epsilon_0

		#print(q)

		# squared polarization
		P_squared = np.dot(self.P, self.P)

		# reaction coefficients
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
		if k_caaa < 0:
			k_caaa = 0
			print("k_caaa < 0")
		if k_aaac < 0:
			k_aaac = 0
			print("k_aaac < 0")
		
		# reaction rates
		r_1 = k_aab * self.n_a**2 - k_baa * self.n_b
		r_2 = k_abc * self.n_a * self.n_b - k_cab * self.n_c
		r_3 = k_aaac * self.n_a**3 - k_caaa * self.n_c

		# delta concentrations from reactions
		dn_a = -2*r_1 - r_2 - 3*r_3
		dn_b = r_1 - r_2
		dn_c = r_2 + r_3
		
		P_eq = 0.0

		# delta polarization (relaxtion of P to E and P_eq if > 0)
		dP_a = (self.n_a*epsilon_ra/epsilon_rges * epsilon_0 * self.E - self.P_a + self.n_a*epsilon_ra/epsilon_rges * P_eq * self.P_a/np.linalg.norm(self.P_a) )/tau_a
		dP_b = (self.n_b*epsilon_rb/epsilon_rges * epsilon_0 * self.E - self.P_b + self.n_b*epsilon_ra/epsilon_rges * P_eq * self.P_b/np.linalg.norm(self.P_b) )/tau_b
		dP_c = (self.n_c*epsilon_rc/epsilon_rges * epsilon_0 * self.E - self.P_c + self.n_c*epsilon_ra/epsilon_rges * P_eq * self.P_c/np.linalg.norm(self.P_c) )/tau_c

		# delta rotation (relaxation to 0 + torque from dipol-dipol interaction)
		dRot_a = -self.rot_a / tau_a + np.cross(self.P_a, self.dipolar_field + self.E)
		dRot_b = -self.rot_b / tau_b + np.cross(self.P_b, self.dipolar_field + self.E)
		dRot_c = -self.rot_c / tau_c + np.cross(self.P_c, self.dipolar_field + self.E)

		# apply updates
		self.n_a += dn_a * dt
		self.n_b += dn_b * dt
		self.n_c += dn_c * dt

		self.P_a += dP_a * dt
		self.P_b += dP_b * dt
		self.P_c += dP_c * dt

		self.rot_a += dRot_a * dt / I_a
		self.rot_b += dRot_b * dt / I_b
		self.rot_c += dRot_c * dt / I_c

		'''
		threshold = 1.00
		if self.rot_a > threshold :
			self.rot_a = threshold
			print("rot_a > 90deg - set to 0")
		if self.rot_b > threshold :
			self.rot_b = threshold
			print("rot_b > 90deg - set to 0")
		if self.rot_c > threshold :
			self.rot_c = threshold
			print("rot_c > 90deg - set to 0")
		if self.rot_a < -threshold :
			self.rot_a = -threshold
			print("rot_a > 90deg - set to 0")
		if self.rot_b < -threshold :
			self.rot_b = -threshold
			print("rot_b > 90deg - set to 0")
		if self.rot_c < -threshold :
			self.rot_c = -threshold
			print("rot_c > 90deg - set to 0")
		'''

		def rotmatrix(rotation):
			if np.linalg.norm(rotation) == 0:
				return np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])	
					

			axis = rotation / np.linalg.norm(rotation)
			theta = np.linalg.norm(rotation)

			costh = np.cos(theta)
			sinth = np.sin(theta)
			
			return np.array([
				[costh + (1-costh)*axis[0]**2, axis[0]*axis[1]*(1-costh) - axis[2]*sinth, axis[0]*axis[2]*(1-costh) + axis[1]*sinth],
				[axis[0]*axis[1]*(1-costh) + axis[2]*sinth, costh + (1-costh)*axis[1]**2, axis[1]*axis[2]*(1-costh) - axis[0]*sinth],
				[axis[0]*axis[2]*(1-costh) + axis[1]*sinth, axis[1]*axis[2]*(1-costh) - axis[0]*sinth, costh + (1-costh)*axis[2]**2]
				])

		#print(rotmatrix([0.1, 0.1, 0.0]))
		#print(rotmatrix([-0.1, -0.1, 0.0]))

		self.P_a = np.tensordot(rotmatrix(self.rot_a), self.P_a, axes=1)
		self.P_b = np.tensordot(rotmatrix(self.rot_b), self.P_b, axes=1)
		self.P_c = np.tensordot(rotmatrix(self.rot_c), self.P_c, axes=1)

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

