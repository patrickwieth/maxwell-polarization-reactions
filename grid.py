import numpy as np
import math
import matplotlib.pyplot as plt

import internal


class grid:
	def __init__(self, constants, parameters, size):
		self.size = size
		self.parameters = parameters
		self.constants = constants

		self.spawn_grid()

		'''
		self.Efield = np.zeros(size)
		self.last_Efield = np.zeros(size)
		self.very_last_Efield = np.zeros(size)
		self.Pfield = np.zeros(size)
		self.last_Pfield = np.zeros(size)
		self.very_last_Pfield = np.zeros(size)
		'''


	def spawn_grid(self):
		self.cells = np.array([internal.cell(self.constants, self.parameters) for x in range(self.size)])


	def set_inital_Efield(self, Efield):
		#new
		for idx, cell in enumerate(self.cells):
			cell.E = math.sin(idx*self.parameters.dx * 2*math.pi / self.parameters.wave_length) #np.array([0.0, math.sin(idx*self.parameters.dx * 2*math.pi / self.parameters.wave_length), 0.0])
			cell.prev2_E = cell.prev_E
			cell.prev_E = cell.E
		'''
		#old
		for pos in range(len(Efield)):
			Efield[pos] = math.sin(pos*self.parameters.dx * 2*math.pi / self.parameters.wave_length)

		self.very_last_Efield[:] = self.last_Efield[:]
		self.last_Efield[:] = self.Efield
		'''

	def calc_first_step(self, Efield):
		#new
		next_E = np.zeros(size)

		for i in range(1, size-1):
			next_E[i] = self.cells[i].E + 0.5*self.parameters.C2*(self.cells[i-1].E - 2*self.cells[i].E + self.cells[(i+1)%size].E)

		for idx, cell in enumerate(self.cells):
			cell.E = next_E[idx]

		self.cells[ 0].E = 0.0 #np.zeros(3); 
		self.cells[-1].E = 0.0 #np.zeros(3); 

		#C2 = (self.parameters.c * self.parameters.dt / self.parameters.dx)**2 # courant number GET THIS THING EVERYWHERE
		# Special formula for first time step
		'''
		#old
		u = np.zeros(size)

		for i in range(1, size-1):
			u[i] = Efield[i] + 0.5*self.parameters.C2*(Efield[i-1] - 2*Efield[i] + Efield[i+1])
			
		u[0] = 0;  u[size-1] = 0

		return u
		'''

	def update_Pfield(self, Efield, Pfield):
		#new
		for cell in self.cells:
			cell.prev2_P = cell.prev_P
			cell.prev_P = cell.P

		for idx, cell in enumerate(self.cells):
			cell.P += (self.parameters.epsilon_0 * cell.E - cell.P)/self.constants.tau_a
		'''
		# tau * dP/dt + P = epsilon_0 * E
		#old
		self.very_last_Pfield[:] = self.last_Pfield
		self.last_Pfield[:] = self.Pfield	

		return Pfield + (self.parameters.epsilon_0 * Efield - Pfield)/self.constants.tau_a
		'''

	def update_Efield(self, prev_Efield, Efield, prev_Pfield, Pfield, next_Pfield):
		#new
		for cell in self.cells:
			cell.prev2_E = cell.prev_E
			cell.prev_E = cell.E

		for idx, cell in enumerate(self.cells):
			cell.E = 0.0

			if(idx > 0 and idx < size-1):

				left_E = self.cells[idx-1].prev_E
				right_E = self.cells[(idx+1)%self.size].prev_E

				P_tt = (cell.prev2_P - 2*cell.prev_P + cell.P)/(self.parameters.dt**2)
				E_xx = (left_E - 2*cell.prev_E + right_E)/(self.parameters.dx**2)
				t_derivative = self.parameters.c**2 * E_xx - 1/self.parameters.epsilon_0 * P_tt

				cell.E = t_derivative * (self.parameters.dt**2) - cell.prev2_E + 2*cell.prev_E

		#old
		# 1 / c^2 d^2u/dt^2 - d^2u/dx^2 + mu0 d^2P/dt^2 = 0
		# d^2u/dt^2 = c^2 d^2u/dx^2 - 1/e0 d^2P/dt^2
		'''
		self.very_last_Efield[:] = self.last_Efield
		self.last_Efield[:] = self.Efield

		next_Efield = np.zeros(size)

		for pos in range(1, size-1):

			P_tt = (prev_Pfield[pos] - 2*Pfield[pos] + next_Pfield[pos])/(self.parameters.dt**2)
			E_xx = (Efield[pos-1] - 2*Efield[pos] + Efield[pos+1])/(self.parameters.dx**2)
			t_derivative = self.parameters.c**2 * E_xx - 1/self.parameters.epsilon_0 * P_tt

			next_Efield[pos] = t_derivative * (self.parameters.dt**2) - prev_Efield[pos] + 2*Efield[pos]

		return next_Efield
		'''

	def set_border_Efield(self, Efield, time):
		#old
		'''
		Efield[0] = 0.0
		Efield[size-1] = 0.0
		'''

		#new
		self.cells[ 0].E = 0.0 #np.zeros(3)
		self.cells[-1].E = 0.0 #np.zeros(3)

	def run(self, steps, saveInterval):

		self.set_inital_Efield(self.Efield)

		self.Efield = self.calc_first_step(self.Efield)

		self.evolution = np.empty([int(math.ceil(steps/saveInterval)), self.size])

		for iteration in range(steps):
			
			#for idx, cell in enumerate(self.cells):
			#	cell.internal_update(self.cells[idx-1].E, self.cells[(idx+1)%self.size].E)
			
			'''self.Pfield[:] = ''' 
			self.update_Pfield(self.Efield, self.Pfield)

			'''self.Efield[:] = ''' 
			self.update_Efield(self.very_last_Efield, self.last_Efield, self.very_last_Pfield, self.last_Pfield, self.Pfield)
			
			self.set_border_Efield(self.Efield, iteration*self.parameters.dt)

			'''
			if(iteration % saveInterval == 0):
				for cell in range(len(self.Efield)):
					self.evolution[int(math.ceil(iteration/saveInterval)), cell] = self.Efield[cell]
			
			'''
			if(iteration % saveInterval == 0):
				for idx, cell in enumerate(self.cells):					
					self.evolution[int(math.ceil(iteration/saveInterval)), idx] = cell.P
			

	def print_evolution(self):
		print(self.evolution)

	def plot_evolution(self):
		fig = plt.figure()
		plt.imshow(np.transpose(self.evolution))
		plt.show()


size = 100

arbitrary_material = internal.material_constants(
					k_aab0 = 0, k_baa0 = 0, k_abc0 = 0, k_cab0 = 0, 
					alpha_aab = 0, alpha_baa = 0, alpha_abc = 0, alpha_cab = 0, beta_aab = 0, beta_abc = 0,
					epsilon_ra = 1, epsilon_rb = 1, epsilon_rc = 1,
					tau_a = 10, tau_b = 1, tau_c = 1)

arbitrary_environment = internal.simulation_parameters(epsilon_0 = 0.1, c = 1.0, dx = 0.5, dt = 0.4, wave_length = (size-0))



test = grid(arbitrary_material, arbitrary_environment, size)
david_gengenbachs_penis = 1000
test.run(david_gengenbachs_penis, 1)

test.print_evolution()
test.plot_evolution()
