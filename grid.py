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


	def spawn_grid(self):
		def init_function(x):
			self.P_a = 0
			self.P_b = 0
			self.P_c = 0
			self.n_a = 1
			self.n_b = 0
			self.n_c = 0
			return self

		self.cells = np.array([internal.cell(self.constants, self.parameters, init_function(x) ) for x in range(self.size)])


	def set_inital_Efield(self):

		for idx, cell in enumerate(self.cells):
			cell.E = math.sin(idx*self.parameters.dx * 2*math.pi / self.parameters.wave_length) #np.array([0.0, math.sin(idx*self.parameters.dx * 2*math.pi / self.parameters.wave_length), 0.0])
			cell.prev2_E = cell.prev_E
			cell.prev_E = cell.E

	def calc_first_step(self):
		#C2 = (self.parameters.c * self.parameters.dt / self.parameters.dx)**2 # courant number GET THIS THING EVERYWHERE
		# Special formula for first time step

		next_E = np.zeros(size)

		for i in range(1, size-1):
			next_E[i] = self.cells[i].E + 0.5*self.parameters.C2*(self.cells[i-1].E - 2*self.cells[i].E + self.cells[(i+1)%size].E)

		for idx, cell in enumerate(self.cells):
			cell.E = next_E[idx]

		self.cells[ 0].E = 0.0 #np.zeros(3); 
		self.cells[-1].E = 0.0 #np.zeros(3); 


	def update_Efield(self):
		# make history
		for cell in self.cells:
			cell.prev2_E = cell.prev_E
			cell.prev_E = cell.E

		for idx, cell in enumerate(self.cells):
			cell.E = 0.0

			# calculate E for all but border cells
			if(idx > 0 and idx < size-1):

				left_E = self.cells[idx-1].prev_E
				right_E = self.cells[(idx+1)%self.size].prev_E

				P_tt = (cell.prev2_P - 2*cell.prev_P + cell.P)
				E_xx = (left_E - 2*cell.prev_E + right_E)
				t_derivative = self.parameters.C2 * E_xx - 1/self.parameters.epsilon_0 * P_tt

				cell.E = t_derivative - cell.prev2_E + 2*cell.prev_E


	def set_border_Efield(self, time):
		self.cells[ 0].E = 0.0 #np.zeros(3)
		self.cells[-1].E = 0.0 #np.zeros(3)

	def run(self, steps, saveInterval):

		self.set_inital_Efield()

		self.Efield = self.calc_first_step()

		self.evolution = np.empty([int(math.ceil(steps/saveInterval)), self.size])

		for iteration in range(steps):
			
			for idx, cell in enumerate(self.cells):
				cell.internal_update(self.cells[idx-1].E, self.cells[(idx+1)%self.size].E)
			

			self.update_Efield()
			
			#self.set_border_Efield(iteration*self.parameters.dt)

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
steps = 1000
saveInterval = 1

arbitrary_material = internal.material_constants(
					k_aab0 = 1.0, k_baa0 = 1.0, k_abc0 = 0.0, k_cab0 = 0.0, 
					alpha_aab = 0, alpha_baa = 0, alpha_abc = 0, alpha_cab = 0, beta_aab = 0, beta_abc = 0,
					epsilon_ra = 1.0, epsilon_rb = 1.0, epsilon_rc = 1.0,
					tau_a = 1, tau_b = 1, tau_c = 1)

arbitrary_environment = internal.simulation_parameters(epsilon_0 = 0.1, c = 1.0, dx = 0.5, dt = 0.4, wave_length = (size-0))



test = grid(arbitrary_material, arbitrary_environment, size)
test.run(steps, saveInterval)

test.print_evolution()
test.plot_evolution()
