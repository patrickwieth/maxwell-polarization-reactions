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
			cell.E =  math.sin(idx*self.parameters.dx * 2*math.pi / self.parameters.wave_length) #np.array([0.0, math.sin(idx*self.parameters.dx * 2*math.pi / self.parameters.wave_length), 0.0])
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


	def update_Efield(self, time):
		# make history
		for cell in self.cells:
			cell.prev2_E = cell.prev_E
			cell.prev_E = cell.E

		for idx, cell in enumerate(self.cells):
			#cell.E = 0.0
			cell.E = math.sin(time*self.parameters.dx * 2*math.pi / self.parameters.wave_length)

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

		self.E_evolution = np.empty([int(math.ceil(steps/saveInterval)), self.size])
		self.P_evolution = np.empty([int(math.ceil(steps/saveInterval)), self.size])
		self.na_evolution = np.empty([int(math.ceil(steps/saveInterval)), self.size])
		self.nb_evolution = np.empty([int(math.ceil(steps/saveInterval)), self.size])
		self.nc_evolution = np.empty([int(math.ceil(steps/saveInterval)), self.size])

		for iteration in range(steps):
			
			for idx, cell in enumerate(self.cells):
				cell.internal_update(self.cells[idx-1].E, self.cells[(idx+1)%self.size].E)
			

			self.update_Efield(iteration)
			
			#self.set_border_Efield(iteration*self.parameters.dt)

			if(iteration % saveInterval == 0):
				for idx, cell in enumerate(self.cells):					
					self.E_evolution[int(math.ceil(iteration/saveInterval)), idx] = cell.E
					self.P_evolution[int(math.ceil(iteration/saveInterval)), idx] = cell.P
					self.na_evolution[int(math.ceil(iteration/saveInterval)), idx] = cell.n_a
					self.nb_evolution[int(math.ceil(iteration/saveInterval)), idx] = cell.n_b
					self.nc_evolution[int(math.ceil(iteration/saveInterval)), idx] = cell.n_c

		self.E_evolution = np.transpose(self.E_evolution)
		self.P_evolution = np.transpose(self.P_evolution)
		self.na_evolution = np.transpose(self.na_evolution)
		self.nb_evolution = np.transpose(self.nb_evolution)
		self.nc_evolution = np.transpose(self.nc_evolution)
			
	def print_evolution(self):
		print(self.E_evolution)

	def plot_evolution(self, name):
		#fig = plt.figure()

		fig, axes = plt.subplots(nrows=5, sharex=True, sharey=True)

		axes[0].set_title('Evolution of E, P, Na, Nb and Nc:')
		axes[0].imshow(self.E_evolution, interpolation=None)
		axes[1].imshow(self.P_evolution, interpolation=None)
		axes[2].imshow(self.na_evolution, interpolation=None)
		axes[3].imshow(self.nb_evolution, interpolation=None)
		axes[4].imshow(self.nc_evolution, interpolation=None)


		fig.subplots_adjust(hspace=0)
		#plt.setp([a.get_xticklabels() for a in axes[:-1]], visible=False)
		[a.set_adjustable('box-forced') for a in axes]


		plt.show()
		#plt.savefig(str(name)+'.png')


#for wave_length in range(100, 200, 10):
size = 100
steps = 6000
saveInterval = 1
wave_length = 300

arbitrary_material = internal.material_constants(
					k_aab0 = 1.0, k_baa0 = 0.8, k_abc0 = 1.0, k_cab0 = 1.0, 
					alpha_aab = 10.0, alpha_baa = -10.0, alpha_abc = 10.0, alpha_cab = -10.0, beta_aab = 100.0, beta_abc = 10.0,
					epsilon_ra = 1.0, epsilon_rb = 1.5, epsilon_rc = 2.0,
					tau_a = 1.0, tau_b = 3.0, tau_c = 10.0)

arbitrary_environment = internal.simulation_parameters(epsilon_0 = 0.1, c = 1.0, dx = 0.5, dt = 0.4, wave_length = wave_length) #wave_length = (size-0)/1)



test = grid(arbitrary_material, arbitrary_environment, size)
test.run(steps, saveInterval)

test.print_evolution()
test.plot_evolution(wave_length)

#alpha_aab = 0.5, alpha_baa = -0.5, alpha_abc = 1.0, alpha_cab = -1.0, beta_aab = 1, beta_abc = 1,
#alpha_aab = 0.0, alpha_baa = -0.0, alpha_abc = 0.0, alpha_cab = -0.0, beta_aab = 0.0, beta_abc = 0.0,