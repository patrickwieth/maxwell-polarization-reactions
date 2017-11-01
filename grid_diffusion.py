import numpy as np
import math, random

import internal


class grid:
	def __init__(self, constants, parameters, size):
		self.size = size
		self.parameters = parameters
		self.constants = constants
		self.current_step = 0
		self.cells = np.empty([size, size])

		self.spawn_grid()

	def iterate(self, fn):
		for idx, line in enumerate(self.cells):
			for idy, cell in enumerate(line):
				fn(cell, idx, idy)

	def spawn_grid(self):
		def init_function(x):
			self.P_a = 0
			self.P_b = 0
			self.P_c = 0
			self.n_a = random.uniform(0, 1)
			self.n_b = 0
			self.n_c = 0
			return self

		self.cells = np.array([
			[internal.cell(self.constants, self.parameters, init_function(x) ) for x in range(self.size)] 
			for y in range(self.size)])

		# neighbors are defined here, this determines topology
		def attach_neighbor(cell, idx, idy):
			if(idx > 0): 
				cell.neighbors.append(self.cells[idx-1, idy])
			if(idy > 0):
				cell.neighbors.append(self.cells[idx, idy-1])
						
			if(idx < self.size-1): 
				cell.neighbors.append(self.cells[idx+1, idy])
			if(idy < self.size-1): 
				cell.neighbors.append(self.cells[idx, idy+1])

		self.iterate(attach_neighbor) 

		self.future_cells = np.copy(self.cells)


	def set_inital_Efield(self):
		def set_init(cell, idx, idy):
			cell.E = 0.0 #math.sin(idx*self.parameters.dx * 2*math.pi / self.parameters.wave_length) #np.array([0.0, math.sin(idx*self.parameters.dx * 2*math.pi / self.parameters.wave_length), 0.0])
			
		self.iterate(set_init)

	def update_Efield(self):
		field = math.sin(self.current_step*self.parameters.dx * 2*math.pi / self.parameters.wave_length)
		#field = 0.0
		#print(field)

		self.cells[ 0, 0].E = field
		self.cells[-1, 0].E = field
		
		# calculate difference of E-Field between Borders
		delta_E = (self.cells[0, 0].E - self.cells[-1, 0].E)/self.size

		def set_E(cell, idx, idy):
			cell.E = self.cells[0, 0].E #+ idx * delta_E

		self.iterate(set_E)


	def apply_diffusion(self):
		def laplace(cell, observable):
			result = sum([getattr(neighbor, observable) for neighbor in cell.neighbors])
			result -= len(cell.neighbors) * getattr(cell, observable)
			return result / self.parameters.dx**2

		observables = ['n_a', 'n_b', 'n_c']
		diffusion_coefficients = np.array([getattr(self.constants, stuff) for stuff in ['D_a', 'D_b', 'D_c']])

		def diffuse(cell, idx, idy):
			laplaces = np.array([laplace(cell, observable) for observable in observables])
			deltas = laplaces * diffusion_coefficients * self.parameters.dt
			
			for obs, delta in zip(observables, deltas):
				setattr(self.future_cells[idx, idy], obs, getattr(cell, obs) + delta)

		self.iterate(diffuse)

	def internal_update(self):
		self.iterate(lambda cell, idx, idy: cell.internal_update() )

		def set_future(cell, idx, idy):
			self.future_cells[idx, idy] = cell

		self.iterate(set_future)

	def make_future_happen(self):
		self.current_step += 1
		self.cells = np.copy(self.future_cells)

	def evolve(self):
		self.set_inital_Efield()
		self.internal_update()
		self.update_Efield()
		self.apply_diffusion()
		self.make_future_happen()
			

	def print_state(self):
		def print_it(cell, idx, idy):
			print(cell.E, cell.P, cell.n_a, cell.n_b, cell.n_c)

		self.iterate(print_it)

	def get_observables(self):
		obs = np.empty([5, self.size, self.size])
		def extract_obs(cell, idx, idy):
			obs[0, idx, idy] = cell.n_a
			obs[1, idx, idy] = 2*cell.n_b
			obs[2, idx, idy] = 3*cell.n_c
			obs[3, idx, idy] = cell.E
			obs[4, idx, idy] = cell.P
		self.iterate(extract_obs)
		return obs

	
	def print_mass(self):
		self.mass = 0
		def add_mass(cell, idx, idy):
			self.mass += cell.n_a + 2*cell.n_b + 3*cell.n_c
		self.iterate(add_mass)
		print(self.mass)

class analyze:
	def __init__(self, grid):
		self.chi = 0
		self.P_T = 0
		self.grid = grid

	def calculate_dielectric_response(self):
		self.chi = 0

		def collect_chi(cell, idx, idy):
			if cell.E != 0:
				self.chi += cell.P / cell.E / self.grid.parameters.epsilon_0

		self.grid.iterate(collect_chi)


		print(self.chi)

	def calculate_total_polarisation():
		self.P_t = 0
