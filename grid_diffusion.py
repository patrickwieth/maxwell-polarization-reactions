import numpy as np
import math, random
import pickle

import internal


class grid:
	def __init__(self, constants, parameters, size, external_fields):
		self.size = size
		self.parameters = parameters
		self.constants = constants
		self.external_fields = external_fields
		self.current_step = 0
		self.cells = np.empty([size, size])

		self.spawn_grid()

		self.set_inital_Efield()

	def iterate(self, fn):
		for idx, line in enumerate(self.cells):
			for idy, cell in enumerate(line):
				fn(cell, idx, idy)

	def spawn_grid(self):
		def init_function(x):
			self.P_a = 0
			self.P_b = 0
			self.P_c = 0
			self.n_a = 1.000 #random.uniform(0, 1)
			self.n_b = 0.000
			self.n_c = 0.000
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

	def apply_external_fields(self):
		def set_E(cell, idx, idy):
			self.future_cells[idx, idy].E = self.external_fields[0].get_strength(self.current_step*self.parameters.dt)
			cell.E = self.external_fields[0].get_strength(self.current_step*self.parameters.dt)

		self.iterate(set_E)

	def apply_diffusion(self):
		def laplace(cell, observable):
			result = sum([getattr(neighbor, observable) for neighbor in cell.neighbors])
			result -= len(cell.neighbors) * getattr(cell, observable)
			return result / self.parameters.dx**2

		concentrations = ['n_a', 'n_b', 'n_c']
		polarizations = ['P_a', 'P_b', 'P_c']
		diffusion_coefficients = np.array([getattr(self.constants, stuff) for stuff in ['D_a', 'D_b', 'D_c']])

		def diffuse(cell, idx, idy):
			n_laplaces = np.array([laplace(cell, observable) for observable in concentrations])
			n_deltas = n_laplaces * diffusion_coefficients * self.parameters.dt
			p_laplaces = np.array([laplace(cell, observable) for observable in polarizations])
			p_deltas = p_laplaces * diffusion_coefficients * self.parameters.dt
			
			for obs, delta in zip(concentrations, n_deltas):
				setattr(self.future_cells[idx, idy], obs, getattr(cell, obs) + delta)

			for obs, delta in zip(polarizations, p_deltas):
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
		self.internal_update()
		self.apply_external_fields()
		self.apply_diffusion()
		self.make_future_happen()
			
	# saves only the state of observables in the cells
	def save(self, filename):
		save_cells = np.array([
				[[
				self.cells[x,y].P_a, 
				self.cells[x,y].P_b, 
				self.cells[x,y].P_c, 
				self.cells[x,y].n_a, 
				self.cells[x,y].n_b, 
				self.cells[x,y].n_c]
				for x in range(self.size)] 
			for y in range(self.size)])

		with open(r""+filename, "wb") as output_file:
			pickle.dump(save_cells, output_file)

	# a valid grid must be constructed when loading
	def load(self, filename):
		with open(r""+filename, "rb") as input_file:
			loaded = pickle.load(input_file)

		def insert_loaded(cell, idx, idy):
			cell.P_a = loaded[idx, idy, 0]
			cell.P_b = loaded[idx, idy, 1]
			cell.P_c = loaded[idx, idy, 2]
			cell.n_a = loaded[idx, idy, 3]
			cell.n_b = loaded[idx, idy, 4]
			cell.n_c = loaded[idx, idy, 5]

		self.iterate(insert_loaded)



class external_field:
	def __init__(self, strength, frequency, wave_length):
		self.strength = strength
		self.frequency = frequency
		self.wave_length = wave_length

	def get_strength(self, t):
		'''
		self.cells[ 0, 0].E = field
		self.cells[-1, 0].E = field
		
		# calculate difference of E-Field between Borders
		delta_E = (self.cells[0, 0].E - self.cells[-1, 0].E)/self.size

		def set_E(cell, idx, idy):
			cell.E = self.cells[0, 0].E #+ idx * delta_E
		'''

		if self.frequency > 0:
			return self.strength * math.sin(t * self.frequency * 2*math.pi) #/ self.wave_length)
		else:
			return self.strength



class analyze:
	def __init__(self, grid):
		self.chi = 0
		self.P_T = 0

		self.chi_accumulated = 0
		
		self.grid = grid

		self.counter = 0

	def calculate_dielectric_response(self, t):
		self.chi = 0
		self.counter = 0

		if(t > 0):
			def collect_chi(cell, idx, idy):
				if cell.E != 0:
					self.chi += cell.P / cell.E / self.grid.parameters.epsilon_0


			self.grid.iterate(collect_chi)
			self.chi = self.chi / (self.grid.size*self.grid.size)

			self.chi_accumulated += self.chi

			print(self.chi_accumulated / (t-0))


	def calculate_total_polarisation(self):
		self.P_t = 0.0

		def collect_P(cell, idx, idy):
			if cell.E != 0:
				self.P_t += cell.P

		self.grid.iterate(collect_P)

		print(self.P_t)

	def print_mass(self):
		self.mass = 0
		def add_mass(cell, idx, idy):
			self.mass += cell.n_a + 2*cell.n_b + 3*cell.n_c

		self.grid.iterate(add_mass)

		print(self.mass)

	def get_observables(self):
		obs = np.empty([5, self.grid.size, self.grid.size])
		def extract_obs(cell, idx, idy):
			obs[0, idx, idy] = cell.n_a
			obs[1, idx, idy] = 2*cell.n_b
			obs[2, idx, idy] = 3*cell.n_c
			obs[3, idx, idy] = cell.E
			obs[4, idx, idy] = cell.P
		self.grid.iterate(extract_obs)
		return obs