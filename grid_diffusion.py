import numpy as np
import math, random
import pickle

import internal

# TODO: relax Diffusion, cooperate Polarization

class grid:
	def __init__(self, constants, parameters, size, external_fields):
		self.size = size
		self.parameters = parameters
		self.constants = constants
		self.external_fields = external_fields
		self.current_step = 0
		self.cells = np.empty([size, size])
		self.dq = 1
		self.f = 0

		self.spawn_grid()

		#self.set_inital_Efield()

	def iterate(self, fn):	
		'''
		# 2D version
		for idx, line in enumerate(self.cells):
			for idy, cell in enumerate(line):
				fn(cell, idx, idy)
		'''
		# 1D version
		for idx, cell in enumerate(self.cells):
			fn(cell, idx)


	def spawn_grid(self):
		def init_function(x):
			self.P_a = 0
			self.P_b = 0
			self.P_c = 0
			self.n_a = 0.3 + random.uniform(-0.1, 0.1)
			self.n_b = 0.2 + random.uniform(-0.1, 0.1)
			self.n_c = 0.1 + random.uniform(-0.1, 0.1)
			return self

		''' # 2D version:
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
		'''

		# 1D version:
		self.cells = np.array([internal.cell(self.constants, self.parameters, init_function(x) ) for x in range(self.size)])

		# neighbors are defined here, this determines topology
		def attach_neighbor(cell, idx):
			if(idx > 0): 
				cell.neighbors.append(self.cells[idx-1])
			if(idx < self.size-1): 
				cell.neighbors.append(self.cells[idx+1])

		self.iterate(attach_neighbor) 

		self.future_cells = np.copy(self.cells)


	def apply_external_fields(self):

		self.f = self.parameters.epsilon_0 * self.external_fields[0].get_strength(self.current_step*self.parameters.dt) - self.cells[0].P

		'''
		def set_E(cell, idx, idy):
			self.future_cells[idx, idy].E = self.external_fields[0].get_strength(self.current_step*self.parameters.dt)
			cell.E = self.external_fields[0].get_strength(self.current_step*self.parameters.dt)

		self.iterate(set_E)
		'''

	def apply_diffusion(self):
		def laplace(cell, observable):
			result = sum([getattr(neighbor, observable) for neighbor in cell.neighbors])
			result -= len(cell.neighbors) * getattr(cell, observable)
			return result / self.parameters.dx**2

		concentrations = ['n_a', 'n_b', 'n_c']
		polarizations = ['P_a', 'P_b', 'P_c']
		diffusion_coefficients = np.array([getattr(self.constants, stuff) for stuff in ['D_na', 'D_nb', 'D_nc']])
		polarization_diffusion_coefficients = np.array([getattr(self.constants, stuff) for stuff in ['D_Pa', 'D_Pb', 'D_Pc']])

		def diffuse(cell, idx):
			n_laplaces = np.array([laplace(cell, observable) for observable in concentrations])
			n_deltas = n_laplaces * diffusion_coefficients * self.parameters.dt
			p_laplaces = np.array([laplace(cell, observable) for observable in polarizations])
			p_deltas = p_laplaces * polarization_diffusion_coefficients * self.parameters.dt
			
			for obs, delta in zip(concentrations, n_deltas):
				setattr(self.future_cells[idx], obs, getattr(cell, obs) + delta)

			for obs, delta in zip(polarizations, p_deltas):
				setattr(self.future_cells[idx], obs, getattr(cell, obs) + delta)

		self.iterate(diffuse)

	def internal_update(self):
		#self.iterate(lambda cell, idx, idy: cell.internal_update(self.q) )

		self.iterate(lambda cell, idx: cell.internal_update(self.f) )

		#def set_future(cell, idx, idy):
		def set_future(cell, idx):
			self.future_cells[idx] = cell

		self.iterate(set_future)

	def make_future_happen(self):
		self.current_step += 1
		self.cells = np.copy(self.future_cells)

	def evolve(self):
		self.internal_update()
		self.apply_external_fields()
		#self.apply_diffusion()
		self.make_future_happen()
			
	# saves only the state of observables in the cells
	def save(self, filename):
		save_cells = np.array([[
				self.cells[x].P_a,
				self.cells[x].P_b,
				self.cells[x].P_c,
				self.cells[x].n_a,
				self.cells[x].n_b,
				self.cells[x].n_c]
				for x in range(self.size)]) 
			

		with open(r""+filename, "wb") as output_file:
			pickle.dump(save_cells, output_file)

		#print("saved state to", filename)

	# a valid grid must be constructed when loading
	def load(self, filename, verbose):
		with open(r""+filename, "rb") as input_file:
			loaded = pickle.load(input_file)

		def insert_loaded(cell, idx):
			cell.P_a = loaded[idx, 0]
			cell.P_b = loaded[idx, 1]
			cell.P_c = loaded[idx, 2]
			cell.n_a = loaded[idx, 3]
			cell.n_b = loaded[idx, 4]
			cell.n_c = loaded[idx, 5]

		self.iterate(insert_loaded)
		if verbose : print("loaded file", filename)


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

		self.chi_accumulated = 0
		
		self.grid = grid

		self.counter = 0

	def calculate_dielectric_response(self):
		self.chi = 0
		self.counter = 0

		def collect_chi(cell, idx):
			self.chi += abs(cell.P)

		self.grid.iterate(collect_chi)
		self.chi = self.chi / self.grid.size

		self.chi_accumulated += self.chi


	def calculate_total_polarisation(self):
		self.P_t = 0.0

		def collect_P(cell, idx):
			if cell.E != 0:
				self.P_t += cell.P


		self.grid.iterate(collect_P)

		print(self.P_t)

	def print_mass(self):
		self.mass = 0
		def add_mass(cell, idx):
			self.mass += cell.n_a + 2*cell.n_b + 3*cell.n_c

		self.grid.iterate(add_mass)

		print(self.mass)

	''' 2D version
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
	'''

	def get_observables(self):
		obs = np.empty([5, self.grid.size])
		def extract_obs(cell, idx):
			obs[0, idx] = cell.n_a
			obs[1, idx] = 2*cell.n_b
			obs[2, idx] = 3*cell.n_c
			obs[3, idx] = cell.E
			obs[4, idx] = cell.P

		self.grid.iterate(extract_obs)
		return obs