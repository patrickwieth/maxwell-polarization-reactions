import numpy as np
import math, random
import pickle

import material2 as material


# TODO: relax Diffusion, cooperate Polarization

class grid:
	def __init__(self, materialc, environment, size):
		self.size = size
		self.parameters = environment
		self.constants = materialc
		
		self.current_step = 0
		#self.cells = np.empty([size, size])
		self.dq = 1
		self.f = 0

		self.spawn_grid()


	def iterate(self, fn):	
		# 1D version
		for idx, cell in enumerate(self.cells):
			fn(cell, idx)


	def spawn_grid(self):
		def init_function(x):
			self.P_a = 0
			self.P_b = 0
			self.P_c = 0
			self.n_a = 0.6 #+ random.uniform(-0.1, 0.1)
			self.n_b = 0.2 #+ random.uniform(-0.1, 0.1)
			self.n_c = 0.1 #+ random.uniform(-0.1, 0.1)
			return self

		# 1D version:
		self.cells = np.array([material.cell(self.constants, self.parameters, init_function(x) ) for x in range(self.size)])

		# neighbors are defined here, this determines topology
		def attach_neighbor(cell, idx):
			if(idx > 0): 
				cell.neighbors.append(self.cells[idx-1])
			if(idx < self.size-1): 
				cell.neighbors.append(self.cells[idx+1])

		self.iterate(attach_neighbor) 

		self.future_cells = np.copy(self.cells)


	def apply_external_fields(self):

		self.f = self.parameters.epsilon_0 * self.parameters.force_fields[0].get_strength(self.current_step*self.parameters.dt) - self.cells[0].P


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
			
	def save(self, filename, verbose):
		save_cells = material.get_data(self)

		with open(r""+filename, "wb") as output_file:
			pickle.dump(save_cells, output_file)

		if verbose : print("saved state to", filename)

	def load(self, filename, verbose):
		with open(r""+filename, "rb") as input_file:
			loaded = pickle.load(input_file)

		insert_loaded = material.load_function(loaded)

		self.iterate(insert_loaded)
		if verbose : print("loaded file", filename)



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