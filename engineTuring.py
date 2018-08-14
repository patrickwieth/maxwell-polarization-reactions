import numpy as np
import math, random
import pickle


class grid2D:
	def __init__(self, material, materialc, environment, size):
		self.size = size
		self.parameters = environment
		self.material = material
		self.constants = materialc

		self.current_step = 0

		self.f = 0

		self.spawn_grid()


	def iterate(self, fn):
		# 2D version
		for idx, line in enumerate(self.cells):
			for idy, cell in enumerate(line):
				fn(cell, idx, idy)

	def spawn_grid(self):
		def init_function(x,y):
			rnd = 0.01
			self.P_a = 0.1
			self.P_b = 0.1
			self.P_c = 0.1
			self.n_a = 0.66 + random.uniform(-rnd, rnd)
			self.n_b = 0.10 + random.uniform(-rnd, rnd)
			self.n_c = 0.04 + random.uniform(-rnd, rnd)
			return self

		#2D
		self.cells = np.array([
			[self.material.cell(self.constants, self.parameters, init_function(x,y) ) for x in range(self.size)]
			for y in range(self.size)])


		# neighbors are needed for diffusion
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



	def apply_external_fields(self):
		strength = self.parameters.epsilon_0 * self.parameters.force_fields[0].get_strength(self.current_step*self.parameters.dt)
		self.f = strength

	def apply_diffusion(self):
		def laplace(cell, observable):
			result = sum([getattr(neighbor, observable) for neighbor in cell.neighbors])
			result -= len(cell.neighbors) * getattr(cell, observable)
			return result / self.parameters.dx**2

		concentrations = ['n_a', 'n_b', 'n_c']
		polarizations = ['P_a', 'P_b', 'P_c']
		diffusion_coefficients = np.array([getattr(self.constants, stuff) for stuff in ['D_na', 'D_nb', 'D_nc']])
		polarization_diffusion_coefficients = np.array([getattr(self.constants, stuff) for stuff in ['D_Pa', 'D_Pb', 'D_Pc']])

		def diffuse(cell, idx, idy):
			n_laplaces = np.array([laplace(cell, observable) for observable in concentrations])
			n_deltas = n_laplaces * diffusion_coefficients * self.parameters.dt
			p_laplaces = np.array([laplace(cell, observable) for observable in polarizations])
			p_deltas = p_laplaces * polarization_diffusion_coefficients * self.parameters.dt

			for obs, delta in zip(concentrations, n_deltas):
				setattr(self.future_cells[idx, idy], obs, getattr(cell, obs) + delta)

			for obs, delta in zip(polarizations, p_deltas):
				setattr(self.future_cells[idx, idy], obs, getattr(cell, obs) + delta)

		self.iterate(diffuse)

	def internal_update(self):
		self.iterate(lambda cell, idx, idy: cell.internal_update(self.f) ) 	#2D

		#def set_future(cell, idx, idy):		#2D
		#	self.future_cells[idx, idy] = cell

		#self.iterate(set_future)

	def make_future_happen(self):
		self.current_step += 1
		self.cells = np.copy(self.future_cells)

	def evolve(self):
		self.internal_update()
		self.apply_external_fields()
		self.apply_diffusion()
		self.make_future_happen()

	def save(self, filename, verbose):
		save_cells = self.material.get_data(self)

		with open(r""+filename, "wb") as output_file:
			pickle.dump(save_cells, output_file)

		if verbose : print("saved state to", filename)

	def load(self, filename, verbose):
		with open(r""+filename, "rb") as input_file:
			loaded = pickle.load(input_file)

		insert_loaded = self.material.load_function(loaded)

		self.iterate(insert_loaded)
		if verbose : print("loaded file", filename)


class analyze:
	def __init__(self, grid):
		self.chi = 0
		self.chi_accumulated = 0
		self.grid = grid
		self.counter = 0

	# assuming the electrical field is in x-direction
	def calculate_dielectric_response(self):
		self.chi = 0
		self.counter = 0

		def collect_chi(cell, idx, idy):
			self.chi += abs(cell.P[1])

		self.grid.iterate(collect_chi)
		self.chi = self.chi / self.grid.size**2

		self.chi_accumulated += self.chi


	def calculate_total_polarisation(self):
		self.P_t = 0.0

		def collect_P(cell, idx):
			self.P_t += cell.P

		self.grid.iterate(collect_P)

		print(self.P_t)

	def print_mass(self):
		self.mass = 0
		def add_mass(cell, idx):
			self.mass += cell.n_a + 2*cell.n_b + 3*cell.n_c

		self.grid.iterate(add_mass)

		print(self.mass)

	# 2D version
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
