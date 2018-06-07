import numpy as np
import math, random
import pickle

import material3D as material

# TODO: relax Diffusion, cooperate Polarization

class grid3D:
	def __init__(self, materialc, environment, size):
		self.size = size
		self.parameters = environment
		self.constants = materialc
		self.current_step = 0
		self.f = np.array([0.0, 0.0, 0.0])

		self.spawn_grid()


	def iterate(self, fn):					# 3D
		for idx, line in enumerate(self.cells):
			for idy, row in enumerate(line):
				for idz, cell in enumerate(row):
					fn(cell, idx, idy, idz)

	def spawn_grid(self):
		def init_function(x,y,z):
			length = 0.003
			angle = random.uniform(-math.pi, math.pi)
			self.P_a = np.array([length*math.cos(angle), length*math.sin(angle), length*math.sin(angle)])
			self.P_b = np.array([length*math.cos(angle), length*math.sin(angle), length*math.sin(angle)])
			self.P_c = np.array([length*math.cos(angle), length*math.sin(angle), length*math.sin(angle)])
			self.n_a = 0.3 #+ random.uniform(-rnd, rnd)
			self.n_b = 0.2 #+ random.uniform(-rnd, rnd)
			self.n_c = 0.1 #+ random.uniform(-rnd, rnd)
			self.rot_a = np.array([0.0, 0.0, 0.0])
			self.rot_b = np.array([0.0, 0.0, 0.0])
			self.rot_c = np.array([0.0, 0.0, 0.0])
			return self

		#3D
		self.cells = np.array([[
					[material.cell(self.constants, self.parameters, init_function(x,y,z) ) 
						for x in range(self.size)] 
							for y in range(self.size)]
								for z in range(self.size)])

		'''
		# neighbors are needed for diffusion, in 3D neighbors were never tested!
		def attach_neighbor(cell, idx, idy, idz):
			if(idx > 0): 
				cell.neighbors.append(self.cells[idx-1, idy, idz])
			if(idy > 0):
				cell.neighbors.append(self.cells[idx, idy-1, idz])
			if(idz > 0):
				cell.neighbors.append(self.cells[idx, idy, idz-1])
									
			if(idx < self.size-1): 
				cell.neighbors.append(self.cells[idx+1, idy, idz])
			if(idy < self.size-1): 
				cell.neighbors.append(self.cells[idx, idy+1, idz])
			if(idz < self.size-1): 
				cell.neighbors.append(self.cells[idx, idy, idz+1])

		self.iterate(attach_neighbor) 
		'''

		self.future_cells = np.copy(self.cells)

	def apply_dipole_field(self):

		def vector_field(p, r):
			r_abs = np.linalg.norm(r)
			E_r = 1/self.parameters.epsilon_0 * (3*np.dot(p, r)*r/r_abs**5 - p/r_abs**3)
			return E_r

		def get_field(cell, idx, idy, idz):		#3D
			field = np.zeros(3)

			for y in range(self.size):
				for x in range(self.size):
					for z in range(self.size):
						field += vector_field(self.cells[x, y, z].P, math.sqrt(x**2 + y**2 + z**2))  if x != 0 or y != 0 or z != 0 else 0
				
			cell.dipolar_field = field

		self.iterate(get_field)


	def apply_external_fields(self):
		strength = self.parameters.epsilon_0 * self.parameters.force_fields[0].get_strength(self.current_step*self.parameters.dt)
		self.f = np.array([0.0, 0.0, strength])


	def internal_update(self):
		self.iterate(lambda cell, idx, idy, idz: cell.internal_update(self.f) ) 	#3D

		#def set_future(cell, idx, idy, idz):		#3D
		#	self.future_cells[idx, idy] = cell 	

		#self.iterate(set_future)

	def make_future_happen(self):
		self.current_step += 1
		self.cells = np.copy(self.future_cells)

	def evolve(self):
		self.internal_update()
		self.apply_external_fields()
		self.apply_dipole_field()
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

	# 2D version
	def get_observables(self):
		obs = np.empty([6, self.grid.size, self.grid.size, self.grid.size])
		def extract_obs(cell, idx, idy, idz):
			obs[0, idx, idy, idz] = cell.n_a
			obs[1, idx, idy, idz] = 2*cell.n_b
			obs[2, idx, idy, idz] = 3*cell.n_c
			obs[3, idx, idy, idz] = cell.P[0]
			obs[4, idx, idy, idz] = cell.P[1]
			obs[5, idx, idy, idz] = cell.P[2]
		self.grid.iterate(extract_obs)
		return obs
	
