import numpy as np
import math, random
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import internal


class grid:
	def __init__(self, constants, parameters, size):
		self.size = size
		self.parameters = parameters
		self.constants = constants

		self.current_step = 0

		self.save_steps = 0

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
			cell.E =  math.sin(idx*self.parameters.dx * 2*math.pi / self.parameters.wave_length) #np.array([0.0, math.sin(idx*self.parameters.dx * 2*math.pi / self.parameters.wave_length), 0.0])
			
		self.iterate(set_init)

	def update_Efield(self):
		field = math.sin(5*self.current_step*self.parameters.dx * 2*math.pi / self.parameters.wave_length)
		print(field)

		self.cells[ 0, 0].E = field
		self.cells[-1, 0].E = field
		
		# calculate difference of E-Field between Borders
		delta_E = (self.cells[0, 0].E - self.cells[-1, 0].E)/self.size

		def set_E(cell, idx, idy):
			cell.E = self.cells[0, 0].E + idx * delta_E

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


	def make_future_happen(self):
		self.current_step += 1
		self.cells = np.copy(self.future_cells)

	def evolve(self):

		self.set_inital_Efield()
		
		self.iterate(lambda cell, idx, idy: cell.internal_update() )

		def set_future(cell, idx, idy):
			self.future_cells[idx, idy] = cell

		self.iterate(set_future)	
			
		self.update_Efield()

		self.apply_diffusion()

		self.make_future_happen()


			

	def print_state(self):
		def print_it(cell, idx, idy):
			print(cell.E, cell.P, cell.n_a, cell.n_b, cell.n_c)

		self.iterate(print_it)

	def get_n_x(self):
		n_x = np.empty([size, size, 3])
		def extract_n(cell, idx, idy):
			n_x[idx, idy, 0] = cell.n_a
			n_x[idx, idy, 1] = 2*cell.n_b
			n_x[idx, idy, 2] = 3*cell.n_c
		self.iterate(extract_n)
		return n_x

	def get_E(self):
		E = np.empty([size, size])
		def extract_E(cell, idx, idy):
			E[idx, idy] = cell.E
			
		self.iterate(extract_E)
		return E



	def print_mass(self):
		self.mass = 0
		def add_mass(cell, idx, idy):
			self.mass += cell.n_a + 2*cell.n_b + 3*cell.n_c
		self.iterate(add_mass)
		print(self.mass)


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
size = 50
steps = 1
saveInterval = 1
wave_length = 50

arbitrary_material = internal.material_constants(
					k_aab0 = 0.5, k_baa0 = 0.1, k_abc0 = 0.2, k_cab0 = 0.1, 
					alpha_aab = 0.5, alpha_baa = -0.5, alpha_abc = 1.0, alpha_cab = -1.0, beta_aab = 1, beta_abc = 1,
					epsilon_ra = 1.0, epsilon_rb = 1.5, epsilon_rc = 2.0,
					tau_a = 1.0, tau_b = 1.0, tau_c = 1.0,
					D_a = 0.2, D_b = 0.1, D_c = 0.05)

arbitrary_environment = internal.simulation_parameters(epsilon_0 = 0.1, c = 1.0, dx = 0.5, dt = 0.4, wave_length = wave_length) #wave_length = (size-0)/1)



test = grid(arbitrary_material, arbitrary_environment, size)
#test.run(steps, saveInterval)


#alpha_aab = 0.0, alpha_baa = -0.0, alpha_abc = 0.0, alpha_cab = -0.0, beta_aab = 0, beta_abc = 0,
#alpha_aab = 0.5, alpha_baa = -0.5, alpha_abc = 1.0, alpha_cab = -1.0, beta_aab = 1, beta_abc = 1,
#alpha_aab = 10.0, alpha_baa = -10.0, alpha_abc = 10.0, alpha_cab = -10.0, beta_aab = 100.0, beta_abc = 10.0,


fig = plt.figure()
#fig, axes = plt.subplots(nrows=2, sharex=True, sharey=True)

def f():
	na = test.get_n_x()
	return na

def g():
	E = test.get_E()
	return E

im = plt.imshow(f(), animated=True)
#axes[0] = plt.imshow(f(), animated=True)
#axes[1] = plt.imshow(g(), animated=True)

def updatefig(*args):
	test.evolve()
	#test.print_mass()
	
	im.set_array(f())
	return im,

	#axes[0].set_array(f())
	#axes[1].set_array(g())
	#return axes,

ani = animation.FuncAnimation(fig, updatefig, interval=1, blit=False)
plt.show()
