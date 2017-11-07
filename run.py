import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import grid_diffusion, internal

#for wave_length in range(100, 200, 10):
size = 50
steps = 1
saveInterval = 1
wave_length = 50

super_simple_material = internal.material_constants(
					k_aab0 = 0, k_baa0 = 0, k_abc0 = 0, k_cab0 = 0, k_caaa0 = 0,
					alpha_aab = 0.0, alpha_baa = 0.0, alpha_abc = 0.0, alpha_cab = -0.0, beta_aab = 0, beta_abc = 0,
					epsilon_ra = 1.0, epsilon_rb = 1.0, epsilon_rc = 1.0, dP = 0.0,
					tau_a = 10.0, tau_b = 15.0, tau_c = 20.0,
					D_a = 0.2, D_b = 0.1, D_c = 0.01)


arbitrary_material = internal.material_constants(
					k_aab0 = 0.1, k_baa0 = 0.02, k_abc0 = 0.1, k_cab0 = 0.01, k_caaa0 = 0.1,
					alpha_aab = 0.5, alpha_baa = -0.005, alpha_abc = 0.001, alpha_cab = -0.001, beta_aab = 1, beta_abc = 1,
					epsilon_ra = 1.0, epsilon_rb = 2.0, epsilon_rc = 3.0, dP = 1.0,
					tau_a = 10.0, tau_b = 15.0, tau_c = 20.0,
					D_a = 0.2, D_b = 0.1, D_c = 0.01)


arbitrary_environment = internal.simulation_parameters(epsilon_0 = 1, c = 1.0, dx = 0.5, dt = 0.4, wave_length = wave_length) #wave_length = (size-0)/1)

E_field = grid_diffusion.external_field(strength=0, frequency=0.0, wave_length=50)

#simulation = grid_diffusion.grid(super_simple_material, arbitrary_environment, size, [E_field])
simulation = grid_diffusion.grid(arbitrary_material, arbitrary_environment, size, [E_field])

observe = grid_diffusion.analyze(simulation)


simulation.load("bla")



subplots = 5

fig, axes, = plt.subplots(1, subplots, sharex='col', sharey='row')

[axes[n].set_adjustable('box-forced') for n in range(subplots)]
[axes[n].set_aspect('equal', 'datalim') for n in range(subplots)]


def get_obs(n):
	return observe.get_observables()[n]


im = [axes[n].imshow(get_obs(n), cmap=plt.get_cmap('viridis'), animated=True, interpolation=None) for n in range(subplots)]
[fig.colorbar(im[n],ax=axes[n]) for n in range(subplots)]





def updatefig(*args):
	simulation.evolve()
	if simulation.current_step == 10:
		simulation.save("bla")

	#observe.calculate_dielectric_response(simulation.current_step)
	observe.calculate_total_polarisation()
	
	#print(simulation.current_step)
	#observe.print_mass()
	

	[im[n].set_array(get_obs(n)) for n in range(subplots)]

	return im

ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=True)

plt.show()
