import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import grid_diffusion, internal

import argparse

parser = argparse.ArgumentParser(description='Run a simulation of reaction-diffusion-polarization equations')

parser.add_argument("--load", help="load a checkpoint")
args = parser.parse_args()


size = 50


super_simple_material = internal.material_constants(
					k_aab0 = 0, k_baa0 = 0, k_abc0 = 0, k_cab0 = 0, k_caaa0 = 0, k_aaac0 = 0,
					alpha_aab = 0.0, alpha_baa = 0.0, alpha_abc = 0.0, alpha_cab = -0.0, beta_aab = 0, beta_abc = 0,
					epsilon_ra = 1.0, epsilon_rb = 1.0, epsilon_rc = 1.0, dP = 0.0,
					tau_a = 10.0, tau_b = 15.0, tau_c = 20.0,
					D_na = 0.1, D_nb = 0.05, D_nc = 0.01, D_Pa = 1, D_Pb = 1, D_Pc = 1)


arbitrary_material = internal.material_constants(
					k_aab0 = 0.01, k_baa0 = 0.01, k_abc0 = 0.01, k_cab0 = 0.005, k_caaa0 = 0.01, k_aaac0 = 0.01,
					alpha_aab = 0.0001, alpha_baa = -0.0005, alpha_abc = 0.0001, alpha_cab = -0.0001, beta_aab = 0.1, beta_abc = 0.1,
					epsilon_ra = 1.0, epsilon_rb = 2.0, epsilon_rc = 3.0, dP = 2.0,
					tau_a = 10.0, tau_b = 15.0, tau_c = 20.0,
					D_na = 0.00028, D_nb = 0.05, D_nc = 0.05, D_Pa = 1, D_Pb = 1, D_Pc = 1)


arbitrary_environment = internal.simulation_parameters(epsilon_0 = 1, c = 1.0, dx = 2.0/size, dt = 0.49*(2.0/size)**2, wave_length = 50) #wave_length = (size-0)/1)

E_field = grid_diffusion.external_field(strength=0.1, frequency=0.1, wave_length=50)

#simulation = grid_diffusion.grid(super_simple_material, arbitrary_environment, size, [E_field])
simulation = grid_diffusion.grid(arbitrary_material, arbitrary_environment, size, [E_field])

observe = grid_diffusion.analyze(simulation)

if not args.load == None:
	simulation.load(args.load)


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
	if simulation.current_step % 1000 == 0:
		simulation.save("checkpoint-"+str(int(simulation.current_step/1000)))
		

	#observe.calculate_dielectric_response(simulation.current_step)
	#observe.calculate_total_polarisation()
	
	#print(simulation.current_step)
	#observe.print_mass()
	

	[im[n].set_array(get_obs(n)) for n in range(subplots)]

	return im

ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=True)

plt.show()
