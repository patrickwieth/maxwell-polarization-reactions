import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import grid_diffusion, internal

import argparse

parser = argparse.ArgumentParser(description='Run a simulation of reaction-diffusion-polarization equations')

parser.add_argument("--load", help="load a checkpoint")
args = parser.parse_args()


size = 1

super_simple_material = internal.material_constants(
					k_aab0 = 0, k_baa0 = 0, k_abc0 = 0, k_cab0 = 0, k_caaa0 = 0, k_aaac0 = 0,
					alpha_aab = 0.0, alpha_baa = 0.0, alpha_abc = 0.0, alpha_cab = -0.0, alpha_aaac = 0, alpha_caaa = 0, 
					epsilon_ra = 1.0, epsilon_rb = 1.0, epsilon_rc = 1.0, dP = 0.0,
					tau_a = 10.0, tau_b = 15.0, tau_c = 20.0,
					D_na = 0.1, D_nb = 0.05, D_nc = 0.01, D_Pa = 0.1, D_Pb = 1, D_Pc = 1)


arbitrary_material = internal.material_constants(
					k_aab0 = 0.01, k_baa0 = 0.02, k_abc0 = 0.05, k_cab0 = 0.05, k_caaa0 = 0.04, k_aaac0 = 0.0001,
					alpha_aab = 0.001, alpha_baa = -0.001, alpha_abc = 0.001, alpha_cab = -0.001, alpha_aaac = 0.01, alpha_caaa = -0.001, 
					epsilon_ra = 1.0, epsilon_rb = 2.0, epsilon_rc = 3.0, dP = 0.0,
					tau_a = 1.0, tau_b = 1.50, tau_c = 2.0,
					D_na = 0.00028, D_nb = 0.05, D_nc = 0.05, D_Pa = 0.00028, D_Pb = 0.05, D_Pc = 0.05)

no_species_C_material = internal.material_constants(
					k_aab0 = 0.00, k_baa0 = 0.00, k_abc0 = 0.00, k_cab0 = 0.00, k_caaa0 = 0.00, k_aaac0 = 0.0000,
					alpha_aab = 0.000, alpha_baa = -0.000, alpha_abc = 0.000, alpha_cab = -0.000, alpha_aaac = 0.00, alpha_caaa = -0.000, 
					epsilon_ra = 1.0, epsilon_rb = 1.0, epsilon_rc = 1.0, dP = 0.0,
					tau_a = 1.0, tau_b = 1.00, tau_c = 1.0,
					D_na = 0.0, D_nb = 0.0, D_nc = 0.0, D_Pa = 0.0, D_Pb = 0.0, D_Pc = 0.0)

#for k, v in vars(no_species_C_material).items():
#	print(k, v)

arbitrary_environment = internal.simulation_parameters(epsilon_0 = 1, c = 1.0, dx = 0.1, dt = 0.01)

E_field = grid_diffusion.external_field(strength=0.0, frequency=0.0, wave_length=50)


simulation = grid_diffusion.grid(no_species_C_material, arbitrary_environment, size, [E_field])


steps = int(10000)


# load start state or equilibrate at beginning
if not args.load == None:
	simulation.load(args.load, False)
else:
	for n in range(steps):
		simulation.evolve()

	simulation.save("equilibrium")


for n in range(500):
	

	freq = 0.01 * 1.1**n
	period_t = 2 * 3.1415 / freq 
	
	arbitrary_environment = internal.simulation_parameters(epsilon_0 = 1, c = 1.0, dx = 0.1, dt = period_t/1000)
	E_field = grid_diffusion.external_field(strength=0.01, frequency=freq, wave_length=50)

	simulation = grid_diffusion.grid(no_species_C_material, arbitrary_environment, size, [E_field])
	simulation.load("equilibrium", False)

	observe = grid_diffusion.analyze(simulation)

	for n in range(steps):
		simulation.evolve()
		observe.calculate_dielectric_response()

	print(observe.chi_accumulated/steps)

'''

if not args.load == None:
	simulation.load(args.load, False)

observe = grid_diffusion.analyze(simulation)

def get_obs(n):
	return observe.get_observables()[n]

fig, ax = plt.subplots()

x = np.arange(0, size)


subplots = 5
#lines = [ax.plot(x, get_obs(n), animated=True) for n in range(subplots)]
line1, = ax.plot(x, get_obs(0), animated=True)
line2, = ax.plot(x, get_obs(1), animated=True)
line3, = ax.plot(x, get_obs(2), animated=True)
line4, = ax.plot(x, get_obs(3), animated=True)
line5, = ax.plot(x, get_obs(4), animated=True)


def animate(i):
	[simulation.evolve() for n in range(10)]

	if simulation.current_step % 10000 == 0:
		simulation.save("checkpoint-"+str(int(simulation.current_step/10000)))

	observe.calculate_total_polarisation()

	#[lines[n].set_ydata(get_obs(n)) for n in range(subplots)]
	#return lines

	line1.set_ydata(get_obs(0))
	line2.set_ydata(get_obs(1))
	line3.set_ydata(get_obs(2))
	line4.set_ydata(get_obs(3))
	line5.set_ydata(get_obs(4))
	return [line1, line2, line3, line4, line5]


# Init only required for blitting to give a clean slate.
def init():
    line1.set_ydata(get_obs(0))
    line2.set_ydata(get_obs(3))
    line3.set_ydata(get_obs(3))
    line4.set_ydata(get_obs(3))
    line5.set_ydata(get_obs(4))
    return [line1, line2, line3, line4, line5]

    #[lines[n].set_ydata(get_obs(n)) for n in range(subplots)]
	#return lines

ani = animation.FuncAnimation(fig, animate, init_func=init, interval=25, blit=True)
plt.show()
'''


''' 2D
subplots = 5

fig, axes, = plt.subplots(1, subplots, sharex='col', sharey='row')

[axes[n].set_adjustable('box-forced') for n in range(subplots)]
[axes[n].set_aspect('equal', 'datalim') for n in range(subplots)]


def get_obs(n):
	#res = observe.get_observables()[n]
	res = np.atleast_2d(observe.get_observables()[n])
	return res


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
'''