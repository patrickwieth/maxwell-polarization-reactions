import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import engineTuring as engine, environment
import material1 as material

import argparse

parser = argparse.ArgumentParser(description='Run a simulation of reaction-diffusion-polarization equations')

parser.add_argument("--load", help="load a checkpoint")
args = parser.parse_args()


size = 100

used_material = material.arbitrary
E_field = environment.external_field(strength=0.1, frequency=0.5, wave_length=50)
turing_env = environment.simulation_parameters(epsilon_0 = 1, dx = 0.1, dt = 0.001, external_fields = [E_field])

simulation = engine.grid2D(material, used_material, turing_env, size)



if not args.load == None:
	simulation.load(args.load, False)

observe = engine.analyze(simulation)

def get_obs(n):
	return observe.get_observables()[n]

#fig, ax = plt.subplots()

x = np.arange(0, size)


'''
subplots = 5

fig, axes, = plt.subplots(1, subplots, sharex='col', sharey='row')

[axes[n].set_adjustable('box-forced') for n in range(subplots)]
[axes[n].set_aspect('equal', 'datalim') for n in range(subplots)]


def get_obs(n):
	res = observe.get_observables()[n]
	#res = np.atleast_2d(observe.get_observables()[n])
	return res


im = [axes[n].imshow(get_obs(n), cmap=plt.get_cmap('viridis'), animated=True, interpolation=None) for n in range(subplots)]
[fig.colorbar(im[n],ax=axes[n]) for n in range(subplots)]
1

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


subplots = 2

fig, axes, = plt.subplots(1, subplots, sharex='col', sharey='row')

plotrange_P = [0.1, 0.3]
plotrange_n_c = [0.110, 0.135]
plotrange = [plotrange_n_c, plotrange_P]

#[axes[n].set_adjustable('box-forced') for n in range(subplots)]
#[axes[n].set_aspect('equal', 'datalim') for n in range(subplots)]

def get_obs(n):
	if(n == 0):
		return observe.get_observables()[2]
	if(n == 1):
		return observe.get_observables()[4]
	else:
		return observe.get_observables()[n]

im = [axes[n].imshow(get_obs(n), vmin=plotrange[n][0], vmax=plotrange[n][1], cmap=plt.get_cmap('viridis'), animated=True, interpolation=None) for n in range(subplots)]
[fig.colorbar(im[n],ax=axes[n]) for n in range(subplots)]


def updatefig(*args):
	simulation.evolve()
	#if simulation.current_step % 1000 == 0:
	#	simulation.save("checkpoint-"+str(int(simulation.current_step/1000)))

	im[0].set_array(get_obs(2))
	im[1].set_array(get_obs(4))

	return im

ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=True)

plt.show()
