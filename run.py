import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import engine, internal, material, environment

import argparse

parser = argparse.ArgumentParser(description='Run a simulation of reaction-diffusion-polarization equations')

parser.add_argument("--load", help="load a checkpoint")
args = parser.parse_args()


size = 1

used_material = material.monoalcohol
E_field = environment.external_field(strength=0.0, frequency=0.0, wave_length=50)

simulation = engine.grid(used_material, environment.arbitrary, size, [E_field])




if not args.load == None:
	simulation.load(args.load, False)

observe = engine.analyze(simulation)

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