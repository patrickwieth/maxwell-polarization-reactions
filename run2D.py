import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import engine2D, environment
import material2D as material

import argparse

parser = argparse.ArgumentParser(description='Run a simulation of reaction-diffusion-polarization equations')

parser.add_argument("--load", help="load a checkpoint")
args = parser.parse_args()

size = 15
freq = 0.1

used_material = material.monoalcohol_Peqd
E_field = environment.external_field(strength=0.1, frequency=freq, wave_length=50)
env = environment.simulation_parameters(epsilon_0 = 1, dx = 1.0, dt = 0.01/freq, external_fields = [E_field])

simulation = engine2D.grid2D(material, used_material, env, size)


#if not args.load == None:
#	simulation.load(args.load, False)

observe = engine2D.analyze(simulation)

def get_obs(n):
	return observe.get_observables()[n]

simulation.evolve()


state = observe.get_observables()

#print(state)



X, Y = np.meshgrid(np.arange(0, size, 1), np.arange(0, size, 1))
U = state[3, X, Y]
V = state[4, X, Y]

fig, ax = plt.subplots(1,1, figsize=(10, 10))

###############################################################################

#plt.figure()
#plt.title('Arrows scale with plot width, not view')
Q = ax.quiver(X, Y, U, V, units='width')
#qk = plt.quiverkey(Q, 0.9, 0.9, 2, r'$2 \frac{m}{s}$', labelpos='E',
#                   coordinates='figure')

#ax.set_xlim(-1, 3)
#ax.set_ylim(-1, 3)

def updatefig(num, Q, X, Y):

	#if simulation.current_step % 10 == 0:
		# save image
		#plt.savefig("dipoles/steps"+str(simulation.current_step), format="svg")
		# save simulationa
		#simulation.save("checkpoint-"+str(int(simulation.current_step/1000)))

	simulation.evolve()

	state = observe.get_observables()

	U = state[3, X, Y]
	V = state[4, X, Y]

	max_c = np.amax(state[2, :, :])
	min_c = np.amin(state[2, :, :])

	if max_c == min_c:
		C = 1
	else:
		C = (state[2, X, Y]-min_c) / (max_c-min_c)

	Q.set_UVC(U,V,C)	# C: yellow = max, purple = min

	#print(state)
	print(simulation.current_step, "min/max:", min_c, max_c)

	return Q,

# you need to set blit=False, or the first set of arrows never gets
# cleared on subsequent frames
anim = animation.FuncAnimation(fig, updatefig, fargs=(Q, X, Y),
                               interval=50, blit=False)
fig.tight_layout()
plt.show()
