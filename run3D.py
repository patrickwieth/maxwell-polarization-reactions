import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import engine3D, environment
import material3D as material

import argparse

parser = argparse.ArgumentParser(description='Run a simulation of reaction-diffusion-polarization equations')

parser.add_argument("--load", help="load a checkpoint")
args = parser.parse_args()

size = 3

used_material = material.monoalcohol
E_field = environment.external_field(strength=0.00, frequency=0.1, wave_length=50)
env = environment.simulation_parameters(epsilon_0 = 1, c = 1.0, dx = 0.1, dt = 0.1, external_fields = [E_field])

simulation = engine3D.grid3D(used_material, env, size)


#if not args.load == None:
#	simulation.load(args.load, False)

observe = engine3D.analyze(simulation)

def get_obs(n):
	return observe.get_observables()[n]

simulation.evolve()


state = observe.get_observables()

#print(state)


X, Y, Z = np.meshgrid(np.arange(0, size, 1), np.arange(0, size, 1), np.arange(0, size, 1))
U = state[3, X, Y, Z]
V = state[4, X, Y, Z]
W = state[5, X, Y, Z]

#fig, ax = plt.subplots(1,1)

fig = plt.figure()
ax = fig.gca(projection='3d')

###############################################################################

#plt.figure()
#plt.title('Arrows scale with plot width, not view')
Q = ax.quiver(X, Y, Z, U, V, W, length=0.5, normalize=True)
#qk = plt.quiverkey(Q, 0.9, 0.9, 2, r'$2 \frac{m}{s}$', labelpos='E',
#                   coordinates='figure')

#ax.set_xlim(-1, 3)
#ax.set_ylim(-1, 3)

simulation.evolve()

def updatefig(num, Q, X, Y, Z):
	simulation.evolve()
	#if simulation.current_step % 1000 == 0:
	#	simulation.save("checkpoint-"+str(int(simulation.current_step/1000)))

	state = observe.get_observables()

	U = state[3, X, Y, Z]
	V = state[4, X, Y, Z]
	W = state[5, X, Y, Z]

	#buttsegs = 

	print(Q)

	Q.set_segments(buttsegs)
	#Q.set_UVC(U,V,W)


	print(state)

	return Q,

# you need to set blit=False, or the first set of arrows never gets
# cleared on subsequent frames
anim = animation.FuncAnimation(fig, updatefig, fargs=(Q, X, Y, Z),
                               interval=50, blit=False)
fig.tight_layout()

plt.show()

