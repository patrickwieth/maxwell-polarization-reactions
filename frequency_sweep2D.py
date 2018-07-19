import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import engine2D, environment
import material2D as material

import argparse

parser = argparse.ArgumentParser(description='Run a frequency sweep of a 2D simulation of reaction-diffusion-polarization equations')
parser.add_argument("--load", help="load a checkpoint")
args = parser.parse_args()


size = 15
steps = int(10000)
used_material = material.monoalcohol

E_field = environment.external_field(strength=1, frequency=0, wave_length=50)
sweep_environment = environment.simulation_parameters(epsilon_0 = 1, dx = 0.1, dt = 0.01, external_fields = [E_field])
simulation = engine2D.grid2D(material, used_material, sweep_environment, size)




# pre-equilibration
for n in range(steps):
	simulation.evolve()

observe = engine2D.analyze(simulation)

print("pre-equilibration done")
#print(observe.get_observables())


for n in range(10):

	#freq = 0.0002 * 1.5**n
	freq = 0.000002 * 1.5**n

	simulation.parameters.dt = 0.001/freq
	simulation.parameters.force_fields[0].frequency = freq

	observe = engine2D.analyze(simulation)

	for i in range(steps):
		simulation.evolve()
		#print(i, simulation.cells[0,0].E[1], simulation.cells[0,0].n_a, simulation.cells[0,0].n_b, simulation.cells[0,0].n_c)
		observe.calculate_dielectric_response()

	print(freq, observe.chi_accumulated/steps)
