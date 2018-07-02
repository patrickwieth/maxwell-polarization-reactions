import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import engine2D, environment
import material2D as material

import argparse

parser = argparse.ArgumentParser(description='Run a frequency sweep of a 2D simulation of reaction-diffusion-polarization equations')
parser.add_argument("--load", help="load a checkpoint")
args = parser.parse_args()


size = 10
steps = int(1000)
used_material = material.monoalcohol

E_field = environment.external_field(strength=1, frequency=0, wave_length=50)
sweep_environment = environment.simulation_parameters(epsilon_0 = 1, dx = 0.1, dt = 0.01, external_fields = [E_field])
simulation = engine2D.grid2D(material, used_material, sweep_environment, size)




# pre-equilibration
for n in range(steps):
	simulation.evolve()

observe = engine.analyze(simulation)

print("pre-equilibration done")
print(observe.get_observables())

	#print("n_a, n_b")
	#print(simulation.cells[0].n_a, simulation.cells[0].n_b)
	#print("yes")


for n in range(150):

	freq = 0.0000002 * 1.1**n
	#period_t = 2 * 3.1415 / freq

	#E_field = environment.external_field(strength=1, frequency=freq, wave_length=50)
	#sweep_environment = environment.simulation_parameters(epsilon_0 = 1, dx = 0.1, dt = period_t/100000000, external_fields = [E_field])


	#simulation = engine2D.grid(material, used_material, sweep_environment, size)
	#simulation.load("equilibrium", False)

	simulation.parameters.force_fields[0].frequency = freq


	for n in range(steps):
		simulation.evolve()
		observe.calculate_dielectric_response()

	print(freq, observe.chi_accumulated/steps)
