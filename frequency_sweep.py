import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import engine, environment
import material1 as material

import argparse

parser = argparse.ArgumentParser(description='Run a frequency sweep of a simulation of reaction-diffusion-polarization equations')

parser.add_argument("--load", help="load a checkpoint")
args = parser.parse_args()


size = 1
steps = int(10000)
used_material = material.simple

simulation = engine.grid(material, used_material, environment.arbitrary, size)


# load start state or equilibrate at beginning
if not args.load == None:
	simulation.load(args.load, False)
else:
	for n in range(steps):
		simulation.evolve()

	simulation.save("equilibrium", False)

	#print("n_a, n_b")
	#print(simulation.cells[0].n_a, simulation.cells[0].n_b)
	#print("yes")


for n in range(150):

	freq = 0.00005 * 1.1**n
	period_t = 2 * 3.1415 / freq

	E_field = environment.external_field(strength=0.1, frequency=freq, wave_length=50)
	sweep_environment = environment.simulation_parameters(epsilon_0 = 1, dx = 0.1, dt = period_t/100000, external_fields = [E_field])

	simulation = engine.grid(material, used_material, sweep_environment, size)
	simulation.load("equilibrium", False)

	observe = engine.analyze(simulation)

	for n in range(steps):
		simulation.evolve()
		observe.calculate_dielectric_response()

	print(freq, observe.chi_accumulated/steps)
