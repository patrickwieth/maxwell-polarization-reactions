import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import engine, environment
import material2 as material

import argparse

parser = argparse.ArgumentParser(description='Run a frequency sweep of a simulation of reaction-diffusion-polarization equations')

parser.add_argument("--load", help="load a checkpoint")
args = parser.parse_args()


size = 1

used_material = material.arbitrary

simulation = engine.grid(used_material, environment.arbitrary, size)


steps = int(10000)


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

	freq = 0.0000001 * 1.1**n
	period_t = 2 * 3.1415 / freq
	
	E_field = environment.external_field(strength=1, frequency=freq, wave_length=50)
	sweep_environment = environment.simulation_parameters(epsilon_0 = 1, c = 1.0, dx = 0.1, dt = period_t/100000000, external_fields = [E_field])
	

	simulation = engine.grid(used_material, sweep_environment, size)
	simulation.load("equilibrium", False)

	observe = engine.analyze(simulation)

	for n in range(steps):
		simulation.evolve()
		observe.calculate_dielectric_response()

	print(observe.chi_accumulated/steps)
