import math

class simulation_parameters:
	def __init__(self, epsilon_0, dx, dt, external_fields):

		#self.wave_length = wave_length
		self.force_fields = external_fields
		self.epsilon_0 = epsilon_0
		self.dt = dt
		self.dx = dx


class external_field:
	def __init__(self, strength, frequency, wave_length):
		self.strength = strength
		self.frequency = frequency
		self.wave_length = wave_length

	def get_strength(self, t):
		'''
		self.cells[ 0, 0].E = field
		self.cells[-1, 0].E = field

		# calculate difference of E-Field between Borders
		delta_E = (self.cells[0, 0].E - self.cells[-1, 0].E)/self.size

		def set_E(cell, idx, idy):
			cell.E = self.cells[0, 0].E #+ idx * delta_E
		'''

		if self.frequency > 0:
			return self.strength * math.sin(t * self.frequency * 2*math.pi) #/ self.wave_length)
		else:
			return self.strength

no_E_field = external_field(strength=0.0, frequency=0.0, wave_length=50)

arbitrary = simulation_parameters(epsilon_0 = 1, dx = 0.1, dt = 0.01, external_fields = [no_E_field])
