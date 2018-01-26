
class constants:
	def __init__(self, k_aab0, k_baa0, k_abc0, k_cab0, k_caaa0, k_aaac0, alpha_aab, alpha_baa, alpha_abc, alpha_cab, alpha_aaac, alpha_caaa, dP, epsilon_ra, epsilon_rb, epsilon_rc, tau_a, tau_b, tau_c, D_na, D_nb, D_nc, D_Pa, D_Pb, D_Pc):
		self.k_aab0 =  k_aab0
		self.k_baa0 = k_baa0
		self.k_abc0 = k_abc0
		self.k_cab0 = k_cab0
		self.k_caaa0 = k_caaa0
		self.k_aaac0 = k_aaac0
		self.alpha_aab = alpha_aab
		self.alpha_baa = alpha_baa
		self.alpha_abc = alpha_abc
		self.alpha_cab = alpha_cab
		self.alpha_aaac = alpha_aaac
		self.alpha_caaa = alpha_caaa
		self.dP = dP
		self.epsilon_ra = epsilon_ra
		self.epsilon_rb = epsilon_rb
		self.epsilon_rc = epsilon_rc
		self.tau_a = tau_a
		self.tau_b = tau_b
		self.tau_c = tau_c
		self.D_na = D_na
		self.D_Pa = D_Pa
		self.D_nb = D_nb
		self.D_Pb = D_Pb
		self.D_nc = D_nc
		self.D_Pc = D_Pc
		
		if self.tau_a == 0:
			print("please don't set relaxation times to 0, if you want instantaneous processes set tau = dt - setting tau_a = 1")
			self.tau_a = 1
		if self.tau_b == 0:
			print("please don't set relaxation times to 0, if you want instantaneous processes set tau = dt - setting tau_b = 1")
			self.tau_b = 1
		if self.tau_c == 0:
			print("please don't set relaxation times to 0, if you want instantaneous processes set tau = dt - setting tau_c = 1")
			self.tau_c = 1
		

simple = constants(
					k_aab0 = 0, k_baa0 = 0, k_abc0 = 0, k_cab0 = 0, k_caaa0 = 0, k_aaac0 = 0,
					alpha_aab = 0.0, alpha_baa = 0.0, alpha_abc = 0.0, alpha_cab = -0.0, alpha_aaac = 0, alpha_caaa = 0, 
					epsilon_ra = 1.0, epsilon_rb = 2.0, epsilon_rc = 3.0, dP = 0.0,
					tau_a = 1.0, tau_b = 1.0, tau_c = 1.0,
					D_na = 0.1, D_nb = 0.1, D_nc = 0.1, D_Pa = 0.1, D_Pb = 0.1, D_Pc = 0.1)

monoalcohol = constants(
					k_aab0 = 0.001, k_baa0 = 0.002, k_abc0 = 0.005, k_cab0 = 0.005, k_caaa0 = 0.004, k_aaac0 = 0.0001,
					alpha_aab = 0.01, alpha_baa = -0.01, alpha_abc = 0.01, alpha_cab = -0.01, alpha_aaac = 0.1, alpha_caaa = -0.01, 
					epsilon_ra = 1.0, epsilon_rb = 1.2, epsilon_rc = 1.3, dP = 0.0,
					tau_a = 1.0, tau_b = 2.0, tau_c = 3.0,
					D_na = 0.00028, D_nb = 0.05, D_nc = 0.05, D_Pa = 0.00028, D_Pb = 0.05, D_Pc = 0.05)

arbitrary = constants(
					k_aab0 = 0.001, k_baa0 = 0.002, k_abc0 = 0.005, k_cab0 = 0.005, k_caaa0 = 0.004, k_aaac0 = 0.0001,
					alpha_aab = 0.01, alpha_baa = -0.01, alpha_abc = 0.01, alpha_cab = -0.01, alpha_aaac = 0.1, alpha_caaa = -0.01, 
					epsilon_ra = 1.0, epsilon_rb = 1.2, epsilon_rc = 1.3, dP = 0.0,
					tau_a = 1.0, tau_b = 2.0, tau_c = 3.0,
					D_na = 0.00028, D_nb = 0.05, D_nc = 0.05, D_Pa = 0.00028, D_Pb = 0.05, D_Pc = 0.05)
