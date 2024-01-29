# TODO: YOUR NAME AND STUDENT NUMBER HERE

import numpy as np

class Contact:
	# p is point of contact
	# n is normal
	# d is penetration depth
	def __init__(self, body1, body2, p, n, d):
		self.body1 = body1
		self.body2 = body2
		self.p = p
		self.n = n
		if (abs(n[0]) < abs(n[1] and abs(n[0]) < abs(n[2]))):
			tmp = np.array([1,0,0])
		elif (abs(n[1]) < abs(n[2])):
			tmp = np.array([0,1,0])	
		else:
			tmp = np.array([0,0,1])
		self.t2 = np.cross(self.n, tmp)
		self.t1 = np.cross(self.t2, self.n)
		self.d = d
		self.lamb = np.zeros(3)

	def compute_jacobian(self):
		# TODO: implement this function
		normal_transpose = self.n.transpose()
		r1 = self.body1.x
		r2 = self.body2.x
		r1_cross_n_transpose = np.cross(r1,self.n).transpose()
		r2_cross_n_transpose = np.cross(r2,self.n).transpose()


		jrow1 = np.array([-normal_transpose, -r1_cross_n_transpose, normal_transpose, r2_cross_n_transpose])
		# u = np.array([v1,omega1,v2,omega2])

		print(jrow1)

		return jrow1
	

	def compute_inv_effective_mass(self):
		# TODO: implement this function

		jacobian = self.compute_jacobian()

		inv_mass_matrix = self.compute_inverse_mass_matrix()

		inv_effective_mass = np.dot(jacobian, np.dot(inv_mass_matrix, jacobian.T))

		return inv_effective_mass
	
	def compute_inverse_mass_matrix(self):
		inv_mass1 = self.body1.mass_inv
		inv_mass2 = self.body2.mass_inv
		inv_inertia1 = np.linalg.inv(self.body1.J)
		inv_inertia2 = np.linalg.inv(self.body2.J)

		# Form the block matrix
		inv_mass_matrix = np.block([
			[np.eye(3) * inv_mass1, np.zeros((3, 3))],
			[np.zeros((3, 3)), inv_inertia1],
			[np.eye(3) * inv_mass2, np.zeros((3, 3))],
			[np.zeros((3, 3)), inv_inertia2]
		])

		return inv_mass_matrix