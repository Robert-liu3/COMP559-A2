# ROBERT LIU 260981372

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
		self.jacobian = np.array([])
		self.lamb = np.zeros(3)
		self.side_b_vector = np.zeros(12)
		self.inverse_mass_matrix = np.zeros((12,12))
		self.inv_effective_mass = np.zeros((3,3))

	def compute_jacobian(self):
		# TODO: implement this function
		normal_transpose = self.n.transpose()
		r1 = self.p - self.body1.x
		r2 = self.p - self.body2.x
		r1_cross_n_transpose = np.cross(r1,self.n).transpose()
		r2_cross_n_transpose = np.cross(r2,self.n).transpose()
		jrow1 = np.concatenate([-normal_transpose, -r1_cross_n_transpose, normal_transpose, r2_cross_n_transpose])
		# u = np.array([v1,omega1,v2,omega2])

		t1 = self.t1
		t2 = self.t2
		t1_transpose = t1.transpose()
		t2_transpose = t2.transpose()
		r1_cross_t1_transpose = np.cross(r1,t1).transpose()
		r2_cross_t1_transpose = np.cross(r2,t1).transpose()
		r1_cross_t2_transpose = np.cross(r1,t2).transpose()
		r2_cross_t2_transpose = np.cross(r2,t2).transpose()

		jrow2 = np.concatenate([-t1_transpose, -r1_cross_t1_transpose, t1_transpose, r2_cross_t1_transpose])
		jrow3 = np.concatenate([-t2_transpose, -r1_cross_t2_transpose, t2_transpose, r2_cross_t2_transpose])

		self.jacobian = np.vstack((jrow1,jrow2,jrow3))
	

	def compute_inv_effective_mass(self):
		# TODO: implement this function

		jacobian = self.jacobian

		self.compute_inverse_mass_matrix()

		inv_effective_mass = np.dot(jacobian, np.dot(self.inverse_mass_matrix, jacobian.transpose()))

		self.inv_effective_mass = inv_effective_mass
	
	def compute_inverse_mass_matrix(self):
		inv_mass1 = self.body1.mass_inv
		inv_mass2 = self.body2.mass_inv
		inv_inertia1 = self.body1.Jinv
		inv_inertia2 = self.body2.Jinv

		# Form the block matrix
		inv_mass_matrix = np.block([
			[np.eye(3) * inv_mass1, np.zeros((3, 3)), np.zeros((3, 3)), np.zeros((3, 3))],
			[np.zeros((3, 3)), inv_inertia1, np.zeros((3, 3)), np.zeros((3, 3))],
			[np.zeros((3, 3)), np.zeros((3, 3)), np.eye(3) * inv_mass2, np.zeros((3, 3))],
			[np.zeros((3, 3)), np.zeros((3, 3)), np.zeros((3,3)), inv_inertia2]
		])
		self.inverse_mass_matrix = inv_mass_matrix