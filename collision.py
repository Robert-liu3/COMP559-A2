# ROBERT LIU 260981372

import igl
import numpy as np
import polyscope as ps
from rigidbody import *
from contact import *

class Collision:
	def __init__(self):
		self.contacts = []
		self.ground_body = RigidBody()

	def reset(self):
		self.contacts = []

	def update_display(self, show_contacts):
		if len(self.contacts) == 0 or not show_contacts:
			self.ps_contacts = ps.register_point_cloud("contacts", np.zeros((0,3)))
		else:
			# can only update the points if they have the same number :(
			pos = np.array([c.p for c in self.contacts])
			depth = np.array([c.d for c in self.contacts])
			normal = np.array([c.n for c in self.contacts])
			t1 = np.array([c.t1 for c in self.contacts])
			t2 = np.array([c.t2 for c in self.contacts])
			force_viz_scale = 2
			force = np.array([force_viz_scale*(c.lamb[0] * c.n + c.lamb[1] * c.t1 + c.lamb[2] * c.t2) for c in self.contacts])
			self.ps_contacts = ps.register_point_cloud("contacts", pos)
			self.ps_contacts.add_scalar_quantity("contact depth", depth, enabled=True)	
			# self.ps_contacts.add_vector_quantity("contact normal", normal, enabled=True, radius=0.01, color=(0,0,1), vectortype='ambient')
			# self.ps_contacts.add_vector_quantity("contact t1", t1, enabled=True, radius=0.01, color=(1,0,0), vectortype='ambient')
			# self.ps_contacts.add_vector_quantity("contact t2", t2, enabled=True, radius=0.01, color=(0,1,0), vectortype='ambient')
			self.ps_contacts.add_vector_quantity("contact force", force, enabled=True, radius=0.01, color=(1,1,0), vectortype='ambient')	
				
	# collision check with ground plane
	def check_ground( self, body ):
		# check if any vertex is below the ground plane
		# if so, compute penetration depth and normal and add to contact list
		vt = body.V @ body.R.T + body.x
		for v in vt:
			if v[1] < 0: # y is up
				self.contacts.append( Contact( self.ground_body, body, v, np.array([0,1,0]), v[1] ) )
		
		# collision check between two bodies
		# for contact in self.contacts:
		# 	print(contact.compute_jacobian())
			
		
	def check_body_pair( self, body1, body2 ):
		# check if any vertex of one body is inside the other
		# NOTE: this is super gross because the signed distance function is expensive
		# thus we check the larger number of vertices with the smaller body
		# but WATCH OUT becaues this appears to be buggy in general.
		# For vericies inside the other body, compute penetration depth and normal
		v1t = body1.V @ body1.R.T + body1.x
		v2t = body2.V @ body2.R.T + body2.x
		if ( v1t.shape[0] > v2t.shape[0] ):
			S,I,C,N = igl.signed_distance( v1t, v2t, body2.F, return_normals=True )
			for i in range(len(S)):
				if S[i] < 0:
					self.contacts.append( Contact( body1, body2, C[i], -N[i], -S[i] ) )
		else:
			S,I,C,N = igl.signed_distance( v2t, v1t, body1.F, return_normals=True )
			for i in range(len(S)):
				if S[i] < 0:
					self.contacts.append( Contact( body2, body1, C[i], -N[i], -S[i] ) )

	def check( self, rigid_body_list ):
		self.contacts = []
		for i in range(len(rigid_body_list)):
			self.check_ground(rigid_body_list[i])
			for j in range(i+1, len(rigid_body_list)):
				self.check_body_pair(rigid_body_list[i], rigid_body_list[j])
		return len(self.contacts) > 0

	def process(self, rigid_body_list, mu, num_iter):
		#TODO: implement this function
		side_b_vector = []
		for contact in self.contacts:
			contact.compute_jacobian()
			v1 = np.concatenate((contact.body1.v,contact.body1.omega))
			v2 = np.concatenate((contact.body2.v, contact.body2.omega)) 
			relative_velocity = np.concatenate((v1,v2))
			contact.side_b_vector = np.dot(contact.jacobian, relative_velocity)
		
		for i in range(0,num_iter):
			for contact in self.contacts:
				for i in range(0,3):
					contact.compute_inv_effective_mass()
					A = contact.inv_effective_mass
					numerator = contact.side_b_vector[i] + np.dot(A[i,:], contact.lamb)
					denominator = A[i,i]
					new_lamb = contact.lamb[i] - numerator/denominator
					old_lambda = contact.lamb[i]

					if (i == 0):
						if (new_lamb < 0):
							new_lamb = 0
						contact.lamb[0] = new_lamb
					else:
						lambda_max = mu*contact.lamb[0]
						lambda_min = -lambda_max
						new_lamb = min(max(lambda_min, new_lamb), lambda_max)
						contact.lamb[i] = new_lamb

					delta_lambda = contact.lamb[i] - old_lambda
					delta_v_inital = np.concatenate((contact.body1.deltav,contact.body2.deltav))
					T = contact.inverse_mass_matrix@contact.jacobian.transpose()
					delta_v = delta_v_inital + T[:,i]*delta_lambda
					contact.body1.deltav = delta_v[0:6]
					contact.body2.deltav = delta_v[6:12]
		
		for rb in rigid_body_list:
			rb.v += rb.deltav[0:3]
			rb.omega += rb.deltav[3:6]
			rb.deltav = np.zeros(6)

			