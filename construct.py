#class created to construct linear system of equations

import numpy as np

class Construct():
	
###############################################################
	
	def __init__(self):
		self
		
###############################################################
		
	def constructSS(self, opts, D, rXS, Gscat, fisXS):

		nGrps = opts.numGroups
		nBins = opts.numBins
		rank = nBins*nGrps
		delta=opts.delta

		self.A = np.zeros((rank, rank))
		self.F = np.zeros((rank, rank))

		for g in range(0, nGrps):
			for x in range(0,nBins):
				i = g*nBins+x
				

				# Calculate the average diffusion theory between current cell and the adjacent cell
				# and account for vacuum boundary conditions
				if x == 0:
					DAdj = 1
					deltaAdj = 0 # 0 for Zero flux, 4 for Marshak
					dLeft = 0 # include for reflective BC only
				else:
					DAdj = D[i-1]
					deltaAdj = delta
					dLeft = (2*D[i]*DAdj)/(deltaAdj*D[i]+delta*DAdj) # include for reflective BC only
				# dLeft = (2*D[i]*DAdj)/(deltaAdj*D[i]+delta*DAdj) # if NOT reflective BC

				if x == (nBins-1):
					DAdj = 1
					deltaAdj = 0 # 0 for Zero flux, 4 for Marshak
					dRight = 0 # include for reflective BC only
				else:
					DAdj = D[i+1]
					deltaAdj = delta
					dRight = (2*D[i]*DAdj)/(deltaAdj*D[i]+delta*DAdj) # include for reflective BC only
				# dRight = (2*D[i]*DAdj)/(deltaAdj*D[i]+delta*DAdj) # if NOT reflective BC

				self.A[i,i] = dLeft+dRight+rXS[i]*delta
				# print self.A

				if x != 0:
					self.A[i,i-1] = -dLeft
				if x != (nBins-1):
					self.A[i,i+1] = -dRight

				# Including downscattering
				if g == 0:
					self.A[i+nBins,i] = -Gscat[x,1]*delta

				# Fission matrix
				if g == 0:
					self.F[i,i] = fisXS[i]*delta
					self.F[i,i+nBins] = fisXS[i+nBins]*delta

		if opts.flux == 'adjoint':
			print("ADJOINT PROBLEM")
			self.A = np.transpose(self.A)
			self.F = np.transpose(self.F)
		

		# print self.A
		# print self.F
				
###############################################################
		
	def constructTD(self, opts, D, PrXS, Gscat, PfisXS):

		nGrps = opts.numGroups
		nBins = opts.numBins
		rank = nBins*nGrps
		delta=opts.delta

		self.A = np.zeros((rank, rank))
		self.F = np.zeros((rank, rank))

		for g in range(0, nGrps):
			for x in range(0,nBins):
				i = g*nBins+x
				

				# Calculate the average diffusion theory between current cell and the adjacent cell
				# and account for vacuum boundary conditions
				if x == 0:
					DAdj = 1
					deltaAdj = 0 # 0 for Zero flux, 4 for Marshak
					dLeft = 0 # include for reflective BC only
				else:
					DAdj = D[i-1]
					deltaAdj = delta
					dLeft = (2*D[i]*DAdj)/(deltaAdj*D[i]+delta*DAdj) # include for reflective BC only
				# dLeft = (2*D[i]*DAdj)/(deltaAdj*D[i]+delta*DAdj) # if NOT reflective BC

				if x == (nBins-1):
					DAdj = 1
					deltaAdj = 0 # 0 for Zero flux, 4 for Marshak
					dRight = 0 # include for reflective BC only
				else:
					DAdj = D[i+1]
					deltaAdj = delta
					dRight = (2*D[i]*DAdj)/(deltaAdj*D[i]+delta*DAdj) # include for reflective BC only
				# dRight = (2*D[i]*DAdj)/(deltaAdj*D[i]+delta*DAdj) # if NOT reflective BC

				if g == 0:
					F1 = PfisXS[i]*delta
					self.A[i,i+nBins] = -PfisXS[i+nBins]*delta
				else:
					F1 = 0
					F2 = 0

				self.A[i,i] = dLeft+dRight+PrXS[i]*delta-F1
				# print self.A

				if x != 0:
					self.A[i,i-1] = -dLeft
				if x != (nBins-1):
					self.A[i,i+1] = -dRight

				# Including downscattering
				if g == 0:
					self.A[i+nBins,i] = -Gscat[x,1]*delta

				# Fission matrix
				if g == 0:
					self.F[i,i] = PfisXS[i]*delta
					self.F[i,i+nBins] = PfisXS[i+nBins]*delta

				# Source vector


		if opts.flux == 'adjoint':
			print("ADJOINT PROBLEM")
			self.A = np.transpose(self.A)
			self.F = np.transpose(self.F)
		

		# print self.A
		# print self.F
				
###############################################################
	def invert(self, A):
		
		self.inv = np.linalg.inv(A)
		return self.inv


#end class 
