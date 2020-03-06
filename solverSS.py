# class to solve by power iteration

from scipy.sparse import coo_matrix
import matplotlib.pyplot as plt
import numpy as np
import copy

class SolveSS:
	
	def __init__(self, opts):
		self.rank = opts.numBins*opts.numGroups
		self.flux = np.zeros(self.rank)
		# self.flux = np.random.rand(self.rank)
		self.B = np.zeros(self.rank)

		# Initial guess 
		self.flux[:] = 1/opts.length
		self.B[:] = 1
		self.k = 1

		# print self.flux
		# plt.plot(self.flux)
		# plt.show()
		# plt.clf()

###############################################################
		
	def NumpySolve(self, opts, Ainv, F, XS):
		

		# Create local array variables.
		RMSflux = np.zeros(self.rank)
		RMSsource = np.zeros(self.rank)

		# Initialize local variables
		self.stored_k = []
		self.stored_ERRsource = []
		self.stored_DR = []
		ERRflux = 1000.
		ERRsource = 1000.
		j = 0

		while ERRflux > opts.FluxConvError or ERRsource > opts.FisConvError:

			# Reset local variables
			j = j+1
			RMSsource[:] = 0
			n = 0

			# Save previous source so that we can compute the
			# residual error for the iteration
			lastFlux = self.flux
			lastB = self.B
			
			self.B = np.dot(F,self.flux)/self.k
			self.flux = np.dot(Ainv,self.B)
			self.k = sum(np.dot(F,self.flux))/sum(self.B) 
			self.stored_k.append(self.k)

			# Perform the source normalization
			# self.flux[:] = self.flux[:]/np.linalg.norm(self.flux)
			# lastB = self.B
			# self.B = np.dot(F,self.flux)/self.k

			

			# Calculate the relative difference in the source between
			# consecutive iterations and take the infinity norm.
			RMSflux[:] = abs((lastFlux[:]-self.flux[:])/self.flux[:])
			ERRflux = np.linalg.norm(RMSflux, 2)/len(self.flux)
			n = 0
			for i in range(0,len(self.B)):
				if self.B[i] != 0:
					n = n+1
					RMSsource[i] = abs((lastB[i]-self.B[i])/self.B[i])
			# 		RMSsource[i] = (lastB[i]/self.B[i]-1.0)**2.0
			# ERRsource = np.sqrt(sum(RMSsource)/n)
			ERRsource = np.linalg.norm(RMSsource, 2)/len(self.B)
			#ERRflux = max(RMSflux)/len(self.flux)
			# print(ERRsource)
			
			self.stored_ERRsource.append(ERRsource)
			if j > 2:
				self.stored_DR.append(self.stored_ERRsource[-1]/self.stored_ERRsource[-2])

			# if j == 1 or j == 10 or j == 20 or j == 30 or j == 40 or j == 50 or j == 60 or j == 70 or j == 80 or j == 90 or j == 100: 
			# 	plt.plot(self.flux)
			# 	plt.savefig('./test')
			

			# Print statement to show eigenvalue convergence by iteration
			# print "Eigenvalue: ", self.k, "(", ERRflux, ";", ERRsource, ") Iteration:", j
			# print self.flux
		# print "Eigenvalue: ", self.k, "(", ERRflux, ";", ERRsource, ") Iteration:", j

		# # Dominance Ratio by direct inspection
		# temp = np.matmul(Ainv,F)
		# ALLeig = np.real(np.sort(np.linalg.eig(temp)[0]))
		# DR = ALLeig[-2]/ALLeig[-1]
		# print('True eig:',ALLeig[-1])
		# print('True DR:',DR)

		self.DR = self.stored_ERRsource[-1]/self.stored_ERRsource[-2]
		self.source_it = j

###############################################################

	def PointJacobi(self, opts, Matrix, XS):

		

		D = np.diag(Matrix.A)
		O = coo_matrix(Matrix.A - np.diag(D))
		# print(O)
		# O = Matrix.A - np.diag(D)

		# Create local array variables
		RMSflux = np.zeros(self.rank)
		RMSsource = np.zeros(self.rank)

		# Initialize local variables
		self.stored_k = []
		self.stored_ERRsource = []
		self.stored_DR = []
		ERRsource = 1000.
		j = 0
		m = 0

		# Fission source iteration (outer loop)
		while ERRsource > opts.FisConvError:
		# for p in range(0,1000):

			# Reset local variables
			j = j+1
			RMSsource[:] = 0
			n = 0
			lastB = self.B
			self.B = np.dot(Matrix.F,self.flux)/self.k

			# Flux iteration (inner loop)
			ERRflux = 1000.
			while ERRflux > opts.FluxConvError:
			# for q in range(0,200):
				m = m + 1
				lastFlux = self.flux
				# for n in range(len(self.flux)):
				# 	select_ind = np.array(n)
				# 	self.flux[n] = self.B[n] - sum(O.tocsr()[select_ind,:]*self.flux[n])
				self.flux = self.B - O.dot(self.flux) #np.dot(O,self.flux) 
				self.flux[:] = self.flux[:]/D[:]
				RMSflux[:] = abs((lastFlux[:]-self.flux[:])/self.flux[:])
				ERRflux = np.linalg.norm(RMSflux, 2)/len(self.flux)
			# print ERRflux

			# print(m)

			# Calculate the fission source in each spatial bin
			self.k = sum(np.dot(Matrix.F,self.flux))/sum(self.B) 
			self.stored_k.append(self.k)

			# Perform the source normalization
			# self.flux[:] = self.flux[:]/np.linalg.norm(self.flux)

			# Calculate the relative difference in the source between
			# consecutive iterations and take the infinity norm.
			# RMSflux[:] = abs((lastFlux[:]-self.flux[:])/self.flux[:])
			for i in range(0,len(self.B)):
				if self.B[i] != 0:
					n = n+1
					RMSsource[i] = abs((lastB[i]-self.B[i])/self.B[i])
			ERRsource = np.linalg.norm(RMSsource, 2)/len(self.B)
			self.stored_ERRsource.append(ERRsource)
			if j > 2:
				self.stored_DR.append(self.stored_ERRsource[-1]/self.stored_ERRsource[-2])
			# print(ERRsource)



			# if j == 1 or j == 10 or j == 20 or j == 30 or j == 40 or j == 50 or j == 60 or j == 70 or j == 80 or j == 90 or j == 100: 
			# 	plt.plot(self.flux)
			# 	plt.savefig('./test')
			

			# Print statement to show eigenvalue convergence by iteration
			# print "Eigenvalue: ", self.k, "(", ERRflux, ";", ERRsource, ") Iteration:", j
			# print self.flux
		# print "Eigenvalue: ", self.k, "(", ERRflux, ";", ERRsource, ") Iteration:", j

		# Dominance Ratio by direct inspection
		# Ainv = np.inverse(Matrix.A)
		# temp = np.matmul(Ainv,F)
		# ALLeig = np.real(np.sort(np.linalg.eig(temp)[0]))
		# self.DR = ALLeig[-2]/ALLeig[-1]

		self.DR = self.stored_ERRsource[-1]/self.stored_ERRsource[-2]
		self.source_it = j
		self.flux_it = m


###############################################################

	def GaussSeidel(self, opts, Matrix, XS):

		D = np.diag(Matrix.A)
		# O = Matrix.A - np.diag(D)
		O = coo_matrix(Matrix.A - np.diag(D))

		# Create local array variables
		RMSflux = np.zeros(self.rank)
		RMSsource = np.zeros(self.rank)

		# Initialize local variables
		self.stored_k = []
		self.stored_ERRsource = []
		self.stored_DR = []
		ERRsource = 1000.
		j = 0
		m = 0

		# Fission source iteration (outer loop)
		while ERRsource > opts.FisConvError:

			# Reset local variables
			j = j+1
			RMSsource[:] = 0
			n = 0
			lastB = copy.copy(self.B)
			self.B = np.dot(Matrix.F,self.flux)/self.k

			# Flux iteration (inner loop)
			ERRflux = 1000.
			while ERRflux > opts.FluxConvError:
				m = m+1
				lastFlux = copy.copy(self.flux)
				self.flux[::2] = (self.B - O.dot(self.flux))[::2]/D[::2]
				self.flux[1::2] = (self.B - O.dot(self.flux))[1::2]/D[1::2]
				# for n in range(len(self.flux)):
					# self.flux[n] = (self.B[n] - np.dot(O[n,:],self.flux))/D[n]
				RMSflux[:] = abs((lastFlux[:]-self.flux[:])/self.flux[:])
				ERRflux = np.linalg.norm(RMSflux, 2)/len(self.flux)
				# print ERRflux

			# Calculate the fission source in each spatial bin
			self.k = sum(np.dot(Matrix.F,self.flux))/sum(self.B) 
			self.stored_k.append(self.k)

			# Perform the source normalization
			# self.flux[:] = self.flux[:]/np.linalg.norm(self.flux)

			# Calculate the relative difference in the source between
			# consecutive iterations and take the infinity norm.
			# RMSflux[:] = abs((lastFlux[:]-self.flux[:])/self.flux[:])
			for i in range(len(self.B)):
				if self.B[i] != 0:
					n = n+1
					RMSsource[i] = abs((lastB[i]-self.B[i])/self.B[i])
			ERRsource = np.linalg.norm(RMSsource, 2)/len(self.B)
			self.stored_ERRsource.append(ERRsource)
			if j > 2:
				self.stored_DR = self.stored_ERRsource[-1]/self.stored_ERRsource[-2]
			print(ERRsource)

			# if j == 1 or j == 10 or j == 20 or j == 30 or j == 40 or j == 50 or j == 60 or j == 70 or j == 80 or j == 90 or j == 100: 
			# 	plt.plot(self.flux)
			# 	plt.savefig('./test')
			

			# Print statement to show eigenvalue convergence by iteration
			# print "Eigenvalue: ", self.k, "(", ERRflux, ";", ERRsource, ") Iteration:", j
			# print self.flux
		# print "Eigenvalue: ", self.k, "(", ERRflux, ";", ERRsource, ") Iteration:", j

		# # Dominance Ratio by direct inspection
		# temp = np.matmul(np.linalg.inv(Matrix.A),Matrix.F)
		# ALLeig = np.real(np.sort(np.linalg.eig(temp)[0]))
		# DR = ALLeig[-2]/ALLeig[-1]
		# print("The true dominance ratio is:", DR)

		self.DR = self.stored_ERRsource[-1]/self.stored_ERRsource[-2]
		self.source_it = j
		self.flux_it = m

