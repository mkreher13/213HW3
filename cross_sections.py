# Class to fill cross section vectors

import numpy as np

class CrossSections:

###############################################################

	def __init__(self, opts, kcrit=1.0):
		self.abs = []
		self.Pfis = []
		self.fis = []
		self.Premoval = []
		self.removal = []
		self.D = []
		self.Gscat = np.zeros((opts.numBins,opts.numGroups*opts.numGroups))

		########################################################
		#     Delayed neutron data
		########################################################
		self.Bi = np.array([0.0218, 0.1023, 0.0605, 0.131, 0.22, 0.06, 0.054, 0.0152]) * 0.01
		self.Btot = sum(self.Bi)
		self.Thalf = [55.6, 24.5, 16.3, 5.21, 2.37, 1.04, 0.424, 0.195]  #half lives [s]
		self.Decayi = np.log(2)/self.Thalf
		self.vel1 = 2200. * 100. * (0.1e4/0.0253)**(0.5)
		self.vel2 = 2200. * 100. * (0.100/0.0253)**(0.5)

		self.Ffactor = (1-self.Btot)/kcrit+sum(
			(self.Bi[:]*self.Decayi[:]*opts.dT)/(1+self.Decayi[:]*opts.dT)/kcrit)
		self.Rfactor = [1/(self.vel1*opts.dT), 1/(self.vel2*opts.dT)]


###############################################################

	def problem1(self, opts, data):

		nBins = opts.numBins
		nGrps = opts.numGroups

		for g in range(1,nGrps+1):
			if g == 1:
				gIn = 2
			else:
				gIn = 1
			for i in range(nBins*(g-1),nBins*g):
				if (i < nBins*(g-1)+30) or (i >=nBins*g-30):
					self.fis.append(data['water']['fisXS'][g])
					self.abs.append(data['water']['absXS'][g])
					self.removal.append((data['water']['absXS'][g]+data['water']['Ex'+str(g)][gIn]))
					self.D.append(data['water']['D'][g])
				elif (i < nBins*(g-1)+200) or (i >=nBins*g-200):
					self.fis.append(data['fuel1']['fisXS'][g])
					self.abs.append(data['fuel1']['absXS'][g])
					self.removal.append((data['fuel1']['absXS'][g]+data['fuel1']['Ex'+str(g)][gIn]))
					self.D.append(data['fuel1']['D'][g])
				else:
					self.fis.append(data['enriched_rod']['fisXS'][g])
					self.abs.append(data['enriched_rod']['absXS'][g])
					self.removal.append((data['enriched_rod']['absXS'][g]+data['enriched_rod']['Ex'+str(g)][gIn]))
					self.D.append(data['enriched_rod']['D'][g])


		# Gscat matrix
		for i in range(0,nBins):
			count = 1
			gcount = 1
			for g in range(0,nGrps*nGrps):
				if i < 30 or i >= nBins-30:
					self.Gscat[i,g] = data['water']['Ex'+str(gcount)][count]
				elif i < 200 or i >= nBins-200:
					self.Gscat[i,g] = data['fuel1']['Ex'+str(gcount)][count]
				else:
					self.Gscat[i,g] = data['enriched_rod']['Ex'+str(gcount)][count]
				count = count+1
				if count == nGrps+1:
					count = 1
					gcount = gcount+1


###############################################################

	def problem2(self, opts, data):

		nBins = opts.numBins
		nGrps = opts.numGroups

		for g in range(1,nGrps+1):
			if g == 1:
				gIn = 2
			else:
				gIn = 1
			for i in range(nBins*(g-1),nBins*g):
				if (i < nBins*(g-1)+30) or (i >=nBins*g-30):
					self.Pfis.append(data['water']['fisXS'][g]*self.Ffactor)
					self.fis.append(data['water']['fisXS'][g])
					self.abs.append(data['water']['absXS'][g])
					self.Premoval.append(data['water']['absXS'][g]+data['water']['Ex'+str(g)][gIn]+self.Rfactor[g-1])
					self.removal.append(data['water']['absXS'][g]+data['water']['Ex'+str(g)][gIn])
					self.D.append(data['fuel1']['D'][g])
				elif (i < nBins*(g-1)+200) or (i >=nBins*g-200):
					self.Pfis.append(data['fuel1']['fisXS'][g]*self.Ffactor)
					self.fis.append(data['fuel1']['fisXS'][g])
					self.abs.append(data['fuel1']['absXS'][g])
					self.Premoval.append(data['fuel1']['absXS'][g]+data['fuel1']['Ex'+str(g)][gIn]+self.Rfactor[g-1])
					self.removal.append(data['fuel1']['absXS'][g]+data['fuel1']['Ex'+str(g)][gIn])
					self.D.append(data['fuel1']['D'][g])
				else:
					self.Pfis.append(data['enriched_rod']['fisXS'][g]*self.Ffactor)
					self.fis.append(data['enriched_rod']['fisXS'][g])
					self.abs.append(data['enriched_rod']['absXS'][g])
					self.Premoval.append(data['enriched_rod']['absXS'][g]+data['enriched_rod']['Ex'+str(g)][gIn]+self.Rfactor[g-1])
					self.removal.append(data['enriched_rod']['absXS'][g]+data['enriched_rod']['Ex'+str(g)][gIn])
					self.D.append(data['enriched_rod']['D'][g])


		# Gscat matrix
		for i in range(0,nBins):
			count = 1
			gcount = 1
			for g in range(0,nGrps*nGrps):
				if i < 30 or i >= nBins-30:
					self.Gscat[i,g] = data['water']['Ex'+str(gcount)][count]
				elif i < 200 or i >= nBins-200:
					self.Gscat[i,g] = data['fuel1']['Ex'+str(gcount)][count]
				else:
					self.Gscat[i,g] = data['enriched_rod']['Ex'+str(gcount)][count]
				count = count+1
				if count == nGrps+1:
					count = 1
					gcount = gcount+1



###############################################################


#end class


