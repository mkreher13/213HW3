# 22.213 Pset3
# This code solves the discrete diffusion equations in a transient
# Last edited by Miriam Kreher on 03/06/2020

import numpy as np 
import copy
from getfilename import *
from diff_opts import *
from material import *
from cross_sections import *
from construct import *
from solverSS import *
from solverTD import *
from plotter import *

fn = FileName()
fn.get_filename()

for f in fn.listfn:
    print f

    options = DiffusionOpts1D()
    options.read(f)
    nBins = options.numBins
    nGrps = options.numGroups

    M = Material()
    M.read()
    # print(M.data)
    # print(M.materialList)

    Matrix=Construct()
    if options.pb_num == 1:
        XS = CrossSections(options)
    	XS.problem1(options, M.data)
        Matrix.constructSS(options, XS.D, XS.removal, XS.Gscat, XS.fis)
        sol = SolveSS(options)
    elif options.pb_num == 2:
        XS = CrossSections(options, kcrit)
    	XS.problem2(options, M.data)
        Matrix.constructTD(options, XS.D, XS.Premoval, XS.Gscat, XS.Pfis)
        sol = SolveTD(options, InitFlux, InitC, kcrit)

    # print(XS.abs)
    # print(XS.fis)
    # print(XS.Gscat)
    # print Matrix.A
    # print Matrix.F
    
    if options.method == 'PointJacobi':
    	sol.PointJacobi(options, Matrix, XS)
    elif options.method == 'GaussSeidel':
    	sol.GaussSeidel(options, Matrix, XS)
    else:
    	Matrix.invert(Matrix.A)
    	sol.NumpySolve(options, Matrix.inv, Matrix.F, XS)
    #print sol.flux

    results = Plotter()
    results.plot_flux(options, sol.flux, Matrix.A, Matrix.F)


    if options.pb_num == 1:
        print "Eigenvalue:", sol.k
        print "Source iterations:", sol.source_it

        # Store and calculate initial values for transient problem
        InitFlux = copy.copy(sol.flux)
        kcrit = copy.copy(sol.k)
        RRfis = np.zeros([nBins])
        for n in range(nBins):
            RRfis[n] = XS.fis[n]*InitFlux[n]+XS.fis[n+nBins]*InitFlux[n+nBins]
        PrecursorFactor = XS.Bi[:]/(XS.Decayi[:]*kcrit)

        InitC = np.zeros([8,nBins])
        for i in range(8):
            InitC[i,:] = PrecursorFactor[i]*RRfis
            plt.ylim([0,0.000005])
            plt.plot(InitC[i])
        plt.savefig('./PrecursorConcentrations'+str(0))
        plt.clf()

    plt.clf()
    plt.plot(InitFlux)
    plt.plot(sol.flux)
    plt.savefig('./CompFlux')
    plt.legend(['Initial','Final'])
    plt.clf()

        # To show the magnitude of the precursors:
        # print(InitC[:,100])

    
    # elif options.pb_num == 2:
    #     Residual = np.zeros(len(sol.flux))
    #     Residual[:] = abs(InitFlux[:]/sol.flux[:])
    #     plt.clf()
    #     plt.plot(Residual)
    #     plt.savefig('./Residual')
    #     plt.clf()



    













