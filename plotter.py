#class to plot

import numpy as np
import matplotlib.pyplot as plt


class Plotter:

####################################################################
    def __init__(self):
    	self

####################################################################

    def plot_DR(self, opts, err, DR, asymDR):
        print("plotting DR")

        # plt.yscale('log')
        plt.plot(DR)
        plt.savefig('./DR'+str(opts.pb_num))
        plt.clf()
        plt.yscale('log')
        plt.plot(err)
        plt.savefig('./ERR'+str(opts.pb_num))


####################################################################

    def plot_flux(self, opts, flux, A, F):

        import matplotlib.pylab as pl
        import copy
        # import scipy.sparse as sps
        # import betterspy


        plt.clf()
    	numBins = opts.numBins
    	numGrps = opts.numGroups
        BinLen = opts.length/numBins
        x = np.linspace(BinLen/2,opts.length+BinLen/2,numBins,endpoint=False)
        pb_num = str(opts.pb_num)

        groups = []
        for g in range(1,numGrps+1):
            l = 0
            x_group = np.zeros(numBins)
            for i in range(numBins*(g-1),numBins*g):
                x_group[l] = flux[i]
                l = l+1

            groups.append('group%i' %g)
            colors = ['r','b','g','o','n','p']
            groups[g-1], = plt.plot(x, x_group, colors[g-1], label=groups[g-1])
            plt.ylabel('Normalized flux')
            plt.xlabel('x [cm]')
            plt.legend(groups)
            plt.savefig('./flux'+pb_num)


        plt.clf()
        pl.spy(A)
        pl.savefig('./spy/spyA'+pb_num)
        pl.clf()
        pl.spy(F)
        pl.savefig('./spy/spyF'+pb_num)
        pl.clf()



####################################################################

    def plot_error(self, opts, eig1, eig2, eig3, eig4,
        norm1, norm2, norm3, norm4, 
        x1, x2, x3, x4):

        import copy
        
        pb_num = str(opts.pb_num)
        pb_type = ['node', 'mesh']

        plt.clf()
        pb_type[0], = plt.plot(x1, eig1, 'r', label='node')
        pb_type[1], = plt.plot(x2, eig2, 'b', label='mesh')
        plt.xlabel('Iteration number')
        plt.ylabel('Eigenvalue')
        plt.legend(pb_type)
        plt.savefig('./eigNOrod')

        plt.clf()
        pb_type[0], = plt.plot(x1, norm1, 'r', label='node')
        pb_type[1], = plt.plot(x2, norm2, 'b', label='mesh')
        plt.xlabel('Iteration number')
        plt.ylabel('L2 norm')
        plt.legend(pb_type)
        plt.savefig('./normNOrod')
        
        plt.clf()
        pb_type[0], = plt.plot(x3, eig3, 'r', label='node')
        pb_type[1], = plt.plot(x4, eig4, 'b', label='mesh')
        plt.xlabel('Iteration number')
        plt.ylabel('Eigenvalue')
        plt.legend(pb_type)
        plt.savefig('./eigWITHrod')

        plt.clf()
        pb_type[0], = plt.plot(x3, norm3, 'r', label='node')
        pb_type[1], = plt.plot(x4, norm4, 'b', label='mesh')
        plt.xlabel('Iteration number')
        plt.ylabel('L2 norm')
        plt.legend(pb_type)
        plt.savefig('./normWITHrod')

        plt.clf()


