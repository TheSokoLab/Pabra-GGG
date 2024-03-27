'''
Created on Jul 14, 2010

@author: nbecker


test interface sampling of a gillespie run.

since direct optimized checking of interface crossing is not 
available as of yet, we try the CrossWrapping to sample interface transitions

this is an all-in-one test including gillespie run, analysis, graphing.

'''


# copy and paste the setup code from IF and gillespie scripts.


import pabra.forester as fstr
import pabra.ifhist as ifht
import pabra.propagate as prop
import numpy as np
import pabra.plot as ppt
import matplotlib.pyplot as plt
import cProfile
import pabra.gillespieCy as gcy
from pabra.gillespieCy import Reaction, Gillespie  #@UnresolvedImport
import itertools as it
import copy as cp
import time
import gc
import os
import h5py as h5
import subprocess


#-------------------------------------------------------------------- utilities

def co(n):
    def cofn(t, x):
        return x[n]
    return cofn 

#------------------------------------------------------------------------ setup

class Setup(object):
    def __init__(self, forestfiles):
        #----------------------------------------------------------------- used
        # strategy
        log_ifs = np.exp(np.linspace(np.log(1.9), np.log(55.5), 9))
        lin_ifs = np.arange(-0., 55, 5)
        manual_ifs = np.cumsum([2, 3, 4, 5, 6, 6, 7]) + .5
        no_ifs = np.array([])
        args1 = {'log_int':1., 'time_bins':np.arange(0., 1020., 20),
        'cross_int':.3, 
        'coord_funs':(co(0),), 
        'interfaces': (manual_ifs,),
        'n_max':5,
        'school_n_max':np.inf,
        'size_target':100}
        self.ifh1 =  ifht.IFHistogram(**args1)
#        # the brute force alternative
        args2 = args1.copy()
        args2['interfaces'] = (no_ifs,)
        args2['cross_int'] = 10000.        
        self.ifh2 = ifht.IFHistogram(**args2)
        # reactions
        self.r1 = Reaction(['mass_action', 1.], [], ['a'], [],  
                               name='0->a')
        self.r2 = Reaction(['mass_action', .05], ['a'], [], ['a'], 
                               name='a->0')
        self.r3 = Reaction(['mass_action', 8.], [], ['b'], [],  
                               name='0->b')
        self.r4 = Reaction(['mass_action', .1], ['b'], [], ['b'], 
                               name='b->0')
        self.r5 = Reaction(['mass_action', .2], ['a'], ['b'], ['a'], 
                               name='a->b')
        self.r6 = Reaction(['mass_action', .05], ['a', 'b'], [], ['a', 'b'], 
                               name='a+b->0')
        # now add a crossing checking wrapper around Gillespie:
        self.XGillespie = prop.cross_wrap(Gillespie)
        self.g12 = self.XGillespie(['a'], [self.r1, self.r2])
        self.g12346 = self.XGillespie(['a', 'b'], [self.r1, self.r2, self.r3,
                                              self.r4, self.r6])
        self.g1to6 = self.XGillespie(['a', 'b'], [self.r1, self.r2, self.r3,
                                            self.r4, self.r5, self.r6])
        # the forester:
        self.fr_12 = fstr.Forester(fstr.Pforest(forestfiles[0]), self.g12, 
                                   self.ifh1, 
                                   initial_gen=np.array([0]), bunch_size=1)
        # the brute forester
        # need an extra instance of the propagator, since that gets initialized
        # in the Forester. bad design?
        self.g12_1 = self.XGillespie(['a'], [self.r1, self.r2])
        self.br_12 = fstr.Forester(fstr.Pforest(forestfiles[1]), self.g12_1, 
                                   self.ifh2, 
                                   initial_gen=np.array([0]), bunch_size=1)
        
    #---------------------------------------------------------------- utilities


    def gen(self, low, high):
        return lambda : np.random.uniform(low=low, high=high, size=1) #1d array


TOTALTIME = 1.5 * 1e3
BF_COMPARE = 0

if __name__ == '__main__':
    print "git rev: {0}".format(subprocess.check_output(
                                ['git', 'rev-parse', 'HEAD']))
    print "numpy {0}; h5py {1}; HDF5 {2}".format(np.version.version, 
                                                 h5.version.version,
                                                 h5.version.hdf5_version)
    # logging
    import pabra
    # what to record 
    pabra.logger.setLevel(pabra.logging.DEBUG)
    # what subset to display
    pabra.console_handler.setLevel(pabra.logging.INFO)
    # write everything to file
    file_handler = pabra.logging.FileHandler('/tmp/tif_gil.log','w')
    file_handler.setLevel(pabra.logging.DEBUG)
    file_handler.setFormatter(pabra.std_fmt)
    pabra.logger.addHandler(file_handler) 
    # setup
    os.chdir("/home/pabra/tmp")
    ffls = ('ff1.h5', 'ff2.h5')
    for i in ffls:
        try: os.unlink(i)
        except OSError: pass
    s1 = Setup([h5.File(i) for i in ffls])
    gcy.therng.seed(3456)  #@UndefinedVariable
    f12 = s1.fr_12
    print 'f12: ', f12.strategy.str_rich()
    cursimtime = 0
    wall_old = time.time()
    while f12.strategy.sim_time < TOTALTIME:
        f12.run()
        print 'total time: ', f12.strategy.sim_time
        print 'f12.strategy.sim_tree_count: ', f12.strategy.sim_tree_count
        print 'f12.forest.height:', f12.forest.height
#    cProfile.run('fabs.run()', sort='time')
    wall_new = time.time()
    print 'forest size: ', f12.forest.size
    print 'forest height: ', f12.forest.height
    print 'valid flux histogram:\n', f12.strategy.valid_n_flx[0,:,:]
    print 'simulated time: ', f12.strategy.sim_time
    print 'elapsed time:', wall_new - wall_old
    #
    # plot something
    #
    #---------------------------------------------------------------- tree plot
    #
    tree_axes = plt.figure().add_subplot(111)
    tree_fg = tree_axes.figure
    tree_axes.set_title('path tree')
    # the first few trees bigger than 20 branches
    def forest_selection(size, number=np.inf, max_size=np.inf):
        nm, cum_size  = 0, 0
        for tr in f12.forest.trees:
            if tr.size > size and nm <= number and cum_size < max_size:
                nm += 1
                cum_size += tr.size
                yield tr
#    import code; code.interact(local=locals())

    for tr in forest_selection(20, 5, 1000):
        tree_axes.plot_time_space(tr, wt_alpha=True)
    for tr in forest_selection(1):
        tree_axes.plot_bp_time_space(tr, zorder=200)
    tree_axes.set_xlabel('time')
    tree_axes.set_ylabel('N')
    print 'tree plot done'
    # unfortunately, the color bar does not work without a ColorMap, 
    # and that does not have an alpha channel.
    #------------------------------------------------------------- scatter plot
    scat_fg, scat_axes  = plt.subplots(2, 1, sharex=1, sharey=1) 
    scat_fg.subplots_adjust(hspace=.3)
    scat_axes[0].set_title('points')
    scat_axes[1].set_title('weighted points')
    scat_axes[1].set_xlabel('time')
    for ax in scat_axes:
        ax.set_ylabel('N')
    gc.disable()   # workaround for a bug in append() which makes it slow
    lgd_t, lgd_x, lgd_wts= [], [], []
    for tr in f12.forest:
        for lgd_pt in tr.trunk:
            lgd_t.append(lgd_pt[0])
            lgd_x.append(lgd_pt[1][0])
            lgd_wts.append(tr.wt)
    gc.enable()
    # we have no automated scatter plot yet! maybe that would be a good idea?
    rgbas = np.array(map(lambda x:[1, 0, 0, 0.0 + 1. * x], 
                         plt.Normalize(vmin=0)(lgd_wts)))
    scat_axes[1].scatter(lgd_t, lgd_x, s=5, c=rgbas,
                    edgecolors='none')
    scat_axes[0].scatter(lgd_t, lgd_x, s=5, c='blue',
                    edgecolors='none')
    print 'len of the logged times:', len(lgd_t)

    #-------------------------------------- a histogram of the occurred weights
    wts = []
    for tr in f12.forest:
        try:
            wts.append(tr.wt)
        except AttributeError:
            print 'no wt!'
    wt_hist_fg, wt_hist = plt.subplots()
    wt_bins = np.exp(np.linspace(*np.log([1e-5,1e3]), num=26))
    wt_hist.hist(wts, bins=wt_bins, log=True)
    wt_hist.set_title('cumulated histogram')
    wt_hist.set_xlabel('weight')
    wt_hist.set_ylabel('sample points')
    wt_hist.set_xscale('log')
    #------------------------------------------------------- density histograms
    #
    hist_fg, (count_hist_2d, wt_hist_2d) = plt.subplots(2, 1, sharex=1,
                                                         sharey=1) 
    hist_fg.subplots_adjust(hspace=.3)
    wt_hist_2d.set_xlabel('time')
    wt_hist_2d.set_ylabel('N')
    count_hist_2d.set_ylabel('N')
    count_hist_2d.set_title('space-time histograms')
    #
    count_density, t_e, x_e = np.histogram2d(lgd_t, lgd_x, 
                                             bins=(f12.strategy.t_bins,) + 
                                             f12.strategy.if_fams, 
                                             normed=False)
    print 'count_density, t_e, x_e shapes', [a.shape for a in (count_density, 
                                                               t_e, x_e)]
    print 'len of the logged times:', len(lgd_t)
    ct_norm = plt.matplotlib.colors.LogNorm()
    ct_norm.autoscale(np.maximum(.1, count_density))
    count_hist = count_hist_2d.pcolor(t_e, x_e, count_density.T, norm=ct_norm)
    #
    wt_density, dmy, dmy = np.histogram2d(lgd_t, lgd_x, 
                                          bins=(f12.strategy.t_bins,) + 
                                          f12.strategy.if_fams,
                                          normed=False,
                                          weights = lgd_wts)
    wt_norm = plt.matplotlib.colors.LogNorm()
    wt_norm.autoscale(np.maximum(.1,wt_density))
    wt_hist = wt_hist_2d.pcolor(t_e, x_e, wt_density.T, norm=wt_norm)
    cbar_ct = hist_fg.colorbar(count_hist, ax=count_hist_2d)
    cbar_ct.set_label('raw counts' )
    cbar_wt = hist_fg.colorbar(wt_hist, ax=wt_hist_2d)
    cbar_wt.set_label('summed weight')
    #--------------------------------------------------- brute force comparison
    if BF_COMPARE:
        b12 = s1.br_12
        print 'b12: ', b12.strategy.str_rich()
        cursimtime = 0
        wall_old_b = time.time()
        while b12.strategy.sim_time < f12.strategy.sim_time:
            b12.run()
            print 'total time: ', b12.strategy.sim_time
        wall_new_b = time.time()
        #cProfile.run('fabs.run()', sort='time')
        print 'brute forest size: ', b12.forest.size
        print 'brute sample points', np.sum(1 for tr in b12.forest
                                            for pt in tr.trunk)
        print 'brute sim time: ', b12.strategy.sim_time
        print 'brute elapsed time', (wall_new_b - wall_old_b)
        print 'branched forest size', f12.forest.size
        print 'branched sample points', np.sum(1 for tr in f12.forest
                                               for pt in tr.trunk)
        print 'branched sim time: ', f12.strategy.sim_time
        print 'branched elapsed time', (wall_new - wall_old)
        print ('branching overhead in per cent of the bf time per sim time:',
               ((wall_new-wall_old)/f12.strategy.sim_time
                - (wall_new_b-wall_old_b)/b12.strategy.sim_time)
                / ((wall_new_b-wall_old_b)/b12.strategy.sim_time) 
                * 100)
        

        bgd_t, bgd_x = [], [] 
        for tr in b12.forest:
            for bgd_pt in tr.trunk:
                bgd_t.append(bgd_pt[0])
                bgd_x.append(bgd_pt[1][0])
    
        print 'length of bf pt list:', len(bgd_t)
        print 'length of branched pt list:', len(lgd_t)
        print 'sum of the branched sampled weights:', np.sum(lgd_wts)
        

        comp_fg, (hax_ct, hax_wt, hax_bf) = plt.subplots(3, 1, sharex=1, 
                                                         sharey=1)
        gridsize_ = (50,20)
        hex_norm = plt.matplotlib.colors.LogNorm(vmin=1, vmax=len(bgd_t)/50.)
        ct_hex = hax_ct.hexbin(lgd_t, lgd_x, 
                               C=np.ones_like(lgd_t), 
                               reduce_C_function=np.sum,
                               gridsize=gridsize_,
                               norm=hex_norm)
        hax_ct.set_title('space-time histograms')
        wt_hex = hax_wt.hexbin(lgd_t, lgd_x, 
                               C=lgd_wts, 
                               reduce_C_function=np.sum,
                               gridsize=gridsize_,
                               norm=hex_norm)
        bf_hex = hax_bf.hexbin(bgd_t, bgd_x, 
                                C=np.ones_like(bgd_t), 
                                reduce_C_function=np.sum,
                               gridsize=gridsize_,
                               norm=hex_norm)
        hax_bf.set_xlabel('time')
        for hex, hax, label in it.izip((ct_hex, wt_hex, bf_hex), 
                                       (hax_ct, hax_wt, hax_bf), 
                            ('raw couts', 'weighted counts', 'brute-force')):
            hax.colorbar_ = comp_fg.colorbar(hex, ax=hax)
            hax.set_ylabel('N')
            hax.colorbar_.set_label(label)
            hax.set_axis_bgcolor(plt.cm.get_cmap()(0)) 
            
        
    #-------------------------------------------------------------------- done:

    plt.show()
    
#==============================================================================
# RESULTS
#==============================================================================

# this test runs fully with:

#git rev: 3553cdd35f9fbdd4c340b41cfa1d828abd80af43
#
#numpy 1.6.1; h5py 2.0.1; HDF5 1.8.8

#simulated time:  26487.4
#elapsed time: 27.6556930542

