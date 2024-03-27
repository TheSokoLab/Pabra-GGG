'''
Created on Jul 7, 2010

@author: nbecker

Test and debug the IFHistogram branching strategy.



'''

import pabra.forester as fstr
import pabra.ifhist as ifht
import pabra.propagate as prop
import numpy as np
import pabra.plot as ppt #@UnusedImport
import matplotlib.pyplot as plt
import cProfile
import copy as cp
import h5py as h5
import os
import subprocess


class Setup1d(object):
    def __init__(self, files):
        #----------------------------------------------------------------- used
        args1 = {'log_int':.02, 'time_bins':np.arange(0.,10.5,.5),
        'cross_int':.005,
        'coord_funs':((lambda t, x: x),), 
        'interfaces': (np.linspace(-0., 5, 21),),
        'n_max':5,
        'school_n_max':150}
        self.ifh1 =  ifht.IFHistogram(**args1)
        self.x_abs_prop = prop.cross_wrap(prop.ItoEulerFwd)(.001, 
                                        force_field=-0.1, D=.5,
                                        abs_cond=(lambda t, x: bool(x < 0.)))
        self.x_abs_prop2 = cp.deepcopy(self.x_abs_prop)
        self.x_abs_ref_prop = prop.cross_wrap(prop.ItoEulerFwd)(.001,
                                 force_field=-.1, D=.5,
                                 abs_cond=(lambda t, x: bool(x < 0.)),
                                 ref_box=[-np.inf, 3.])
        print 'x_abs_prop:', self.x_abs_prop
        # the forester:
        self.fr_abs = fstr.Forester(fstr.Pforest(files[0]), self.x_abs_prop,
                                    self.ifh1, initial_gen=1., bunch_size=1)
        self.fr_abs_ref = fstr.Forester(fstr.Pforest(files[1]), 
                                        self.x_abs_ref_prop,
                                    self.ifh1, initial_gen=1., bunch_size=1)
        
        #------------------------------------------------------------- not used
        self.dummyprop = prop.SystemIntegrator()
        unibins = tuple(np.arange(0., 10, .5)) + (10.,)
        nonunibins = (list(np.arange(0., 4, .2)) + list(np.arange(4, 10, .6))
                       + [10.])
        x_ifs = np.array([-2, -1, 0, 1, 2])
        x_ifs = np.arange(-2, 2, .5)
        self.ifh0 = ifht.IFHistogram(log_int=.02, 
                                      cross_int=.04,
                                      time_bins=unibins,
                                      coord_funs=((lambda t, x: x),),
                                      interfaces=(x_ifs,),
                                      flux_directions=[1,1],
                                      n_max=5, school_n_max=150)
        self.x_abs_double_prop = prop.cross_wrap(prop.ItoEulerFwd)(.001,
                                 force_field=np.zeros(1), D=.5,
                                 abs_cond=(lambda t, x: bool(x < 0 or x > 3)))

        print 'x_abs_prop.abs_cond(2,-2) ', self.x_abs_prop.abs_cond(2,-2)
        print ('fr_abs.propagator.abs_cond(2,-2) ', 
                self.fr_abs.propagator.abs_cond(2,-2) )
        def sinesweeper(t, x):
            return  -1.2 * (x - 2. * np.sin(1.2 * t))
        self.x_ito_sin = prop.cross_wrap(prop.ItoEulerFwd)(.001,
                                                    force_field=sinesweeper)
        self.bf = fstr.BruteForce(log_int=.01, t_final=10.)
        self.fr_sin = fstr.Forester(fstr.Pforest(files[2]), self.x_ito_sin,
                                      self.ifh0, self.gen(-3, 3),
                                       bunch_size=1)
        self.fr_bf_abs = fstr.Forester(fstr.Pforest(files[3]),
                                         self.x_abs_prop2, self.bf,
                                         self.gen(.5, .5), bunch_size=1)
        #------------------------------------------------------------------ end
        # combination to try:

    def gen(self, low, high):
        return lambda : np.random.uniform(low=low, high=high, size=1) #1d array

def co(n):
    def cofn(t, pt):
        return pt[n]
    return cofn



global TOTALTIME
TOTALTIME = 4.15
FR_TYPE = 'ref'

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
    file_handler = pabra.logging.FileHandler('/tmp/tifhist.log','w')
    file_handler.setLevel(pabra.logging.DEBUG)
    file_handler.setFormatter(pabra.std_fmt)
    pabra.logger.addHandler(file_handler)
    # setup
    os.chdir("/home/pabra/tmp")
    ffls = ['ff'+i+'.h5' for i in '1234']
    for i in ffls:
        try: os.unlink(i)
        except OSError: pass
    s1 = Setup1d([h5.File(i) for i in ffls])    
    fore = {'abs':s1.fr_abs, 'ref':s1.fr_abs_ref}[FR_TYPE]
    print 'the strategy:\n', fore.strategy.str_rich()
    fore.bunch_size=1
    fore.strategy.school_n_max = 84
    cursimtime = 0
    while fore.strategy.sim_time < TOTALTIME:
        fore.run()
        print 'total time: ', fore.strategy.sim_time
        print 'forest time: ', fore.forest.length
        print 'fore.strategy.sim_tree_count: ', fore.strategy.sim_tree_count
        print 'forest population: ', fore.forest.population 
    cProfile.run('for i in range(5): fore.run()', sort='cumulative')

    print 'forest size: ', fore.forest.size
    print 'forest height: ', fore.forest.height
#    print 'flux histogram:\n', fore.strategy.fluxes[0,:,:]
    #
    # plot something
    #
    axes = plt.figure().add_subplot(111)
    
#    sel = fore.forest.trees[187]
#    print 'length, height:', sel.length, ',', sel.height
#    print 'total no of branches:', sel.size
#    mo1 = sum(tr.branch_no -1 for tr in sel if tr.branch_no > 1)
#    print 'no of extra child branches:', mo1
#    sel.dump(elems=('base', 'tip'), short=True)
#
    forest_selection = fore.forest.trees[-20:]
#    for tr in forest_selection:
#        tr.dump(elems=('base', 'tip'), short=True)
    for tr in forest_selection:
        axes.plot_time_space(tr, wt_alpha=True, max_size=200)
    for tr in forest_selection:
        axes.plot_bp_time_space(tr)
    print 'axes', axes
    axes.set_xlabel('time')
    axes.set_ylabel('x')
    print 'forest size', fore.forest.size
    print 'total simulated time', fore.strategy.sim_time
    plt.show()

#==============================================================================
# RESULTS
#==============================================================================

# this test runs fully with:

#git rev: 3553cdd35f9fbdd4c340b41cfa1d828abd80af43
#
#numpy 1.6.1; h5py 2.0.1; HDF5 1.8.8
#


#         362973 function calls (362953 primitive calls) in 0.769 seconds
#
#   Ordered by: cumulative time
#
#   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
#        1    0.000    0.000    0.769    0.769 <string>:1(<module>)
#        5    0.000    0.000    0.769    0.154 forester.py:184(run)
#        5    0.000    0.000    0.645    0.129 forester.py:169(cont_run)
#        5    0.007    0.001    0.645    0.129 ifhist.py:641(propagate_school)
#      167    0.003    0.000    0.303    0.002 ifhist.py:595(_branch_on_crossing)
#      167    0.002    0.000    0.165    0.001 ptree_h5.py:916(branch)
#      168    0.005    0.000    0.157    0.001 ptree_h5.py:1006(grow_trunk)
#      395    0.000    0.000    0.137    0.000 propagate.py:485(int)
#      395    0.005    0.000    0.137    0.000 propagate.py:463(_prop_abs_cross)
#        5    0.000    0.000    0.112    0.022 forester.py:174(finish_run)
#        5    0.001    0.000    0.108    0.022 ifhist.py:718(update_from_school)
#      173    0.010    0.000    0.105    0.001 ptree_h5.py:371(flush)
#     1409    0.005    0.000    0.101    0.000 base.py:212(get)