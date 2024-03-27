'''
Created on Jan 13, 2012

@author: nbecker
'''

import os
import os.path as pp
import pabra.prop_ext as pe
import numpy as np
import shutil
import glob
import pabra.ifhist as ihst
import pabra.forester as fstr
import pabra.ptree_h5 as p5
import cProfile


testdir = pp.join(pp.dirname(__file__), '../../sandbox/prop_ext/')

print testdir

global TOTALTIME
TOTALTIME = 1430


class Setup(object):
    def __init__(self):
        for f in glob.glob(testdir + '/safe/*'):
            shutil.copy(f, testdir)
        args1 = {'log_int':np.inf, 
                 'time_bins':np.arange(0.,200.,20.),
                 'cross_int':1., # not meaningful here; defined in inputAI
                 'coord_funs':(lambda t,x: x[-1],),
                 'interfaces': (np.arange(.45, .75, .01),), 
                 'n_max':5,
                 'weight_limits':(1e-5, 5.),
                 'flux_directions': [1,0] # 
                 }
        self.ifh1 =  ihst.IFHistogram(**args1)
        # required by gillespie.X
        os.chdir(testdir)
        print "changed working dir to {0}".format(os.getcwd())
        self.prop = pe.TIntegrator(testdir + 'runfile.10.3.10.14', 
                                   working_dir=testdir)
#        self.prop.init_from_ifhist(self.ifh1)
        self.fr = fstr.Forester(p5.PforestExt(  #@UndefinedVariable
                                        shadow_group=p5.h5.File( 
                                        pp.join(testdir,'../tpropext.h5')), 
                                        root_dir=testdir), 
                                self.prop, 
                                self.ifh1, 
                                initial_gen=[1.,1.,1.,1.,1.,1.], 
                                bunch_size=1)

def cleanupDir():
    total_files = os.listdir(testdir)
    for f in total_files:
        if not f.startswith('safe'):
            shutil.rmtree(pp.join(testdir,f), ignore_errors=1)
            # and again for the files , d'oh.
            try: os.remove(pp.join(testdir,f))
            except OSError: pass
    try: os.remove(pp.join(testdir, '../tpropext.h5'))
    except OSError: pass

if __name__ == '__main__':
    cleanupDir()
    # logging
    import pabra
    # what to record 
    pabra.logger.setLevel(pabra.logging.DEBUG)
    # what subset to display
    pabra.console_handler.setLevel(pabra.logging.INFO)
    # write everything to file
    file_handler = pabra.logging.FileHandler('/tmp/tpropext.log','w')
    file_handler.setLevel(pabra.logging.DEBUG)
    file_handler.setFormatter(pabra.std_fmt)
    pabra.logger.addHandler(file_handler)
    # setup
    s = Setup()
    print 'the strategy:\n', s.fr.strategy.str_rich()
    s.fr.strategy.school_n_max = 84
    cursimtime = 0
    while s.fr.strategy.sim_time < TOTALTIME:
        s.prop.current_seed = np.random.randint(1001)
        s.fr.run()
        print 'total time: ', s.fr.strategy.sim_time
        print 'forest time: ', s.fr.forest.length
        print 's.fr.strategy.sim_tree_count: ', s.fr.strategy.sim_tree_count
        print 'forest population: ', s.fr.forest.population 
    cProfile.run('for i in range(5): s.fr.run()', sort='cumulative')

    print 'forest size: ', s.fr.forest.size
    print 'forest height: ', s.fr.forest.height

#    #
#    # plot something
#    #
#    axes = plt.figure().add_subplot(111)
#    
##    sel = fore.forest.trees[187]
##    print 'length, height:', sel.length, ',', sel.height
##    print 'total no of branches:', sel.size
##    mo1 = sum(tr.branch_no -1 for tr in sel if tr.branch_no > 1)
##    print 'no of extra child branches:', mo1
##    sel.dump(elems=('base', 'tip'), short=True)
##
#    forest_selection = fore.forest.trees[-20:]
##    for tr in forest_selection:
##        tr.dump(elems=('base', 'tip'), short=True)
#    for tr in forest_selection:
#        axes.plot_time_space(tr, wt_alpha=True, max_size=200)
#    for tr in forest_selection:
#        axes.plot_bp_time_space(tr)
#    print 'axes', axes
#    axes.set_xlabel('time')
#    axes.set_ylabel('x')
#    print 'forest size', fore.forest.size
#    print 'total simulated time', fore.strategy.sim_time
#    plt.show()

    cleanupDir()
    