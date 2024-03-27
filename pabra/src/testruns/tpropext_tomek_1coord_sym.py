# -*- coding: utf-8 -*-
'''
Created 2012-2020

@author: nbecker, tsokolowski
'''

import os
import os.path as pp
import sys
import string
import pabra.prop_ext as pe
import numpy as np
import shutil
import glob
import pabra.ifhist as ihst
import pabra.forester as fstr
import pabra.ptree_h5 as p5
import cProfile



class Setup(object):
    def __init__(self):
        #for f in glob.glob(pp.join(inputdir, './*')):
        #    shutil.copy(f, outputdir)
        args1 = {'log_int':300.0, 
                 'time_bins':np.arange(0.,TOTALTIME/N_TRAJ+BINTIME,BINTIME),
                 'uniform':True,
                 'crossings_per_bin':0.25,
                 'cross_int':1., # not meaningful here; defined in inputAI
                 'coord_funs': (lambda t,x: (max(x[0], x[2]) + max(x[1],x[3]))/x[4],), 
                 'interfaces': ( np.concatenate( (np.arange(.30, .60, .05), np.arange(.60, .96, .02), np.arange(.96, 1.00, .004) )), ), 
                 'n_max':5,
                 'size_limits':(1., 32.),
                 'weight_limits':(1e-4, 5.),
                 'flux_directions': [1,1],
                 'traverse': 'breadth' # 
                 }
        self.ifh1 =  ihst.IFHistogram(**args1)
        # required by gillespie.X
        os.chdir(outputdir)
        print "changed working dir to {0}".format(os.getcwd())
        self.prop = pe.TIntegrator(outputdir + runfilename, 
                                   working_dir=outputdir)
#        self.prop.init_from_ifhist(self.ifh1)
        self.fr = fstr.Forester(p5.PforestExt(  #@UndefinedVariable
                                        shadow_group=p5.h5.File( 
                                        pp.join(outputdir,'./tpropext.h5')), 
                                        root_dir=outputdir), 
                                        self.prop,
                                        self.ifh1,
                                        initial_gen=np.ones(12+40, dtype=np.dtype('float64')), 
                                        bunch_size=1)


def set_variables(input_basedir, output_basedir, pars, totaltime, n_traj, bintime, comment):
    global TOTALTIME
    global N_TRAJ
    global BINTIME
    global runfilename
    global outputdir
    global inputdir

    TOTALTIME = string.atof(totaltime)
    N_TRAJ    = string.atof(n_traj)
    BINTIME   = string.atof(bintime)
    runfilename = 'runfile.' + pars
    inputdir    = input_basedir  #  + 'input.' + pars + '/'
    outputdir   = output_basedir # + 'test_1coord_' + comment + '_' + pars +'/'
    print('Input directory  = ' + inputdir)
    print('Output directory = ' + outputdir)


def cleanupDir():
    total_files = os.listdir(outputdir)
    for f in total_files:
        if not f.startswith('safe'):
            shutil.rmtree(f, ignore_errors=1)
            # and again for the files , d'oh.
            try: os.remove(f)
            except: pass    


if __name__ == '__main__':
    # logging
    import pabra
    # what to record 
    pabra.logger.setLevel(pabra.logging.DEBUG)
    # what subset to display
    pabra.console_handler.setLevel(pabra.logging.INFO)

    # Set output paths and total time
    input_basedir=sys.argv[1]
    output_basedir=sys.argv[2]
    pars=sys.argv[3]
    totaltime=sys.argv[4]
    n_traj=sys.argv[5]
    bintime=sys.argv[6]
    comment=sys.argv[7]
    set_variables(input_basedir, output_basedir, pars, totaltime, n_traj, bintime, comment)
    
    # Start action
    #os.mkdir(outputdir)
    # write everything to file
    file_handler = pabra.logging.FileHandler(pp.join(outputdir,'./tpropext.log'),'w')
    file_handler.setLevel(pabra.logging.DEBUG)
    file_handler.setFormatter(pabra.std_fmt)
    pabra.logger.addHandler(file_handler)

    # setup
    s = Setup()
    print 'the strategy:\n', s.fr.strategy.str_rich()
    s.fr.strategy.school_n_max = 64
    cursimtime = 0
    while s.fr.strategy.sim_time < TOTALTIME:
        s.prop.current_seed = np.random.randint(1001)
        s.fr.run()
        print 'total time: ', s.fr.strategy.sim_time
        print 'forest time: ', s.fr.forest.length
        print 's.fr.strategy.sim_tree_count: ', s.fr.strategy.sim_tree_count
        print 'forest population: ', s.fr.forest.population 
    #cProfile.run('for i in range(5): s.fr.run()', sort='cumulative')

    #print 'forest size: ', s.fr.forest.size
    #print 'forest height: ', s.fr.forest.height

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

    #cleanupDir()
    
