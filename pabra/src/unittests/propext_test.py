'''
Created on Nov 18, 2011

@author: nbecker
'''
import unittest as ut
import os
import os.path as pp
import pabra.ptree_h5 as p5
import pabra.prop_ext as pe
import numpy as np
import shutil
import glob
from pabra.propagate import Crossing
import pabra.ifhist as ihst
import pabra.forester as fstr


testdir = pp.join(pp.dirname(__file__), '../../sandbox/prop_ext')


class Test(ut.TestCase):
    def setUp(self):
        os.chdir(testdir)  # this is mandatory since gillespie.X requires it.
        for f in glob.glob('./safe/*'):
            shutil.copy(f, '.')
        self.tfr = p5.PforestExt()  #@UndefinedVariable
        self.tfr.tx_dtype = p5.construct_dtype((0.,0.))
        self.tfr.log_int = np.inf
        self.ts = self.tfr.add_school()
        self.ttr = self.ts.add_child()
        self.ttr.base = (0., 0.)
        self.ttr.tip = (0., 0.)
        self.flxdirs =  np.zeros((1, 10-1, 10, 2))  # dtjijb
        self.prop = pe.TIntegrator('runfile.10.3.10.14',
                                   coord_funs = (lambda t,x:x,),
                                   interfaces = (np.arange(-10, 10, 2),),
                                   time_bins = np.arange(0,10),
                                   flux_directions = self.flxdirs)
        self.prop2 = pe.TIntegrator('runfile.10.3.10.14', 
                                    coord_funs = (lambda t, x:x,), 
                                    interfaces = (np.arange(0, 1, .1),), 
                                    time_bins = np.arange(0, 10),
                                    flux_directions = self.flxdirs)
        args1 = {'log_int':np.inf, 
                 'time_bins':np.arange(0.,10.,1.),
                 'cross_int':1., # not meaningful here.
                 'coord_funs':(lambda t,x:x,), 
                 'interfaces': (np.arange(-10, 10, 2),), 
                 'n_max':5,
                 'school_n_max':150,
                 'flux_directions': [0,0] # no crossings valid
                 }
        self.ifh1 =  ihst.IFHistogram(**args1)
        self.fr = fstr.Forester(
                        p5.PforestExt(root_dir=testdir), #@UndefinedVariable
                        self.prop, self.ifh1, initial_gen=1., bunch_size=1)



    def tearDown(self):
        total_files = os.listdir(testdir)
        for f in total_files:
            if not f.startswith('safe'):
                shutil.rmtree(f, ignore_errors=1)
                # and again for the files , d'oh.
                try: os.remove(f)
                except: pass


    def testSchool(self):
        self.assertEqual(self.ts.root_dir, self.tfr.node_dir)
        self.assertTrue(pp.isfile(pp.join(self.ts.node_dir,'inputAI.end')) )

    def testTree(self):
        self.assertEqual(self.ttr.parent.node_dir, self.ts.node_dir)
        self.assertTrue(pp.isfile(pp.join(self.ttr.node_dir,'inputAI')) )

    def testProp(self):
        self.prop.current_dir = self.ttr.node_dir
        self.assertAlmostEqual(self.prop._prop_or_abs(0)[0], 1., 2)

    def testGrowth(self):
        self.ttr.grow_trunk(self.prop, 2.)
        self.assertAlmostEqual(self.ttr.tip[0], 2., 2)

    def testCross(self):
        self.assertRaises(Crossing, self.ttr.grow_trunk, *(self.prop2, 2.))

    def testForest(self):
        self.fr.strategy.school_n_max = 84
        totsimtime = 10
        while self.fr.strategy.sim_time < totsimtime:
            self.fr.run()
            print 'total time: ', self.fr.strategy.sim_time
            print 'forest time: ', self.fr.forest.length
        self.assertGreater(self.fr.strategy.sim_time, totsimtime)
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    ut.main()
