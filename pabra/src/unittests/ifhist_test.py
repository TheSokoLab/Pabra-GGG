'''
Created on May 7, 2010

@author: nbecker

Collect tests of the forester module: Pforests and branching strategies
'''


import unittest as ut
from pabra import ptree_h5 as p5
from pabra.propagate import *
from pabra.forester import *
from pabra.ifhist import *
import pabra
import cProfile

pabra.logger.setLevel(pabra.logging.WARNING)
pabra.console_handler.setLevel(pabra.logging.WARNING)

import time

p5.DEFAULT_FILE_IN_MEM = 0

class TestForester(ut.TestCase):
    """test branching strategies"""

    def setUp(self):
        self.dtype_1 = p5.construct_dtype((1., 1.))
        self.base_f1 = Pforest(tx_dtype=self.dtype_1)
        self.onegen = self.base_f1.add_child()
        self.onegen.add_child(0.1, (0., 1.))  # logwt and base
        self.onegen.add_child(0.4, (0., -1))
        self.base_f1.shadow.copy(self.onegen.shadow, self.base_f1.shadow,
                                 name='t1')
        self.onegen_dup = p5.Ptree(self.base_f1.shadow['/t1'])
        self.absprop = SystemIntegrator((lambda t, x: bool(t > 1.)))
        self.ito_no_f = ItoEulerFwd(.01, 0.)
        self.ito_sin_sweep = ItoEulerFwd(.01, force_field=self.sinesweeper)
        self.ito_harmonic_2d = ItoEulerFwd(.01, force_field=self.harmonic_2d)
        self.dtype_2 = p5.construct_dtype((1., [1., 1.]))
        self.base_f2 = Pforest(tx_dtype=self.dtype_2, log_int=.1)
        self.onegen2d = self.base_f2.add_child(base=(0., np.zeros((2,))))
        for wt in - 1., 1.:
            self.onegen2d.\
            add_child(wt, base=(0., np.zeros((2,))))


    # class-level methods; these are also available
    def x_coord(self, t, x):
        return x
    def sinesweeper(self, t, x):
        return - 1.2 * (x - 2. * math.sin(1.2 * t))
    def harmonic_2d(self, t, pt):
        return - pt
    def co(self, n):
        return lambda t, x : x[n]
    def init_gen(self, n):
        return lambda: np.random.uniform(size=(n,))

#------------------------------------------------------------- the actual tests

    @ut.skipIf(0, 'no reason')
    def test_if_setup(self):
        ifh = IFHistogram(log_int=.1,
                         time_bins=np.arange(0., 10, .5),
                         coord_funs=(self.x_coord,),
                         interfaces=(np.linspace(-5., 5, 11),),
                         flux_directions=None,
                         uniform=False,
                         n_max=5)
        self.assertTrue(ifh.rs_betas is rs_betas_default)
        self.assertTrue( ifh.n_flx.shape == (1, 19, 11, 2))
        self.assertTrue( ifh.p_flx.shape == (1, 19, 11, 2))
        self.assertFalse( ifh.n_flx.any())
        self.assertFalse( ifh.p_flx.any())
        self.assertTrue(  ifh.school_counter <= 0 )
        # this is not checked anymore: c_funs are a thing of the propagator
#        self.assertRaises(BinningError, IFHistogram,
#                      **{'log_int':.1,
#                         'time_bins':(0, 1, 2),
#                         'coord_funs':(self.x_coord,),
#                         'interfaces':((0, 1, 2), (1, 2))})
        # the volume of the patch with direction 0, time-bin 2, interface no. 3
        self.assertAlmostEqual(ifh._patch_vol(0, 2, 3), .5 , 3)
        print ifh._patch_vol(0, 2, 2)
        self.assertRaises(Exception, ifh._patch_vol, *(1, 2, 2))
#        assert_raises(Exception, IFHistogram, 
#                      **{'log_int':.1,
#                         'coord_funs':(self.x_coord,), 
#                         'interfaces':(np.linspace(-5.,5,11),), 
#                         'flux_directions':[3,4]} )


    @ut.skipIf(0, 'no reason')
    def test_if_run(self):
        self.base_f1.dump()
        f1 = self.base_f1
        h5file = f1.shadow
        tlist = (self.onegen, self.onegen_dup)
        args = {'log_int':.1,
                'time_bins':np.arange(0., 10.5, .5),
                'coord_funs':(self.x_coord,),
                'interfaces': (np.linspace(-5, 5, 11),),
                'cross_int':.01,
                'n_max':5}
        ifh = IFHistogram(**args)
#        # now run it a bit.
#        # put a few trees in school
        # make a new school and register it
        f1.log_int = .1
        sc = f1.add_school()
        ifh._init_school(sc)
        self.assertTrue(ifh.school_counter == 0)
        for tr in tlist:
            self.assertIsInstance(tr, p5.PtreeBase)
            h5file.copy(tr.shadow, sc.shadow)
        sc.dump()
        print 'sdf'
        self.assertEqual(len(f1.children[0].trees), 2)
        print 'leaves', list(sc.leaves())
        therng.seed(12)
        for tr in f1.children[0].trees:
            tr.grow_tips(self.ito_sin_sweep, 7.)
        #print [lf.tip for lf in  fh.school.leaves()]
        for lf in f1.children[0].leaves():
            self.assertAlmostEqual(lf.tip[0], 7, 2)
        self.assertEqual( sum(1 for dmy in f1.children[0].leaves()), 4)
        abranch = tuple(f1.children[0].leaves())[0]
        print list(abranch._group_path_from_root())
        print list(gr.get('trunk') for gr in abranch._group_path_from_root())
        print 'trunk from root length', len(abranch.trunk_from_root)
        print 'tip', abranch.tip
#        print 'trunk', abranch.trunk
        abranch.dump()
        print 'interface families:', ifh.if_fams
        np.random.seed(123)
        for tr in f1.children[0].trees:
            tr.grow_tips(self.ito_sin_sweep, 9.)
        print 'length', sc.length
        # no crossing yet
        self.assertFalse(np.any(sc._hist_get('n_flx')))
        sc.cap_tips()
        ifh.update_from_school(sc)
        # another school.
        sc = f1.add_school()
        ifh._init_school(sc)
        print 'new size', f1.children[1].size
        f1.children[1].dump()
        self.assertEqual(f1.children[1].size, 1)
        h5file.copy(self.onegen.shadow, f1.children[1].shadow)
        # dangerous: now tr cannot be flushed!
        tr = f1.children[1].trees[0]
        np.random.seed(1230)
        tr.grow_tips(self.ito_sin_sweep, .4)
        self.assertTrue(tr.alive)
        self.assertEqual(tr.leaf_no(), 2)
#        print tr.size
        np.random.seed(23)
        self.assertTrue(tr.alive)
        ItoCross = cross_wrap(ItoEulerFwd)
        self.assertEqual(ItoCross.__name__, 'CrossItoEulerFwd')
        self.assertTrue(issubclass(ItoCross, ItoEulerFwd))
        ssweep_cross = ItoCross(delta=.1, force_field=self.sinesweeper)
        ssweep_cross.init_from_ifhist(ifh)
        self.assertEqual(ssweep_cross.cross_int, .01)
        self.assertEqual(ssweep_cross.c_funs, ifh.c_funs)
        self.assertTrue(f1.children[1].alive)
        print f1.children[1].leaves().next().tip
        # no danger, since this is a forest
        sc1 = f1.children[1]
        ifh._init_school(sc1)
        ifh.propagate_school(sc1, ssweep_cross, school_n_max=13)
        self.assertFalse(tr.leaf_no())
        self.assertFalse(list(tr.postorder_leaves()))
        #print tr.size
        self.assertFalse(tr.alive)
        #assert_almost_equal(tr.height, 10, 3)
        self.assertFalse(list(sc1.leaves()))
        ifh.update_from_school(sc1)
        # not: one more from increasing the new_tree_ct somewhere above.
        self.assertEqual(ifh.sim_tree_count, 1)
        st = ifh.sim_time
        ifh.propagate_school(sc, ssweep_cross, school_n_max=134)
        self.assertEqual(st, ifh.sim_time)
        ifh.update_from_school(sc)
        self.assertGreater(ifh.sim_time,st)
    
    @ut.skipIf(0, 'no reason')
    def test_2d(self):
        t0 = time.time()
        tforest = self.base_f2
        print 'dtype:', tforest.tx_dtype
        self.base_f2.shadow.copy(self.onegen2d.shadow, self.base_f2.shadow,
                                 name='1')
        print 'dtype:', tforest.tx_dtype
        args = {'log_int':.1, 'time_bins':np.arange(0., 10.5, .5),
                'cross_int':.02,
                'coord_funs':(self.co(0),),
                'interfaces': (np.linspace(-5, 5, 11),),
                'n_max':5,
                'school_n_max':19}
        ifh = IFHistogram(**args)
        harmonic_cross = cross_wrap(ItoEulerFwd)(.01,
                                                 force_field=self.harmonic_2d)
        frst1 = Forester(tforest, harmonic_cross, ifh, self.init_gen(2))
        print 'dtype', tforest.tx_dtype
        np.random.seed(12)
        repeats = 3
        print 'init gen', frst1.init_gen()
        for dmy in range(repeats):
            frst1.run()
#        print frst1.strategy.fluxes[0,:,:]
        args2 = {'log_int':.2,
                 'time_bins':np.arange(0., 4.5, .5),
                 'cross_int':.1,
                 'coord_funs':(self.co(0),
                              lambda t, pt: .01 * t + pt[0] - pt[1]),
                 'interfaces': (np.linspace(-5, 5, 11), np.linspace(-3, 3, 7)),
                 'n_max':5,
                 'school_n_max':9,
                 'flux_directions': [True, False]}
        ifh2 = IFHistogram(**args2)
        h_c2 = cross_wrap(ItoEulerFwd)(.01, force_field=self.harmonic_2d)
        frst2 = Forester(tforest, h_c2, ifh2, self.init_gen(2))
        np.random.seed(12)
        repeats = 4
        for dmy in xrange(repeats):
            frst2.run()
        self.assertEqual(frst2.strategy.n_flx.shape, (2, 8, 11, 7, 2))
        print 'time elapsed:', time.time() - t0
#        print '2-d, 2-coord shape:', frst2.strategy.fluxes.shape
#        print 'fluxes in the 0 direction'
#        print frst2.strategy.fluxes[0].sum(axis=2)
#        print 'fluxes in the 1 direction'
#        print frst2.strategy.fluxes[1].sum(axis=1)
#        print 'weighted flux in the 0 direction'
#        print frst2.strategy.wts[0].sum(axis=2)


    @ut.skipIf(0, 'no reason')
    def test_forester(self):
        tforest = self.base_f1
        args = {'log_int':.1,
                'time_bins':np.arange(0., 10.5, .5),
                'cross_int':.05,
                'coord_funs':(self.x_coord,),
                'interfaces': (np.linspace(-5, 5, 11),),
                'n_max':5,
                'school_n_max':19}
        ifh = IFHistogram(**args)
        def ig():
            return - .4 + 8 * np.random.random()
        ssweep_cross = cross_wrap(ItoEulerFwd)(.01,
                                               force_field=self.sinesweeper)
        frstr = Forester(tforest, ssweep_cross, ifh, ig, bunch_size=3)
        self.assertEqual(ssweep_cross.c_funs, ifh.c_funs)
        self.assertEqual(ssweep_cross.cross_int, .05)
        self.assertEqual(frstr.strategy.log_int, .1)
        self.assertEqual(frstr.propagator, ssweep_cross)
        np.random.seed(12)
        print 'population before', frstr.forest.population
        repeats = 2
        for dmy in range(repeats):
            frstr.run()
#        print frstr.strategy.fluxes[0,:,:]
        print frstr.forest.branch_no
        self.assertEqual(frstr.forest.branch_no, 2 + repeats)
        self.assertEqual( repeats * frstr.bunch_size,  sum(1
                                                for tr in frstr.forest.trees
                                                for dmy_ch in tr.children
                                                if not tr == self.onegen
                                                and not tr == self.onegen_dup))
        print 'sum of sim starts', sum(1 for tr in frstr.forest.trees
                                       for dmy_ch in tr.children
                                       if not tr == self.onegen)
#        print '++++++ END OF TESTS ++++++'
#        frstr.forest.dump()

if __name__ == '__main__':
    from pstats import Stats
    from tempfile import NamedTemporaryFile
    tf1 = TestForester('test_if_run')
    tf2 = TestForester('test_2d')
    tf1.run()
    with NamedTemporaryFile() as fl:
        cProfile.run("tf2.run()", filename=fl.name)
        stats = Stats(fl.name)
        stats.strip_dirs().sort_stats('time')  # ('cumulative')
        stats.print_stats(20)

