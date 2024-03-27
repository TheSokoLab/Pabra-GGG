'''
Created on Apr 13, 2010

@author: nbecker

collect test for propagation
'''

from pabra.propagate import *
from pabra import ptree_h5 as p5
import numpy as np
import unittest as ut
from unittests.fixtures import H5TestCase

class TestPropagation(H5TestCase):
    """test the structural operations on the p5 class"""

    def setUp(self):
        self.setuptrees()
        self.setupito()

    def tearDown(self):
        self.teardowntrees()
        
    def test_dummyprop(self):
        tr = self.onegen
        # lets see
        def oc(): 
            return self.onegen.children
        tr.log_int = 1
        oc()[0].grow_trunk(self.dummyprop, .5)
        oc()[0].dump(['base', 'tip', 'status'])
        self.assertFalse(self.dummyprop.abs_cond(9, 'sdf'))
        self.assertEqual(tr.leaf_no(.1), 1)
        self.assertEqual(tr.leaf_no(.6), 0)
        #trying to grow back
        self.assertRaises(p5.GrowthError, 
                          oc()[0].grow_trunk, self.dummyprop, .4)
        for t in np.arange(0, 3, .15):
            oc()[1].grow_trunk(self.absprop, t)
        self.assertTrue(oc()[1].absorbed)
        self.assertEqual(oc()[1].height, 1.05)
        tr.log_int = .1
        self.assertEqual(oc()[0].log_int, .1)
#        print oc()[0].tip
#        print oc()[0].log_int
#        print tr.log_int
#        print oc()[0].parent is tr
#        print repr(tr)
        # with this logging interval the absprop is checked every .1:
        oc()[0].grow_trunk(self.absprop, 2.)
#        print oc()[0].tip
        self.assertEqual(oc()[0].height, 1.1)

    @ut.skipIf(0, '')
    def test_growtips(self):
        def oc(): 
            return self.onegen.children
        tr = self.onegen
        tr.log_int = .1
        oc()[0].grow_trunk(self.absprop, 2)
        assert oc()[0].absorbed
        tr.grow_tips(self.dummyprop, 1.5)
        assert tr.alive
        tr.grow_tips(self.absprop2, 3)
        self.assertEqual(tr.height, 2.1)

    @ut.skipIf(0, '')
    def test_ItoEuler(self):
        def oc(): 
            return self.onegen.children
        tr = self.onegen
        tr.log_int = .1
        # see if random growing works ok
        np.random.seed(112)
        #print 'no force:' , self.ito_no_f.force_field(3,'sd')
        tr.grow_tips(self.ito_no_f, 2.)
        #print oc()[1].trunk
#        print oc()[1].tip
#        print oc()[1].base
        self.assertNotAlmostEqual(oc()[1].tip[1], oc()[1].base[1], 1)
        self.assertEqual(len(oc()[1].trunk_from_root), 20)
        # now with a potential
        np.random.seed(123)
        tr.grow_tips(self.ito_sin_sweep, 13.)
        self.assertEqual(len(oc()[1].trunk_from_root), 130)
        # one more for the case without an explicit log_int
        # should default to inf!
        tr_no_log = self.onegen_dup
        self.assertTrue(tr_no_log.log_int == np.inf)
        np.random.seed(412)
        tr_no_log.grow_tips(self.ito_no_f, 2.)
        tr_no_log.cap_tips()
        self.assertFalse(tr_no_log.children[0].trunk_from_root)
        self.assertAlmostEqual(tr_no_log.height, 2, 1)

    @ut.skipIf(0, '')
    def test_Ito_2d(self):
        tr = self.onegen2d
        tr.dump()  # no specifications for Pforest
        # fine. now propagate that.
        np.random.seed(23)
        tr.grow_tips(self.ito_harmonic_2d, 1.)
        tr.dump()
        for t in tr.leaves():
            self.assertFalse(np.allclose(t.tip[1], t.base[1]))
            self.assertAlmostEqual(t.tip[0], 1., 2)

    @ut.skipIf(0, '')
    def test_Ito_2d_aniso(self):
        tr = self.onegen2d
        tr.dump()
        # fine. now propagate that.
        np.random.seed(23)
        tr.grow_tips(self.ito_harmonic_aniso, 1.)
        tr.dump()
        for t in tr.leaves():
            self.assertFalse(np.allclose(t.tip[1], t.base[1]))
            self.assertAlmostEqual(t.tip[0], 1., 2)
            # one direction not moving:
            self.assertAlmostEqual(t.tip[1][0], 0., 3)


if __name__ == '__main__':
    print 'running manually...'
    tp = TestPropagation()
    tp.setup()
    tp.test_Ito_2d_aniso()


#'''
#Created on Apr 14, 2010
#
#@author: nbecker
#'''

