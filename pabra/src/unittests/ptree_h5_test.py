'''
Created on Apr 13, 2010

@author: nbecker

collect test for the ptree_h5 module; slight variation of ptree_test
'''

import unittest as ut
from pabra.ptree_h5 import *
from numpy.testing.decorators import skipif
import tempfile
import numpy as np
import h5py as h5
from unittests.fixtures import H5TestCase

## pass by failing
#@raises(Exception)
#def nose1_test():
#    raise Exception

class TestTreeStructure(ut.TestCase):
    """test the structural operations on the ptree class"""

    def setUp(self):
        self.tf = tempfile.NamedTemporaryFile()
        self.h5rt = h5.File(self.tf.name)
        self.tx_dt = construct_dtype((1., 1.))
        self.p5rt = Pforest(self.h5rt)
        self.p5rt.tx_dtype = self.tx_dt
        self.emptytree = self.p5rt.add_child()
        # now we need to make one nontrivial one
        self.onegen = self.p5rt.add_child()
        self.onegen.add_child(logwt=.1, base=(0., 1.))
        self.onegen.add_child(logwt=.4, base=(0., -1))
        self.another = self.p5rt.add_child()
        # a nested tree
        # attention not to run into an infinite recursion
        self.ntree = self.p5rt.add_child()
        for dmy in range(3):
            subl = list(self.ntree.postorder())  # freeze
            for subtr in subl:
                if not subtr.children:
                    for dmy in '12':
                        subtr.add_child()

    def tearDown(self):
        del self.emptytree, self.onegen, self.ntree, self.h5rt, self.p5rt
        self.tf.close()
        del self.tf

    def test_empty(self):
        tr = self.emptytree
        self.assertIsInstance(tr, Ptree)
        for ass in zip([tr.size, tr.height, tr.is_leaf,
                        tr.survived, tr.branch_no, tr.absorbed, tr.alive,
                        tr.logwt, tr.wt],
                        [1, np.inf, True, False, False, False, True,
                          0., 1.]):
            self.assertEqual(*ass)
        self.assertEqual(tr.__str__(), '<base_t  0, tip_t  0,  growing>')

    @ut.skipIf(0, '')
    def test_children(self):
        print 'test_children'
        tr = self.onegen
        self.assertEqual(self.p5rt.generations, 4)
        self.assertEqual(tr.generations, 1)
        self.assertEqual(tr.size, 3)
        for ch in tr.children:
            self.assertEqual(ch.parent, tr)
            self.assertEqual(ch.depth, 2) # there is a Pforest underneath tr
            self.assertEqual(ch.size, 1)
            self.assertEqual(ch.generations, 0)
        self.assertEqual(tr.total_weight(), np.exp(.1) + np.exp(.4))
        tr.children[0].wt = 3
        self.assertEqual(tr.total_weight(), 3. + np.exp(.4))

    @ut.skipIf(0, '')
    def test_branches_and_tips(self):
        tr = self.onegen
        tr0 = lambda : tr.children[0]
        tr0().tip = (2., 4)
        tr0().branch(0)
        self.assertFalse(tr0().is_leaf)
        self.assertTrue(tr.alive)
        tr1 = lambda: tr.children[1]
        tr1().absorb((3., 1))
        for t, no in [(.2, 2), (2.2, 1), (4, 0)]:
            self.assertEqual(tr.leaf_no(t), no)
        self.assertFalse(tr.alive)
        self.assertEqual(tr.height, 3.)
        self.assertEqual(tr.total_weight(2.2), np.exp(.4))
        #tr.dump(short=1)
        trn = self.ntree
        self.assertTrue(trn.alive)
        self.assertEqual(trn.size, 2 ** 3 + 2 ** 2 + 2 ** 1 + 2 ** 0)
        trn.dump()
        self.assertEqual(trn.leaf_no(), 2 ** 3)
        trn.cap_tips()
        trn.dump()
        print 'leaves:', trn.leaf_no()
        print [tr.is_leaf for tr in trn.leaves()]
        self.assertEqual(trn.leaf_no(), 0)
        self.assertFalse(trn.alive)



    @ut.skipIf(0, '')
    def test_branching(self):
        tr = self.onegen
        tr0 = lambda : tr.children[0]
        tr1 = lambda : tr.children[1]
        tr0().tip = (2., 4)
        tr0().branch(3, point=tr0().tip)
        self.assertRaises(BranchingError, tr0().branch, 2)
        self.assertRaises(BranchingError, tr0().absorb)
        tr1().absorb((3., 1.))
        self.assertFalse(tr0().is_leaf)
        self.assertFalse(tr0().pruned)
        for i, str in enumerate(tr0().children):
            str.tip = (i + 3, 3 * i)
            del str     # THIS IS ESSENTIAL FOR FLUSHING !!!
        self.assertEqual(tr1().height, 3.)
        self.assertEqual(tr0().height, np.inf)
        self.assertEqual(tr.height, np.inf)
        self.assertTrue(tr.alive)
        tr.cap_tips()
        self.assertEqual(tr0().height, 5)
        self.assertEqual(tr.height, 5)
        self.assertFalse(tr.alive)

    @ut.skipIf(0, '')
    def test_traversal(self):
        tr = self.onegen
        tr0 = lambda : tr.children[0]
        tr1 = lambda : tr.children[1]
        tr0().tip = (2., 4)
        tr0().branch(0)
        self.assertEqual(tr0().preorder_t(1.6).next(), tr0())
#        print tr0().preorder_t(2.6).next()
        self.assertRaises(StopIteration, tr0().preorder_t(2.6).next)
        tr1().absorb((3., 1))
        tmpit = tr.preorder_t(2.6)
        tmpit.next()
        self.assertEqual(tmpit.next(), tr1())
        self.assertFalse(tr.alive)
        print list(tr.preorder_t())
        self.assertRaises(StopIteration, tr.preorder_t().next)
        trn = self.ntree
        trn.cap_tips()
        self.assertRaises(StopIteration, trn.preorder_t().next)
        self.assertRaises(BranchingError, trn.cap)
        self.assertEqual(tr1().trunk_from_root, [])

    @ut.skipIf(0, '')
    def test_grafting(self):
        tr = self.onegen
        tr0 = lambda : tr.children[0]
        tr1 = lambda : tr.children[1]
        tr0().tip = (2., 4)
        tr0().branch(0)
        tr1().absorb((3., 1))
        self.assertFalse(tr.alive)
        tw = tr.total_weight(1.)
        tr.graft(2)
        self.assertTrue(tr.alive)
        self.assertEqual(tr.size, 5)
        tr.children[2].tip = (1.2, 3.)
        tr.children[3].tip = (1.2, 3.5)
        self.assertEqual(tr.total_weight(1.), tw)
        self.assertRaises(BranchingError, tr1().graft)

    @ut.skipIf(0, '')
    def test_forest_structure(self):
        tforest = Pforest(tx_dtype=self.tx_dt)
        self.h5rt.copy(self.onegen.shadow, tforest.shadow)
        tforest.dump()
        onegen_copy = tforest.children[0]
        self.assertEqual(onegen_copy.shadow['/'], tforest.shadow['/'])
        self.assertEqual(tforest.population, 1)
        self.assertEqual(tforest.size, 1 + self.onegen.size)
        self.assertFalse(tforest.is_leaf)
        self.assertTrue(tforest.alive)
        self.assertTrue(tforest.branch_no)
        self.assertEqual(tforest.height, np.Inf)
        self.assertEqual(tforest._name_depth('/'), 0)
        self.assertEqual(tforest._name_depth('/s1'), 1)
        self.assertEqual(tforest._name_depth('/s1/t1'), 2)
        self.assertEqual(tforest.depth, self.onegen.depth - 1)
        self.assertEqual(onegen_copy.depth, self.onegen.depth)
        self.assertRaises(BranchingError, self.onegen.absorb)
        self.assertEqual(tforest.log_int, np.inf)
        tforest.log_int = .1
        self.assertEqual(onegen_copy.log_int, .1)
        self.assertEqual(onegen_copy.children[0].trunk_from_root, \
                self.onegen.children[0].trunk_from_root)

    @ut.skipIf(0, '')
    def test_essentials(self):
        self.assertTrue(laba_daba())


'''
Created on Apr 14, 2010

@author: nbecker
'''
