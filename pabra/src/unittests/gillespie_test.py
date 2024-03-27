'''
Created on Jul 14, 2010

@author: nbecker
'''

from pabra.gillespieCy import Reaction, Gillespie  #@UnresolvedImport
import numpy as np

import pabra
from unittests.fixtures import H5TestCase
from pabra.propagate import Absorption
pabra.logger.setLevel(pabra.logging.INFO)
import pabra.ptree_h5 as p5



class TestPropagation(H5TestCase):
    """test the structural operations on the ptree class"""

    def setUp(self):
        self.setupemptytrees()
        self.r1 = Reaction(['mass_action', 10.], [], ['a'], [],  name='0->a')
        self.r2 = Reaction(['mass_action', .1], ['a'], [], ['a'], name='a->0')
        self.r3 = Reaction(['mass_action', 8.], [], ['b'], [],  name='0->b')
        self.r4 = Reaction(['mass_action', .1], ['b'], [], ['b'], name='b->0')
        self.r5 = Reaction(['mass_action', .2], ['a'], ['b'], ['a'], 
                           name='a->b')
        self.r6 = Reaction(['mass_action', .05], ['a', 'b'], [], ['a', 'b'], 
                           name='a+b->0')
        self.g12 = Gillespie(['a'], [self.r1, self.r2])
        self.g12346 = Gillespie(['a', 'b'], [self.r1, self.r2, self.r3,
                                              self.r4, self.r6])
        self.g1to6 = Gillespie(['a', 'b'], [self.r1, self.r2, self.r3,
                                            self.r4, self.r5, self.r6])
        self.g_stop = Gillespie(['a'], [self.r2])
        
    def test_propagation(self):
        tf, st_f = self.g12.step(0, [0])
        print 'time', tf, 'state', st_f
        assert tf > 0
        assert st_f[0] > 0
        # standstill
        ti, st_i = tf, st_f
        tf, st_f = self.g12.int(tf, ti, st_i)
        assert tf == ti
        tf, st_f = self.g12.int(2 * ti, ti, st_i)
        assert tf == 2 * ti
        self.assertRaises(Absorption, self.g_stop.int, 2325.3, 0, [1000])
        # another one
        try:
            self.g_stop.int(2345.3, 0., [1000])
        except Absorption as abs_:
            self.assertTrue(abs_.point[0] < 2000)
        else:
            raise Exception('Absorption was expected!')
        
    def test_growth(self):
        self.onegen.log_int = 1.
        self.onegen.tx_dtype = p5.construct_dtype((0.,[0]))
        # this is not allowed: when another reference holds on to
        # a tree leaf, it is not flushed to disk on deletion of the 
        # temporary Ptree instance!
        # ch = self.onegen.children[0]
        print 'uninitialized t, state:', self.onegen.children[0].tip
        print 'log int', self.onegen.log_int
        for lf in self.onegen.leaves():
            lf.tip=(0,[0])
        self.onegen.grow_tips(self.g12, 200.)
        # safe after propagation:
        ch = self.onegen.children[0]
        print ch
#        self.assertTrue(ch.tip[0] == 200.)
        print ch.trunk_from_root[-20:]
        self.onegen_dup.log_int = 10.
        self.onegen_dup.tx_dtype = p5.construct_dtype((0.,[0,0]))        
        for tr in self.onegen_dup.children:
            tr.tip = 0., [23, 3456]
        self.onegen_dup.grow_tips(self.g12346, 2000.)
        ch_dup = self.onegen_dup.children[0]
        print ch_dup.trunk[-20:]        
        self.assertTrue(len(ch_dup.trunk_from_root) == 200)

    def tearDown(self):
        self.teardowntrees()        
    
