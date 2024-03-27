'''
Created on Apr 26, 2010

@author: nbecker
'''
import unittest as ut
from pabra.propagate import *
import numpy as np
from pabra.branching import *
from unittests.fixtures import H5TestCase

class TestPropagation(H5TestCase):
    """test the structural operations on the ptree class"""    
    
    def setUp(self):
        self.setuptrees()
        self.setupito()
                    
    def tearDown(self):
        self.teardowntrees()
        
        
    def test_eq(self):
        betas = [1,2,1,1]
        self.assert_aae(5./4 * np.array([np.Inf,1,1./2,1./3])[1:],
                eq_weights(betas)[1:])
        self.assert_aae(5./4 * np.array([np.Inf,1,1./2,1./3])[1:],
                        eq_weights(betas)[1:])
        betas = [0,0,0,4]
        self.assert_aae(eq_weights(betas)[1:], 
                        np.array([np.inf,1.,1./2,1./3])[1:])
    
    def test_betas(self):
        self.assert_aae(betas_top_0(1.7), np.array([.15,0,1./2*1.7]))
        self.assert_aae(betas_top_0(4., 3 ), betas_top_0(3., 3 ))
        self.assert_aae(betas_top_two(2.8, 3),[0,0,.2,.8])
        
    def test_bn(self):
        bts1 = betas_top_two(2.8, 3)
        for dmy_i in xrange(100):
            assert (2 <= branch_number(bts1) <= 3)
        bts2 = betas_top_0(2, 4)
        for dmy_i in xrange(100):
            assert branch_number(bts2) in set((0,4)) 
        
    def test_wt(self):
        l = []
        for dmy_i in xrange(1000):
            l.append(exp_wait_time(rate=3))
        assert 1./4 < np.array(l).mean() < 1./2
        
    def test_now(self):
        l = []
        for dmy_i in xrange(10000):
            l.append(branch_trigger())
        self.assert_aae(np.array(l).astype(int).mean(),171./271, 1)
        
        
