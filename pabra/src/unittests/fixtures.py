'''
Created on Jan 17, 2012

@author: nbecker

shared setup and teardown
'''

import unittest as ut
import os
import h5py as h5
import pabra.ptree_h5 as p5
from pabra.propagate import SystemIntegrator, ItoEulerFwd
import numpy as np

TESTDIR = "/tmp/pabra"
try: os.makedirs(TESTDIR)
except OSError: pass

class H5TestCase(ut.TestCase):
    """set up a few Pforests in the TESTDIR"""
    
    treenames = ('emptytree', 'onegen', 'onegen_dup', 'onegen2d')

    def setuptrees(self):
        self.setupemptytrees()
        for pf in (self.onegen,self.onegen_dup):
            pf.tx_dtype=p5.construct_dtype((0.,0.))
            pf.add_child(0.1, base=(0., 1.))
            pf.add_child(0.4, base=(0., -1))
        self.onegen2d.tx_dtype = p5.construct_dtype((0.,(0.,0.)))
        for lw in (-1.,1.):
            self.onegen2d.add_child(lw)
    
    def setupemptytrees(self):
        self.teardowntrees()
        for name in self.treenames:
            setattr(self, name, p5.Pforest(h5.File(name + '.h5')))
        for pf in (self.onegen,self.onegen_dup):
            pf.add_child()
            pf.add_child()
        
    
    def setupito(self):
        self.dummyprop = SystemIntegrator()
        self.absprop = SystemIntegrator((lambda t, x: bool(t > 1.)))
        self.absprop2 = SystemIntegrator((lambda t, x: bool(t > 2.)))
        self.ito_no_f = ItoEulerFwd(.01, (lambda t, x:0.), D=1.)
        def sinesweeper(t, x):
            def mean(tt):
                return np.sin(tt)
            return - (x - mean(t))
        self.ito_sin_sweep = ItoEulerFwd(.01, force_field=sinesweeper)
        def harmonic_2d(t, pt):
            return - pt + .5 * ([1,-1]*pt[::-1])
        def harmonic_2d_diag(t, pt):
            return - pt
        self.ito_harmonic_2d = ItoEulerFwd(.01, force_field=harmonic_2d, D=2.)
        #
        self.ito_harmonic_aniso = ItoEulerFwd(.01, 
                                              force_field=harmonic_2d_diag,
                                              M=[1, 1.], D=np.array([0., 1.]))
        
        
    def teardowntrees(self):
        os.chdir(TESTDIR)
        for name in self.treenames:
            try: os.unlink(name + ".h5")
            except OSError: pass
        
    def assert_aae(self, *args):
        return np.testing.assert_allclose(*args)

