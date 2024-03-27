'''
Created on Apr 20, 2010

@author: nbecker
'''


import copy
from pabra.propagate import *
from pabra import ptree_h5 as ptree
import h5py as h5
import matplotlib.pyplot as plt
import os
import subprocess


def setup():
    global emptytree, onegen, dummyprop, absprop, absprop2, harmonic_2d
    global ito_sin_sweep, ito_harmonic_2d, ito_no_f, onegen2d 
    os.chdir('/home/pabra/tmp')
    for i in '123':
        try: os.unlink('tmp'+i+'.h5')
        except OSError: pass
    emptytree = ptree.Ptree(h5.File('tmp1.h5',driver='core'))
    emptytree.tx_dtype=ptree.construct_dtype((0.,0.))
    onegen = ptree.Ptree(h5.File('tmp2.h5',driver='core'))
    onegen.tx_dtype=ptree.construct_dtype((0.,0.))
    onegen.add_child(0.1)
    onegen.add_child(0.4)
    dummyprop = SystemIntegrator()
    absprop = SystemIntegrator(lambda t, x: bool(t > 1.))
    absprop2 = SystemIntegrator(lambda t, x: bool(t > 2.))
    ito_no_f = ItoEulerFwd(.001, force_field=0.)
    def sinesweeper(t,x):
        def mean(tt):
            return np.sin(tt)
        return -(x - mean(t))
    ito_sin_sweep = ItoEulerFwd(.001, force_field = sinesweeper, kT=.3)
    def harmonic_2d(t, pt):
        return (- pt +  .2 * pt[::-1] * [1,-1]) 
    ito_harmonic_2d = ItoEulerFwd(.01, force_field=harmonic_2d)
    onegen2d = ptree.Ptree(h5.File('tmp3.h5',driver='core'))
    onegen2d.tx_dtype=ptree.construct_dtype((0.,np.zeros((2,))))
    onegen2d.log_int=.1
    onegen2d.logwt = 0.
    for lwt in (-1.,1.):
        onegen2d.add_child(lwt)       
    
def fig_setup():
    global fig, ax, fig2, ax2, t
    fig = plt.figure()
    ax = fig.add_subplot(111)
    t = ax.set_title('visualization tests')
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)

#==============================================================================
# the script
#==============================================================================

# better would be to set up as much as possible and then use interactive 
# drawing commands combined with logging. can then pick out what was
# successful.

if __name__ == '__main__':
    print "git rev: {0}".format(subprocess.check_output(
                                ['git', 'rev-parse', 'HEAD']))
    print "numpy {0}; h5py {1}; HDF5 {2}".format(np.version.version, 
                                                 h5.version.version,
                                                 h5.version.hdf5_version)
    setup()
    fig_setup()
    #onegen.dump()
    tr = onegen
    print onegen.base
    tr.log_int = .05
    # see if random growing works ok
    np.random.seed(12)
    tr.grow_tips(ito_no_f, 12.)
    print tr.children[0].trunk
    for tree in tr:
        if tree.trunk:
            t, x = np.transpose(tree.trunk)
            ax.plot(t, x, color = 'blue')
            # np.clip(np.exp(tree.logwt),0,1)
    # now with a potential
    np.random.seed(123)
    #tr.grow_tips(ito_sin_sweep, 43.)
    for tree in tr:
        if tree.trunk:
            t, x = np.transpose(tree.trunk)
            ax2.plot(t, x, color = (0,0,.2,np.clip(np.exp(tr.logwt),0,1) ))
    #now 2d
    tr2d = onegen2d
#    print tr2d.log_int
#    print ito_harmonic_2d.__dict__
#    print ito_harmonic_2d.deti(0,np.ones(2),.1)
#    print ito_harmonic_2d.D(0,np.ones((2,)))
#    print ito_harmonic_2d.delta
#    print ito_harmonic_2d.sigma(0,np.ones(2),.1)
#    for i in range(4):
#        print ito_harmonic_2d.incr(0,np.ones(2))
    tr2d.grow_tips(ito_harmonic_2d, to_t=12.,log_int=.01)
    #tr2d.children[0].grow_trunk(ito_harmonic_2d,2.)
    #print tr2d.children[0].trunk
    for st in tr2d:
        if st.trunk:
            #print st.trunk
            x, y = np.array([txy[1] for txy in st.trunk]).T
            ax.plot(x,y)
    plt.show()
    
#==============================================================================
# RESULT
#==============================================================================
    
# runs fully in 

#git rev: 3553cdd35f9fbdd4c340b41cfa1d828abd80af43