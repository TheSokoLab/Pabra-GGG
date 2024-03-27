'''
test the layout of all trunks in a single dataset vs in individual groups
Created on Apr 3, 2012

@author: nbecker
'''


import h5py as h5
import os
import itertools as it
import time
import numpy as np
import subprocess
from pabra.ptree_h5 import construct_dtype, PtreeBase
import os.path as pp

print "git rev: {0}".format(subprocess.check_output(
                            ['git', 'rev-parse', 'HEAD']))
print "numpy {0}; h5py {1}; HDF5 {2}".format(np.version.version, 
                                             h5.version.version,
                                             h5.version.hdf5_version)

os.chdir('/tmp/')

_MB = 1024*1024
# make a new file
h5_file = h5.File('p5layout.h5', mode='w', libver='latest', driver='sec2')

rt = PtreeBase(h5_file)

dt = construct_dtype((0.,[0.,0.]))
constant_dataset = np.zeros(23, dtype=dt)
constant_dataset[:] = (1.,[2.,3.])

# the central repo
ba = h5_file.create_dataset ('bigarray', dtype=dt, chunks=True,
                            shape=(5,), maxshape=(None,),
                            compression=None)

def make_2(grp):
    for n in '01':
        grp.create_group(n)    

# fill with datasets
def fill_direct(grp):
    for n in '01':
        grp[n].create_dataset('d', data=constant_dataset)

# create a linked tree

class fill_linked(object):
    def __init__(self, counter):
        self.counter = counter
    def __call__(self, grp):
        ln = constant_dataset.shape[0]
        nextcounter = self.counter + 2 * ln
        if nextcounter >= ba.shape[0]:
            ba.resize((ba.shape[0] + 100 * ln,))
#            print 'new shape', ba.shape
        for i, n in zip((0,1),'01'):
            r = ba.regionref[self.counter+i*ln:self.counter+(i+1)*ln]
            # r is only usable as address after flushing (!)
            ba[self.counter+i*ln:self.counter+(i+1)*ln] = constant_dataset
            grp[n].attrs['ref']=r
        self.counter = nextcounter


# first 2 children up to I
I = 2
# then 1 child up to J
J = 12

print "I {0}; J {1}".format(I,J)

t0 = time.time()
# need always two copies of the single use iterators...
cur1, cur2 = it.tee([h5_file])

linkmaker = fill_linked(counter=0)

for i in xrange(I):
    for nd in cur1:
        make_2(nd)
#        fill_direct(nd)
        linkmaker(nd)
    h5_file.flush()
    cur1, cur2 = it.tee([g[i] for g in cur2 for i in g if i in '01'])

def prt(grp):
    print grp
    return None

h5_file.visit(prt)

exec_str = """
for i in xrange(J):
    for nd in cur1:
        make_2(nd)
        fill_direct(nd)
        linkmaker(nd)
    h5_file.flush()
    cur1, cur2 = it.tee([g[i] for g in cur2 for i in g if i in '01'])
"""
# J extra levels!
import cProfile
cProfile.run(exec_str, sort=1)    

#Pforest(h5_file).dump()


h5_file.close()    
    
    