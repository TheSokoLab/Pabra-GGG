'''
This is a test to find out if our ptree implementation is good
in terms of performance and memory handling.


Created on Nov 30, 2010

@author: nbecker
'''

import h5py as h5
import os
import itertools as it
import time
import numpy as np
from pabra import ptree_h5 as p5
import subprocess


print "git rev: {0}".format(subprocess.check_output(
                            ['git', 'rev-parse', 'HEAD']))
print "numpy {0}; h5py {1}; HDF5 {2}".format(np.version.version, 
                                             h5.version.version,
                                             h5.version.hdf5_version)

os.chdir('/tmp/')

_MB = 1024*1024
# make a new file
h5_file = h5.File('p5.h5', mode='w', libver='latest', driver='sec2')
# make a PtreeBase of the file
#rt = p5.PtreeBase(h5_file)
rt = p5.PtreeBase(h5_file)

# create an empty tree
def make_2(tr, i):
    # leaving away the names because these are generated automatically
    tr.add_child()
    tr.add_child()
        
def make_1(tr, i):
    c = tr.add_child()
    c._set_attr('anattr', np.random.random((1234,)))


# first 2 children up to I
I = 10
# then 1 child up to J
J = 40

print "I {0}; J {1}".format(I,J)

t0 = time.time()
# need always two copies of the single use iterators...
cur1, cur2 = it.tee([rt], 2)
for i in xrange(I):
    for nd in cur1:
#        print 'nd inner loop', nd
        make_2(nd, i)
#        nd._g_flushGroup()
    h5_file.flush()
    cur1, cur2 = it.tee(it.chain(*[nd.children for nd in cur2]), 2)

print 'done proliferating, time', time.time() - t0
#print 'no of nodes', sum(1 for nd in h5_file.walkGroups())
def incr(dmy):
    global i
    i+=1
i=1
h5_file.visit(incr)
print 'no of nodes', i
time.sleep(2)

assert 0


t0 = time.time()
for j in xrange(I,J):
    for nd in cur1:
        make_1(nd,j)
    h5_file.flush()
    cur1, cur2 = it.tee(it.chain(*[nd.children for nd in cur2]), 2)
        
h5_file.flush()
print 'done extending', time.time() - t0
time.sleep(2)
i=1
t0 = time.time()
# either:
h5_file.visit(incr)
print 'tot no of nodes', i
print 'counting time', time.time() - t0
#print h5_file

i=1
t0 = time.time()
# or:
itr = rt.preorder(); itr.next()
for nd_ in itr:
    i += 1
print 'tot no of nodes', i
print 'counting time 2', time.time() - t0

i=1
t0 = time.time()
# or:
itr = rt._group_preorder(); itr.next()
for nd_ in itr:
    i += 1
print 'tot no of nodes', i
print 'counting time 3', time.time() - t0

import cProfile
exec_str = """
for k in xrange(J, J + 3):
    for nd in cur1:
        make_1(nd,j)
    cur1, cur2 = it.tee(it.chain(*[nd.children_iter() for nd in cur2]), 2)
"""
cProfile.run(exec_str, sort=1)    

h5_file.close()

#os.remove('h5_file.h5')

#==============================================================================
# RESULTS
#==============================================================================

# core is _slower_ in writing!!!

# also: visit() is much better for broad trees than for tall trees
# for tall trees, the self-made iterators are a lot better

# extending (ie making new nodes is _linear_ in the total size, in the range
# tested. this is independent of broad or tall trees!

#-------------------------------------------------------------- #driver='core:'

#git rev: 3553cdd35f9fbdd4c340b41cfa1d828abd80af43
#
#numpy 1.6.1; h5py 2.0.1; HDF5 1.8.8
#I 10; J 40
#done proliferating, time 0.535041093826
#no of nodes 2047
#done extending 24.3482589722
#tot no of nodes 32767
#counting time 1.47802305222
#tot no of nodes 32767
#counting time 2 4.60841107368
#tot no of nodes 32767
#counting time 3 3.12693095207
#         432135 function calls in 1.321 seconds
#
#   Ordered by: internal time
#
#   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
#     7168    0.150    0.000    0.182    0.000 {h5py.h5o.open}
#     3072    0.129    0.000    0.153    0.000 {h5py.h5g.create}
#     7168    0.100    0.000    0.111    0.000 group.py:252(__iter__)
#    10240    0.085    0.000    0.086    0.000 base.py:49(_shared)
#     3072    0.074    0.000    0.096    0.000 {h5py.h5a.create}
#     3072    0.060    0.000    0.060    0.000 {method 'random_sample' of 'mtrand.RandomState' objects}
#    20480    0.057    0.000    0.057    0.000 weakref.py:62(__contains__)
#     3072    0.052    0.000    0.072    0.000 {method 'write' of 'h5py.h5a.AttrID' objects}
#        1    0.036    0.036    1.321    1.321 <string>:2(<module>)
#     3072    0.033    0.000    0.526    0.000 ptree_h5.py:248(add_child)


#driver='sec2:'


#git rev: 3553cdd35f9fbdd4c340b41cfa1d828abd80af43
#
#numpy 1.6.1; h5py 2.0.1; HDF5 1.8.8
#I 10; J 40
#done proliferating, time 0.517347097397
#no of nodes 2047
#done extending 19.9993340969
#tot no of nodes 32767
#counting time 1.5415160656
#tot no of nodes 32767
#counting time 2 4.52957105637
#tot no of nodes 32767
#counting time 3 3.20304608345
#         432135 function calls in 1.292 seconds
#
#   Ordered by: internal time
#
#   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
#     7168    0.149    0.000    0.179    0.000 {h5py.h5o.open}
#     3072    0.139    0.000    0.163    0.000 {h5py.h5g.create}
#     7168    0.101    0.000    0.112    0.000 group.py:252(__iter__)
#     3072    0.081    0.000    0.103    0.000 {h5py.h5a.create}
#    10240    0.080    0.000    0.082    0.000 base.py:49(_shared)
#     3072    0.056    0.000    0.056    0.000 {method 'random_sample' of 'mtrand.RandomState' objects}
#    20480    0.054    0.000    0.054    0.000 weakref.py:62(__contains__)
#     3072    0.051    0.000    0.069    0.000 {method 'write' of 'h5py.h5a.AttrID' objects}
#        1    0.036    0.036    1.292    1.292 <string>:2(<module>)
#     3072    0.031    0.000    0.520    0.000 ptree_h5.py:248(add_child)
