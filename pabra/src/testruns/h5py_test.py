'''
This is a test to find out if the h5py approach can give a good result
in terms of performance and memory handling.


Created on Nov 30, 2010

@author: nbecker
'''

import h5py as h5
import os
import itertools as it
import time
import numpy as np


os.chdir('/home/pabra/tmp/')

_MB = 1024*1024
# make a new file
h5f = h5.File('h5pyf.h5', mode='w')
# here, root==file
rt = h5f

# tried this: the slow traversal is not in the 
# iteration itself but in getting the children which are
# loaded lazily. The following works but is not faster.
def preord(rt_nd):
    stack = [rt_nd]
    while stack:
        nd_ = stack.pop()
        yield nd_
        stack.extend(nd_.itervalues())    

# create an empty tree
def make_2(grp, i):
    grp.create_group('c'+str(i)+'a')
    grp.create_group('c'+str(i)+'b')
    
def make_1(grp, i):
    grp.create_group('c'+str(i))
    grp.attrs['an_attribute'] = np.random.random((1234,))

# first 2 children up to I
I = 12
# then 1 child up to J
J = 13

t0 = time.time()
# need always two copies of the single use iterators...
cur1, cur2 = it.tee([rt],2)
for i in xrange(I):
    for nd in cur1:
#        print 'nd inner loop', nd
        make_2(nd, i)
#        nd._g_flushGroup()
    h5f.flush()
    cur1, cur2 = it.tee(it.chain(*[nd.itervalues() for nd in cur2]), 2)

print 'done proliferating, time', time.time() - t0
#print 'no of nodes', sum(1 for nd in h5f.walkGroups())
def incr(dmy):
    global i
    i+=1
i=1
h5f.visit(incr)
print 'no of nodes', i
time.sleep(2)

t0 = time.time()
for j in xrange(I,J):
    for nd in cur1:
        make_1(nd,j)
    h5f.flush()
    cur1, cur2 = it.tee(it.chain(*[nd.itervalues() for nd in cur2]), 2)
        
h5f.flush()
print 'done extending', time.time() - t0
time.sleep(2)
i=1
t0 = time.time()
# either:
h5f.visit(incr)
# or:
#itr = preord(rt); itr.next()
#for nd_ in itr:
#    i += 1
print 'tot no of nodes', i
print 'counting time', time.time() - t0
#print h5f

import cProfile
exec_str = """
for k in xrange(J, J + 3):
    for nd in cur1:
        make_1(nd,j)
    cur1, cur2 = it.tee(it.chain(*[nd.itervalues() for nd in cur2]), 2)
"""
cProfile.run(exec_str, sort=1)    

h5f.close()

#os.remove('h5f.h5')