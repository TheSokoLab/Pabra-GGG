'''
Created on Jun 14, 2010

@author: nbecker
'''

from pabra.gillespieCy import Reaction, Gillespie  #@UnresolvedImport
import subprocess
import numpy as np

if __name__ == '__main__':
    import numpy
    print "git rev: {0}".format(subprocess.check_output(
                            ['git', 'rev-parse', 'HEAD']))
    print 'numpy ', numpy.version.version
    print "numpy {0}".format(np.version.version)

    re1 = Reaction(['mass_action', .1], ['a'], ['b'], ['a'], name='a->b')
    re2 = Reaction(['mass_action', 12.], [], ['a'], [],  name='0->a')
    re3 = Reaction(['mass_action', .0002], ['a', 'b'], [], ['a', 'b'],
                   name='a+b->0')
    re4 = Reaction(['mass_action', 3.], ['b'], [], ['b'],
                   name='b->0')
    for (i, re) in enumerate([re1, re2, re3]):
        print 'reaction ' + str(i) + ', rate: ', re.rate
    g1 = Gillespie(['a', 'b'], [re1, re2, re3])
    print 'simulator:\ro    n', g1.__dict__ 
#    print [re.prop_fn(*g1.state[g1.prop_inputs[i]]) for (i, re) 
#                                in enumerate(g1.reactions)]
#    print [re.prop_fn(*g1.state[re.prop_inputs]) for re 
#                                in  g1.reactions]
     
#    print 'testcls', testcls().amethod(9)

#                   init_state=np.array([23, 12], dtype=np.int_)
    exec_str = ('for dmy in xrange(1234): tf, state_f = '
                'g1.step(0,np.array([23, 12],'
                     'dtype=np.int_), N=10000)')
    print exec_str
    import cProfile
    cProfile.run(exec_str, sort='time')
    print 'final ...'
    print 'time: ', tf  #@UndefinedVariable
    print 'state: ', state_f  #@UndefinedVariable
    
#==============================================================================
# RESULTS
#==============================================================================

# with cython 0.15.1

#git rev: 3553cdd35f9fbdd4c340b41cfa1d828abd80af43
#
#numpy  1.6.1
#numpy 1.6.1
#reaction 0, rate:  0.1
#reaction 1, rate:  12.0
#reaction 2, rate:  0.0002
#simulator:
#o    n {'reactions': [MAReaction a->b, MAReaction 0->a, MAReaction a+b->0], 'prop_inputs': array([[0, 0],
#       [0, 0],
#       [0, 1]], dtype=uint64), 'species_names': ['a', 'b'], 'input_lengths': array([1, 0, 2], dtype=uint64), 'affected_lengths': array([3, 3, 3], dtype=uint64), 'abs_cond': <function <lambda> at 0xedcb18>, 'k_reactions': 3, 'n_species': 2, 'prop_affected': array([[0, 1, 2],
#       [0, 1, 2],
#       [0, 1, 2]], dtype=uint64), 'changes': array([[-1,  1],
#       [ 1,  0],
#       [-1, -1]]), 'step_chunk': 1000}
#for dmy in xrange(1234): tf, state_f = g1.step(0,np.array([23, 12],dtype=np.int_), N=10000)
#         14810 function calls in 0.244 seconds
#
#   Ordered by: internal time
#
#   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
#     1234    0.180    0.000    0.188    0.000 gillespieCy.pyx:220(steps)
#     3702    0.029    0.000    0.029    0.000 gillespieCy.pyx:396(prop)
#     1234    0.014    0.000    0.230    0.000 gillespieCy.pyx:137(step)
#     1234    0.010    0.000    0.010    0.000 {numpy.core.multiarray.array}
#        1    0.004    0.004    0.244    0.244 <string>:1(<module>)
#     2468    0.003    0.000    0.008    0.000 numeric.py:65(zeros_like)
#     2468    0.003    0.000    0.003    0.000 {numpy.core.multiarray.empty_like}
#     2468    0.002    0.000    0.002    0.000 {method 'fill' of 'numpy.ndarray' objects}
#        1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
#
#
#final ...
#time:  44.5573026215
#state:  [ 77 271]

