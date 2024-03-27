# encoding: utf-8
# cython: profile=True

'''
Created on Jun 11, 2010

@author: nbecker

cythonized and optimized 
Gillespie direct method simulator  

'''


#----------------------------------------------------------------- cython stuff

import cython
import numpy as np
cimport cython
cimport numpy as np 
from cgillespieCy cimport *

# type declarations


ctypedef np.uint_t index_t 
ctypedef np.float64_t float_t
# need signed int because of error checking for neg. counts
ctypedef np.int64_t int_t

#------------------------------------------------------------ pabra integration

from pabra.propagate import SystemIntegrator, Absorption

#---------------------------------------------------------------------- logging

import logging
log_ = logging.getLogger(name=__name__)

#--------------------------------------------------------------- random numbers

# seed used for each new Simulator generator
rand_seed = 123

# the rng to be used:
cdef gsl_rng* the_rng
the_rng = gsl_rng_alloc(gsl_rng_default)
#gsl_rng_set(the_rng, rand_seed)

# wrapping class to get the same Python API
cdef class Rng(object):
    cdef gsl_rng * rng
    cdef int a
    def __init__(self, seed=1):
        self.rng = the_rng
        self.seed(seed)
    def set_a(self, int b):
        self.a = b
    def seed(self, int seed):
        gsl_rng_set(self.rng, seed)
    def _print_state(self):
        gsl_rng_print_state(self.rng)


# now the wrapped Python rng instance; note the missing underscore
therng = Rng(rand_seed)
therng.seed(123)

class Gillespie(SystemIntegrator):
    """Gillespie direct method;
    
    :species_names: a list of string names for the species; species are
    internally referred to by their index in this list
    :reactions: a list of Reactions
    :init_state: the initial state vector, of the same length as species_names
    :step_chunk:  a max.step number if no crossings occur
    :**sysint_kwargs: init arguments of SystemIntegrator. Currently this is
    only the absorption condition.
    
    This is a SystemIntegrator instance which remains in pure python; it calls
    the optimized helper functions.
    
    """
    def __init__(self, species_names, reactions, step_chunk=1000, 
                 **sysint_kwargs):
        # parent initialization; this is mainly the absorption condition
        SystemIntegrator.__init__(self, **sysint_kwargs)
        # attributes
        self.species_names = species_names
        self.n_species = len(species_names)
        self.reactions = reactions
        self.k_reactions = len(self.reactions)
        self.step_chunk = step_chunk  # a max step number
        # gsl rng. declaration is done at top-level, outside 
        global the_rng
        
        # now the arrays at the core of the simulation
        # construct an array of increments and decrements of mol numbers
        self.changes = np.array([[re.products.count(s) - re.reactands.count(s)
                                  for s in self.species_names]
                                 for re in self.reactions], dtype=np.int64)
        # a k-tuple of lists which contain the indices in the state
        # vector of the input species to the propensity of the l-th reaction
        self.prop_inputs = [[self.species_names.index(s)
                             for s in re.prop_species]
                            for re in self.reactions]
        # array-ify by filling up the zeros:
        self.input_lengths, self.prop_inputs = self.arrayify(self.prop_inputs)
        # Coupling between reactions:
        # list each reaction which can be affected by the change in counts
        # due to the l-th reaction.
        self.prop_affected = [# make a sorted list of unique elements
                                sorted(list(set(
                                # of: reactions r
                                [r for r in range(self.k_reactions)
                                # such that species i is changed by reaction l
                                for (i, s) in enumerate(self.changes[l]) if s
                                # and r depends on i;
                                if i in self.prop_inputs[r] ] )))
                                # listed over all reactions l 
                                for l in range(self.k_reactions)]
        # arrayify.
        (self.affected_lengths, self.prop_affected) = self.arrayify(
                                                            self.prop_affected)
        # initialize the prop fns
        for re, inp in zip(self.reactions, self.prop_inputs):
            re.init_prop_fn(inp)


    # for convenience
    @property
    def reac_dict(self):
        return dict((r.name, r) for r in self.reactions)
    
    # tested: the overhead of calling self.step is 
    # negligible from step number N=1000 on. At N=100, about 20%

    def step(self, tinit, init_state, N=None):
        # only time needs extra assignment
        if not N: 
            N = self.step_chunk
        state = np.array(init_state)
        time = tinit
        props = np.array([r.prop(state) for r in self.reactions])
        time = steps(self, state, props, time)
        return time, state
    
    def _prop_or_abs(self, tfinal, tinit, init_state):
        """Propagate to the final time if possible
        
        An Absorption() may occur before tfinal; this has to be 
        handled by the caller. 
        If the chunk size is insufficient, propagation is looped and
        a warning is logged. No protection against infinite loops.
        
        If tinit and xinit are given, the inital time and state are reset
        to these values before starting.
        """
        state = np.array(init_state)
        time = tinit
        props = np.array([r.prop(state) for r in self.reactions])
        if self.abs_cond(time, state):   # Absorption only @ start of run
            raise Absorption(time, state)
        while time < tfinal:
            time = steps(self, state, props, time, tfinal)
            if time < tfinal and log_.isEnabledFor(logging.INFO):
                log_.info('final t={0} not reached within'
                             ' {1} steps; continuing.'
                             .format(tfinal, self.step_chunk))
        return (time, state)
        
    # utility for making enveloping arrays
    # these are necessary to avoid the slower looping over python lists. 
    @classmethod   # staticmethod does not work in cython...
    def arrayify(cls, lilist):
        lens = np.array([len(s) for s in lilist], dtype=np.uint)
        vals = np.zeros((len(lens), max(lens)), dtype=np.uint)
        for l in range(len(lens)):
            vals[l, :lens[l]] = lilist[l]
        return lens, vals


#----------------------------- optimized functions for cum sum and index search

# cdef implementation.
# need to pass arrays for output since
# these cannot be allocated easily inside the cdef f'n.

@cython.boundscheck(False)
@cython.profile(False)   # set to False, otherwise nonsense profiling result
# passing pointers to float_t; i.e. we need to get the array data 
# out of the ndarrays on call of this function; as in steps_2
cdef int cum_sum(index_t k, float_t *parray, float_t *carray):
    cdef float_t cs_run = 0
    cdef index_t i
    for i in range(k):
        cs_run += parray[i]
        carray[i] = cs_run
    return 0


# cdef implementation.
# returning index_t (which is unsigned); error handling with unlikely value.
@cython.boundscheck(False)
@cython.profile(False)
# again, we pass a pointer to a C array.
cdef index_t search_up(index_t k, float_t *carray, 
                         double c_val) except? 1000:
    cdef index_t i
    for i in range(k):
        if c_val <= carray[i]:
            return i
    # this should only happen if an out-of-bounds r.n. was given
    raise RateError

#------------------------------------------------ fast typed stepping function

# since this is called from python anyway, no need to cdef it.
@cython.boundscheck(False)
def steps(self, 
          np.ndarray[int_t, ndim=1, mode='c'] state,
          np.ndarray[float_t, ndim=1, mode='c'] props,
          float_t time, 
          float_t up_to=-1.):
    """Make N steps, except for Exceptions
    
    """
    # import rng
    global the_rng
    # import stuff from self
    cdef list reacs = self.reactions
    cdef index_t N = self.step_chunk
    cdef np.ndarray[int_t, ndim=2, mode='c'] changes = self.changes
    cdef np.ndarray[index_t, ndim=1, 
                    mode='c'] affected_lengths = self.affected_lengths
    cdef np.ndarray[index_t, ndim=2, 
                    mode='c'] prop_affected = self.prop_affected
    # hot loop locals
    cdef np.ndarray[float_t, ndim=1, mode='c'] cum_props = np.zeros_like(props)
    cdef np.ndarray[int_t, ndim=1, mode='c'] temp_state = np.zeros_like(state)
    cdef float_t prop_sum
    cdef index_t dmy_n, i, l, m, j
    cdef index_t k_reactions = props.shape[0]
    # stop condition
    # hot loop
    for dmy_n in range(N):
        # calling the cdef function with proper 'safe' data access.
        cum_sum(k_reactions,
                  <float_t *>props.data, 
                  <float_t *>cum_props.data
                  )
        prop_sum = cum_props[k_reactions - 1]
        if prop_sum <= 0:
            if log_.isEnabledFor(logging.DEBUG):
                log_.debug('Absorbed: no more reactions possible')
            raise Absorption(time, state)
        # add the first reaction time
        # this is a little faster than gsl_ran_exponential
        time -= logf(gsl_rng_uniform_pos(the_rng)) / prop_sum 
        if 0. < up_to < time:
            return up_to
        # which reaction?
        l = search_up(k_reactions, <float_t *>cum_props.data,
                        prop_sum * gsl_rng_uniform(the_rng))
        # update
        for i in range(state.shape[0]):
            state[i] += changes[l,i]
        # update affected propensities
        for m in range(affected_lengths[l]):
            j = prop_affected[l, m]
            try:
                # a try with pointing to functions. most flexible option!
                # fastest version: no temp variable;
                # first paren: the list item as PropFn, which is a pointer
                # type (although this is not visible explicitly)
                # eval part: direct access to the cdef class method
                props[j] = (<ReactionBase>(reacs[j])).prop_eval(
                                <int_t*>state.data, <int_t*>temp_state.data)
            except RateError, e:
                print 'last l', l
                print 'state ', state
                print 'cum', cum_props
                print 'props', props
                raise
    return time

#--------------------------------------------------------- reaction definitions

# defining outside the class.. workaround for a cython limitation
def _new_switch(cls, kinetics, *args, **kwargs):
    if kinetics[0].startswith('mass'):
        return MAReaction(kinetics[1], *args, **kwargs)
    else: 
        raise ReactionInstantiationError('unknown reaction kinetics type')

class Reaction(object):
    """ a reaction. 
    
    This makes a specialized reaction type based on the first argument.
    :kinetics: gives the type of registered kinetics.
    
    """
    __new__ = _new_switch

cdef class ReactionBase(object):
    """reaction base class
    
    :reactands: input molecule list
    :products: output molecule list
    :prop_species: list of species on which the propensity fn depends
    """
    description = 'Generic Reaction'
    cdef public object reactands, products, prop_species, name
    def __init__(self, reactands, products, prop_species=None, name='noname'):
        self.reactands = reactands
        self.products = products
        self.prop_species = prop_species or self.reactands
        self.name = name
    def __repr__(self):
        return self.__class__.description + ' ' + self.name
    
    # signature for declaring in Gillespie
    cdef float_t prop_eval(self, int_t* state, int_t* temp_state) except -1.:
        raise ReactionInstantiationError
    
#    def __reduce__(self):
#        return type(self), (),  self.__getstate__()

    def __getstate__(self):
        return (self.reactands, self.products, self.prop_species, self.name)
    
    def __setstate__(self, state):
        (self.reactands, self.products, self.prop_species, self.name) = state
#    
# the known types of reaction kinetics

cdef class MAReaction(ReactionBase):
    """A mass action kinetics reaction.

    extra attribute:
    :rate: the rate constant multiplying the molecule numbers
    
    this makes an appropriate propensity function of class MaPropFn
    """
    description = 'MAReaction'

    cdef public float_t rate
    cdef index_t length
    cdef index_t * inputs
    cdef object input_inds

    def __init__(self, rate, *args, **kwargs):
        ReactionBase.__init__(self, *args, **kwargs)
        self.rate = rate
        self.length = len(self.reactands)
        self.input_inds = None
        
    def init_prop_fn(self, np.ndarray[ndim=1, dtype=np.uint_t, mode='c'] 
                     input_indices=None):
        if input_indices == None: 
            if self.input_inds:
                input_indices = self.input_inds
        elif input_indices.shape[0] < self.length:
            raise ReactionInstantiationError('bad input indices list')
        self.input_inds = input_indices
        self.inputs = <index_t *>input_indices.data
        
    @cython.boundscheck(False)
    @cython.profile(False)
    cdef float_t prop_eval(self, int_t* state, int_t* temp_state) except -1.:
        """propensity for a single mass action reaction; 
        
        takes rate, length of input, a 1d array of input species indices,
        and the current state"""
        cdef index_t i
        cdef float_t res = self.rate
        # copy state over
        for i in range(self.length):
            temp_state[self.inputs[i]] = state[self.inputs[i]]
        for i in range(self.length):
            # this tests integers and should be safe
            if temp_state[self.inputs[i]] <= 0:
                return 0.
            # build up the product
            res *= temp_state[self.inputs[i]]
            # take one molecule out of the pool that participates
            # only relevant if a species occurs more than once in the reactands
            temp_state[self.inputs[i]] -= 1
        if res < 0: 
            raise RateError('rate product', res, '< 0')
        return res
    
    # slow evaluation for calls from python space
    def prop(self, np.ndarray state):
        res = self.rate
        temp_state = state.copy()
        for i in range(self.length):
            if temp_state[self.inputs[i]] == 0:
                return 0.
            res *= temp_state[self.inputs[i]]
            temp_state[self.inputs[i]] -= 1
        return res

    def __getstate__(self):
        state = ReactionBase.__getstate__(self)
        return state + (self.rate, self.length, self.input_inds)

    def __setstate__(self, state):
        self.rate, self.length, self.input_inds = state[-3:]
        ReactionBase.__setstate__(self, state[:-3])
        self.init_prop_fn()
        
        
    
# TODO: untested!
cdef class MMReaction(ReactionBase):
    """A Michaelis-Menten kinetics reaction.
    
    Has to be called with 1 reactand (the substrate)
    The enzyme concentration is irrelevant; this is a coarse-
    graining of the scheme
    E + S <-> ES -> E + P
    with a fast equilibrium and a constant total enzyme concentration;
    propensity = vmax [S] / (K_M + [S]).
    
    need as positional arguments:    
    :K_M: the Michaelis-Menten constant K_M
    :vmax: the maximum reaction speed
    here, vmax is prop. to the total enzyme concentration and K_M
    is a combination of the rate constants.
    
    """
    description = "MMReaction"

    cdef public float_t K_M, vmax
    cdef index_t input_i
    cdef object input_ind

    def __init__(self, K_M, vmax, *args, **kwargs):
        if len(args[0]) != 1:
            raise ReactionInstantiationError('only 1 substrate allowed')
        ReactionBase.__init__(self, *args, **kwargs)
        self.K_M = K_M
        self.vmax = vmax

    def init_prop_fn(self, input_index=None):
        if input_index is None:
            input_index = self.input_ind
        self.input_i = input_index
    
    cdef float_t prop_eval(self, int_t* state, int_t* temp_state) except -1.:
        """propensity for a single mass action reaction; 
        
        takes rate, length of input, and a 1d array of input species indices,
        and the current state"""
        cdef float_t sub = state[self.input_i]
        try:
            sub = self.vmax * sub / (self.K_M + sub)
        except:
            return -1.0
        return sub
    # slow evaluation for calls from python space
    def prop(self, np.ndarray state, np.ndarray temp_state):
        sub = state[self.input_i]
        return self.vmax * sub / (self.K_M + sub)
    
    def __getstate__(self):
        state = ReactionBase.__getstate__(self)
        return state + (self.K_M, self.vmax, self.input_ind)

    def __setstate__(self, state):
        self.K_M, self.vmax, self.input_ind = state[-3:]
        ReactionBase.__getstate__(self, state[:-3])
        self.init_prop_fn()


#----------------------------------------------------------------------- errors


class ReactionInstantiationError(Exception):
    pass

class RateError(Exception):
    pass

