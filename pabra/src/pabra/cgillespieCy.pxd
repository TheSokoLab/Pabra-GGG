#@PydevCodeAnalysisIgnore
"""
definition file augmenting gillespieCy.py
"""

#cimport numpy as np
cimport cpython.list as pyl

# test of external imports; not actually used

cdef extern from "math.h": 
    double sin(double)
    double logf(double)
    
cdef extern from "gsl/gsl_rng.h" :
    # random number generator type (mersenne, etc) 's type
    ctypedef struct gsl_rng_type:
        pass
    # random number generator's type
    ctypedef struct gsl_rng:
        pass
    # give the name of a generator
    char* gsl_rng_name(gsl_rng * r)
    # default generator
    gsl_rng_type* gsl_rng_default
    # initialize
    gsl_rng* gsl_rng_alloc (gsl_rng_type* T)
    # seed
    void gsl_rng_set (gsl_rng* r, unsigned long int s)
    # mainly for debugging
    void gsl_rng_print_state (gsl_rng* r)
    # sample a uniform distribution
    double gsl_rng_uniform (gsl_rng* r)
    # exclude 0.0:
    double gsl_rng_uniform_pos (gsl_rng * r)
    # deallocate
    void gsl_rng_free (gsl_rng * r)

cdef extern from "gsl/gsl_randist.h":
    # sample an exponential distribution
    double gsl_ran_exponential (gsl_rng * r, double mu)

