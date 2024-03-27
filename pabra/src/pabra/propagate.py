'''
Created on Apr 1, 2010

@author: nbecker
This module contains the system propagators. First, there are toy ones,
but later wrappers around fast simulators should follow.
'''
from __future__ import division

import numpy as np
from numpy.random import RandomState
import copy as cp
# logging; overall level and format controlled in __init__.py
import logging
log_ = logging.getLogger(name=__name__)
import math
from pabra import NIError

# random generator management
therng = RandomState()
therng.seed(123)



class SystemIntegrator(object): 
    """Propagate the system dynamics forward in time.

    When an absorption condition is met, Absorption is thrown.
    :abs_cond: absorption condition, either static or a function;
    if it's a function, it must take time and state as input
    can include absorption by generating True with certain probability!

    """
    #TODO: this is a major candidate for optimization.
    #we need to get a good idea of how to handle the
    #events: branching, pruning, absorption efficiently,
    #and how often to check for absorption and branching.
    #
    #ideally, the whole integrator should be a cython 
    #cdef class and the abs_cond as well as the actual
    #dynamics should be compiled C-functions.

    
    # dummy for testing; does nothing except yelling if absorbed
    def __init__(self, abs_cond=False):
        if not callable(abs_cond):
            self.abs_cond = (lambda t, x: abs_cond)
        else:
            self.abs_cond = abs_cond

    def _prop_or_abs(self, tfinal, tinit, xinit):
        """ just checks for absorption at the initial time! """
        if self.abs_cond(tinit, xinit):
            raise Absorption(tinit, xinit)
        return (tfinal, xinit)

    def int(self, tfinal, tinit, xinit):
        """This may throw an Absorption event that needs to be handled!"""
        return self._prop_or_abs(tfinal, tinit, xinit)


class ItoEulerFwd(SystemIntegrator):
    '''Euler-Maruyama forward integration

    This integrates the Ito Langevin equation; 
    everything is evaluated pre-point.
    
    dx = M(t,x).f(t,x) dt + sqrt(2) B(t,x).dw(t) 

    Initialized with 
    :delta: timestep
    :force field: scalar or 1xd ndarray; constant or function of t, x
        do not forget to include the friction force here!!!!!!
    :D: diffusion constant. 
        as scalar: constant or fn of t, x; 
        as  1d ndarray: only constant supported. Interpreted as diagonal.
        Then the noise strength B := sqrt(D) (but see above formula)
    :kT: thermal energy scale, defaults to 1.
    :M: mobility. defaults to D / kT. 
        If given, must obey the same restrictions as D
    :abs_cond: condition of absorption: returns bool; constant or fn of t, x
    
    These arguments default to overdamped Langevin dynamics, but allow also
    Langevin dynamics with mass, e.g. by letting x = p,q (momentum and pos.)
    In that case, since D and M are restricted to be diagonal, have to put
    the relation \dot q = p/m into the force vector.  Do not forget to include 
    the friction force here!!!!!!

    '''
    def __init__(self, delta, force_field, D=1., kT=1., 
                 M=None, abs_cond=False, ref_box=None, inc_type='gaussian'):
        # define the local step update functions.
        self.delta = delta  # timestep
        # here the if statement is evaluated only once, on assignment.
        self.force_field = (force_field if callable(force_field) 
                            else lambda t, x: force_field)
        if callable(D):
            self.B = lambda t, x: self.sqrt(D(t, x))  # only scalar D
        else:
            self.D_ = np.array(D)
            # handle possible nonrandom directions (0-entries in D) 
            self.D_inds = self.D_.nonzero()
            if self.D_inds and len(self.D_inds[0]) < self.D_.size:
                # reduced dimension for B:
                self.red = True
                # the square root of D for inside the increment
                self.B_red = np.sqrt(self.D_[self.D_inds]) 
                self.B = lambda t, x: self.B_red
            else:
                self.red = False
                self.B_ = np.sqrt(self.D_)
                self.B = lambda t, x: self.B_
        if M is None:
            self.M_ = self.D_ / kT
            self.M = lambda t, x: self.M_
        else:
            if callable(M): 
                self.M = M
            else:
                self.M_ = np.array(M)
                self.M = lambda t, x: self.M_
        # small speedups: localizing function calls in tight loops
        # need to take care of the arrays in sqrt!
        self.sqrt = math.sqrt
        self.sqr3 = self.sqrt(3)
        # SPEED: specifically for the random numbers:
        # could generate a large np array and use them one by one.
        self.normal = therng.normal
        self.uniform = therng.uniform
        self.switch = lambda: -1 + 2 * therng.randint(1)
        self._incr = {'gaussian':self.incr_gauss,
                         'uniform':self.incr_uni, 
                         'switch':self.incr_switch}[inc_type]
        # SPEED: hot!
#        self.deti = (lambda t, x, dt: 
#                     dt * self.M(t, x) * self.force_field(t, x) )
        # sdev of stochastic step; no anisotropy implemented
        # SPEED: hot!
        self.sigma = lambda t, x, dt: self.sqrt(2. * dt) * self.B(t, x)
        # reflection
        if ref_box:
            self.ref_box = np.array(ref_box)
            self._prop_or_abs = self._prop_abs_ref
            if self.ref_box.shape == (2,):
                self.reflect = self.reflect_scalar
            else:
                self.reflect = self.reflect_multi
        else:
            self.ref_box = None
            self._prop_or_abs = self._prop_abs_noref
        
        # parent init for the absorption part
        SystemIntegrator.__init__(self, abs_cond)
        
    def incr(self, t, x, dt):
        return self._incr(t, x, dt)
    
    # inlining did not help here; leave as is.
    def deti(self, t, x, dt):
        return dt * self.M(t, x) * self.force_field(t, x) 

    # SPEED: this is hot. optimization would be nice.
    # if we actually do something with the langevin stuff, cythonify here!
    # (would probably want to make an external Ito Cy module, like gillespie)
    def incr_gauss(self, t, x, dt):   # works n-dim; dim is set by x input
        if dt <= 0.:
            return np.zeros_like(x)
        inc = self.deti(t, x, dt)
        ran = self.normal(0., self.sigma(t, x, dt))
        if self.red:
            inc[self.D_inds] += ran
        else:
            inc += ran
        return inc
        
    def incr_uni(self, t, x, dt):   # works n-dim; dim is set by x input
        if dt <= 0.:
            return np.zeros_like(x)
        inc = self.deti(t, x, dt)
        half_width = self.sqr3 * self.sigma(t, x, dt)
        ran = self.uniform(- half_width, half_width)
        if self.red:
            inc[self.D_inds] += ran
        else:
            inc += ran
        return inc

    def incr_switch(self, t, x, dt):   # works n-dim; dim is set by x input
        if dt <= 0.:
            return np.zeros_like(x)
        inc = self.deti(t, x, dt)
        size = self.sigma(t, x, dt)
        ran = self.switch() * size
        if self.red:
            inc[self.D_inds] += ran
        else:
            inc += ran
        return inc

    def reflect_scalar(self, x, mi, ma):
        '''reflect the particle back to the box
        
        mi and ma are the lower and upper bounds, respectively
        (which may be inf). only scalars are supported
        '''
        #return 2 * np.clip(x, *box) - x # this is slow.
        if x < mi:
            return 2 * mi - x
        if x > ma:
            return 2 * ma - x
        return x

    def reflect_multi(self, x, mis, mas):
        '''reflect the particle back to the box
        
        mis, mas are the lower and upper bounds of shape (k,) respectively.
        where k = len(x)
        All coordinates are reflected; possibly some entries of box
        have to be set to +-inf. 
        '''
        xr = np.clip(x, mis, mas)
        xr *= 2
        xr -= x
        return xr
    
    def _prop_abs_noref(self, tfinal, tinit, xinit):
        # checked: localizing incr it is not helpful
        incr = self.incr  # Gaussian random numbers by default
        delta = self.delta
        xcur = cp.deepcopy(xinit)  # xinit is mutable!! prevents modification
        tcur = tinit
        if tfinal == tinit:
            return (tfinal, xcur)
        while tcur < tfinal - delta:
            xcur += incr(tcur, xcur, dt=delta)
            if self.abs_cond(tcur, xcur):   # always
                raise Absorption(tcur, xcur)
            tcur += delta
        xcur += incr(tcur, xcur, dt=(tfinal - tcur))  # final incr
        if self.abs_cond(tcur, xcur):   # always
            raise Absorption(tcur, xcur)
        return (tfinal, xcur)

    def _prop_abs_ref(self, tfinal, tinit, xinit):
        # checked: localizing incr it is not helpful
        incr = self.incr  # Gaussian random numbers by default
        delta = self.delta
        xcur = cp.deepcopy(xinit)  # xinit is mutable!! prevents modification
        tcur = tinit
        if tfinal == tinit:
            return (tfinal, xcur)
        while tcur < tfinal - delta:
            xcur += incr(tcur, xcur, dt=delta)
            if self.abs_cond(tcur, xcur):   # Absorption always: needed now!
                raise Absorption(tcur, xcur)
            xcur = self.reflect(xcur, *self.ref_box)
            tcur += delta
        xcur += incr(tcur, xcur, dt=(tfinal - tcur))
        if self.abs_cond(tcur, xcur):   # last step: check again.
            raise Absorption(tcur, xcur)
        xcur = self.reflect(xcur, *self.ref_box)
        return (tfinal, xcur)
  

class ItoEulerImplicit(ItoEulerFwd):
    """
    An extension to the ItoEulerFwd propagator; 
    allows to enlarge the state x to contain trajectory functionals which are 
    updated online, based on increments of the actual state variables.
    
    The update goes in two steps. First we update based on the pre-point data:
    dx_1 = M(t,x).f(t,x) dt + sqrt(2) B(t,x).dw(t)
    then we update F:
    dx_2 = ff(t, x, dt, dx_1)
    in total, we have
    dx = dx_1 + dx_2
    
    The typical application is where the new state vector x = (x', F). 
    Here, e.g. F = \int ff dx(t). Then dx_1 = dx' updates only x, and 
    dx_2 = dF updates only F, and dx_2 is independent of dt.
    If on the other hand F = \int ff(x(t)) dt, then dx_2 is independent of 
    dx_1.
    In particular, the force and D will have 0 in the F components.

    :fun_increment: a callable.
        fun_increment(t, x, dx_1) where x = (x', F) and dx_1 is the first-step
        increment, returns dx_2
    """
    def __init__(self, delta, force_field, fun_increment, **kwargs):
        ItoEulerFwd.__init__(self, delta, force_field, **kwargs)
        if not callable(fun_increment):
            raise ValueError("Need an F increment function ff(t, x, dx)")
        if not self.red:
            log_.warning("The updated functional F has its own noise?!")
        # reshuffle the increment function
        self.funi = fun_increment
        
    def _incr1(self, t, x, dt):
        return ItoEulerFwd.incr(self, t, x, dt)
    
    def incr(self, t, x, dt):
        dx_1 = self._incr1(t, x, dt)  # regular Ito increment
        dx_2 = self.funi(t, x, dt, dx_1) # the second step
        return dx_1 + dx_2


# for deterministic dynamics; no noise, but chaos could be enough
class EulerFwd(SystemIntegrator):  # pragma: no cover
    def __init__(self):
        raise NIError
    
# better for deterministic dynamics; no noise
class RungeKuttaFwd(SystemIntegrator):  # pragma: no cover
    def __init__(self):
        raise NIError

#-------------------------------- adaptation of integrators for crossing events

class CrossingIntegrator(SystemIntegrator):
    """Something which throws Crossing exceptions.

    This happens on crossing events specified by the input coordinate functions
    and interface locations. In this base class, just a do-nothing propagator
    is implemented; it checks for crossing and absorption at the initial time.
    
    :cross_int: interval for checking a crossing
    :coord_funs: k-tuple of 'interfaced' coordinates as f'ns of t, x
    :interfaces: k-tuple of the corresponding interface values
    :flux_directions:  accepted flux direction specification.
    :time_bins: bin edges for time.
    
    """
    def __init__(self,
                 cross_int,
                 coord_funs=(),
                 interfaces=(),
                 flux_directions=None,
                 time_bins=None):
        self.c_funs = coord_funs
        self.if_fams = tuple(np.array(iflist) for iflist in interfaces)
        self.if_nos = tuple(if_fam.shape[0] for if_fam in self.if_fams)
        self.k_fluxes = len(self.c_funs)
        # need sane input in the following: (see IFHistogram)
        self.flux_dirs = flux_directions
        # the time binning for generating complete crossing events
        self.t_bins = time_bins
        # crossing check interval
        self.cross_int = cross_int
        
    # another copy of this function, when __init__ is overridden
    init_cross = __init__
        
    def init_from_ifhist(self, ifh):
        if hasattr(self, 'initialized_from_ifhist'):
            raise BinningError('attempt to re-initialize' 
                               ' CrossingIntegrator. Make a copy instead.')
        if ifh.c_funs == ():
            raise BinningError('no coordinate function info in {0}...'
                               .format(ifh))
        self.init_cross(coord_funs=ifh.c_funs,
                        interfaces=ifh.if_fams,
                        flux_directions=ifh.flux_dirs,
                        time_bins=ifh.t_bins,
                        cross_int=ifh.cross_int)
        self.initialized_from_ifhist = True
        
    def _cross(self, pcur, iprev):
        """Check if crossings have occured; then throw an exception.

        Called with the current pt and previous index to avoid 
        computing indices twice.
        """
        # current indices
        icur = self._get_inds([c_f(*pcur) for c_f in self.c_funs])
        try:
            # only t_ind_right works here; otherwise loads of complaints
            # at the final time.
            t_i_cur = self._get_t_ind_right(pcur[0])  
        except BinningError as e:
            # do not raise a Crossing if time binning fails.(e.g. out of range)
            if log_.isEnabledFor(logging.DEBUG):
                log_.debug(e)
            return icur
        # the bin changes as an array.
        idelta = icur - iprev
        if not idelta.any():   # nothing to do:
            return icur
        # this reports leaps over two or more barriers
        # attention: we should properly handle traversal by two bins 
        # on the next if.
        jump = (np.abs(idelta) > 1).any()
        if jump:
            if log_.isEnabledFor(logging.DEBUG): 
                log_.debug('A multiple interface jump has occurred')
        # detect crossings and collect Crossing data
        # all crossings are counted; no direction testing here.
        cr = None
        try:
            for dir, delta in enumerate(idelta):
                # only on an actual crossing anything happens:
                if delta: # note that for multi crossings, |delta| > 1!
                    # bwd crossing?
                    bwd = delta < 0
                    # the crossed interface number (True==0 etc.)
                    n_i = iprev[dir] + delta + bwd
                    jij = icur.copy()
                    jij[dir] = n_i
                    for djr, ind in enumerate(jij):
                        # we should cancel the Crossing...
                        if djr != dir:
                            if 0 > ind or ind >= self.if_nos[djr] - 1:
                                #...if the other djr are out of binned range
                                raise BinningError
                    #  we have a crossing within range. is it valid?
                    dtjijb = (dir, t_i_cur) + tuple(jij) + (int(bwd),)
                    valid = self.flux_dirs[dtjijb]
                    # only the lowest crossing direction dir is detected...
                    cr = Crossing(pcur, dtjijb, valid, multi=jump)
                    if valid:  # directly raise on the first valid crossing
                        #log_.warning("CROSSING!")
                        raise cr
        except BinningError:
            pass  # getting out of the inner loop
        if not cr is None:
            # there were only invalid crossings. 
            # raise the last one (arbitrary)
            raise cr
        # no crossing attempt was successful:
        return icur

    @staticmethod
    def _bin_search(array_, value):
        """Binary bin search.
        
        Returns index i for the half-open interval [ array[i], array[i+1] )
        Meaningful bin indices are 0 to len(array) - 2. 
        Returns -1 and len(array) - 1 for out-of-range left and right values.
        """
        return np.searchsorted(array_, value, side='right') - 1 
    
    @staticmethod
    def _bin_search_right(array_, value):
        """Binary bin search.
        
        Returns index i for the half-open interval ( array[i], array[i+1] ]
        Meaningful bin indices are 0 to len(array) - 2.
        Returns -1 and len(array) - 1 for out-of-range left and right values.
        """
        return np.searchsorted(array_, value, side='left') - 1 

    # get the bin numbers; if outside range
    def _get_inds(self, coords):
        return np.array([self._bin_search(fam, co) for fam, co in
                         zip(self.if_fams, coords)], dtype=int)
        
    def _get_inds_right(self, coords):
        return np.array([self._bin_search_right(fam, co) for fam, co in
                         zip(self.if_fams, coords)], dtype=int)

    def _get_t_ind(self, t):
        if not self.t_bins[0] <= t < self.t_bins[-1]:
            raise BinningError('time {0} out of t_bins'.format(t))
        return self._bin_search(self.t_bins, t)


    def _get_t_ind_right(self, t):
        if not self.t_bins[0] < t <= self.t_bins[-1]:
            raise BinningError('time {0} out of t_bins'.format(t))
        return self._bin_search_right(self.t_bins, t)


    # augmenting _prop_or_abs with crossing event generation.
    # SPEED: this is potentially slow; introduces a new 
    # time scale cross_int.
    def _prop_abs_cross(self, tfinal, tinit, xinit, cross_int=None):
        cross_int = cross_int or self.cross_int
        # loop init
        # @todo: the copy operations here are probably unneccessary since 
        # the _prop_or_abs functions copy their input to prevent modification.
        tprev, xprev = tinit, cp.deepcopy(xinit)  # xinit is mutable.
        # this should catch the marginal case of crossing at the initial time
        # has to be the same as the one in _cross! beware!
        iprev = self._get_inds([c_f(tprev, xprev) for c_f in self.c_funs])
        # loop: propagate and check
        while tprev < tfinal - cross_int:
            tcur, xcur = self._prop_or_abs(tprev + cross_int, tprev, xprev)
            iprev = self._cross((tcur, xcur), iprev)  #  raise?
            tprev, xprev = tcur, cp.deepcopy(xcur)
        # propagate to final t and x.
        tcur, xcur = self._prop_or_abs(tfinal, tprev, xprev)
        # final check to prevent overspill
        if tcur < tfinal:
            self._cross((tcur, xcur), iprev)
        return tcur, xcur

    # override SystemPropagator.int to include crossings
    def int(self, tfinal, tinit, xinit):
        """This may throw an Absorption or a Crossing!"""
        return self._prop_abs_cross(tfinal, tinit, xinit)

#-------------------------------------------------- crossing propagator factory

def cross_wrap(SysIntCls):
    """
    Factory function to generate a CrossingIntegrator subclass from a
    standard SystemIntegrator class.

    Takes the class of SystemIntegrator and construction arguments
    of CrossingIntegrator.

    """
    # this is tricky multiple inheritance. since sys_int_cls is left,
    # the methods of SystemIntegrator that it overrides, should
    # stay overridden even though CrossingIntegrator also inherits them
    class CrossWrap(SysIntCls, CrossingIntegrator):
        def __init__(self, *int_args, **int_kwargs):
            SysIntCls.__init__(self, *int_args, **int_kwargs)
#        def init_cross(self, *cr_args, **cr_kwargs):
#            CrossingIntegrator.__init__(self, *cr_args, **cr_kwargs)
    # note: the class definition of CrossWrap is local to this
    # function; the class CrossWrap is not accessible outside.
    # a later call of cross_wrap with different arguments
    # does not change the methods of the class returned here!
    CrossWrap.__name__ = 'Cross' + SysIntCls.__name__
    return CrossWrap  # the class!


#------------------------------------------------------------------- exceptions

class Absorption(Exception):
    def __init__(self, t, x):
        self.point = (t, x)
        if log_.isEnabledFor(logging.DEBUG):
            log_.debug("An absorption has occured")

class Crossing(Exception):
    """
    Interface crossing; contains direction and number
    """
    def __init__(self, point, full_inds, valid=True, **flags):
        # put state info
        self.point = point
        # flux direction, time bin index
        self.dir, self.t_i = full_inds[:2]
        # lower bin inds, index of the crossed interface, higher bin inds
        self.jij = full_inds[2:-1]
        # fwd is True if the crossing was forward, False if backward
        self.fwd = not bool(full_inds[-1])
        # a branch-generating crossing?
        self.valid = valid
        # full indices appropriate for histograms
        self.full_indices = full_inds
        # other flags 
        self.__dict__.update(flags)
        if log_.isEnabledFor(logging.DEBUG):
            log_.debug('Crossing event.'
                       '\n  direction:  {0}'
                       '\n forward: {1}'
                       '\n  jij coord bins: {2}'
                       '\n valid: {3}' 
                       '\n other flags set: {4}'
                       .format(self.dir, self.fwd, self.jij, self.valid,
                               flags.keys()))

class BinningError(Exception):
    pass







