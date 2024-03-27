'''
Created on Apr 1, 2010

@author: nbecker
This module contains the branching decision handling based on
the propensities which may depend on position, time and forest statistics.
It also contains the reweighting.
'''

from __future__ import division

import numpy as np
# needed for numpy > 1.3 :
np.seterr(divide='ignore')

# logging; overall level and format controlled in __init__.py
import logging
log_ = logging.getLogger(name=__name__)

# random ng management

therng = np.random.RandomState()
therng.seed(123)

def eq_weights(betas):
    '''Simple reweighting formula.

    r(n) = 1/n 1/sum(beta(n),n >= 1) where beta are the
    branch number weights, normalized incl n=0.
    attention: returns Inf for the wt at n=0
    
    input need not be normalized
    '''
    betas = np.array(betas).astype(float)
    c = betas.sum() / betas[1:].sum()
    return c / np.arange(0,len(betas))
    # dummy weight Inf for 0. Should lead to errors if used! 
    # speed: could store computed betas (there is a technique with a decorator)


def inv_mean_weights(betas):
    '''Even simpler reweighting formula.
    
    Reweighting does not depend on the number of children;
    it's just the inverse of the mean child number.
    
    input need not be normalized.
    '''
    betas = np.array(betas).astype(float)
    betas /= betas.sum()
    mean_c_no = np.arange(betas.shape(0)) * betas 
    return np.ones_like(betas) / mean_c_no


def betas_top_0(navg, nmax=2):
    # TODO speed: could store computed betas
    # (there is a technique with a decorator)
    '''Branch number weights for outcomes 0 and nmax with prescribed mean

    The avg child number is navg but clipped to (0,nmax).
    '''
    if nmax < 1:
        raise BranchingError('need nmax at least 1 in betas')
    bn = np.clip(float(navg), 0, nmax) / nmax
    wlist = np.zeros(nmax + 1)
    wlist[0] = 1 - bn
    wlist[-1] = bn
    return wlist


def betas_top_two(navg, nmax=2):
    '''Branch number weights for outcomes nmax and nmax-1 with prescribed mean

    Must have nmax-1 <= navg <= nmax, otherwise navg is clipped to that range
    '''
    if nmax < 1:
        raise BranchingError('need nmax at least 1 in betas')
    navclp = np.clip(float(navg), nmax - 1, nmax)
    wlist = np.zeros(nmax + 1)
    bn = navclp + 1 - nmax
    wlist[-2:] = [1 - bn, bn]
    return wlist


def branch_number(betas):
    '''Give the random outcome based on a list of branch weights

    '''
    betacum = np.cumsum(betas, dtype=float)
    betacum /= betacum[-1]
#    old and probably slower
#    r = therng.uniform()
#    for i, B in enumerate(betacum):
#        if r < B:
#            return i
    return np.searchsorted(betacum, therng.uniform())


# the different Prellberg-like reweighting choices

# default:


def rs_betas_1(ratio, n_max, n_min=None):
    """
    Gives branching propensities as a list, a variation of the balanced
    reweighting scheme in flatPERM.

    """
    n_min = 1. / n_max if n_min is None else float(n_min)
    # mkI:
    # slight deviation from prellberg: admit two outcomes in the
    # enriching case, nmax and nmax - 1
    # this does include pruning for r < 1.
    # the weights then have to be taken as 1/#children or weight_est
    # in the r > 1 and r < 1 cases, respectively. this comes out
    # automatically if the eq_weights are used.
    #======================================================================
    navg = np.clip(ratio, n_min, n_max)
    betas = betas_top_two(navg, np.ceil(navg))
    #======================================================================
    return eq_weights(betas), betas

# variation with a uniform weight
# actually, this is just the same as rs_betas_target except for the 
# non-occurring child numbers!

def rs_betas_2(ratio, n_max, n_min=None):
    """
    like rs_betas_1 but the reweighting is done with the same factor no matter
    which of the 2 possible results comes out. 
    
    this does not respect per-event conservation of weight even in the 
    enriching case.

    """
    n_min = 1. / n_max if n_min is None else float(n_min)
    # slight deviation from prellberg: admit two outcomes in the
    # enriching case, nmax and nmax - 1
    # this does include pruning for r < 1.
    #======================================================================
    navg = np.clip(ratio, n_min, n_max)
    betas = betas_top_two(navg, np.ceil(navg))
    #======================================================================
    same_wt_for_all = np.ones_like(betas) / navg
    return same_wt_for_all, betas

# exact copy of prellberg scheme.

def rs_betas_prell(ratio, n_max, n_min=None):
    """
    Gives branching propensities as a list; this is exactly
    flatPERM; not reweighting if arriving wt is between 1 and 2 time target

    """
    n_min = 1. / n_max if n_min is None else float(n_min)
    # the version rs_betas_1 always branches even if
    # the ratio is already almost 1; this may be inefficient.
    #
    # mk II:
    # the following is exactly like prellberg; take the floor:
    #======================================================================
    ratio = np.clip(ratio, n_min, n_max)
    n_children = ratio if ratio <= 1 else min(np.floor(ratio), n_max)
    betas = betas_top_0(n_children, np.ceil(n_children))
    #======================================================================
    # this solution does not produce branches for weights
    # between 1 and 2 times the target wt.
    # advantage: multiple branching is self-regulating; however this
    # works only if count updates are instantaneous. (now implemented)
    #
    # disadvantage: no reduction of variance below the factor 2 limit.
    # can get 'relaxation oscillations' in the counts
    return eq_weights(betas), betas


# note: this is essentially a duplicate of rs_betas_2

def rs_betas_target(ratio, n_max, n_min=None):
    """
    Gives branching propensities as a list; resulting wt is always exactly
    on target

    """
    n_min = 1. / n_max if n_min is None else float(n_min)
    # mkIII:
    # this produces _exactly_ the estimated weight in all
    # surviving outcomes. that way, repeated branching events have no
    # effect; no need to take care that only one event occurs per bin.
    # for this, need to take
    ratio = np.clip(ratio, n_min, n_max)
    betas = betas_top_two(ratio, np.ceil(ratio))
    # to fulfill the no-bias condition, need to set the outgoing weight
    # to the target weight, or when hitting n_max, to 1/n_max * incoming wt
    # this is accomplished by 1/ratio in both cases.
    rs = np.zeros_like(betas)
    #rs = np.array((0.,) * (len(betas)-2) + (1. / ratio,) * 2)
    rs[0] = np.inf
    rs[-2:] = (1. / ratio,) * 2
    return rs, betas


def exp_wait_time(rate=1.):
    '''Give an exponential waiting time.

    Use with a propensity function. Note: changes in propensity while waiting
    go unnoticed.

    '''
    return therng.exponential(1. / rate)


def branch_trigger(rate=1., dt=1.):
    '''Branching yes/no in time step dt
    '''
    return bool(exp_wait_time(rate) < dt)


# more general question: what happens when we need to collect branch
# statistics in spatial bins? does this require a revision of the
# program strategy?

class BranchingError(Exception):
    pass
