'''
Created on Jul 1, 2010

@author: nbecker


The extended FFS sampling strategy. The basic idea is to branch only
at interfaces defined in one or a few coordinates; branching happens on
crossing the interface, so it is event-based. Branching decisions are
based on the accumulated histograms on the interfaces, in bins of the
remaining dimensions.
'''

#==============================================================================
#Guidelines:
#
#The target weight at each crossing event is set based on the following
#considerations:
#
#1. The plane families control the component of the flux normal to the plane
#direction.
#2. Equalization of flux between i-planes makes the i component of the flux,
#integrated over the other directions, uniform along the i direction
#3. Time plays a special role since the integrated flux in the time direction
#is always +1 at each time, _except_ for absorption and creation.
#4. Plane intersections define bins on both sets of planes;
#5. Time is binned, not sectioned, because of its special role.
#
#so we end up with:
#d space dimensions, 1 time dimension,
#k <= d families of interfaces which are not normal to time
#bin edges on each interface wherever it intersects another one
#extra time bin edges defined separately

#We then have a target weight rule for crossing a 'patch' on plane family i,
#which is a bin with edges in time and in the other plane family dimensions j.
#
#there are several possible rules for enforcing interface crossing :
#
#alternative 1: 'one pass across each interface of each family
#per start, over the whole time course'. this is then distributed
#over the patches according to their relative area (density) or just
#patch counts (counts)
#
#alternative 2: one could take: 'one pass across one of the interfaces
#of each (or any) family per time bin.'
#then again distribute according to area or counts.
#
#alternative 3: one could adjust the global target weight such that at the
#injection at time=0, the target weight of the next barrier crossing is 1
#this avoids frequent abort at injection and a resulting high starting
# weight of the survivors.
#
#in density mode, we get a flat flux density of flux component i
# across each interface in the t and j directions
#
#In the i direction itself this is neither possible nor necessary.
#==============================================================================

from __future__ import division

import numpy as np

import branching
from forester import BranchingStrategy

from propagate import Crossing, BinningError

import logging
import math
log_ = logging.getLogger(name=__name__)

import ptree_h5 as p5
Pforest = p5.Pforest

from collections import deque

# choice of prellberg-like reweighting
def rs_betas_default(ratio, n_max):
    # this turns off the minimum ratio
    return branching.rs_betas_1(ratio, n_max, n_min=0.)
    # this turns on the minimum ratio
    #return branching.rs_betas_1(ratio, n_max)


class IFHistogram(BranchingStrategy):
    """Strategy to iteratively generate equalized flux in specified directions

    :log_int: logging interval for recording trajectories.
    :time_bins: sequence of time bin boundaries. Minimum length is 2.
    :t_initial, t_final: optional bounds for the simulated duration;
        otherwise the outer time bin borders are taken.
    :coord_fns: A k-tuple of functions of state (t, x)
        (combinations of coords, possibly time-dependent). optional.
    :interfaces: A k-tuple of arrays containing the locations
        of interfaces in the coordinate boundaries.
    :cross_int: the interval for interface crossing event checking.
        needed if using a cross_wrapped SystemIntegrator; not needed
        if the SystemIntegrator has built-in crossing event generation.
    :flux_dirs: a specification flagging those crossing directions which
        trigger branching. default: True ( => [True, True] for fwd and bwd).
        This is upcasted to the shape of the histograms in the standard
        numpy way. E.g. for giving counting dir spec per if. family if there
        are two lambda coords, of which only the first one triggers:
            flux_dirs=np.array([[1, 1], [0, 0]])[:, None, None, None],
        where the 'None's are for time and the two lambda indices, respectively
    :crossings_per_bin: a global proliferation factor which defaults to 1. 
        this sets how many crossings happen per bin per started tree in an
        equilibrated simulation.
    :n_max: An upper bound on the number of branches in each crossing event
    :min_counts: the minimum number of prior crossings in a bin 
        before branching is switched on.
    :school_n_max: An upper bound on the number of new branches per school run
    :uniform: if set to True, weighted fluxes across patches are counted; 
        otherwise the flux is weighted with the patch area to get normalized
        flux density. This makes a difference only for nonuniform bins.
    :flat_number: if True, branching is biased by an extra number dependent 
        factor to generate a flatter histogram
    :size_target: None, or the target tree size for explosion containment, or
        a function of the current tree size, giving the wt modifier
        the default is lambda x: x/size_target
    :weigth_limits: None, or the limits for absolute branch weight. Useful
        for containing a simulation in the initial phase.
    :rs_betas: The rule of choosing child numbers based on weight ratio
    :traverse: Depth first or breadth first tree traversal

    Events and Timescales:
    Histograms on the patches are augmented at each crossing event.
    Absorption is checked at an interval set by the SystemIntegrator.
    Logging of points along the trunk is set by log_int.
    Time bins for histogramming may be uniform or non-uniform.

    """
    def __init__(self,
                 log_int,
                 time_bins=None,
                 t_initial=None, t_final=None,
                 coord_funs=(),
                 interfaces=(),
                 cross_int=None,
                 flux_directions=None,
                 crossings_per_bin=1.,
                 decorrelate=None,
                 rs_betas=None,
                 flat_number=None,
                 n_max=5,
                 min_counts=0,
                 uniform=True,
                 school_n_max=np.inf,
                 size_target=None,
                 size_limits=None,
                 weight_limits=None,
                 traverse='depth'):
        # parent init
        BranchingStrategy.__init__(self, log_int)
        #-------------------------------------------------------------- binning
        # time bounds

        try:
            if time_bins is None:
                self.t_bins = np.linspace(t_initial, t_final, 11)
            else:
                self.t_bins = np.array(time_bins, dtype=np.float)
        except TypeError:
            raise BinningError('Need to give at least time bins or '
                               'initial and final times')
        self.t_initial = t_initial or self.t_bins[0]
        self.t_final = t_final or self.t_bins[-1]
        # number of samples
        self.t_samples = (self.t_bins[-1] - self.t_bins[0]) / self.log_int
        # a simulated-tree count array with the same shape as the time bins
        # sim_tree_counts[ti] = 'started trees which could reach t_bins[ti]'
        self.sim_tree_counts = np.zeros(self.t_bins.shape[0] - 1, dtype=int)
        # set up binning and slicing parameters
        # a k-tuple of interface location arrays
        self.if_fams = tuple(np.array(iflist) for iflist in interfaces)
        # coord functions
        self.c_funs = coord_funs
        if len(self.c_funs) == 0:
            log_.warn("no coord functions given; propagator needs them!")
        # convenience
        self.k_fluxes = len(self.if_fams)  # how many flux directions?
        self.if_nos = tuple(if_fam.shape[0] for if_fam in self.if_fams)  # #ifs
        self.patch_nos = tuple(np.prod([len(fam) - 1 for fam in
                                         _skip(dir, self.if_fams)])
                                            for dir in range(self.k_fluxes))
        # the master histograms. interpretation of the array:
        # The first index is the flux component (direction). If the first
        # dimension has length I, then the whole array is 2+I+1 dimensional;
        # The second index is the time bin.
        # If the first index is i then the 2 + ith index
        # is filled to the end, this is the interface index;
        # all other indices including time have a 0 as last entry.
        # they represent bins on the i-th interface.
        # the last index indicates fwd / bwd crossings.
            # can use np.swapaxes etc to manipulate...
        self.hist_shape = ((self.k_fluxes,)  #  flx components
                           + (self.t_bins.shape[0] - 1,)  # time bins
                           + self.if_nos  # the actual interface dims tuple
                           + (2,))  # last dimension for fwd-bwd
        log_.debug('histogram shape:\n{0}'.format(self.hist_shape))
        # number flux
        self.n_flx = np.zeros(self.hist_shape, dtype=np.int)
        # 'uncorrelated' number flux
        self.u_flx = np.zeros(self.hist_shape, dtype=np.float)
        # probability flux
        self.p_flx = np.zeros(self.hist_shape, dtype=np.float)
        # number flux after branching
        self.a_flx = np.zeros(self.hist_shape, dtype=np.int)
        # weight flux after branching
        self.q_flx = np.zeros(self.hist_shape, dtype=np.float)
        # which events trigger branching:
        # [forward?, backward?]
        self.flux_dirs = ([True, True] if flux_directions is None
                          else flux_directions)  # see property below
        log_.debug('flux direction array:\n {0}'.format(self.flux_dirs))
        # global factor
        self.cpb = crossings_per_bin
        #------------------------------------------------------- cross-sections
        # d - 1 dimensional coord volumes of each interface family
        self.if_c_vols = [np.prod([other_if_f[-1] - other_if_f[0]
                                   for other_if_f in self.if_fams
                                   if other_if_f is not if_f])
                          for if_f in self.if_fams]
        # d - dimensional time - coord volumes of each interface family
        self.if_vols = [(self.t_bins[-1] - self.t_bins[0]) * c_vol
                        for c_vol in self.if_c_vols]
        #-------------------------------------------------------------- options
        # crossing check interval; default: 1/10 of the time bin interval
        self.cross_int = (cross_int or
                          (self.t_bins[-1] - self.t_bins[0])
                          * .1
                          / (self.t_bins.shape[0] - 1))
        # crossing statistics
        self.crossings = 0
        self.multi_crossings = 0
        # branching child number choice
        self.rs_betas = rs_betas or rs_betas_default
        # flattening options
        self.uniform = uniform
        # extra number flattening factor
        self.flat_number = flat_number
        # correlation handling
        self.deco = decorrelate
        # initialization 
        # TODO: this could possibly be handled with the weight hist alone also
        self.min_counts = min_counts 
        # size control
        self.n_max = n_max  # max child no in each branching event
        self.school_n_max = school_n_max  # max child no per school
        if school_n_max < np.inf:
            log_.warning('limiting branching can bias branching '
                         'to favor early times!')
        self.size_target = size_target
        self.size_limits = size_limits
        self.weight_limits = weight_limits
        if self.weight_limits:
            self.weight_limits = np.array(self.weight_limits)
            if not self.weight_limits.shape == (2,):
                raise BinningError('need two weight limits') 
        # traversal
        self.traversal = {'pre':0,"dep":0, 'pos':1, 'bre':1}[traverse[:3]]
        #------------------------------------------------------------------ log
        log_.info('Making new ' + self.str_rich())


    # flux dir formatting handled
    def _blow_up_dirs(self, spec):
        return np.array(np.broadcast_arrays(self.n_flx, spec).pop()
                        )  #  no casting to bool for use in plotting 
    @property
    def flux_dirs(self):
        return self._bcast_flux_dirs
    @flux_dirs.setter
    def flux_dirs(self, spec):
        try:  # blow this up to the histogram shape
            self._bcast_flux_dirs = self._blow_up_dirs(spec).astype(bool)
        except ValueError:
            raise BinningError('Incompatible flux dir spec.')

    # the fwd / bwd fluxes which cause branching, summed
    @property
    def valid_n_flx(self):
        return (self.n_flx * self.flux_dirs).sum(axis= -1)
    @property
    def valid_p_flx(self):
        return (self.p_flx * self.flux_dirs).sum(axis= -1)

    def str_rich(self):
        description_str = ('IFHistogram strategy.\n' +
                           '\n'.join('  ' + k + ': ' +
                                     str(self.__getattribute__(k))
                                     for k in ['t_initial',
                                               't_final',
                                               't_bins',
                                               'c_funs',
                                               'k_fluxes',
                                               'if_nos',
                                               'if_fams',
                                               'flux_dirs',
                                               'hist_shape',
                                               'cross_int',
                                               'uniform',
                                               'deco',
                                               'n_max',
                                               'size_target',
                                               'school_n_max',
                                               'log_int',
                                               't_samples']))
        return description_str

    # the total sim tree count
    @property
    def sim_tree_count(self):
        return self.sim_tree_counts[0]
    @sim_tree_count.setter
    def sim_tree_count(self, dmy_value):
        # this should be called only by parent class __init__
        log_.debug('sim_tree count cannot be set directly; doing nothing')

    @property
    def sim_tree_counts_full(self):
        '''The simulated tree counts, but with the hist_shape

        This is a slow convenience function'''
        full = np.empty(self.hist_shape, dtype=int)
        # trick to broadcast to a sandwiched dimension
        full.swapaxes(1, -1)[:] = self.sim_tree_counts
        return full

    # augment parent: need the school-level histograms
    # the shadow_group needs to be made first; it should be at the first 
    # level below the root of the file, and be named starting with 'school'
    def _init_school(self, new_school, t_final_ind=None):
        # first part:
        BranchingStrategy._init_school(self, new_school)
        # now put the necessary data:
        # histograms
        for l in 'nupaq':
            for no in ['new', 'orig']:
                new_school.__setattr__(no + '_' + l + '_flx',
                               self.__getattribute__(l + '_flx').copy())
        # number of started trees = total injected weight
        # this is needed for the target wt determination.
        # note: these are now arrays of the size of time_bins, 
        # to allow for short time initialization. thus it is possible
        # to combine shorter and longer run length together.
        for no in ['new', 'orig']:
            new_school.__setattr__(no + '_tree_cts', 
                                   self.sim_tree_counts.copy())
        # handle the new branch counting inside the school
        # these are temporary Python attributes!
        # also new_tree_ct (without the 's') is a Python attribute.
        new_school.new_branches = 0
        new_school.no_branch = False
        # started tree counting
        new_school.new_tree_ct += new_school.population  # one bunch of new trees.
        new_school.new_branches += new_school.leaf_no() # starting branches.
        # now according to the final time: t_final_ind = None case ok!
        new_school.new_tree_cts[:t_final_ind] += new_school.population  
        new_school.crossings = 0
        new_school.multi_crossings = 0        

    def _patch_vol(self, dir, t_ind, *jij_inds):
        """
        The k-1+1 dim volume of a patch in the transverse directions and time

        Call with an unpacked element of the return of _get_patches
        """
        if not 0 <= dir < self.k_fluxes:
            raise ValueError("direction index out of range")
        # time
        vol = self.t_bins[t_ind + 1] - self.t_bins[t_ind]
        # the other dimensions; one 1 at position i stays.
        try:
            for djr, fam in _skip(dir, enumerate(self.if_fams)):
                vol *= fam[jij_inds[djr] + 1] - fam[jij_inds[djr]]
        except IndexError:
            print 'djr: {0}, jij_inds: {1}, fam: {2}'.format(djr,
                                                             jij_inds, fam)
            raise
        # ugly, but this seems to be it;
        return vol

    def _patch_ratio(self, school, arriving_weight, dtjij):
        """
        Gives the ratio of incoming weight and the target weight

        :dtjij: the patch dtjij: dir, time_ind, j_bin_inds, if_number,
            j_bin_inds as returned by a list element of _get_patches

        """
        # inherit options
        uniform = self.uniform
        deco = self.deco
        flat_number = self.flat_number
        size_t = self.size_target
        size_limits = self.size_limits
        if not size_limits is None:
            size_limits = np.array(size_limits, ndmin=1)
        weight_limits = self.weight_limits
        # make a fn out of it
        tgt_mod = size_t if callable(size_t) else (lambda x: x / size_t)
        # index handling
        (dir, t_ind), jij_inds = dtjij[:2], dtjij[2:]
        # new trees are already counted before the call, so:
        started_trees = school.new_tree_cts[t_ind]
        # the total crossing weight so far in the school, in the counted direc.
        total_patch_weight = np.dot(school.new_p_flx[dtjij], 
                                    self.flux_dirs[dtjij])
        # new initial-pass handling: update all histograms before this call
        # in propagate_school
        #
        # for decorrelation:
        if deco:
            total_deco_counts = np.dot(school.new_u_flx[dtjij],
                                       self.flux_dirs[dtjij])
        # now compute the ratio.
        # there is an issue with the global proliferation propensity.
        # the no-branching on first crossing approach works only if there is
        # no other global factor. i.e. we are asking for
        # "one crossing per time bin per start"
        # then the question is: should there be an extra modulation by
        # the area fraction of each bin?

        # uniform case
        weight_target = total_patch_weight / started_trees
        # flat density case: weighting with patch surface
        if not uniform :
            # transverse fractional volume
            patch_frac_tj_vol = self._patch_vol(*dtjij) / self.if_vols[dir]
            # get back the same average propensity
            modulation = patch_frac_tj_vol * self.patch_nos[dir]
            # strategy: one crossing of any patch on avg, weighted with
            # frac volume
            weight_target *= modulation
            
        # explicit local flattening of counts
        if flat_number:
            total_patch_number = np.dot(school.new_n_flx[dtjij],
                                        self.flux_dirs[dtjij])
            # the goal is one crossing per patch per start;
            # note that the global count can be adjusted by the tree size 
            # factor
            weight_target *= total_patch_number / started_trees

        # if decorrelated counts are on, there is an extra factor
        # counting the number of decorrelated hits per bin
        # what exactly decorrelated means, is set in self.uc_leaf_wt
        if deco:
            # each decorrelated count is smaller than 1.
            weight_target *= total_deco_counts / started_trees

        # decorate with size control factor:
        if size_t:
            # one of the nastiest bugs ever: integer division!
#            weight_target *= float(max(1., school.new_branches)) / size_target
            # max to repair -1 and 0
            weight_target *= tgt_mod(float(max(1., school.new_branches)))
            
        # global proliferation factor
        weight_target /= self.cpb

        # finally, compute the ratio with this target:
        ratio = arriving_weight / weight_target
        # but not if branching is disabled because of size limits
        if not size_limits is None:
            fb = float(max(1., school.new_branches))
            try:
                if fb > size_limits[1]:
                    ratio = min(ratio, 1.)
                if fb < size_limits[0]:
                    ratio = max(ratio, 1.)
            except IndexError:
                raise branching.BranchingError('invalid size_limits')
            
        # or if branching is disabled because of weight limits
        if not weight_limits is None:
            wmin, wmax = weight_limits
            old_ratio = ratio
            if arriving_weight < wmin:
                ratio = min(ratio, 1.)
            elif arriving_weight > wmax:
                ratio = max(ratio, 1.)
            if log_.isEnabledFor(logging.DEBUG):
                if not old_ratio == ratio: 
                    log_.debug('ratio clipped from {0} to {1} due to'
                               'weight limiting'.format(old_ratio, ratio))
        # logging.
        # see where the extreme cases are localized
        if log_.isEnabledFor(logging.INFO):
            if not ratio > 1e-10:
                log_.info('low ratio.\n'
                             ' ratio = {1}\n'
                             ' wt = {2}\n'
                             ' at dtjij={0}'
                             .format(dtjij, ratio, arriving_weight))
        # log
        if log_.isEnabledFor(logging.DEBUG):
            log_.debug("patch ratio computed.\n"
                       "  flux direction: {0}\n"
                       "  interface: {1}\n"
                       "  time index: {2}\n"
                       "  patch jij indices: {3}\n"
                       "  started trees: {4}\n"
                       "  incoming wt: {5}\n"
                       "  total wt in bin: {6} (incl. incoming)\n"
                       "  target wt: {8}\n"
                       "  ratio inc. weight / target: {9}\n"
                       "  target size: {10}\n"
                       "  current branch size: {11}\n"
                       "  final ratio: {12}\n"
                       .format(
                        dir, dtjij[dir + 2], t_ind, jij_inds,
                        started_trees, arriving_weight, total_patch_weight,
                        uniform, weight_target, arriving_weight/weight_target,
                        size_t, school.new_branches, ratio))
        return ratio

    def uc_leaf_wt(self, lf):
        if not self.deco:
            return 1
        try:
            if self.deco.startswith('rel'):
                return self._uc_relative_wt(lf)
            if self.deco.startswith('tot'):
                return self._uc_total_fraction_wt(lf)
        except AttributeError:
            return self._uc_exp_wt(self.deco, lf)

    def _uc_relative_wt(self, lf):
        # the Prellberg convention: use the fraction of the path length
        tt, tb, ti = lf.tip[0], lf.base[0], self.t_initial
        try:
            return (tt - tb) / (tt - ti)
        except ZeroDivisionError:
            log_.error('branching at zero length?')
            return 1.

    def _uc_total_fraction_wt(self, lf):
        # Since we have a scale for the max length, we can use that
        # to get a scale for the correlation time
        tt, tb, ti, tf = lf.tip[0], lf.base[0], self.t_initial, self.t_final
        return (tt - tb) / (tf - ti)

    def _uc_exp_wt(self, corr_time, lf):
        # scale with an exponentially decaying correlation
        tb, tt = lf.base[0], lf.tip[0]
        if tb <= self.t_initial:
            return 1.
        rescaled_time = (tt - tb) / float(corr_time)
        return 1 - math.exp(-rescaled_time)
    
    def _count_crossing(self, school, lf, crossing):
        """
        Histograms are updated on the fly. Call before actually branching

        :lf: the leaf which is crossing
        :crossing: a Crossing event, containing indices and direction.

        """
        # where?
        cr_ind = crossing.full_indices  # last index cast to int
        # direct update of the histograms
        # update now with the pre-branch data: 
        # the initial case of empty bins is then taken care of automatically.
        # counts
        school.new_n_flx[cr_ind] += 1
        # system weights
        school.new_p_flx[cr_ind] += lf.wt
        # for the uncorrelated count, add the uncorrelated fraction
        uc_wt = self.uc_leaf_wt(lf)
        school.new_u_flx[cr_ind] += uc_wt
        if log_.isEnabledFor(logging.DEBUG):
            log_.debug('n flux counts added: {0}; u flux wt added: {1}; '
                       'p flux wt added: {1}'.format(1, uc_wt, lf.wt))

    def _count_exit(self, school, lvs, crossing):
        """
        Exit count histogram, updated on the fly. Call after branching.
        This is not done for the weights since that should be the same 
        as the entering weights

        :lf: the leaf which is crossing
        :crossing: a Crossing event, containing indices and direction.

        """
        # where?
        cr_ind = crossing.full_indices  # last index cast to int
        # how many?
        no = len(lvs)
        # how much?
        twt = sum(lf.wt for lf in lvs)
        # update with the exiting number: 
        school.new_a_flx[cr_ind] += no
        # after-branching weight
        school.new_q_flx[cr_ind] += twt
        if log_.isEnabledFor(logging.DEBUG):
            log_.debug('a flux counts added: {0}; q flux added: {1}'.
                       format(no, twt))


    def _branch_on_crossing(self, school, lf, crossing, rs_betas=None):
        """
        The main branching function.

        :lf: the leaf which is crossing
        :crossing: a Crossing event, containing indices and direction.

        """
        # where?
        cr_ind = crossing.full_indices
        dtjij = cr_ind[:-1]
        if school.no_branch:
            if log_.isEnabledFor(logging.DEBUG):
                log_.debug('1-child branch. (school threshold reached)')
            c_no, rewt = 1, 1.
        else:
            # propensities generator
            rs_betas = rs_betas or self.rs_betas
            p_ratio = self._patch_ratio(school, lf.wt, dtjij)
            # inhibit proliferation only, if below count limit:
            if (self.min_counts and 
                p_ratio > 1 and 
                np.dot(school.orig_n_flx[dtjij], self.flux_dirs[dtjij]) 
                    < self.min_counts):
                c_no, rewt = 1, 1.
                if log_.isEnabledFor(logging.DEBUG):
                    log_.debug('1-child branch. (not yet {0} counts)'
                               .format(self.min_counts))
            else: 
                rs, betas = rs_betas(self._patch_ratio(school, lf.wt, dtjij),
                                  self.n_max)
                c_no = branching.branch_number(betas)
                rewt = rs[c_no]
                if log_.isEnabledFor(logging.DEBUG):
                    log_.debug('_branch_on_crossing call.\n'
                               '  branching probabilities: {0}'.format(betas))
        # branch always; branching decision outside this function.
        if log_.isEnabledFor(logging.DEBUG):
            log_.debug('child no.: {0}'.format(c_no))
        # _including_ the n=1 branching events:
        c_list = lf.branch(c_no, rewt, point=crossing.point)
        # finalize: return all children;
        # (the new_branch count is handled in simulate_school)
        return c_list


    def propagate_school(self, school, cross_propa, t_final=None,
                        traversal=None, rs_betas=None, school_n_max=None):
        # initialize / inherit
        if school_n_max is None:  # allow for setting school_n_max = 0
            school_n_max = self.school_n_max
        if t_final is None:
            t_final = self.t_final
            to_end = True
        else:
            to_end = False
        if traversal is None:
            traversal = self.traversal
        # branching propensities
        rs_b = rs_betas or self.rs_betas        
        lvs_dq = deque(school.leaves())  # for access from both sides
        while lvs_dq:  # this still uses the height update!
            lf = lvs_dq.popleft()
            try:
                lf.grow_trunk(cross_propa, t_final, log_int=self.log_int)
                # here, we have made it to t_final or to absorption.
                # if we simulate to the end and do depth first,
                # take out the finished branches since they cannot grow.
                if lf.absorbed or (to_end and not traversal):  
                    school.new_branches -= 1
                lf.t_data.flush()  # 
                # nothing to be appended to lvs_new
            except Crossing as cr:
                school.crossings += 1
                school.multi_crossings += cr.multi
                # hard switch on max
                if (not school_n_max is np.inf and not school.no_branch and
                    school.new_branches > school_n_max):
                    school.no_branch = True
                    log_.warning('setting no_branch flag after {0} new'
                                 ' branches'.format(school.new_branches))

                # the actual branching;
                # this works since lf is still local in exception handling
                # first, the counts:
                self._count_crossing(school, lf, cr)
                # now, the branching
                if not cr.valid:
                    if log_.isEnabledFor(logging.DEBUG):
                        log_.debug('Not branching on inactive crossing type.')
                    # count one exiting branch
                    self._count_exit(school, [lf], cr)
                    # push the old branch back where it came from
                    lvs_dq.appendleft(lf) 
                else:
                    new_brs = self._branch_on_crossing(school, lf, cr, rs_b)
                    # count exiting children
                    self._count_exit(school, new_brs, cr)                    
                    # append the new branches to lvs_new:
                    if new_brs:  # save fn calls for pruning
                        if traversal:
                            lvs_dq.extend(new_brs)  # 'rolling' -> postorder
                        else: 
                            lvs_dq.extendleft(new_brs) # 'stack' -> preorder
                    # new_branches update: #children-1. This should give
                    # the number of alive tips at each time.
                    # TODO: rename?
                    school.new_branches += len(new_brs) - 1
                    lf.t_data.flush()
        pass
                
    def finalize_school(self, school):
        school.cap_tips()
        try:
            multi_ratio = school.multi_crossings / school.crossings
        except ZeroDivisionError:
            log_.warning("no crossing detected in school no. {0}!".
                         format(school.counter))
            multi_ratio = np.nan
        log_.info('school no. {0} got a total of {1} new branches;'
                     '\n multi-crossing ratio: {2}'.
                  format(school.counter, school.new_branches, multi_ratio))       

    def update_from_school(self, school):
        # dump school data
        for l in 'nupaq':
            school._hist_set(l+'_flx', 
                             school.__getattribute__('new_'+l+'_flx') - 
                             school.__getattribute__('orig_'+l+'_flx'))
#        
#        # we first write the school data persistent:
#        school._hist_set('n_flx', school.new_n_flx - school.orig_n_flx)
#        school._hist_set('u_flx', school.new_u_flx - school.orig_u_flx)
#        school._hist_set('p_flx', school.new_p_flx - school.orig_p_flx)
#        school._hist_set('a_flx', school.new_a_flx - school.orig_a_flx)
#        school._hist_set('q_flx', school.new_q_flx - school.orig_q_flx)

        # now store that back to the strategy
        for l in 'nupaq':
            a_ = self.__getattribute__(l+'_flx')
            a_ += school._hist_get(l+'_flx')
#        self.n_flx += school._hist_get('n_flx')
#        self.p_flx += school._hist_get('p_flx')
#        self.u_flx += school._hist_get('u_flx')
#        self.a_flx += school._hist_get('a_flx')
        # time-bin-wise tree counts
        school._hist_set('extra_tree_cts', 
                         school.new_tree_cts - school.orig_tree_cts)
        self.sim_tree_counts += school._hist_get('extra_tree_cts')
        extralen = school.length  # expensive!
        self.sim_time += extralen
        self.crossings += school.crossings
        self.multi_crossings += school.multi_crossings
        log_.info('updating from school no. {2}'
                  '\n  new trees added: {0}'
                  '\n  sim time used: {1}'.
                  format(school._hist_get('extra_tree_cts')[0], 
                         extralen, school.counter))

#-------------------------------------------------------------------- utilities

# skippy generator!
def _skip(i, iterable):
    for j, elem in enumerate(iterable):
        if j != i:
            yield elem

