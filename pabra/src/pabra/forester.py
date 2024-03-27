'''
Created on Mar 29, 2010

@author: nbecker

This is the forester class, which steers simulations. it calls 
BranchingStrategies, which are imported from separate modules;
only brute force integration is kept in here.
'''
from __future__ import division

from ptree_h5 import construct_dtype

from ptree_h5 import Pforest  #@UnusedImport; for import in scripts 


# logging; overall level and format controlled in __init__.py
import logging
from pabra import NIError
log_ = logging.getLogger(name=__name__)
#


class BranchingStrategy(object):
    """A generic strategy.

    Takes a forest and gives back a set of branching propensities
    as a function of weight, time and state, applicable for a new tree in the
    forest. In this way, forest statistics can be incorporated iteratively.

    """
    def __init__(self, log_int=None):
        self.log_int = log_int
#        self.schools = {}
        # strategy level counters
        self.sim_tree_count = 0
        self.sim_time = 0
        self.school_counter = -1

    def _init_school(self, new_school):
        # note: each school is an extra forest, not rooted at the
        # main forest root!
        self.school_counter += 1
        log_.info('Making new tree school no. {1}, log_int {0}'
                  .format(self.log_int, self.school_counter))
#        self.schools[self.school_counter] = new_school
        # the new school should know about its index:
        new_school.counter = self.school_counter
        # number of started trees = total injected weight
        (new_school.new_tree_ct,
         new_school.orig_tree_ct) = self.sim_tree_count, self.sim_tree_count
        # do not return the newly created school

    def simulate_school(self):
        raise NIError()

    def update_from_school(self):
        """Re-integrate the new collected statistics"""
        raise NIError()



class BruteForce(BranchingStrategy):
    """Strategy to do brute-force simulation and forest update

    :forest: The Pforest that the strategy is applied to.

    """
    def __init__(self, log_int, t_final, t_initial=0):
        self.t_final = t_final
        self.t_initial = t_initial
        BranchingStrategy.__init__(self, log_int)
        log_.info('Making new ' + self.str_rich())

    def str_rich(self):
        description_str = ('BruteForce strategy.\n' +
                           '\n'.join('  ' + k + ': ' + str(self.__dict__[k])
                                     for k in
                                       ['t_initial',
                                        't_final',
                                        'log_int']))
        return description_str

    # note: this is the old api. not used anymore in IFHistogram!
    def simulate_school(self, school, propa, t_final=None, t_initial=None):
        """Propagate one school, using the propagate.SystemIntegrator propa.

        No branching!

        """
        if t_final is None:
            t_final = self.t_final
        if t_initial is None:
            t_initial = self.t_initial
        school.grow_tips(propa, t_final)
        school.cap_tips()
        log_.debug('just advancing time to {0}'.format(t_final))
        # a python attribute of the school
        school.new_tree_ct += school.population  # one bunch of new trees done.
        log_.info('school no. {0} got a total of {1} new trees'.
                  format(school.counter,
                         school.new_tree_ct - school.orig_tree_ct))

    def update_from_school(self, school):
        self.sim_tree_count += school.new_tree_ct - school.orig_tree_ct
        self.sim_time += sum(tr.tip[0] - tr.base[0] for tr in school
                             if tr.tip and tr.base)



class Forester(object):
    """The class in charge of building up the simulation.

    forest: a Pforest; this should be at a file root.
    propagator: a propagate.SystemIntegrator
    strategy: a BranchingStrategy
    intitial_gen: a function which returns an initial point;
        initial_gen() -> x which may be deterministic or drawn from some
        time-0 starting equilibrium distribution.

    """

    #TODO. first, write the actual updating method to see where it belongs.
    def __init__(self, forest, propagator, strategy, initial_gen,
                 bunch_size=1):
        self.forest = forest
        self.propagator = propagator
        self.strategy = strategy
        # set the log_int at the root of the forest
        self.forest.log_int = self.strategy.log_int
        # automatic initialization for crossing propagators.
        # this is only necessary if the strategy is IFHist
        try:
            self.propagator.init_from_ifhist(self.strategy)
        except AttributeError:
            log_.error('Cross wrapping not applicable for propagator '
                      + str(self.propagator) +
                      '\n  and/or strategy ' + str(self.strategy), exc_info=0)
            pass
        self.init_gen = (initial_gen if callable(initial_gen)
                         else lambda: initial_gen)
        # construct and assign the data type 
        # (assuming the log int is the same type as the time itself)
        self.tx_dtype = construct_dtype((self.strategy.log_int, 
                                         self.init_gen()))
        self.forest.tx_dtype = self.tx_dtype
        # how many concurrent trees?
        self.bunch_size = bunch_size
        # no need to set the log_int, since
        # the actual growing happens in the school anyway.
        # self.forest.log_int = self.strategy.log_int

    #==========================================================================
    # Now the steering function. What it needs to do is:
    # add (a) new tree(s) to school
    # for each new tree:
    # grow the tree while logging the path
    # in some branching interval or on some other trigger, initiate branching.
    #==========================================================================

    def init_run(self):
        sc = self.forest.add_school()
        for dmy in range(self.bunch_size):
            sc.add_tree(logwt=0.,
                    base=(self.strategy.t_initial, self.init_gen()) )
        self.strategy._init_school(sc)
        log_.info('now simulating school no. {0}, {1} trees'\
              .format(sc.counter, sc.population))
        return sc
    
    def cont_run(self, sc, **kwargs):
        # ok, now we go:
        self.strategy.propagate_school(sc, self.propagator, **kwargs)    
        return sc
    
    def finish_run(self, sc):
        # cap
        self.strategy.finalize_school(sc)
        # update
        self.strategy.update_from_school(sc)
        # log
        log_.info('added one school with {0} trees; max height {1}'.
                  format(sc.population, sc.height))
        return sc

    def run(self, **kwargs):
        """Add a school of new trees, propagate and branch each one
        according to the strategy. Finally update the strategy data and
        the simulation trajectories and delete the school.

        """
        sc = self.init_run()
        sc = self.cont_run(sc, **kwargs)
        return self.finish_run(sc)


class Glossary:
    """A class which holds a dictionary of the use of words to describe the
    elements of a Ptree etc.

    """
    # TODO
    pass


class UpdateError(Exception):
    pass

# remove
class ForesterInfo(Exception):
    pass

