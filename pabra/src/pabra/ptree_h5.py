'''
Path trees to hold the simulation results and steer the simulation itself.

There is a 'shadow-structure': Ptree is a wrapper class which has all the 
methods for tree inspection, propagation initiation and branching, 
but all data are written into an underlying h5py.Group instance.

The Groups make up the whole tree as an h5py file, which contains the tree
hierarchy. The Ptrees have no references to children or parents, and are 
therefore garbage collectable. They are created on demand when children are 
requested.

The Ptree instances are transient. They get deleted upon completion of the
path segment. Whenever they are deleted, their temporary data get flushed
to the HDF5 file. 

Important: interference with the deletion _will_ lead to runtime errors or
nonsense. In particular, no references from python space to Ptree instances
that are being propagated, are allowed! This is very fragile, but fast.
A safe way (not tested) would be to use weak references:
tr = weakref.proxy(pforest.children[0]) where pforest is a Pforest.


Created on Dec 4, 2010

@author: nbecker
'''
#---------------------------------------------------------------------- imports

from __future__ import division

import numpy as np
import sys
import itertools as it
import h5py as h5
import math
import posixpath as pp

import cPickle as pkl  # for storing dtypes

from collections import deque

from propagate import Absorption, Crossing

# logging; overall level and format controlled in __init__.py
import logging
from weakref import ref as Wref
log_ = logging.getLogger(name=__name__)

#--------------------------------------------------------- for default h5 files

from uuid import uuid1
from tempfile import NamedTemporaryFile
import os
DEFAULT_FILE_IN_MEM = False
TEMP_FILES = []
TEMP_DIR = '/tmp/pabra'
try:
    os.makedirs(TEMP_DIR)
except OSError:
    pass

def _default_file():
    global DEFAULT_FILE_IN_MEM, TEMP_FILES, TEMP_DIR
    if DEFAULT_FILE_IN_MEM:
        return h5.File(str(uuid1()), driver='core', backing_store=False)
    else:
        a = NamedTemporaryFile(suffix='.h5', prefix='for', dir=TEMP_DIR,
                               delete=True)
        TEMP_FILES.append(a)  # to delete later!
        return h5.File(name=a.name, mode='r+')

def purge_temp_files():
    """Call only at the end when Pforests do not need the temp files anymore!
    
    """
    global TEMP_FILES
    for a in TEMP_FILES:
        a.close()  # this should delete them as well.

#---------------------------------------------------------------------- support


# a public helper function to make a dtype out of a list of t and x;
# can be called to make the global dtype
def construct_dtype(tx):
    t, x = map(np.array, tx)
    return np.dtype([('t', t.dtype), ('x', x.dtype, x.shape)])

#--------------------------------------- trunk storage and compression settings

# compression for individual datasets
# 'szip' does not work, the others are typically useless
# since the file structure is not compressed
TRUNK_COMPRESSION = None  # this cannot be changed anymore for non chunked trk.

#------------------------------------------------------------------ the classes



class PtHelpers(object):

    __slots__ = ()  # otherwise a dict is passed down!

    """
    Only a namespace for helper functions 
    
    since these are used to in the evaluation of the Ptree class definitions,
    they need to be outside of that namespace.
    """

    #------------------------------- bootstrap functions for defining PtreeBase
    # these need to be outside of PtreeBase class definition to work.

    @staticmethod
    def _t_pt_gs(name):
        def getter(self):
            try:
                return getattr(self.t_data, name)
            except AttributeError:
                zeros = tuple(np.zeros(shape=(), dtype=self.tx_dtype)[()])
                setattr(self.t_data, name, zeros)
                return zeros
        def setter(self, val):
            self.t_data.dirty = True
            setattr(self.t_data, name, val) 
        return getter, setter
    
    @staticmethod
    def _finalize(ref):
        """Internal weakref callback to flush temp data"""
        if not ref.dirty:
            return
        try:
            ref.flush()
        except Exception:
            log_.error("Exception running {0}:".format(ref.flush), 
                       exc_info=1)
            raise


    #----------------------------------------------------- other static methods

    # depth in the tree based on name parsing
    @staticmethod
    def _name_depth(name):
        if name == '/':
            return 0
        return name.count('/')  # '/s1' -> 1

    # shortcut based on survival at least to t: 
    # if the height attribute is already lower, don't bother
    @staticmethod
    def _shortcut(t):
        if t is None:
            t = np.inf
        a = np.zeros((), float)
        def cond(grp):
            id = grp.id
            if not h5.h5a.exists(id, 'height'):
                return True
            else:
                h5.h5a.open(id, 'height').read(a)
                return a >= t
        return cond

    # running tree ht update

    class _RunningMax(object):
        __slots__ = 'val'
        def __init__(self):
            self.val = -1.
        def update(self, sib_ht):
            if sib_ht > self.val:
                self.val = sib_ht

    @staticmethod
    def _tree_code(name):
        lastpart = name.rsplit('/', 1)[-1]
        return 'r' if lastpart == '' else lastpart[0]
    
    # pruned/absorbed/survived management
    
    @staticmethod
    def _pas_gs(p_a_s):
        def getter(self):
            pas = self.shadow.attrs.get('pas', None)
            return pas == p_a_s  # '==' is crucial, 'is' fails because of dtype
        def setter(self, val):
            self.shadow.attrs['pas'] = p_a_s if val else 'n'
        return getter, setter

#------------------------------- ABCs for organizing the contents thematically.


class _PtreeStructure(PtHelpers):
    """ The basic tree structure incl. traversal and HDF5 attribute access. 
    
    naming convention: subtree names start with 't' or 's' (for schools)
    """
    __slots__ = ('shadow',)
    #---------------------------------------------------------- special methods

    # ATTENTION: PtreeBase nodes are initialized all the time when traversing
    # trees purely for data retrieval. Never call __init__ with extra 
    # attrs in these circumstances: data gets overwritten!

    def __init__(self, shadow_group):
        # the underlying HDF5 data store; must be an h5py.Group instance
        self.shadow = shadow_group


    def __iter__(self):
        return self.preorder()

    # def __sizeof__(self): missing

    def __eq__(self, other):
        return self.shadow == other.shadow

    # def __le__(self, other): etc.: missing; not implemented in h5py


    #----------------------------------------------------------- helper methods

    def _group_keys(self, group=None):
        group = group or self.shadow
        return it.ifilter(lambda k:k.startswith(('t', 's')), group.iterkeys())

    # subtree names start with 't' or 's' (for schools) 
    def _group_items(self, group=None):
        group = group or self.shadow
        for k in self._group_keys(group):
                yield k, group.get(k)

    # this is only to allow the subclasses to create their own instances
    def _make_from(self, group):
        return self.__class__(group)


    #------------------------------------------------------ basic tree handling
    # this is a basic implementation which always gives Nodes of the same type. 
    @property
    def parent(self):
        sh_p = self.shadow.parent
        return self._make_from(sh_p)

    def children_iter(self):
        return (self._make_from(v) for dmy_k, v in self._group_items())

    @property
    def children(self):
        return list(self.children_iter())

    # set what identity children have; overridden in Pforest
    @classmethod
    def _child_class(cls):
        return cls

    def add_child(self):
        # child names have the form 't' + str(int)
        try:
            # throw away the leading group indicator 't' or 's'
            max_ch_n = max(int(k[1:]) for k in self._group_keys())
            new_k = ''.join(('t', str(max_ch_n + 1)))
        except ValueError:
            new_k = 't0'
        sh_ch = self.shadow.create_group(new_k)
        # making also the PtreeBase object; this makes handling the 
        # initialization easier in what follows
        return self._child_class()(sh_ch)

    # def addChildrenFromList(self): not needed.


    #------------------------------------------ basic and conditional traversal

    def _name_path_to_root(self):
        cur_name = self.shadow.name
        while cur_name > '/':
            yield cur_name
            cur_name = pp.dirname(cur_name)
        yield cur_name

    def _group_path_to_root(self):
        shd = self.shadow
        while shd.name > '/':
            yield shd
            shd = shd.parent
        yield shd

    def _group_path_from_root(self):
        l = list(self._group_path_to_root())
        for trb in reversed(l):
            yield trb

    def path_to_root(self):
        return it.imap(self._make_from, self._group_path_to_root())

    def path_from_root(self):
        return it.imap(self._make_from, self._group_path_from_root())

    def get_root(self):
        #check if needed
        return self._make_from(self.shadow.file)

    def _group_preorder(self, *func):
        """preorder traversal of the underlying Groups
        
        Optional shortcut conditions as functions func which take 
        h5.Group arguments (!)
        Whenever conditions are not fulfilled, the entire subtree is skipped.
        This is different from filtering the unconditioned iterator
        """
        # TODO: try if using ids or full names is quicker!
        pass_ = func == ()
        stack = [self.shadow]
        while stack:
            shd = stack.pop()
            if pass_ or all(f(shd) for f in func):
                yield shd
                stack.extend(v for dmy_k, v in self._group_items(shd))

    def _group_postorder(self, *func):
        """as _group_preorder, but postorder
        
        """
        pass_ = func == ()
        stack = deque([self.shadow])
        while stack:
            shd = stack.popleft()
            yield shd
            if pass_ or all(f(shd) for f in func):
                stack.extend(v for dmy_k, v in self._group_items(shd))

    def preorder(self, *func):
        """preorder traversal of the tree; optional conditions
        
        Here, the conditions take PtreeBase nodes as input
        all PtreeBase instances are made in any case, filtering is afterward"""
        it_ = it.imap(self._make_from, self._group_preorder())
        if func:
            it_ = it.ifilter(lambda ptb: all(f(ptb) for f in func), it_)
        return it_

    def postorder(self, *func):
        """like preorder, but in postorder"""
        it_ = it.imap(self._make_from, self._group_postorder())
        if func:
            it_ = it.ifilter(lambda ptb: all(f(ptb) for f in func), it_)
        return it_

class _PtreeTempData(Wref):
    # Ptree will create one such weakref on initialization;
    # this allows temp data to be flushed when the Ptree is gc'd
    
    # ATTENTION: dangerous when more than one Ptree open the same node!
    # no precautions are taken for that case; the respective temp data
    # WILL get mangled!
    
    __slots__ = ('shadow', 't_base', 't_tip', 't_trunk', 'dirty')
    
    def __new__(cls, ptb, shadow, callback):  # __init__ not enough!
        obj = Wref.__new__(cls, ptb, callback)
        return obj

    def __init__(self, ptb, shadow, callback):
        self.dirty = False   # set True by write access to temp data
        self.shadow = shadow
#        if 'lock' in self.shadow.attrs:
#            if not ptb.__class__ is Pforest:
#                raise ValueError("multiple tree node access attempt")
#        self.shadow.attrs['lock'] = 1
        try:
            btt = shadow['d_btt'][:]
            self.t_base, self.t_tip = map(tuple, btt[:2])
            self.t_trunk = map(tuple, list(btt[2:]))
        except KeyError:
            pass
        
    #-------------------------------------------------- flushing temporary data

    def flush(self):
        """flushing the temp data"""
#        try:
#            del self.shadow.attrs['lock']
#        except KeyError: pass
        self.dirty = 0
        dtype_ = pkl.loads(self.shadow.file.attrs['tx_dtype'])
        ttr = getattr(self, 't_trunk', [])
        btt = np.zeros((2 + len(ttr),), order='c', dtype=dtype_)
        btt[:1] = getattr(self, 't_base', np.zeros(1,dtype=dtype_))
        btt[1:2] = getattr(self, 't_tip', np.zeros(1,dtype=dtype_))
        btt[2:] = ttr
        if 'd_btt' in self.shadow:
            # this could potentially be optimized; direct overwrite: no.
            del self.shadow["d_btt"]
        self.shadow["d_btt"] = btt
    

class _PtreeData(PtHelpers):
    """All the basic Ptree Attributes and Datasets needed.
    
    Pforest and Ptree specific ones are in the respective classes
    
    naming convention: Dataset nodes have names starting with 'd_'
    naming convention: strategy related metadata are stored in an
                        extra subgroup node 'h/'
    
    """
    __slots__ = ()  # with t_data, it does not work. -> abstract class only!

    #----------------------------------------------------------------- Datasets

    # naming convention: Dataset nodes have names starting with 'd_'
    # base and tip are stored as one Dataset 'd_ends'
    # this is a 1-dim record array with a 't' and an 'x' field
    # the [0] entry is the base, [1] is the tip

    # actually setting the properties base and tip
    
    base = property(*PtHelpers._t_pt_gs('t_base'))
    tip = property(*PtHelpers._t_pt_gs('t_tip'))

    @property
    def trunk(self):
        try:
            return self.t_data.t_trunk
        except AttributeError:
            emptytrunk = []
            self.t_data.t_trunk = emptytrunk
            return emptytrunk
    @trunk.setter
    def trunk(self, val_l):
        self.t_data.dirty = 1
        self.t_data.t_trunk = val_l

    # it makes more sense to actually count from the base to the tip
    @property
    def trunk_length(self):
        return self.tip[0] - self.base[0]

    # a convenience function:
    @property
    def btrunkt(self):
        res = [self.base]
        res.extend(self.trunk)
        res.append(self.tip)
        return res
    
    #----------------------------------------------------------- h5.Group attrs

    # direct lookup without constructing the dict
    def _get_attr(self, at, default=None):
        return self.shadow.attrs.get(at, default)

    def _set_attr(self, at, val):
        self.shadow.attrs[at] = val

    # the shadow attributes
    @property
    def attrs(self):
        return dict(self.shadow.attrs)

    #--------------------------------------------------- attributes of the tree

    @property
    # the dtype is stored at the root: fixed for the tree!
    def tx_dtype(self):
        return pkl.loads(self.shadow.file.attrs['tx_dtype'])
    @tx_dtype.setter
    def tx_dtype(self, dtype):
        if not self.shadow.name == '/':
            raise ValueError("The dtype can be set only at the root.")
        self.shadow.file.attrs['tx_dtype'] = pkl.dumps(dtype)

    @property
    def log_int(self):
        # looking this up upward
        v = self.shadow.file.attrs.get('log_int', np.inf)
        if v:
            return v
        raise ValueError('something went horribly wrong with log_int')
    @log_int.setter
    def log_int(self, log_int):
        self.shadow.file.attrs['log_int'] = log_int

    # child number retrieval
    @property
    def branch_no(self):
        return sum(1 for dmy in self._group_items())

    @property
    def height(self):
        return self._get_attr('height', np.inf)

    @property
    def alive(self):
        return bool(self.height is np.inf)

    # this may look slow, but a test with h5.Group.visit instead showed
    # no improvement in speed, probably because of python call overhead
    @property
    def size(self):
        cts = 0
        for dmy in self._group_preorder():
            cts += 1
        return cts

    @property
    def depth(self):
        """how deep under the root"""
        return self._name_depth(self.shadow.name)

    @property
    def generations(self):
        """how many generations of descendants; potentially slow!"""
        return max(sub.depth for sub in self) - self.depth


    @property
    def length(self):
        # danger: only works if subtrees have been flushed to disk!
        l = 0.
        for grp in self._group_preorder():
            try:
                bt = grp["d_btt"][:2]
                l += bt[1][0] - bt[0][0]
            except KeyError:
                pass
        return l

    # a variant as an iterator; nicer
    def trunk_from_root_iter(self):
        # if trunk is empty
        def grp_get_trunk(grp):
            try: # no empty selections
                return grp.get('d_btt', [])[2:]
            except ValueError:
                return []
        return it.chain(*(grp_get_trunk(grp) for grp in 
                          self._group_path_from_root()))

    @property
    def trunk_from_root(self):
        return list(self.trunk_from_root_iter())

    # setting a default logwt of 0.
    @property
    def logwt(self):
        return self._get_attr('logwt', 0.)
    @logwt.setter
    def logwt(self, val):
        self._set_attr('logwt', val)

    @property
    def wt(self):
        return math.exp(self.logwt)
    @wt.setter
    def wt(self, val):
        self.logwt = math.log(val)

    def total_weight(self, t=None):
        return sum(l.wt for l in self.leaves(t))


    #----------------------------------------- management of running attributes

    @staticmethod
    def _upd_max(grp, runningmax):
        def fn(sub_name):
            if sub_name.startswith(('r', 's', 't')):
                sib_ht = grp.get(sub_name).attrs.get('height', np.inf)
                runningmax.update(sib_ht)
        return fn


    def _update_height(self, ht):
        # this could streamlined by not making a Group but accessing the
        # attributes directly.
        pathit = it.imap(self.shadow.get, self._name_path_to_root())
        if ht is np.inf:
            # we do not set the height if it was and is np.inf
            selfgrp = pathit.next()  # the starting group
            if selfgrp.attrs.get('height', None):
                del selfgrp.attrs['height']
            for grp in pathit:
                if grp.attrs.get('height', np.inf) == np.inf:  # parent ht
                    return
                try:
                    del grp.attrs['height']
                except KeyError:
                    pass
        else:  # finite ht case: here we have to check the siblings
            pathit.next().attrs['height'] = ht   # the starting group
            for grp in pathit:
                running_max = self._RunningMax()
                h5.h5g.iterate(grp.id, self._upd_max(grp, running_max))
                if running_max.val == np.inf:
                    try:
                        del grp.attrs['height']
                    except KeyError:
                        pass
                    return
                grp_ht = grp.attrs.get('height', np.inf)
                if running_max.val == grp_ht:  # ie. the grp was up to date
                    return
                grp.attrs['height'] = running_max.val  # tolerate np.inf 



class PtreeBase(_PtreeData, _PtreeStructure):
    """
    A class made for inheritance by Ptree and by Pforest
    
    this class inherits the basic structure and the common attributes.

    """

    __slots__ = ('t_data', '__weakref__')
    
    def __init__(self, shadow_group):
        _PtreeStructure.__init__(self, shadow_group)
        # a weakref to self as a member attribute.
        # when self gets collected, weakref's callback 
        # flushes temp data
        self.t_data = _PtreeTempData(self, shadow_group, PtHelpers._finalize)
    #---------------------------------------------- height attribute management

    # now adding the attribute update

    def add_child(self, logwt=0., base=None, tip=None, height=None,
                   **new_attrs):
        # this is a self.__child_class instance.
        ch = _PtreeStructure.add_child(self)
        if logwt:  # extra treatment for easier correspondence to old vsn. 
            ch.logwt = logwt
        for k, v in new_attrs.iteritems():
            ch.__setattr__(k, v)
        # extra treatment for height attr
        if height:
            ch._update_height(height)
        # extra treatment for the Dataset parts, tip and base
        if not base is None:
            ch.base = base
            ch.tip = base if tip is None else tip
        elif not tip is None:
            ch.tip = tip
        # could be useful:
        return ch


    #----------------------------------------------------- time-based traversal

    @property
    def is_leaf(self):
        return self._group_is_leaf()

    # overrides parent class. defined here since this has to be in Pforest too!
    def _group_is_leaf(self, group=None):
        """The branch has a growing tip; should be False for all branches in 
        a finished tree."""
        # this is convoluted but  about 2x as fast as using 
        # h5.Group.attrs['pas'] directly
        id = group.id if not group is None else self.shadow.id
        val = np.empty((), np.str)
        try:
            h5.h5a.open(id, 'pas').read(val) # pruned, survived, or absorbed
            if val[()] in 'pas':  # indexing required; membership uses '=='
                return False
        except KeyError:
            pass
        try:
            self._group_keys(group).next()
            return False
        except StopIteration:
            return True


    # this needs both the height attribute and the basic traversal.

    def _group_preorder_t(self, t=None):
        return self._group_preorder(self._shortcut(t))

    def _group_postorder_t(self, t=None):
        return self._group_postorder(self._shortcut(t))

    def _group_t_leaf(self, t):
        if t is None:
            # this should still work after overriding _group_is_leaf
            def cond(grp):
                return self._group_is_leaf(grp)
        else:
        # note that for not-yet flushed groups, this will yield False
            def cond(grp):
                btt = grp.get('d_btt', [])[:2]
                try:
                    return btt[0]['t'] < t <= btt[1]['t']  # inbetween
                except IndexError:
                    return False
        return cond

    def preorder_t(self, t=None):
        return it.imap(self._make_from, self._group_preorder_t(t))

    def postorder_t(self, t=None):
        return it.imap(self._make_from, self._group_postorder_t(t))

    def _group_leaves(self, t=None):
        return it.ifilter(self._group_t_leaf(t), self._group_preorder_t(t))

    def _group_postorder_leaves(self, t=None):
        return it.ifilter(self._group_t_leaf(t), self._group_postorder_t(t))

    def leaves(self, t=None):
        return it.imap(self._make_from, self._group_leaves(t))

    def postorder_leaves(self, t=None):
        return it.imap(self._make_from, self._group_postorder_leaves(t))

    def leaf_no(self, t=None):
        return sum(1 for dmy in self._group_leaves(t))

    #------------------------------- growth methods common to Ptree and Pforest


    def grow_tips(self, propa, to_t, from_t=None, log_int=None):
        '''Grow all alive leaves up to time to_t.

        propa gives the dynamics, see ptree.propagate.
        giving the start time from_t is supposed to speed up the tip retrieval;
        this helps only if called on the root of a larger tree.

        '''
        try:
            l_i = log_int or self.log_int # ValueError if not set for a parent
        except ValueError as e:  # no problem, write-out is just turned off.
            log_.debug(e)
            l_i = None
        for subptr in self.leaves(from_t):
            subptr.grow_trunk(propa, to_t, l_i)

    def cap_tips(self):
        '''End the growth phase by capping all remaining leaves

        '''
        for ptr in self.leaves():
            ptr.cap()


    #----------------------------------------------------------------- printing

    def dump(self, outf=sys.stdout):
        for i in self.preorder():
            print >> outf, "  "*(i.depth), i


class Pforest(PtreeBase):

    # no slots: no memory issues here.
    # BranchingStrategies may use the __dict__ to store school data.

    # ATTENTION: when referring to a Pforest as a parent, 
    # an empty Pforest instance is made; the dict is not recovered, since
    # not saved on disk

    #--------------------------------------------------------------------- init

    def __init__(self, shadow_group=None, **attrs):
        # no shadow_group: make a temporary h5 file
        if shadow_group is None:
            shadow_group = _default_file()
        PtreeBase.__init__(self, shadow_group)
        for k, v in attrs.iteritems():
            self.__setattr__(k, v)
        # we add a Group node for storing the strategy-dependent metadata
        shadow_group.require_group('h')

    #----------------------------------- attrs and data in the strategy subnode

    # persistent strategy metadata are in a subgroup 'h/'
    # do not use for frequent updates: _slow_ !
    
    def _hist_get(self, flxname, default=[]):
        return self.shadow.get('h/' + flxname, default)[:]
#        return self.shadow['/'.join(('h', flxname))]

    def _hist_set(self, flxname, val_ar):
        self.shadow.require_group('h')
        self.shadow.require_dataset('h/' + flxname,
                                    val_ar.shape, val_ar.dtype)
        self.shadow['/'.join(('h', flxname))][:] = val_ar

    def _hist_get_attr(self, at, default=None):
        hist = self.shadow.get('h', default)
        if hist is default:
            return default
        return hist.attrs.get(at, default)

    def _hist_set_attr(self, at, val):
        self.shadow.require_group('h')
        self.shadow['h'].attrs[at] = val

    #----------------------------------------- making new tree schools Pforests

    def add_school(self, **new_attrs):        
        # iterating over school names takes too long: O(n^2) in the 
        # population
        try:
            self.last_tree_index += 1
        except AttributeError:
            self.last_tree_index = 1
        new_k = 's' + str(self.last_tree_index)
        sh_ch = self.shadow.create_group(new_k)
        # making the school Pforest object; 
        sc = self.__class__(sh_ch)
        for k, v in new_attrs.iteritems():
            sc.__setattr__(k, v)
        return sc
    
    def add_school_old(self, **new_attrs):        
        # school names are 's' + strings of integers
        try:
            max_school_n = max(int(k[1:])
                           for k in self._group_keys()
                           if k.startswith('s'))
            new_k = 's' + str(max_school_n + 1)
        except ValueError:
            new_k = 's0'
        sh_ch = self.shadow.create_group(new_k)
        # making the Pforest object; 
        sc = self.__(sh_ch)
        for k, v in new_attrs.iteritems():
            sc.__setattr__(k, v)
        return sc

    #--------------------------------------------------------------- attributes

    # compatibility
    @property
    def pruned(self):
        return False

    #----------------------------------------------------- children are Ptrees

    @classmethod
    def _child_class(cls):
        return Ptree

    # this is used for parents or children: has to be general
    # 
    def _make_from(self, group):
        if PtHelpers._tree_code(group.name) == u't':
            return self._child_class()(group)
        if group == self.shadow:
                return self
        return self.__class__(group)

    #----------------------------------------------------------------- renaming

    population = PtreeBase.branch_no
    add_tree = PtreeBase.add_child

    # here, there is a subtlety: h5py returns links by default in 
    # alphabetical order, not in creation order. this does not 
    # work when skipping the first trees to improve statistics.
    # since no way of globally redefining iteration order seems to
    # exist, we have to manually sort the tree list.
    @property
    def trees(self):
        return [self._make_from(self.shadow[k]) for k in self.tree_keys()]
    
    def tree_keys(self):
        num = lambda key: int(pp.basename(key)[1:])
        ntgrps = sorted((num(k), k) for k in self._group_keys())
        return (k for dmy_n, k in ntgrps)
        

    #----------------------------------------------------------------- printing

    def str_rich(self):
        d = {'log_int': 'logging interval {0:2g}', 'population': '{0} trees'}
        report_str = []
        for k, v in d.iteritems():
            try:
                report_str.append(v.format(self.__getattribute__(k)))
            except ValueError: pass
        return '<< ' + ','.join(report_str) + ' >>'

    def __str__(self):
        return self.str_rich()



class Ptree(PtreeBase):

    # temporary data to save write operations.
    __slots__ = ()

    # inherit __init__ from base.

    #------------------ extra properties which only make sense when propagating

    # default values for absorbed, pruned, survived implemented in the getters
    # we lump these three into a compound attribute p_a_s to speed up access.

    absorbed = property(*PtHelpers._pas_gs('a'))
    pruned = property(*PtHelpers._pas_gs('p'))
    survived = property(*PtHelpers._pas_gs('s'))


    # this is for the Ptree / Pforest combination
    # The 0-depth level '/' in a file exclusively is a Pforest.
    @property
    def parent(self):
        sh_p = self.shadow.parent
        # forests are always at the file root, or named 'school'
        is_forest = PtHelpers._tree_code(sh_p.name) in 'rs'
        return  (Pforest(sh_p) if is_forest else Ptree(sh_p))

    # children of trees are trees; inherit children property.

    #--------------------------------------------------- structure manipulation

    # branching.
    def branch(self, nfold, rewt=1., point=None, **child_attrs):
        if not self.is_leaf:
            raise BranchingError('Trying to branch a non-leaf')
        if point:
            temp_tip = point  # tip -> branching point (get?set is costly!)
        else:
            temp_tip = self.tip
        if not nfold:
            self.pruned = True
            self._update_height(temp_tip[0])
            self.tip = temp_tip
            return []
        # no extra handling of the single child case: re-weighting
        # may be necessary, so this really is a branching event.
        clogwt = self.logwt + math.log(rewt)
        clist = []
        for dmy in xrange(nfold):
            clist.append(self.add_child(base=temp_tip,
                                        logwt=clogwt, **child_attrs))
        self.tip = temp_tip
        return clist  # the children are made anyway.


    # absorption. need to call when appropriate for the system at hand
    def absorb(self, point=None):
        '''Be absorbed here.

        Part of the system dynamics; not pruning!

        '''
        try:
            if not self.is_leaf:
                raise BranchingError('Trying to absorb a non-tip')
            if point:
                self.tip = point  # tip -> branching point
            self.absorbed = True
            self._update_height(self.tip[0])
        except BranchingError as e:
            log_.error(e)
            raise

    # survival. called after running to the final time
    def cap(self, point=None):
        '''Cap a still growing tip here.

        do this at the end of a simulation

        '''
        try:
            if not self.is_leaf:
                raise BranchingError('Trying to cap a non-tip')
            if point:
                self.tip = point  # tip -> branching point
            self.survived = True
            self._update_height(self.tip[0])
        except BranchingError as e:
            log_.error(e)
            raise

    # graft new children at an existing node.

    def graft(self, nextra=1, rewt=None):
        '''Graft new branches onto an existing tree.

        add nextra new zero-length children to the current branch; 
        all children are reweighted to conserve total wt.

        '''
        try:
            if not self.branch_no:
                raise BranchingError('Grafting only works on fertile branches')
            cbase = self.children[0].base
            owt = sum([math.exp(c.logwt) for c in self.children])
            cwt = (owt / self.branch_no)  # new wt is avg. old wt
            if not rewt:
                rewt = owt / (nextra * cwt + owt)  # conserve children's weight
            for dmy in range(nextra):
                self.add_child(base=cbase, logwt=math.log(cwt))
            self._update_height(np.infty)
            for i in self.children:
                i.logwt += math.log(rewt)
        except BranchingError as e:
            log_.error(e, exc_info=1)
            raise


    #-------------------------------------------------------------- propagation

    # growth, possibly interrupted by absorption.

    def grow_trunk(self, propa, upto_t, log_int=None):
        '''Grow the branch up to time upto_t and log at intervals

        propa is a propagate.SystemIntegrator().
        Logging appends the time and state
        to the Ptree.trunk at multiples of the log_int.

        '''
        if not self.is_leaf:  # silently exit for internal branches
            return
        temp_tip = self.tip[:]
        from_t = temp_tip[0]
        if from_t > upto_t:
            raise GrowthError('trying to grow backwards')
        try:
            if log_int is None:
                log_int = self.log_int
            # for eqtimes:
            # arange cuts the last element off.
            # trunk writing convention: [from_t, to_t)
            # open interval on the right. only relevant if
            # logging times coincide with start or stop times
            eqtimes = log_int * np.arange(np.ceil(from_t / log_int),
                                          np.ceil(upto_t / log_int))
            # alternative:
#            eqtimes = log_int * np.arange(- (-from_t // log_int),
#                                          - (-upto_t // log_int))

        except ValueError as e:
            log_.warning(e, exc_info=1)
            eqtimes = []
        # fine, now go on with the actual propagation.
        # temp data for speedup and for appending.
        temp_trunk = list(self.trunk)
        try:
            for t in eqtimes:  # integrate up to the next logging interval
                temp_tip = propa.int(t, *temp_tip)
                temp_trunk.append(temp_tip[:])
            self.tip = propa.int(upto_t, *temp_tip)
        except Absorption as abs:
            self.tip = temp_tip
            self.absorb(abs.point)  # handles tip and trunk etc.
        except Crossing as cr:
            self.tip = cr.point
            raise  # for branching outside
        finally:
            self.trunk = temp_trunk


    #----------------------------------------------------------------- printing

    def str_rich(self, elems=None, short=False):
        """give a customized representation of the node.

        the elements can be combined from
        ('base_t','tip_t', 'base', 'tip', 'status')
        the short switch abbreviates.
        """
        if not elems:
            elems = ('base_t', 'tip_t', 'status')
        # only one of these options should be true:
        try:
            status = {self.is_leaf:'growing',
                          self.absorbed:'absorbed',
                     bool(self.branch_no):str(self.branch_no) + ' children',
                          self.pruned:'pruned',
                          self.survived:'survived'
                          }[True]
        except KeyError:  # pragma: no cover
            raise PrintingError('undefined status')
        d = {'base_t': ['base_t', '[0]:2g', self.base],
             'base': ['base', '', self.base],
             'tip_t': ['tip_t', '[0]:2g', self.tip],
             'tip': ['tip', '', self.tip],
             'status': ['', '', status] }
        if short:  # the abbreviations
            for v in d.itervalues():
                v[0] = v[0][:1]
                v[-1] = v[-1][:2]
        substitutions = [d[el][-1] for el in elems]
        format_str = ('<' + ', '.join(d[el][0] + ' ' + 
                                      '{' + str(i) + d[el][1] + '}'
                                for i, el in enumerate(elems))
                       + '>')
        return format_str.format(*substitutions)
        # return '<base {0[0]:3g}, {2}>'.format(self.base, self.tip, status)

    def __str__(self):
        return self.str_rich(('base_t', 'tip_t', 'status'), short=False)

    # override base class: extra kwargs
    def dump(self, elems=None, short=False, outf=sys.stdout):
        """Dump a formatted representation of this tree to the specified file
            descriptor.

            :outf Output file descriptor.
        """
        for i in self:
            print >> outf, "  "*(i.depth), i.str_rich(elems, short)



#------------------------------------------------------------------ error types

class BranchingError(Exception):
    """An error making a branching attempt fail"""
    pass

class PruningError(Exception):
    """An error when printing a tree object"""
    pass

class GrowthError(Exception):
    """An error that occurs during growth"""
    pass

class PrintingError(Exception):
    """Error during printing"""
    pass


#------------------------------------------------------------------- essentials

def laba_daba():
    return 'doo-wap!'

