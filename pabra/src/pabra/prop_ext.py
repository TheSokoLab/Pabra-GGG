'''
Created on Nov 10, 2011

@author: nbecker

Here we make a wrapper class around an external propagator module.
It implements either the SystemIntegrator interface or the 
CrossingIntegrator interface (the latter if interface crossing detection is
built in.
'''
from pabra.propagate import CrossingIntegrator

import subprocess
import os.path as pp
import os
import re
import pabra.ptree_h5 as p5
import copy as cp
import shutil
from h5py.h5g import get_objinfo
import numpy as np

# logging; overall level and format controlled in __init__.py
import logging
log_ = logging.getLogger(name=__name__)

#==============================================================================
# Utilities
#==============================================================================

class CallError(Exception):
    """unsuccessful calls to external programs"""
    pass


def oscall(command) :
    try:
        retcode = subprocess.call(command, shell=True)
        if retcode <> 0:
            raise CallError("Call returned signal {0}".format(-retcode))
    except OSError as e:
        raise CallError("Call failed with: {0}".format(e))

def oscall_with_output(command, args):
    try: 
        # Python 2.7 only
        # output = subprocess.check_output(command+" "+args, shell=True) 
        output = subprocess.Popen([command, args],
                                  stdout=subprocess.PIPE).communicate()[0]
        return output
    except OSError as e:
        raise CallError("Call failed with: {0}".format(e))
    
#----------------------------------------- module constants: filenames and such

# no parameterID in these files, in order to handle them from PtreeExt
_infile_name_ = "inputAI"
_outfile_name_ = "inputAI.end"
_tempfile_name_ = 'temp'
_exe_name_ = 'gillespie.X'
_runfile_name_ = '_runfile'

#==============================================================================
# Ptree structure adaptation to external file system objects
#==============================================================================

#In order to link a Ptree with files for each node we need to set some extra 
#attributes. also we need to call the propagator (below) with these attributes.

class _PtreeExtras(object):
    '''a mixin class with extra machinery to maintain a tree on disk.
    at the moment, the directory structure is flat'''
    @staticmethod
    def _mkndir(ndir):
        try:
            os.makedirs(ndir)
        except OSError as e:
            if not e.errno == 17: raise  # 'dir already exists'

    @property
    def node_dir(self):
        return self.shadow.attrs['node_dir']
    @node_dir.setter
    def node_dir(self, dirname):
        self._mkndir(pp.join(self.root_dir, dirname))
        self.shadow.attrs['node_dir'] = dirname

    @property
    def root_dir(self):
        # note that the root dir of the PforestExt is independent of the 
        # dir in which the HDF5 file resides.
        try: return self.shadow['/'].attrs['node_dir']
        except KeyError: return os.getcwd()


class _PforestExt(_PtreeExtras, p5.Pforest):
    def __init__(self, shadow_group=None, 
                 node_dir=None, root_dir=None, **attrs):
        p5.Pforest.__init__(self, shadow_group, **attrs)
        # PForestExt node dirs are absolute.
        if not root_dir is None:
            self.shadow['/'].attrs['node_dir'] = pp.abspath(root_dir)
        if not node_dir is None:
            self.node_dir = pp.relpath(node_dir, self.root_dir)
        else:
            try: self.node_dir
            except KeyError: self.node_dir = './' #  ie. = root dir

    def _child_class(self):
        return p5.PtreeExt  #@UndefinedVariable

    def add_school(self, **new_attrs):
        # this returns PforestExt (by child class mechanism)
        sc = p5.Pforest.add_school(self, **new_attrs)
        # this sets and makes the correct node dir
        sc.node_dir = pp.join(self.node_dir, 's' + 
                        'x'.join(map(str, get_objinfo(sc.shadow.id).objno)))
        # copy input file
        os.symlink(pp.join(self.root_dir, self.node_dir, _outfile_name_),
                pp.join(self.root_dir, sc.node_dir, _outfile_name_))
        return sc
    
    def add_child(self, **attrs):
        # child class instance, initialized!
        ch = p5.PtreeBase.add_child(self, **attrs)
        # copy state data
        os.symlink(pp.join(self.root_dir, self.node_dir, _outfile_name_),
                pp.join(self.root_dir, ch.node_dir, _infile_name_))
        return ch
    
    def add_tree(self, **attrs):
        return self.add_child(**attrs)


class _PtreeExt(_PtreeExtras, p5.Ptree):
    def __init__(self, shadow_group):
        p5.Ptree.__init__(self, shadow_group)
        # determine the directory of the parent school
        schoolname = self.shadow.name[:self.shadow.name.find('/',1)]
        # make a directory and link it with the shadow node.
        self.node_dir = pp.join(self.shadow[schoolname].attrs['node_dir'],
                        'x'.join(map(str, get_objinfo(self.shadow.id).objno)))

    def add_child(self, **attrs):
        # child class instance, initialized!
        ch = p5.PtreeBase.add_child(self, **attrs)
        # copy state data
        os.link(pp.join(self.root_dir, self.node_dir, _outfile_name_),
                pp.join(self.root_dir, ch.node_dir, _infile_name_))
        return ch

    @property
    def parent(self):
        sh_p = self.shadow.parent
        # forests are always at the file root, or named 'school'
        is_forest = p5.PtHelpers._tree_code(sh_p.name) in 'rs'
        return  (p5.PforestExt(sh_p)#@UndefinedVariable
                 if is_forest 
                 else p5.PtreeExt(sh_p))  #@UndefinedVariable

    def grow_trunk(self, propa, upto_t, log_int=None):
        propa.current_dir = pp.join(self.root_dir, self.node_dir)
        return p5.Ptree.grow_trunk(self, propa, upto_t, log_int)

    def _cleanup(self):
        err = []
        for n in (_outfile_name_, _infile_name_,
                  _runfile_name_, _tempfile_name_):
            try:
                os.unlink(pp.join(self.root_dir, self.node_dir, n))
            except OSError as e:
                err.append("{0}".format(e).split(':')[-1].strip())
        try:
            os.removedirs(pp.join(self.root_dir, self.node_dir))
        except OSError as e:
            err.append("{0}".format(e).split(':')[-1].strip())
        if err:
            log_.warning('runfile cleanup failed at time {0}. files:\n'
                         .format(self.tip[0]) + 
                         '\n'.join(err))

    def branch(self, nfold, rewt=1., point=None, **child_attrs):
        clist = p5.Ptree.branch(self, nfold, rewt, point, **child_attrs)
        self._cleanup()
        return clist

    def absorb(self, point=None):
        p5.Ptree.absorb(self, point)
        self._cleanup()

    def cap(self, point=None):
        p5.Ptree.cap(self, point)
        self._cleanup()



p5.PforestExt = _PforestExt
p5.PtreeExt = _PtreeExt

#==============================================================================
# System propagation using external calls
#==============================================================================

    
#TODO here: what happens to the state x? there is a lot of 
#state handling in the ptree and in propagate...
#in principle, nobody needs to know what x was inside the python code.
#we could present a dummy or boiled down state x to the pabra system, from the 
#external propagator. that would allow everything in ptree to remain the same.
#the dead simplest variant is actually to report the progress coord itself 
#as x...


class TIntegrator(CrossingIntegrator):
    """A SystemIntegrator which runs Tomeks Gillespie code through OS calls
    initargs:
    :runfile_name: self explanatory
    :working_dir: defaults to current wd
    :ifhist: if given, no other kwargs are used and the interfaces are read
    from the given IFHistogram instance
    :other: if ifhist is not given, kwargs of CrossingIntegrator."""
    #a class constant
    _cur = re.compile('CURDIR')
    # init
    # problem here: the flux direction spec is quite complicated; 
    # would be nice to initialize from a given ifhist like in crosswrap.
    def __init__(self, runfile_name, working_dir=None, ifhist=None, 
                 **cikwargs):
        if not ifhist is None:
            self.init_from_ifhist(ifhist)
            self._cross_int = None  # this is in the external input file.
        else:
            CrossingIntegrator.__init__(self, cross_int=None, **cikwargs)
        # now the external stuff.
        if not working_dir is None:
            self.wd = pp.abspath(working_dir)
        else:
            self.wd = os.getcwd()
        # set up the file names
        self.runfile_name = _runfile_name_
        self.template_runfile = pp.join(self.wd, runfile_name)
        with open(pp.join(self.wd, _outfile_name_), 'r') as parfile:
            self._cross_int = float(parfile.readlines()[4])
        # mutable state:
        self.current_dir = self.wd  # this should be a full path
        self.current_seed = 123
    
    # readonly: set in the parameter file of the external propagator.
    @property
    def cross_int(self):
        return self._cross_int
    @cross_int.setter
    def cross_int(self, val):
        if not val is None:
            log_.warning("trying to set cross_int in external integrator;"
                         " not set.")

    def _init_segment(self):
        # the state file must be already copied to the new subdir.
        #
        # two options. either this is a new branch, then no output for 
        # this branch exists.
        cur_dir = pp.abspath(self.current_dir)
        seg_out = pp.join(cur_dir, _outfile_name_)
        seg_in = pp.join(cur_dir, _infile_name_)
        if pp.isfile(seg_out):
            shutil.move(seg_out, seg_in)
        else:  # need to do actual init
            current_runfile = pp.join(cur_dir, self.runfile_name)
            new_tempfile = pp.join(cur_dir, _tempfile_name_)
            self.current_seed += 1
            self.current_seed = self.current_seed % 1000000
            with open(new_tempfile, 'w') as ntf:
                ntf.write('0\n{0}\n'.format(self.current_seed))
            with open(current_runfile, 'w') as nrf:
                with open(self.template_runfile, 'r') as trf:
                    nrf.write(re.sub(self._cur, 
                                     pp.relpath(cur_dir, self.wd), trf.read()))

    def _sample_segment(self):
        try:
            assert os.getcwd() == self.wd
        except AssertionError:
            print "C propagator must be run from its own directory!"
            raise
        start_command = pp.join(self.wd, _exe_name_)
        start_arg = pp.join(pp.abspath(self.current_dir), self.runfile_name)
        # Start it and capture the output
        L = oscall_with_output(start_command, start_arg).split()  # list of str
        # the following returns a 0d array if only one coord is returned
#        li = L.__iter__()
#        t = float(li.next())
#        x = li.next()  # still str
#        try:
#            x = [x, li.next()]
#            while 1:
#                x.append(li.next())
#        except StopIteration:
#            return t, np.array(x, dtype=np.float)  # float conversion.
        return  float(L[0]), np.array(L[1:], dtype=np.float)

    # this only knows how to propagate by cross_int; 
    # no time difference argument. this means that cross_int is the finest 
    # time scale division.
    # xcur is contained in the state file...
    def _prop_or_abs(self, tcur):
        self._init_segment()
        dt, newx = self._sample_segment()
        return dt + tcur, newx


    def _prop_abs_cross(self, tfinal, tinit, xinit):
        # here, propagation proceeds always in cross_int steps.
        cross_int = self.cross_int
        tcur, xcur = tinit, cp.deepcopy(xinit)  # xinit is mutable.
        icur = self._get_inds([c_f(tcur, xcur) for c_f in self.c_funs])
        while tcur < tfinal - .5 * cross_int:  # to-last
            tcur, xcur = self._prop_or_abs(tcur)
            icur = self._cross((tcur, xcur), icur)  #  raise?
            #log_.warning("lambda = " + str(xcur)+" , icur = 
            # "+str(icur)+", tcur = "+str(tcur))
        # last propagation to roughly the final t and x.
        # tcur, xcur = self._prop_or_abs(tcur)
        # final check to prevent overspill
        if tcur < tfinal:
            self._cross((tcur, xcur), icur)
        return tcur, xcur



class ExternalCrossingPropagator(CrossingIntegrator):
    pass

