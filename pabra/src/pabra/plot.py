'''
Created on Apr 27, 2010

@author: nbecker
'''
import matplotlib.pyplot as plt
import numpy as np
# logging; overall level and format controlled in __init__.py
import logging
log_ = logging.getLogger(name=__name__)
import matplotlib
import functools as fut

#plt.rcParams['backend'] = 'Qt4Agg'

def time_space(ptr, coord=0):
    """Extract time and space from a tree trunk in a form 
    suitable for plotting with plot. Returns a shape (2,n) array
    
    """
    extract_ = ((lambda t, x: np.array([t, coord(t, x)])) if callable(coord)
                # this is needed for the case that x is just a number: 
                else lambda t, x: np.array([t, np.array(x, ndmin=1)[coord]]))
    if not ptr.trunk and ptr.base[0]==ptr.tip[0]:
        raise EmptyTrunk(ptr)
    tru = ptr.btrunkt
    return np.array([extract_(*pt) for pt in tru]).T

def space_space(ptr, coords=[0, 1]):
    """Extract two space coords from a tree trunk in a form 
    suitable for plotting with plot. Returns a shape (2,n) array.
    If no coords are given, take the first two
    
    """
    extract_ = (coords if callable(coords)
                # x must be an array: 
                else lambda t, x: x[list(coords)])
    if not ptr.trunk and ptr.base[0]==ptr.tip[0]:
        raise EmptyTrunk(ptr)
    tru = ptr.btrunkt
    return np.array([extract_(*pt) for pt in tru]).T

alpha_fn = (lambda wt:np.clip(.2 + wt * .6, 0, 1))


def _plot_tree(ptree, to_plot, attr_string=None,
               axes_obj=None, color='black', wt_alpha=None, max_size=np.inf,
               order='pre', **plot_args):
    """Plot the whole ptree into the axes_obj.
    
    what coords to plot is specified by to_plot.
    """
    if wt_alpha is None:
        al_fn = lambda wt: 1
    else: 
        al_fn = wt_alpha if callable(wt_alpha) else alpha_fn
    if not axes_obj:  # not used since bad: isinstance...
        new_fig = plt.figure()
        axes_obj = new_fig.add_subplot(1, 1, 1)
    iter = ptree.preorder() if order.startswith('pre') else ptree.postorder()
    for i, sub in enumerate(iter):
        if i <= max_size:
            try:
                x, y = to_plot(sub)
                lines = axes_obj.plot(x, y, color=color, alpha=al_fn(sub.wt),
                                      **plot_args)[0]
                if attr_string:
                    sub.__setattr__(attr_string, lines)
            except EmptyTrunk as which_tree:
                log_.debug('not plotting: ' + repr(which_tree))
        else:
            break
    axes_obj.figure.canvas.draw()
    return axes_obj


def plot_time_space(ptree, coord=0, **kwargs):
    """Do a space-time plot of the whole ptree into the axes_obj.
    
    """
    to_plot = fut.partial(time_space, coord=coord)
    # a full functools version is possible, but then ptree would become
    # keyword-only. therefore, we do not do it.
    return _plot_tree(ptree, to_plot, **kwargs)

def plot_space_space(ptree, coords=[0, 1], **kwargs):
    """Do a trajectory plot of the whole ptree into the axes_obj.
    
    """
    to_plot = fut.partial(space_space, coords=coords)
    return _plot_tree(ptree, to_plot, **kwargs)

# branching points plots

def bpts_time_space(ptr, coord=0):
    extract_ = (coord if callable(coord)
                # this is needed for the case that x is just a number: 
                else lambda t, x: np.array(x, ndmin=1)[coord])
    if not ptr.alive:
        c_val = extract_(*ptr.tip)
        # also put the number of children
        return np.array((ptr.tip[0], c_val, ptr.branch_no))

def bpts_space_space(ptr, coords=[0, 1]):
    extract_ = (coords if callable(coords)
                # x must be an array: 
                else lambda t, x: np.array(x, ndmin=1)[coords])
    if not ptr.alive:
        c_val = extract_(*ptr.tip)            
        # also the number of children
        l = list(c_val)
        l.append(ptr.branch_no)
        return np.array(l)

# not used at the moment:
def n_color(n):
    basecol = np.array([.7, .7, 1])
    cd = dict(zip(range(6), (.2 * basecol * i for i in range(6))))
    try:
        c = cd[int(n)]
    except KeyError:
        c = basecol
    return c

_custom_clist = matplotlib.colors.ColorConverter().to_rgba_array(
            ['none', 'black', 'orange', 'blue', 'darkgreen', 'green'])

# this does not support alpha (!)
#custom_cmap = matplotlib.colors.ListedColormap(_custom_clist, name='custom')

def custom_cfun(n):
    try:
        return _custom_clist[n]
    except IndexError:
        return (0.,0,0,.5)

# set the size of branching and pruning points: 1 child case is smaller.
def n_size(n):
    try:
        return {1:10}[n]
    except KeyError:
        return 10

# set the edge style depending on progeny

def n_edge(n):
    return 'black' if n == 0 else 'none'


def _plot_bp(ptr, to_plot, axes_obj=None, size_fn=n_size, 
             vmax=4, edge_fn=n_edge, max_size=np.inf, order='pre', **kwargs):
    """Plot the branch points for visualization"""
    if not axes_obj:
        new_fig = plt.figure()
        axes_obj = new_fig.add_subplot(1, 1, 1)
    iter = ptr.preorder() if order.startswith('pre') else ptr.postorder()
    try:
        acc = [to_plot(sub) for i, sub in enumerate(iter) if i <= max_size]
        x, y, n = np.array(acc).T
    except TypeError:  # happens when to_plot gives empty list
        log_.info('No branching points to plot for {0}'.format(ptr))
        return
    return axes_obj.scatter(x, y, s=map(size_fn, n), 
                            c=map(custom_cfun, n), 
                            vmin=0, vmax=vmax,
                            marker='o',
                            edgecolors=map(edge_fn, n),
                            **kwargs)

def plot_bp_time_space(ptree, coord=0, **kwargs):
    """Do a space-time plot of the ptree branching points into the axes_obj.
    
    """
    to_plot = fut.partial(bpts_time_space, coord=coord)
    return _plot_bp(ptree, to_plot, **kwargs)

def plot_bp_space_space(ptree, coords=[0, 1], **kwargs):
    """Do a space-space plot of the ptree branching points into the axes_obj.
    
    """
    to_plot = fut.partial(bpts_space_space, coords=coords)
    return _plot_bp(ptree, to_plot, **kwargs)


#---------------------------------- adding this stuff as methods to matplotlib:


def _make_method(plot_fun):
    def method(self, *args, **kwargs):
        return plot_fun(*args, axes_obj=self, **kwargs)
    return method

matplotlib.axes.Axes.plot_time_space = _make_method(plot_time_space)
matplotlib.axes.Axes.plot_space_space = _make_method(plot_space_space)
matplotlib.axes.Axes.plot_bp_time_space = _make_method(plot_bp_time_space)
matplotlib.axes.Axes.plot_bp_space_space = _make_method(plot_bp_space_space)


#def method_plot_time_space(self, ptree, **kwargs):
#    return plot_time_space(ptree, axes_obj=self, **kwargs)
#
#matplotlib.axes.Axes.plot_time_space = method_plot_time_space

#------------------------------------------------------------------------- done



class EmptyTrunk(Exception):
        pass

