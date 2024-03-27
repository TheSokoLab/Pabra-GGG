'''
Created on May 6, 2010

@author: nbecker

Collect tests for the plot facility
'''

from pabra.propagate import *  #@UnusedWildImport
from pabra.plot import *  #@UnusedWildImport
from unittests.fixtures import H5TestCase


class TestPlotting(H5TestCase):
    """test the structural operations on the p5 class"""    
    
    def setUp(self):
        self.setuptrees()
        self.setupito()
        self.onegen.log_int = .05
        # see if random growing works ok
        np.random.seed(412)
        self.onegen.grow_tips(self.ito_no_f, 2.)
        np.random.seed(23)
        self.onegen2d.log_int = .05
        self.onegen2d.grow_tips(self.ito_harmonic_2d, 10.) 
                     
    def tearDown(self):
        self.teardowntrees()
        
    def test_timespace(self):
        ax = plot_time_space(self.onegen)
        print ax.figure
        self.onegen.cap_tips()
        #print ax.lines[0].get_data()
        self.assertAlmostEqual(ax.lines[0].get_data()[0][-1], 
                               self.onegen.height, 2)
        ax2 = plot_time_space(self.onegen2d, coord=1, color='b', wt_alpha=True)
        print ax2
        self.assertAlmostEqual(ax2.lines[0].get_data()[1][0], .00115, 2)
        #print ax2.lines[0].get_data()
        plot_time_space(self.onegen2d.children[0], 
                        coord=1, axes_obj=ax, wt_alpha=1)
        ax.lines[-1].set(lw=2., marker='o')
        for tr in self.onegen2d:
            print tr.wt


    def test_spacespace(self):
        axs = plot_space_space(self.onegen2d)
        print axs.figure
        #print axs.lines[0].get_data()
        # the last point of the drawn line is the tip:
        self.assertAlmostEqual(axs.lines[0].get_data()[0][-1],
                            self.onegen2d.children[1].tip[1][0], 2) 
        # the first point of the drawn line is the base:
        self.assertAlmostEqual(axs.lines[0].get_data()[0][0],
                            self.onegen2d.children[1].base[1][0], 2) 
        #this is just to cover the case with the wt alpha 
        dmy_axs2 = plot_space_space(self.onegen2d, coords=(0,1), 
                                    color='b', wt_alpha=True)
#        plt.show()
        pass
