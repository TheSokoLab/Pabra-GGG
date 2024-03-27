# -*- coding: utf-8 -*-
'''
Created on Feb 29, 2012

@author: sokolowski
'''

import sys
import os
import numpy
import pabra.ptree_h5 as p5


class plot():

    name = "bogus"
    command = "plot"
    rows = "1:2"
    xlabel = "x"
    ylabel = "y"
    xrange = "[0:*]"
    yrange = "[0:*]"
    xtics  = ""
    ytics  = ""
    xformat = ""
    yformat = ""
    definitions = ""


##################
# PATH EXTRACTOR #
##################
class PathExtractor(object):

    def __init__(self, filename, teh, t_total, dt_log, dt_min):
        # Make the forest accessible        
        self.pf=p5.Pforest(p5.h5.File(filename,'r'))

        # The first layer of branches is the schools
        self.schools = self.pf.children

        # Define how many schools should be extracted
        self.total_extract_height = float(teh)
        self.min_generations = 10

        # The logging time of the FFS sim.
        self.t_total = float(t_total)   # the maximal simulation time of one trajectory
        self.dt_log  = float(dt_log)    # the logging interval of the FFS sim.
        self.dt_min  = float(dt_min)    # the minimal update interval of the FFS sim.; typically 60 s

        # The row of the RC array
        # at which the concentration data starts      
        self.conc_start = 12
        # and the no. of nuclei recorded
        self.no_z_nucl = 32

        # Set some default values
        self.paths = []
        self.paths_longest = []
        self.outputdir = './paths'

        # The height of the longest (filtered) school in the tree
        self.longest_filtered_height = 0
        self.longest_filtered_no     = -1
        
        # Parameters of the snapshot acquisition
        self.basejump = 3
        self.t_tol = self.dt_min  # time tolerance for snapshots; the loop will make a snapshot at t' if t'>t-t_tol/2 and t'<t+t_tol/2
        self.t_snap = 1800.0      # time between two snapshots; make sure this divides t_total
        self.N_snap = int(self.t_total/self.t_snap) + 1 # no. of snapshots (including first and last)                                
        self.t_extra = self.dt_log # half of the time interval within which the scheme will
                                   # make additional snapshots around the main snapshot defined
                                   # Set this to zero if you want only the N_snap snapshots
                                   # defined above

        # The dictionary containing the colors for each path
        self.colors = {}

        print("# Created "+self.__class__.__name__+" with filename = "+filename)
        print "# t_snap=%s, N_snap=%s, t_tol=%s,  t_extra=%s" % (self.t_snap, self.N_snap, self.t_tol, self.t_extra)
        

    def extract(self):
        """ Constructs an enumerated list of all paths, i.e. a list of 
            tuples (path_no, path).
        """
        print("# Extracting paths...")
        print("# total_extract_height = " + str(self.total_extract_height))
        
        # First filter the schools that have more than 1 generation
        # to get rid of immediately pruned schools
        schools_filt = self.filter_schools()

        print("# Building paths list...")       

        for s_no, s in schools_filt :

            new_paths = [(s_no, node.btrunkt) for node in s if not node.base[0] == node.tip[0] ] # discard nodes without trunk
            self.paths.extend(new_paths)
            if s_no == self.longest_filtered_no :
                
                leaves = [node for node in s if not node.children and node.tip[0]>=0.99*self.longest_filtered_height] # only take long sequences
                # Sort the leaves by reaction coordinate values (symmetric RC); tip[1] is the RC vector
                leaves_sorted = sorted(leaves, key=lambda n: (max(n.tip[1][0], n.tip[1][2])+max(n.tip[1][1], n.tip[1][3])) / n.tip[1][4] )
                # Construct the path list for the branch with the highest final RC value (destruction path)
                new_paths = [(s_no, node.btrunkt) for node in leaves_sorted[-1].path_from_root() if not node.base[0] == node.tip[0] ]
                self.paths_longest.extend(new_paths)

        self.paths         = list(enumerate(self.paths))
        self.paths_longest = list(enumerate(self.paths_longest)) # to ensure equal format, even though this should have only one element

        ll = len(self.paths)
        lll = len(self.paths_longest)
        print("#   Extracted "+str(ll)+" paths.")
        print("#   Extracted "+str(lll)+" paths from school with largest height.")

    def filter_schools(self):

        print("# Filtering schools ...")
        schools_filt = []
        ll = len(self.schools)

        i=0
        h=0.0
        while i < ll and h < self.total_extract_height :
            i = i+1
            s = self.schools[-i]
            if(s.generations >= self.min_generations):
                schools_filt.append( (i, s) )
                h = h + s.height
                if(s.height > self.longest_filtered_height):
                    self.longest_filtered_height = s.height
                    self.longest_filtered_no = i
                print("#   Appending a school with school.generations = "+str(s.generations)+" and school.height = "+str(s.height))

        print("#   Longest filtered school no. = "+str(self.longest_filtered_no)+" with height = "+str(self.longest_filtered_height))

        return schools_filt
        
    def write(self):
        """ Creates a folder containing all paths and write them
            into separate files into this folder.
        """
        assert(self.paths != [])
        assert(self.paths_longest != [])

        # Create the folder
        if not os.path.exists(self.outputdir):
            print("# Creating directory " + self.outputdir)
            os.mkdir(self.outputdir)
        else:
            print("# Output directory exists, omitting creation...")

        # OUTPUT ITERATION 1
        # Now go through the path list and output each path
        # into an own file in the created folder
        for pid, path_with_numbers in self.paths_longest: # in self.paths or self.paths_longest
            
          filename = self.pathfilename(pid)
          print "# Writing " + filename
          fout = open(filename, "w")

          # Get school number and btrunkt
          s_no    = path_with_numbers[0]
          btrunkt = path_with_numbers[1]

          # Set path colour based on school number
          self.colors[str(pid)] = s_no
    
          # Go through the trunk and append to output for this path
          for t, xvector in btrunkt :
              if t > 0 : # discard initial values
                  output = str(t)
                  i=0
                  for x in xvector:
                      if i<self.conc_start :
                          output = output + " " + str(x)
                      i+=1

                  output = output + " " + str(self.colors[str(pid)]) + "\n"
                  fout.write(output)

          fout.close()

        # OUTPUT ITERATION 2 - "SNAPSHOTS"
        # As a second step, output the concentrations as a function of
        # the axis coordinate z
        for pid, path_with_numbers in self.paths_longest: # in self.paths or self.paths_longest

          # Get school number and btrunkt
          s_no    = path_with_numbers[0]
          btrunkt = path_with_numbers[1]
    
          # Go through the trunk and append to output for this path
          for t, xvector in btrunkt :
 
              lfh = self.longest_filtered_height

              def have_snapshot_at(t):
            
                  #lfh = self.longest_filtered_height
                  have = False                  

                  for t_off in [-self.t_extra, 0, +self.t_extra]:

                      # first and last frame
                      if (t + t_off > self.t_extra - self.t_tol/2 and t + t_off < self.t_extra + self.t_tol/2) or \
                         (t + t_off > (self.t_total-self.t_extra) - self.t_tol/2 and t + t_off < (self.t_total-self.t_extra) + self.t_tol/2) :
                              
                              have = True

                      # The very first snapshot has always to be treated separately because for t=0 the
                      # lambda value is yet undefined => take t=dt_min instead
                      if (t==0.0):
                              have = False

                      if (t > self.dt_min - self.t_tol/2 and t < self.dt_min + self.t_tol/2):
                              have = True

                      # in-between frames
                      for n in range(1, self.N_snap-1):                                            
                      
                          if (t + t_off > 1.0*n*self.t_snap - self.t_tol/2) and (t + t_off <= 1.0*n*self.t_snap + self.t_tol/2) :
                          
                              have = True

                  return have


              if have_snapshot_at(t):

                  # Tag output filename with reaction coordinate value and time
                  # l   = xvector[5] # asymmetric lambda
                  l = (max(xvector[0],xvector[2])+max(xvector[1],xvector[3])) / xvector[4] # symmetric lambda

                  # Construct filename and open the file
                  filename = self.outputdir + "/snapshot_s" + str(s_no) + "_t" + str(int(round(t,-1))).zfill(6) + "_l" + str(round(l,2))
                  print "# Writing " + filename
                  fout = open(filename, "w")

                  i=0
                  for x in xvector:

                      if i>=self.conc_start :
                      # Only consider the relevant data in the RC array
                      # Concentration data starts at index self.conc_start

                          output = str(t) + " " + str(i-self.conc_start)

                          # Extract the separate gene copy numbers from
                          # the long integer word
                          value = numpy.zeros([5])
                          xrest = x
                          for j in [3, 2, 1, 0] :
                              base = 10**(j * self.basejump)
                              value[j] = int(xrest / base)
                              xrest -= value[j]*base

                          # Append them to the output string
                          # in reverse order, i.e. starting from Hb
                          for j in [0, 1, 2, 3] :
                              output += " " + str(int(value[j]))

                          # As a control add the whole string and the school no.
                          output += " " + str(int(x)) + " " + str(s_no) + " " + str(l) + "\n"
                          fout.write(output)

                      i+=1

                  fout.close()

    def create_gnuplot_file(self, plot):
        """ Creates a rudimentary gnuplot script for path files
        """
        assert(self.paths != [])       

        # Determine filename
        filename = self.outputdir + '/plot_' + plot.name + '.gnu'
        print "# Writing " + filename
        fout = open(filename, "w")          

        # Write a preamble with some settings
        fout.write('set autoscale\n')
        fout.write('unset key\n')
        fout.write('set xrange '+plot.xrange+'\n')
        fout.write('set yrange '+plot.yrange+'\n')
        fout.write('set xlabel \''+plot.xlabel+'\'\n')
        fout.write('set ylabel \''+plot.ylabel+'\'\n')
        fout.write('set xtics '+plot.xtics+'\n')
        fout.write('set ytics '+plot.ytics+'\n')
        fout.write('set format x '+plot.xformat+'\n')
        fout.write('set format y '+plot.yformat+'\n')
        fout.write('set terminal postscript color enhanced solid \"Arial\" 20\n')
        fout.write('set output \''+plot.name+'.ps\'\n')
        fout.write(plot.definitions+'\n')

        # Filtering of input files
        #filtercommand='< tail --lines=+3 '
        filtercommand=''

        # Construct the plot command
        plots = ''
        for pid, path in self.paths:
          
          # Name of the file to plot
          plotname = self.pathfilename(pid)
          # The plotstyle
          color = self.colors[str(pid)]
          style = 'u '+plot.rows+' w l lw 3 lc '+str(color)

          if plots:
              plots = plots + ',\\\n' + ' \''+filtercommand+plotname+'\' ' + style
          else:
              plots = ' \''+filtercommand+plotname+'\' ' + style
        
        # Write the constructed plot command
        plotcommand = plot.command + ' ' + plots + '\n'
        fout.write(plotcommand)
        fout.close()

    def pathfilename(self, pid):

      return self.outputdir + '/path' + str(pid).zfill(6)


#####################
# SEGMENT EXTRACTOR #
#####################
class SegmentExtractor(object):

    def __init__(self, filename, teh):
        # Make the forest accessible
        self.pf=p5.Pforest(p5.h5.File(filename,'r'))

        # The first layer of branches is the schools
        self.schools = self.pf.children

        # Define how many schools should be extracted
        self.total_extract_height = float(teh)

        # Set some default values
        self.segments = []
        self.max_wt = 1.0
        self.outputdir = './segments'

        print("# Created "+self.__class__.__name__+" with filename = "+filename)

    def extract(self):
        """ Constructs a list of all segments
        """
        print("# Extracting segments...")
        print("# total_extract_height = " + str(self.total_extract_height))

        # First filter the schools that have more than 1 generation
        # get rid of immediately pruned schools
        schools_filt = self.filter_schools()

        for s in schools_filt:

            self.segments.extend([tr for tr in s if not tr.base[0] == 0.0])
            # discarding segments with time = 0 initially

        ll = len(self.segments)
        print("#   Extracted "+str(ll)+" segments.")

    def filter_schools(self):

        print("# Filtering schools ...")
        schools_filt = []
        ll = len(self.schools)

        i=0
        h=0.0
        while i < ll and h <= self.total_extract_height :
            i = i+1
            s = self.schools[-i]
            if(s.generations > 1):
                schools_filt.append(s)
                h = h + s.height
                print("#   Appending a school with school.generations = "+str(s.generations)+" and school.height = "+str(s.height))

        return schools_filt

    def write(self):
        """ Creates a folder containing all paths and write them
            into separate files into this folder.
        """
        assert(self.segments != [])

        # Create the folder
        if not os.path.exists(self.outputdir):
            print("# Creating directory " + self.outputdir)
            os.mkdir(self.outputdir)
        else:
            print("# Output directory exists, omitting creation...")

        for sid, seg in list(enumerate(self.segments)):
            
          filename = self.segmentfilename(sid)
          print "# Writing " + filename
          fout = open(filename, "w")

          t, xvector = seg.base
          output = str(t)
          for x in xvector:
              output = output + " " + str(x)
          output = output + " " + str(seg.wt) + "\n"
          fout.write(output)

          t, xvector = seg.tip
          output = str(t)
          for x in xvector:
              output = output + " " + str(x)
          output = output + " " + str(seg.wt) + "\n"
          fout.write(output)

          fout.close()

    def create_gnuplot_file(self, plot):
        """ Creates a rudimentary gnuplot script for segment files
        """
        assert(self.segments != [])

        # Determine maximal weight for color adjustment
        max_wt = self.get_max_weight(self.segments)

        # Determine filename
        filename = self.outputdir + '/plot_' + plot.name + '.gnu'
        print "# Writing " + filename
        fout = open(filename, "w")          

        # Write a preamble with some settings
        fout.write('set autoscale\n')
        fout.write('unset key\n')
        fout.write('set palette defined (0 "grey", 1 "blue")\n')
        fout.write('set cbrange [0:1]\n')
        fout.write('set xrange '+plot.xrange+'\n')
        fout.write('set yrange '+plot.yrange+'\n')
        fout.write('set xlabel \''+plot.xlabel+'\'\n')
        fout.write('set ylabel \''+plot.ylabel+'\'\n')
        fout.write('set xtics '+plot.xtics+'\n')
        fout.write('set ytics '+plot.ytics+'\n')
        fout.write('set format x '+plot.xformat+'\n')
        fout.write('set format y '+plot.yformat+'\n')
        fout.write('set terminal postscript color enhanced solid \"Arial\" 20\n')
        fout.write('set output \''+plot.name+'.ps\'\n')
        fout.write(plot.definitions+'\n')

        # Filtering of input files
        filtercommand=''

        # Construct the plot command
        plots = ''
        for sid, segment in list(enumerate(self.segments)):      
          
          # Name of the file to plot
          plotname = self.segmentfilename(sid)          
          # Get the weight of the path
          wt = self.wt_mapping(segment.wt / max_wt)
          assert(wt >= 0.0 and wt <= 1.0)
          # The plotstyle
          style = 'u '+plot.rows+' w l lw 3 lc palette frac '+str(wt)+' '      

          if plots:
              plots = plots + ',\\\n' + ' \''+filtercommand+plotname+'\' ' + style
          else:
              plots = ' \''+filtercommand+plotname+'\' ' + style

          
        plotcommand = plot.command + ' ' + plots + '\n'
        fout.write(plotcommand)
        fout.close()

    def get_max_weight(self, segments):
        """ Gets the weights for a list of segments and
            returns the maximal weight in the set.
            The list must be a ptree list.
        """
        print("# Extracting weights...")        

        max_wt = max([tr.wt  for tr in segments])
        print("max_wt = "+str(max_wt))
        self.max_wt = max_wt

        return max_wt

    def wt_mapping(self, wt):
        """ Defines a function which maps the weight on
            the interval [0:1]
        """
        base=1.0e5
        return numpy.log(1.0+(base-1.0)*wt)/numpy.log(base)

    def segmentfilename(self, sid):

      return self.outputdir + '/segment' + str(sid).zfill(6)


###################
# TRUNK EXTRACTOR #
###################
class Basin:

    def __init__(self, timebins, rectangular=False, circular=False, lo=0, hi=0, lt=0, rt=0, ctr=[0,0], rad=0):

        assert( timebins > 0 )
        assert( rectangular != circular )

        if rectangular :
            assert( lo <= hi )
            assert( lt <= rt )

        if circular :
            assert( rad >= 0 )
            assert( len(ctr)==2 )

        self.wt = numpy.zeros([timebins])
        self.timebins = timebins

        self.is_rectangular = rectangular
        self.lo = lo
        self.hi = hi
        self.lt = lt
        self.rt = rt

        self.is_circular = circular
        self.ctr = ctr
        self.rad = rad

class TrunkExtractor(object):

    def __init__(self, filename, nPSbins, nTbins, psmax, histtime, rc_projection_mode='default', extra_tag=None):

        # Make the forest accessible        
        self.pf=p5.Pforest(p5.h5.File(filename,'r'))

        # The first layer of branches is the schools
        self.schools = self.pf.children

        # Set some default lists empty
        self.trunks = []        
        self.displacements = []
        self.displacements_shifted = []
        self.displacements_symmetric = []

        # Convert to floats to be sure
        histtime_f = float(histtime)
        psmax_f    = float(psmax)

        # Initialize histogram arrays
        # H = histogram, E = bin edges
        self.nPSbins = int(nPSbins)
        self.nTbins = int(nTbins)
        self.normed = True
        self.rc_projection_mode = rc_projection_mode
          # set to 'default', 'promoters', 'Hb-Kni' or 'Kr-Gt' (see below for more detail)
        self.rc_velocities_min_t = 3600.0   # minimial time from which on trajectories are
                                            # considered for the velocity averaging

        # Extra tag for output directory
        if extra_tag != None:
            extra_tag = '_' + str(extra_tag)
        else:
            extra_tag = ''
        # Construct output directory name
        if self.rc_projection_mode == 'default':
            self.outputdir = './trunkdata' + extra_tag + '/'
        else:
            self.outputdir = './trunkdata_' + str(self.rc_projection_mode) + extra_tag + '/'

        # Common histogram ranges
        # Abbreviations:   rc = reaction coordinate
        #                  ps = phase space
        self.range_rc  = [[0,histtime_f],[0.0,1.0]]
        self.range_ps  = [[0,psmax_f],[0,psmax_f]]
        self.tbinsize  = histtime_f / float(self.nTbins)
        self.psbinsize = psmax_f / float(self.nPSbins)

        # weighted histograms
        # reaction coord. vs. time
        self.H_rc = []  # bins
        self.E_rc = []  # edges
        # phasespace
        self.H_ps = []  # bins
        self.E_ps = []  # edges
        # phasespace density snapshots
        self.ps_snapshot_frame  = 300.0  # times t_ps_snapsh +- ps_snapshot_frame will be binned into one histogram
        self.time_tolerance     = 1.01
        t_ps_snapsh_min         = 3600.0 # time of first snapshot
        t_ps_snapsh_max         = histtime_f - self.ps_snapshot_frame # time of first snapshot        
        # list of snapshot times
        self.t_ps_snapsh = [ t_ps_snapsh_min, (histtime_f - t_ps_snapsh_min)/2.0, t_ps_snapsh_max ]
        # we will use dictionaries of numpy histogram objects; the snapshot time acts as a key
        self.H_ps_snapsh = {}
        self.E_ps_snapsh = {}
        self.H_ps_snapsh_Nbin = {} # to count binned datapoints
        # construct dictionaries of empty lists        
        for ts in self.t_ps_snapsh:
            self.H_ps_snapsh[ts] = []
            self.E_ps_snapsh[ts] = []
            self.H_ps_snapsh_Nbin[ts] = 0

        # unweighted histograms
        # reaction coord. vs. time
        self.H_rc_unw = []  # bins
        self.E_rc_unw = []  # edges
        # phasespace
        self.H_ps_unw = []  # bins
        self.E_ps_unw = []  # edges

        # unweighted, unnormed histograms (counts)
        # reaction coord. vs. time
        self.H_rc_cnt = []  # bins
        self.E_rc_cnt = []  # edges
        # phasespace
        self.H_ps_cnt = []  # bins
        self.E_ps_cnt = []  # edges
        
        # Arrays of local displacement statistics
        # We will construct 3 (= N_av_modes) arrays: one in which we average
        # the displacements with respect to the [lx(t),ly(t)] point in phase
        # space (index 0), one in which we assume that the displacement
        # actually corresponds to the interpolated point 
        # [1/2*(lx(t+dt)+lx(t)),1/2*(ly(t+dt)+ly(t))] (index 1, DEPRECATED),
        # and a symmetrized version, in which we average for each timepoint t 
        # the "incoming" and "outgoing" displacements (index 2).
        N_av_modes = 3        
        # Avg. phasespace x-displacements
        self.A_ps_dX  = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])
        # Avg. phasespace y-displacements
        self.A_ps_dY  = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])
        # Avg. phasespace squared displacements
        self.A_ps_dR2 = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])        
        # Avg. correlation of the incoming and outgoing trajectories, x-component
        self.A_ps_cX  = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])
        # Avg. correlation of the incoming and outgoing trajectories, y-component
        self.A_ps_cY  = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])
        # Auxiliary arrays for the variances of incoming and outgoing trajectories
        # used to normalized the above two quantities
        self.A_ps_dX2_in  = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])
        self.A_ps_dX2_out = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])
        self.A_ps_dY2_in  = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])
        self.A_ps_dY2_out = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])
        self.A_ps_dX1_in  = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])
        self.A_ps_dX1_out = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])
        self.A_ps_dY1_in  = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])
        self.A_ps_dY1_out = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])
        # Avg. probability flux (time integral), x-component
        self.A_ps_Jx  = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])
        # Avg. probability flux (time integral), y-component
        self.A_ps_Jy  = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])        
        # The total accumulated weight for the bin
        self.A_ps_wt  = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])
        # The total sample count for the bin
        self.A_ps_cnt = numpy.zeros([N_av_modes, self.nPSbins+1, self.nPSbins+1])

        # Arrays of rc moments as a function of time
        self.m_lambda  = numpy.zeros([5, self.nTbins+10])
        self.m_psX     = numpy.zeros([5, self.nTbins+10])
        self.m_psY     = numpy.zeros([5, self.nTbins+10])
        self.m_cnt     = numpy.zeros([5, self.nTbins+10])

        # Functions of basin occupancies as a function of time
        self.basin_numbers = [0, 1, 2, 3, 4, 5, 6]
        self.basin_wt = numpy.zeros([len(self.basin_numbers)+1, self.nTbins+1])
        self.total_wt = numpy.zeros([self.nTbins+1])

        # Construct basins list
        self.basin = []
        self.basin.append( Basin(circular=True,    ctr=[0.3, 0.3], rad=0.15, timebins=self.nTbins+1) )           # stable basin
        self.basin.append( Basin(rectangular=True, lo=0.20, hi=0.43, lt=0.20, rt=0.45, timebins=self.nTbins+1) ) # stable basin v2
        self.basin.append( Basin(rectangular=True, lo=0.00, hi=0.43, lt=0.45, rt=1.00, timebins=self.nTbins+1) ) # right destroyed basin
        self.basin.append( Basin(rectangular=True, lo=0.43, hi=1.00, lt=0.00, rt=0.45, timebins=self.nTbins+1) ) # upper destroyed basin
        self.basin.append( Basin(rectangular=True, lo=0.43, hi=1.00, lt=0.45, rt=1.00, timebins=self.nTbins+1) ) # upper right (final) destroyed basin

        self.basin.append( Basin(rectangular=True, lo=0.00, hi=0.05, lt=0.40, rt=1.00, timebins=self.nTbins+1) ) # Hb-Kni coordinates: Kni destroyed
        self.basin.append( Basin(rectangular=True, lo=0.40, hi=1.00, lt=0.00, rt=0.05, timebins=self.nTbins+1) ) # Hb-Kni coordinates: Hb destroyed

        assert len(self.basin) == len(self.basin_numbers)

        self.rc_cutoff = 1.0

        print("# Created "+self.__class__.__name__+" with filename = "+filename)
        

    def extract(self):
        """ Constructs an enumerated list of all paths, i.e. a list of 
            tuples (path no, path data).
        """
        print("# Extracting trunks...")
        #end_segments = [tr for tr in self.pf if not tr.branch_no]
        self.trunks = [ (tr.trunk, tr.wt) for tr in self.pf.preorder() if not tr.trunk == [] ]


    def process_reaction_coordinates(self, rclist, mode='default'):

        assert mode in ['default', 'promoters', 'Hb-Kni', 'Kr-Gt', 'differences']

        if mode != 'promoters':
        # standard reaction coordinate set
        # the two components are the max total copy number 
        # of the two antagonistic genes
            rc_Hb       = rclist[0]
            rc_Kr       = rclist[1]
            rc_Kni      = rclist[2]
            rc_Gt       = rclist[3]

            rc_maxHbKni = max(rc_Hb, rc_Kni)
            rc_maxKrGt  = max(rc_Kr, rc_Gt)
            rc_total    = rclist[4]

        elif mode == 'promoters':
        # Promoter-based RCs
            rc_Hb       = rclist[6]
            rc_Kr       = rclist[7]
            rc_Kni      = rclist[8]
            rc_Gt       = rclist[9]

            rc_maxHbKni = max(rc_Hb, rc_Kni)
            rc_maxKrGt  = max(rc_Kr, rc_Gt)
            rc_total    = rclist[10]
        
        # Set the coordinates
        # The principal one used for biasing is always the same
        rc_lambda = (rc_maxHbKni + rc_maxKrGt) / rc_total

        # The phase space components differ from case to case
        if mode == 'default' or mode == 'promoters':
        # Plotting the max of two antagonists vs. the max of the other pair
        # in phase space (both for the standard and the promoter based case)
            rc_psX = rc_maxHbKni / rc_total
            rc_psY = rc_maxKrGt  / rc_total

        elif mode == 'Hb-Kni':
        # Plotting Kni vs. Hb in phase space (normalized by total protein)    
            rc_psX = rc_Hb  / rc_total
            rc_psY = rc_Kni / rc_total

        elif mode == 'Kr-Gt':
        # Plotting Kni vs. Hb in phase space (normalized by total protein)    

            rc_psX = rc_Kr / rc_total
            rc_psY = rc_Gt / rc_total

        elif mode == 'differences':
        # Plotting (normed) [Hb]-[Kni] vs [Kr]-[Gt]
        # Values should be between 0 and 1 and positive, so we shift and
        # renormalize a bit
        
            rc_psX = 0.5 + 0.5 * (rc_Hb - rc_Kni) / rc_total
            rc_psY = 0.5 + 0.5 * (rc_Kr - rc_Gt)  / rc_total


        return rc_lambda, rc_psX, rc_psY


    def accumulate_data(self, histtime):
        
        assert(self.trunks != [])

        print("# Accumulating data, histtime = "+str(histtime)+", tbinsize = "+str(self.tbinsize)+" ...")
        
        # Initialize lists for data point accumulation
        # These later will be binned into numpy dd-histograms
        list_H_rc = []
        list_H_ps = []
        list_wt = []
        # Same for the density snapshot histograms
        # Here we use dictionaries with snapshot times as keys
        list_H_ps_snapsh = {}
        list_H_ps_snapsh_wt = {}
        for ts in self.t_ps_snapsh:
            list_H_ps_snapsh[ts] = []
            list_H_ps_snapsh_wt[ts] = []

        #for twt_id, twt_data in self.trunks:

        for twt in self.trunks:

            trunk = twt[0]      # trunk with weight
            wt    = twt[1]      # weight
           
            have_previous_rc_values = 0
            pre_last_time           = -1
            last_time               = -1
            last_timediff           = -1
            last_wt                 = -1
            last_rc_psX             = -1
            last_rc_psY             = -1
            this_timediff           = -1

            for data_tuple in trunk:

                time = data_tuple[0]                
                
                if time>0 and time <= histtime :

                    # Check whether timeline is correct
                    assert time > last_time
                    assert time > pre_last_time
                    if last_time > -1:
                        assert last_time > pre_last_time
                                        
                    if have_previous_rc_values > 0:

                        this_timediff = time - last_time

                        if this_timediff != last_timediff and last_timediff > -1:
                            print "# Warning: %s==this_timediff != last_timediff==%s" % (this_timediff, last_timediff)

                        # will set last_timediff = this_timediff only below, because we need last_timediff


                    # Extract the relevant reaction coordinates
                    # Get the reaction coordinate list
                    rclist = data_tuple[1]
                    
                    # Filter out the reaction coordinates                    
                    rc_lambda, rc_psX, rc_psY = self.process_reaction_coordinates(rclist, mode=self.rc_projection_mode)
                    # Coordinates in phasespace:
                    #   psX = phasespace X
                    #   psY = phasespace Y

                    # Calculate the displacements in phasespace
                    # as compared to the values at t-dt
                    if have_previous_rc_values > 1:
                          
                        rc_dX  = rc_psX - last_rc_psX
                        rc_dY  = rc_psY - last_rc_psY
                        # Also calculate squared displacment(s)
                        rc_dX2 = rc_dX * rc_dX
                        rc_dY2 = rc_dY * rc_dY
                        rc_dR2 = rc_dX2 + rc_dY2
                        
                        # In the following two cases, the correlations between incoming and outgoing trajectories
                        # are force-set to zero because in these cases they are meaningless (compare to 3rd case below!)
                        rc_Corr_dX = 0.0
                        rc_Corr_dY = 0.0
                        # Same for the variances: quantities undefined
                        rc_dX2_in  = 0.0
                        rc_dX2_out = 0.0
                        rc_dY2_in  = 0.0
                        rc_dY2_out = 0.0
                        # and first moments
                        rc_dX1_in  = 0.0
                        rc_dX1_out = 0.0
                        rc_dY1_in  = 0.0
                        rc_dY1_out = 0.0
                        
                        # Store the data, which will be averaged later in another subroutine
                        self.displacements.append( (last_wt, last_time, last_rc_psX, last_rc_psY, rc_dX, rc_dY, rc_dR2, \
                                                    rc_Corr_dX, rc_Corr_dY, rc_dX2_in, rc_dX2_out, rc_dY2_in, rc_dY2_out, \
                                                    rc_dX1_in, rc_dX1_out, rc_dY1_in, rc_dY1_out) )

                        self.displacements_shifted.append( (0.5*(wt+last_wt), 0.5*(time+last_time), \
                                                            0.5*(rc_psX+last_rc_psX), 0.5*(rc_psY+last_rc_psY), rc_dX, rc_dY, rc_dR2, \
                                                            rc_Corr_dX, rc_Corr_dY, rc_dX2_in, rc_dX2_out, rc_dY2_in, rc_dY2_out, \
                                                            rc_dX1_in, rc_dX1_out, rc_dY1_in, rc_dY1_out) )
                                                            # DEPRECATED!
                                                
                        # For the symmetric definition of the displacements, we average the "outgoing"
                        # and the "incoming displacement" with respect to the timepoint with respect to last_time
                        # i.e. we calculate dX = 0.5*(dX_21 - dX_10) = 0.5*[(X2-X1) - (X1-X0)]
                        # where t2 = (current) time, t1 = last_time and t0 = pre_last_time
                        rc_dX_21  = rc_psX - last_rc_psX
                        rc_dX_10  = last_rc_psX - pre_last_rc_psX
                        rc_dY_21  = rc_psY - last_rc_psY
                        rc_dY_10  = last_rc_psY - pre_last_rc_psY
                        
                        rc_dX     = 0.5 * (rc_dX_21 + rc_dX_10) # = 0.5 * (rc_psX - pre_last_rc_psX)
                        rc_dY     = 0.5 * (rc_dY_21 + rc_dY_10) # = 0.5 * (rc_psY - pre_last_rc_psY)
                        # Second moment / variance calculated as above
                        rc_dX2    = rc_dX * rc_dX
                        rc_dY2    = rc_dY * rc_dY
                        rc_dR2    = rc_dX*rc_dX + rc_dY*rc_dY
                        # Here we can also measure the correlation between the incoming and outgoing displacement
                        # as a test for the Markovian assumption
                        rc_Corr_dX = rc_dX_21 * rc_dX_10
                        rc_Corr_dY = rc_dY_21 * rc_dY_10                        
                        # For this we also want the separate second moments (separate variances) of the incoming
                        # vs. the outgoing trajectories, for both coordinates X and Y
                        rc_dX2_in  = rc_dX_10 * rc_dX_10
                        rc_dX2_out = rc_dX_21 * rc_dX_21
                        rc_dY2_in  = rc_dY_10 * rc_dY_10
                        rc_dY2_out = rc_dY_21 * rc_dY_21
                        # First moments as well (just a renaming)
                        rc_dX1_in  = rc_dX_10
                        rc_dX1_out = rc_dX_21
                        rc_dY1_in  = rc_dY_10
                        rc_dY1_out = rc_dY_21
                        # Add the data to the list
                        self.displacements_symmetric.append( (last_wt, last_time, last_rc_psX, last_rc_psY, rc_dX, rc_dY, rc_dR2, \
                                                              rc_Corr_dX, rc_Corr_dY, rc_dX2_in, rc_dX2_out, rc_dY2_in, rc_dY2_out, \
                                                              rc_dX1_in, rc_dX1_out, rc_dY1_in, rc_dY1_out) )

                    # Record phasespace x- and y-values for next
                    # displacement calculation                    
                    pre_last_time   = last_time
                    pre_last_wt     = last_wt
                    pre_last_rc_psX = last_rc_psX
                    pre_last_rc_psY = last_rc_psY
                    last_time       = time
                    last_wt         = wt
                    last_rc_psX     = rc_psX
                    last_rc_psY     = rc_psY
                    last_timediff   = this_timediff
                    have_previous_rc_values += 1

                    # Start binning and moment accumulation
                    if rc_lambda <= self.rc_cutoff:

                        # Build the data points lists to bin
                        list_H_rc.append([time,  rc_lambda])
                        list_H_ps.append([rc_psX, rc_psY])
                        # Don't forget the weights; order must be the same as for data!
                        list_wt.append(wt)

                        # Now the same for the density snapshots
                        for ts in self.t_ps_snapsh:
                            if time >= ts - self.ps_snapshot_frame * self.time_tolerance and \
                               time <= ts + self.ps_snapshot_frame * self.time_tolerance :

                                  list_H_ps_snapsh[ts].append([rc_psX, rc_psY])
                                  list_H_ps_snapsh_wt[ts].append(wt)
                                  self.H_ps_snapsh_Nbin[ts] += 1

                        # Also calculate the running averages in time
                        tindex = int(round(time/self.tbinsize, 0))
                        self.m_lambda[1][tindex] += wt*rc_lambda
                        self.m_lambda[2][tindex] += wt*rc_lambda*rc_lambda
                        self.m_lambda[3][tindex] += wt*rc_lambda*rc_lambda*rc_lambda
                        self.m_psX[1][tindex]    += wt*rc_psX
                        self.m_psX[2][tindex]    += wt*rc_psX*rc_psX
                        self.m_psX[3][tindex]    += wt*rc_psX*rc_psX*rc_psX
                        self.m_psY[1][tindex]    += wt*rc_psY
                        self.m_psY[2][tindex]    += wt*rc_psY*rc_psY
                        self.m_psY[3][tindex]    += wt*rc_psY*rc_psY*rc_psY
                        self.m_cnt[1][tindex]    += wt
                        # ... and the summed weight in the separate basins as a function of time
                        for bn in self.basin_numbers: 

                            if ( self.basin[bn].is_rectangular and \
                                 (rc_psX >= self.basin[bn].lt and rc_psX < self.basin[bn].rt) and \
                                 (rc_psY >= self.basin[bn].lo and rc_psY < self.basin[bn].hi)        )\
                               or \
                               ( self.basin[bn].is_circular and \
                                 (rc_psX - self.basin[bn].ctr[0])**2 + (rc_psY - self.basin[bn].ctr[1])**2 <= self.basin[bn].rad**2 ) :

                                  self.basin[bn].wt[tindex] = self.basin[bn].wt[tindex] + wt

                        # Don't forget the total wt at that time
                        self.total_wt[tindex] = self.total_wt[tindex] + wt

        # Convert lists to arrays, as required by numpy.histogramdd
        array_H_rc = numpy.array(list_H_rc)
        array_H_ps = numpy.array(list_H_ps)
        array_wt   = numpy.array(list_wt)
        
        # weighted histograms
        self.H_rc, self.E_rc = numpy.histogramdd(array_H_rc, bins=self.nTbins, range=self.range_rc, weights=array_wt, normed=self.normed)
        self.H_ps, self.E_ps = numpy.histogramdd(array_H_ps, bins=self.nPSbins, range=self.range_ps, weights=array_wt, normed=self.normed)
        # unweighted histograms
        self.H_rc_unw, self.E_rc_unw = numpy.histogramdd(array_H_rc, bins=self.nTbins, range=self.range_rc, normed=self.normed)
        self.H_ps_unw, self.E_ps_unw = numpy.histogramdd(array_H_ps, bins=self.nPSbins, range=self.range_ps, normed=self.normed)
        # unweighted, unnormed histograms
        self.H_rc_cnt, self.E_rc_cnt = numpy.histogramdd(array_H_rc, bins=self.nTbins, range=self.range_rc, normed=False)
        self.H_ps_cnt, self.E_ps_cnt = numpy.histogramdd(array_H_ps, bins=self.nPSbins, range=self.range_ps, normed=False)

        # density snapshot histograms
        for ts in self.t_ps_snapsh:
            # first again convert lists to arrays
            array_H_ps = numpy.array(list_H_ps_snapsh[ts])
            array_wt   = numpy.array(list_H_ps_snapsh_wt[ts])
            # then, bin it
            self.H_ps_snapsh[ts], self.E_ps_snapsh[ts] = numpy.histogramdd(array_H_ps, bins=self.nPSbins, range=self.range_ps, weights=array_wt, normed=self.normed)


    def calculate_displacement_statistics(self):
        
        # Calculate the local displacement sums for each bin
        # Weight by the weight of the trajectories

        # We construct three sets of average arrays:
        #  - one in which we assume that a displacement belongs to 
        #    its starting point (self.displacements)
        #  - one in which it belongs to the interpolated point between 
        #    start and end point (self.displacements_shifted) [DEPRECATED]
        #  - one in which at a given timepoint t we average the "incoming"
        #    trajectory/displacement (x(t)-x(t-)) with the "outgoing" one
        #    (x(t+)-x(t)) (self.displacements_symmetric)
        
        disp_lists = [self.displacements, self.displacements_shifted, self.displacements_symmetric]

        for li in [0,1,2]:

            total_disp_wt = 0.0

            for (wt, time, rc_psX, rc_psY, rc_dX, rc_dY, rc_dR2, \
                 rc_Corr_dX, rc_Corr_dY, rc_dX2_in, rc_dX2_out, rc_dY2_in, rc_dY2_out, \
                 rc_dX1_in, rc_dX1_out, rc_dY1_in, rc_dY1_out ) in disp_lists[li]:

                if time >= self.rc_velocities_min_t:

                    xi = (int)(rc_psX / self.psbinsize)
                    yi = (int)(rc_psY / self.psbinsize)

                    # Note that the normalization (sum of weights) is
                    # the same for all three calculated quantities
                    self.A_ps_dX[li][xi][yi]  += wt * rc_dX
                    self.A_ps_dY[li][xi][yi]  += wt * rc_dY
                    self.A_ps_dR2[li][xi][yi] += wt * rc_dR2
                    
                    # Correlations between in- and outgoing trajectories
                    # This really only is meaningful for the 3rd case (li=2)
                    self.A_ps_cX[li][xi][yi]      += wt * rc_Corr_dX
                    self.A_ps_cY[li][xi][yi]      += wt * rc_Corr_dY
                    self.A_ps_dX2_in[li][xi][yi]  += wt * rc_dX2_in
                    self.A_ps_dX2_out[li][xi][yi] += wt * rc_dX2_out
                    self.A_ps_dY2_in[li][xi][yi]  += wt * rc_dY2_in
                    self.A_ps_dY2_out[li][xi][yi] += wt * rc_dY2_out
                    self.A_ps_dX1_in[li][xi][yi]  += wt * rc_dX1_in
                    self.A_ps_dX1_out[li][xi][yi] += wt * rc_dX1_out
                    self.A_ps_dY1_in[li][xi][yi]  += wt * rc_dY1_in
                    self.A_ps_dY1_out[li][xi][yi] += wt * rc_dY1_out
                    
                    self.A_ps_wt[li][xi][yi]  += wt

                    # Calculate the avg. flux components
                    self.A_ps_Jx[li][xi][yi]  += wt * rc_dX
                    self.A_ps_Jy[li][xi][yi]  += wt * rc_dY
                    total_disp_wt             += wt            

                    # Samples count
                    self.A_ps_cnt[li][xi][yi]  += 1

            # Normalize by the respective weight of the bin,
            # or the total weight for the flux quantities
            assert total_disp_wt > 0.0
            
            # CF = "Centering flag"
            # Set this to 1.0 to use centered moments in the computation of correlations
            # between incoming and outgoing trajectories below;
            # set to 0.0 for using non-centered moments.
            CF = 0.0;

            for xi in range(self.nPSbins):
                for yi in range(self.nPSbins):

                    if self.A_ps_wt[li][xi][yi]:
                    # avoid division by zero
                        self.A_ps_dX[li][xi][yi]  /= self.A_ps_wt[li][xi][yi]
                        self.A_ps_dY[li][xi][yi]  /= self.A_ps_wt[li][xi][yi]
                        self.A_ps_dR2[li][xi][yi] /= self.A_ps_wt[li][xi][yi]
                        
                        # stuff needed for the correlation coefficient btw. incoming and outgoing trajectories/displacements
                        self.A_ps_cX[li][xi][yi]      /= self.A_ps_wt[li][xi][yi]
                        self.A_ps_cY[li][xi][yi]      /= self.A_ps_wt[li][xi][yi]                                           
                        # first moments of in- and outgoing trajectories, component-wise
                        self.A_ps_dX1_in[li][xi][yi]  /= self.A_ps_wt[li][xi][yi]
                        self.A_ps_dX1_out[li][xi][yi] /= self.A_ps_wt[li][xi][yi]
                        self.A_ps_dY1_in[li][xi][yi]  /= self.A_ps_wt[li][xi][yi]
                        self.A_ps_dY1_out[li][xi][yi] /= self.A_ps_wt[li][xi][yi]
                        # second moments, either centered or not, depending on the value of CF (see above)
                        self.A_ps_dX2_in[li][xi][yi]  = self.A_ps_dX2_in[li][xi][yi]  / self.A_ps_wt[li][xi][yi] - CF*self.A_ps_dX1_in[li][xi][yi]*self.A_ps_dX1_in[li][xi][yi]
                        self.A_ps_dX2_out[li][xi][yi] = self.A_ps_dX2_out[li][xi][yi] / self.A_ps_wt[li][xi][yi] - CF*self.A_ps_dX1_out[li][xi][yi]*self.A_ps_dX1_out[li][xi][yi]
                        self.A_ps_dY2_in[li][xi][yi]  = self.A_ps_dY2_in[li][xi][yi]  / self.A_ps_wt[li][xi][yi] - CF*self.A_ps_dY1_in[li][xi][yi]*self.A_ps_dY1_in[li][xi][yi]
                        self.A_ps_dY2_out[li][xi][yi] = self.A_ps_dY2_out[li][xi][yi] / self.A_ps_wt[li][xi][yi] - CF*self.A_ps_dY1_out[li][xi][yi]*self.A_ps_dY1_out[li][xi][yi]                                                        

                    # Replace the correlations between in- and outgoing trajectories
                    # in the bin by the correlation coefficient, by normalizing with
                    # the respective (bin-wise) variances (if nonzero)
                    # Depending on the value of CF (1 or 0), we will calculate the centered or non-centered variant
                    # (1) x-component of the given obseravables
                    if self.A_ps_dX2_in[li][xi][yi]<>0.0 and self.A_ps_dX2_out[li][xi][yi]<>0.0:
                        self.A_ps_cX[li][xi][yi] = (self.A_ps_cX[li][xi][yi] - CF*self.A_ps_dX1_in[li][xi][yi]*self.A_ps_dX1_out[li][xi][yi]) \
                                                      / numpy.sqrt(numpy.abs(self.A_ps_dX2_in[li][xi][yi])) / numpy.sqrt(numpy.abs(self.A_ps_dX2_out[li][xi][yi]))                                                      
                    # (2) y-component of the given obseravables
                    if self.A_ps_dY2_in[li][xi][yi]<>0.0 and self.A_ps_dY2_out[li][xi][yi]<>0.0:
                        self.A_ps_cY[li][xi][yi] = (self.A_ps_cY[li][xi][yi] - CF*self.A_ps_dY1_in[li][xi][yi]*self.A_ps_dY1_out[li][xi][yi]) \
                                                      / numpy.sqrt(numpy.abs(self.A_ps_dY2_in[li][xi][yi])) / numpy.sqrt(numpy.abs(self.A_ps_dY2_out[li][xi][yi]))
                    
                    # Normalize the fluxes
                    self.A_ps_Jx[li][xi][yi]  /= total_disp_wt
                    self.A_ps_Jy[li][xi][yi]  /= total_disp_wt


    def write_one_histogram(self, H, E, nbins, name, comment=''):

        # TODO Assert that H and E have the right format,
        # i.e. are output by numpy.histogramdd

        filename = self.outputdir + '/' + name
        print "# Writing " + filename
        fout = open(filename, "w")

        if comment != '':
            fout.write('# ' + str(comment) + '\n')

        dx = (E[0][1] - E[0][0])/2.0
        dy = (E[1][1] - E[1][0])/2.0

        for i in range(nbins):

            this_x = E[0][i] + dx

            for j in range(nbins):

                  this_y = E[1][j] + dy
                  value =  H[i][j]

                  fout.write('%s %s %s\n' % (this_x, this_y, value) )

            fout.write('\n')

        fout.close

    def write_histograms(self):

        # Create the folder
        if not os.path.exists(self.outputdir):
            print("# Creating directory " + self.outputdir)
            os.mkdir(self.outputdir)
        else:
            print("# Output directory exists, omitting creation...")
        
        # Histogram (time, lambda)
        self.write_one_histogram(self.H_rc, self.E_rc, self.nTbins, 'hist_rc.dat')
        # Histogram (lambda_x, lambda_y)
        self.write_one_histogram(self.H_ps, self.E_ps, self.nPSbins, 'hist_ps.dat')

        # Unweighted Histogram (time, lambda)
        self.write_one_histogram(self.H_rc_unw, self.E_rc_unw, self.nTbins, 'hist_rc_unw.dat')
        # Unweighted Histogram (lambda_x, lambda_y)
        self.write_one_histogram(self.H_ps_unw, self.E_ps_unw, self.nPSbins, 'hist_ps_unw.dat')

        # Unweighted, unnormed Histogram (time, lambda)
        self.write_one_histogram(self.H_rc_cnt, self.E_rc_cnt, self.nTbins, 'hist_rc_cnt.dat')
        # Unweighted, unnormed Histogram (lambda_x, lambda_y)
        self.write_one_histogram(self.H_ps_cnt, self.E_ps_cnt, self.nPSbins, 'hist_ps_cnt.dat')

        # Phase space density snapshots (weighted)
        s = 1
        for ts in self.t_ps_snapsh:
            
            comm = 'snapshot time = ' + str(ts) + ', total binned (x,y) points = ' + str(self.H_ps_snapsh_Nbin[ts])
            self.write_one_histogram(self.H_ps_snapsh[ts], self.E_ps_snapsh[ts], self.nPSbins, 'hist_ps_snapsh_' + str(s) + '.dat', comment = comm)
            s = s+1


    def write_displacement_statistics(self):

        # Create the folder
        if not os.path.exists(self.outputdir):
            print("# Creating directory " + self.outputdir)
            os.mkdir(self.outputdir)
        else:
            print("# Output directory exists, omitting creation...")        

        filenames  = ['displacement_stats.dat','displacement_stats_shifted.dat','displacement_stats_sym.dat']

        for li in [0,1,2]:

            # Open the output file
            name = filenames[li]
            filename = self.outputdir + '/' + name
            print "# Writing " + filename
            fout = open(filename, "w")

            # Go through the average arrays and output them
            for xi in range(self.nPSbins):
                for yi in range(self.nPSbins):

                      this_x = (xi+0.5) * self.psbinsize
                      this_y = (yi+0.5) * self.psbinsize

                      avg_dX  = self.A_ps_dX[li][xi][yi]
                      avg_dY  = self.A_ps_dY[li][xi][yi]
                      avg_dR2 = self.A_ps_dR2[li][xi][yi]

                      avg_Jx  = self.A_ps_Jx[li][xi][yi]
                      avg_Jy  = self.A_ps_Jy[li][xi][yi]
                      
                      avg_cX  = self.A_ps_cX[li][xi][yi]
                      avg_cY  = self.A_ps_cY[li][xi][yi]

                      cnt     = self.A_ps_cnt[li][xi][yi]

                      fout.write('%s %s %s %s %s %s %s %s %s %s\n' \
                                  % (this_x, this_y, avg_dX, avg_dY, avg_dR2, avg_Jx, avg_Jy, cnt, avg_cX, avg_cY) )

                fout.write('\n')

            fout.close


    def write_running_averages(self):        

        filename = self.outputdir + '/' + 'running_avg.dat'
        print "# Writing " + filename
        fout = open(filename, "w")
        
        # Reminder:
        # psX = phasespace X coord. = typically either [Hb] or max([Hb],[Kni])
        # psY = phasespace Y coord. = typically max([Kr],[Gt])
        for tindex in range(self.nTbins):

            time          = tindex * self.tbinsize
            total_cnt     = self.m_cnt[1][tindex]
            mean_lambda   = self.m_lambda[1][tindex] / total_cnt
            mean_psX      = self.m_psX[1][tindex]    / total_cnt
            mean_psY      = self.m_psY[1][tindex]    / total_cnt

            var_lambda    = self.m_lambda[2][tindex] / total_cnt - mean_lambda*mean_lambda
            var_psX       = self.m_psX[2][tindex]    / total_cnt - mean_psX*mean_psX
            var_psY       = self.m_psY[2][tindex]    / total_cnt - mean_psY*mean_psY

            skew_lambda   = (self.m_lambda[3][tindex]/total_cnt  - 3*mean_lambda*var_lambda   - mean_lambda**3)  / (var_lambda**(3.0/2.0))
            skew_psX      = (self.m_psX[3][tindex]/total_cnt     - 3*mean_psX*var_psX         - mean_psX**3)     / (var_psX**(3.0/2.0))
            skew_psY      = (self.m_psY[3][tindex]/total_cnt     - 3*mean_psY*var_psY         - mean_psY**3)     / (var_psY**(3.0/2.0))

            fout.write('%s %s %s %s %s %s %s %s %s %s %s %s\n' % (tindex, time, total_cnt, mean_lambda, var_lambda, skew_lambda, mean_psX, var_psX, skew_psX, mean_psY, var_psY, skew_psY) )

        fout.close

    def write_basin_occupancies(self):        

        filename = self.outputdir + '/' + 'basin_occupancies.dat'
        print "# Writing " + filename
        fout = open(filename, "w")
        
        fout.write('# basin boundaries:\n')
        for bn in self.basin_numbers :
            if self.basin[bn].is_rectangular:
                fout.write('# basin no. %s (rectangular) : lo = %s / hi = %s / lt = %s / rt = %s \n' % (bn, self.basin[bn].lo, self.basin[bn].hi, self.basin[bn].lt, self.basin[bn].rt) )
            else:
                fout.write('# basin no. %s (circular) : ctr = (%s, %s) / rad = %s \n' % (bn, self.basin[bn].ctr[0], self.basin[bn].ctr[1], self.basin[bn].rad) )
        fout.write('# FORMAT: time_index, time, total_wt, d_total_wt/d_time, wt_in_basin[0], d_wt_in_basin[0]/d_time, wt_in_basin[1], [..]\n#\n#\n')

        for tindex in range(self.nTbins):                        

            time = tindex * self.tbinsize

            # Output total weight and time derivative of the total weight
            if tindex: # calculate the derivative
                der = (self.total_wt[tindex]-self.total_wt[tindex-1])/self.tbinsize
            else:
                der = None

            fout.write('%s %s %s %s ' % (tindex, time, self.total_wt[tindex], der ))

            # Now do the same for each individual basin
            for bn in self.basin_numbers :

                occ = self.basin[bn].wt[tindex]

                if tindex: # calculate the derivative
                    der = (self.basin[bn].wt[tindex] - self.basin[bn].wt[tindex-1]) / self.tbinsize
                else: 
                    der = None

                fout.write('%s %s ' % (occ, der) )
            
            fout.write('\n')

        fout.close

    def write(self): # TODO
        """ Creates a folder containing all paths and write them
            into separate files into this folder.
        """
        assert(self.paths != [])

        # Create the folder
        if not os.path.exists(self.outputdir):
            print("# Creating directory " + self.outputdir)
            os.mkdir(self.outputdir)
        else:
            print("# Output directory exists, omitting creation...")

        # Now go through the path list and output each path
        # into an own file in the created folder
        for pid, path_with_names in self.paths:
            
          #path, pname = path_and_name
          filename = self.pathfilename(pid)
          print "# Writing " + filename
          fout = open(filename, "w")

          for pwn in path_with_names:
              t, xvector = pwn[0]

              # determine the tree number
              # from the name property
              name = pwn[1]
              split_name = name.split('/')
              tree_id = split_name[1]
              self.colors[str(pid)] = tree_id[1:]
        
              output = str(t)
              for x in xvector:
                  output = output + " " + str(x)
              output = output + " " + str(self.colors[str(pid)]) + "\n"
              fout.write(output)

          fout.close()

    def create_gnuplot_file(self, plot): # TODO
        """ Creates a rudimentary gnuplot script for path files
        """

        # Determine filename
        filename = self.outputdir + '/plot_' + plot.name + '.gnu'
        print "# Writing " + filename
        fout = open(filename, "w")          

        # Write a preamble with some settings
        fout.write('set autoscale\n')
        fout.write('unset key\n')
        fout.write('set xrange '+plot.xrange+'\n')
        fout.write('set yrange '+plot.yrange+'\n')
        fout.write('set xlabel \''+plot.xlabel+'\'\n')
        fout.write('set ylabel \''+plot.ylabel+'\'\n')
        fout.write('set xtics '+plot.xtics+'\n')
        fout.write('set ytics '+plot.ytics+'\n')
        fout.write('set format x '+plot.xformat+'\n')
        fout.write('set format y '+plot.yformat+'\n')
        fout.write('set terminal postscript color enhanced solid \"Arial\" 20\n')
        fout.write('set output \''+plot.name+'.ps\'\n')

        # Write the constructed plot command
        plotdata = '\'' + self.outputdir + plot.name +'.dat\''

        lstyle = " w l lw 5 "
        # TODO The if statement here is not very elegant...
        if plot.name == "running_avg" :

            av_rows = plot.rows.replace('X', "4:(sqrt($5))")            
            av_style = "w yerrorbars lw 5 "
            av_title = " t \"average ${/Symbol l}\" "

            sk_rows = plot.rows.replace('X', "6")
            sk_style = lstyle
            sk_title = " t \"skewness\" "

            plotcommand = plot.command + ' ' + plotdata + ' u ' + av_rows + " axes x1y1 " + av_style + av_title + ',' \
                                             + plotdata + ' u ' + sk_rows + " axes x1y2 " + sk_style + sk_title + '\n'

        elif plot.name == "basin_occupancies" :

            style = lstyle
            plots=''
            for bn in self.basin_numbers :

                row_no = 5 + 2*bn # basin data starts at row 5; each basin has two rows (second one is the derivative)
                rows = plot.rows.replace('X', str(row_no))
                title = ' t \"basin %s\" ' % bn
                if plots:
                    plots = plots + ',\\\n' + plotdata + ' u ' + rows + style + title
                else :
                    plots = ' ' + plotdata + ' u ' + rows + style + title

            plotcommand = plot.command + ' ' + plots + '\n'

        else : # histograms, i.e. 2D plots (splot and pm3d)

            assert(plot.command == "splot")
            plotcommand = plot.command + ' ' + plotdata + ' u ' + plot.rows + ' w pm3d\n'

        fout.write(plot.definitions+'\n')
        fout.write(plotcommand)
        fout.close()

    def trunkfilename(self, pid):

      return self.outputdir + '/trunk' + str(pid).zfill(6)


##############################
# AUTOCORRELATION CALCULATOR #
##############################
class AutocorrelationCalculator(object):

    def __init__(self, filename, tpoints, dt):
        # Make the forest accessible        
        self.pf=p5.Pforest(p5.h5.File(filename,'r'))

        # The first layer of branches is the schools
        self.schools = self.pf.children

        # Set some default values
        self.outputdir = './corr_func/'
        self.datafilename = 'corr_func.dat'
        self.trunks = []

        self.dt = dt
        self.tpoints = tpoints + 1
        self.tmin = 0.25*dt*tpoints
        self.tmax = 1.00*dt*tpoints
        self.tdiff_max = 0.50*dt*tpoints

        self.rc_index = 5
        self.rc_cutoff = 0.7

        self.detrend = 1

        # Initialize statistical counters / arrays
        self.c = numpy.zeros([self.tpoints])
        self.w = numpy.zeros([self.tpoints])
        self.n = numpy.zeros([self.tpoints])

        self.rc_wt_sum  = numpy.zeros([self.tpoints+10])
        self.rc2_wt_sum = numpy.zeros([self.tpoints+10])
        self.wt_sum     = numpy.zeros([self.tpoints+10])
        self.wt_cnt     = numpy.zeros([self.tpoints+10])
        # Index for global averages (stored in the same arrays)
        self.g  = self.tpoints + 3       # global mean
        self.gd = self.tpoints + 4       # global mean, detrended

        print("# Created "+self.__class__.__name__+" with filename = "+filename)
        

    def time_index(self, t) :
        """ Just a helper function
        """

        tindex = int(round(t / self.dt, 0))

        assert (tindex < self.tpoints)

        if tindex == 0:
            assert (t == 0.0)

        return tindex


    def within_ranges(self, t, rc) :
        """ Checks whether a time and reaction coordinate pair is within
            the specified ranges
        """

        return (t > self.tmin) and (t <= self.tmax) and (rc <= self.rc_cutoff)


    def iterate(self):
        """ Iterates through all nodes of the tree and calculates the correlations
            between the reaction coordinate values at different t.
        """

        tree_iterator = self.pf
        # TODO Maybe different iteration type is better?

        # First calculate the running and global average
        for node in tree_iterator :

            if node.trunk :

                for T, RC_array in node.trunk :

                    RC_val = RC_array[self.rc_index]
                    wt     = node.wt

                    if self.within_ranges(T, RC_val) :
                    # Only ever do sth. if the time and RC is in the specified range

                        Ti = self.time_index(T)

                        # Running average
                        self.rc_wt_sum[Ti]  += (wt * RC_val)
                        self.rc2_wt_sum[Ti] += (wt * RC_val * RC_val)                        
                        self.wt_sum[Ti]     +=  wt
                        self.wt_cnt[Ti]     +=  1
                        # Global average
                        self.rc_wt_sum[self.g]  += (wt * RC_val)
                        self.rc2_wt_sum[self.g] += (wt * RC_val * RC_val)                        
                        self.wt_sum[self.g]     +=  wt
                        self.wt_cnt[self.g]     +=  1

        
        # Now do the double iteration and calculate the correlations
        for node in tree_iterator :

            if node.trunk :

                for T, RC_array in node.trunk :

                    RC_val = RC_array[self.rc_index]
                    wt     = node.wt

                    if self.within_ranges(T, RC_val) :
                    # Only ever do sth. if the time and RC is in the specified range                        

                        if self.detrend :

                            # Calculate global mean and variance of detrended variables
                            # Make sure RC_val remains unchanged here, for it is still used later!
                            Ti     = self.time_index(T)
                            RC_detrend = self.rc_wt_sum[Ti]/self.wt_sum[Ti]
                            RC_val_d = RC_val - RC_detrend
                            self.rc_wt_sum[self.gd]  += (wt * RC_val_d)
                            self.rc2_wt_sum[self.gd] += (wt * RC_val_d * RC_val_d)
                            self.wt_sum[self.gd]     +=  wt
                            self.wt_cnt[self.gd]     +=  1

                        # Now do the second loop to calculate correlations with
                        # rc values at distant timepoints
                        for t, rc_array in node.trunk_from_root :

                            rc_val = rc_array[self.rc_index]
                            # RC_val is already calculated above

                            if (t > self.tmin) and (t <= T) and (rc_val <= self.rc_cutoff) :
                            # t > 0.0:  rc-values for t=0 are artificial (not measured)
                            # t <= T:   only correlate with earlier time points

                                ti     = self.time_index(t)

                                tdiff  = T - t
                                tdiffi = self.time_index(tdiff)

                                if (tdiff <= self.tdiff_max) and (tdiffi < self.tpoints) :
                                # if tdiff is smaller than running frame size and histogram size

                                    # Detrending
                                    if self.detrend :
                                        rc_detrend = self.rc_wt_sum[ti]/self.wt_sum[ti]
                                        RC_detrend = self.rc_wt_sum[Ti]/self.wt_sum[Ti]
                                    else :
                                        rc_detrend = 0.0
                                        RC_detrend = 0.0

                                    # Accumulate the correlations and weights
                                    self.c[tdiffi] += (wt * (rc_val-rc_detrend)*(RC_val-RC_detrend))
                                    self.w[tdiffi] +=  wt
                                    self.n[tdiffi] +=  1

    def write(self):
        """ Creates the output directory and writes the correlations file
        """

        # TODO Check whether arrays contain data
        print "Mean recorded weight = " + str(self.wt_sum[self.g]/self.wt_cnt[self.g])

        # Create the folder
        if not os.path.exists(self.outputdir):
            print("# Creating directory " + self.outputdir)
            os.mkdir(self.outputdir)
        else:
            print("# Output directory exists, omitting creation...")

        
        filename = self.outputdir+'/'+self.datafilename
        print "# Writing " + filename
        fout = open(filename, "w")

        # Choose the right global mean index
        if self.detrend :
            g = self.gd
        else :
            g = self.g
          
        glob_mean = self.rc_wt_sum[g]  / self.wt_sum[g]
        glob_var  = self.rc2_wt_sum[g] / self.wt_sum[g] - (glob_mean*glob_mean)

        for t in range(0, self.tpoints) :                            
        
            time = t * self.dt
            mean = self.rc_wt_sum[t]  / self.wt_sum[t]
            var  = self.rc2_wt_sum[t] / self.wt_sum[t] - (mean*mean)

            output = str(t)+" "+str(time)+" "+str(self.c[t])+" "+str(self.w[t])+" "+str(self.n[t])\
                           +" "+str(mean)+" "+str(var)+" "+str(glob_mean)+" "+str(glob_var)+"\n"
            fout.write(output)

        fout.close()

    def create_gnuplot_file(self, plot):
        """ Creates a rudimentary gnuplot script for path files
        """

        # Determine filename
        filename = self.outputdir + '/plot_' + plot.name + '.gnu'
        print "# Writing " + filename
        fout = open(filename, "w")          

        # Write a preamble with some settings
        fout.write('set autoscale\n')
        fout.write('set xrange '+plot.xrange+'\n')
        fout.write('set yrange '+plot.yrange+'\n')
        fout.write('set xlabel \''+plot.xlabel+'\'\n')
        fout.write('set ylabel \''+plot.ylabel+'\'\n')
        fout.write('set xtics '+plot.xtics+'\n')
        fout.write('set ytics '+plot.ytics+'\n')
        fout.write('set format x '+plot.xformat+'\n')
        fout.write('set format y '+plot.yformat+'\n')
        fout.write('set terminal postscript color enhanced solid \"Arial\" 20\n')
        fout.write('set output \''+plot.name+'.ps\'\n')
        fout.write(plot.definitions+'\n')

        # Write the constructed plot command
        plotdata = '\'' + self.outputdir+'/'+self.datafilename+'\''
        plotcommand = plot.command + ' ' + plotdata + ' u ' + plot.rows + ' w l lw 5\n'
        fout.write(plotcommand)
        fout.close()

    def trunkfilename(self, pid):

      return self.outputdir + '/trunk' + str(pid).zfill(6)



################
# MAIN ROUTINE #
################
if __name__ == '__main__':

    thisname=sys.argv[0]
    filename=sys.argv[1]
    taskname=sys.argv[2]
    histtime=sys.argv[3]
    nTbins=sys.argv[4]  # time histogram bins; normally sth. like 50.
                        # Make sure that the size of one simulation timebin is a multiple of histtime/nTbins.
    nPSbins=sys.argv[5] # phase space histogram bins
    dt_log=sys.argv[6]  # the logging time step of the FFS simulation (typically 300 s)
    dt_min=sys.argv[7]  # the std. time update in the FFS sim.; this sets the minimal time difference
                        # between two recorded timepoints

    plot_only=bool(int(sys.argv[8]))
    

    total_extract_height = 9*float(histtime)

    format_standard   = "\"%.1f\""
    format_standard2  = "\"%.2f\""
    format_integer    = "\"%.0f\""
    format_scientific = "\"%.1t*10^{%T}\""

    print("\n### TASK = "+str(taskname)+" ###")
    if plot_only:
        print("plot_only = true")

    ## Paths
    if str(taskname) == "paths" :

        PE = PathExtractor(filename, total_extract_height, histtime, dt_log, dt_min)

        if not plot_only:          
          PE.extract()
          PE.write()

        p = plot()
        plot.name = "paths"
        plot.rows = "($1/3600.0):7"
        plot.xlabel = "t [h]"
        plot.ylabel = "{/Symbol l}"
        plot.xrange = "[0:*]"
        plot.yrange = "[0:1]"
        plot.xtics  = ""
        plot.ytics  = "0, 0.2"
        plot.xformat = format_standard
        PE.create_gnuplot_file(p)

    # Segments
    if str(taskname) == "segments" :

        SE = SegmentExtractor(filename, total_extract_height)

        if not plot_only:          
          SE.extract()
          SE.write()

        p = plot()
        plot.name = "segments"
        plot.rows = "($1/3600.0):7"
        plot.xlabel = "t [h]"
        plot.ylabel = "{/Symbol l}"
        plot.xrange = "[0:*]"
        plot.yrange = "[0:1]"
        plot.xtics  = ""
        plot.ytics  = "0, 0.2"
        if float(histtime) / 3600.0 < 1000:
            plot.xformat = format_standard
        else:
            plot.xformat = format_scientific
        SE.create_gnuplot_file(p)

        plot.name = "phasespace"
        plot.definitions="max(x,y)=x>y ? x : y"
        plot.rows = "2:(max($3,$5))"
        plot.xlabel = "[Hb]"
        plot.ylabel = "max([Kr],[Gt])"
        plot.xrange = "[0:*]"
        plot.yrange = "[0:*]"
        plot.xtics  = "0, 2000"
        plot.ytics  = "0, 2000"
        plot.xformat = format_integer
        SE.create_gnuplot_file(p)

    # Trunks and histogram
    if str(taskname) == "histograms" :
        psmax=1.0 # the max. extension of the phase space histogram; 5000 for CN=15 if not normalized, scale up with CN accordingly

        for mode in [ 'default', 'Hb-Kni', 'Kr-Gt', 'differences' ]: # 'default', 'Hb-Kni', 'Kr-Gt', 'differences', 'promoters'

            TE = TrunkExtractor(filename, nPSbins, nTbins, psmax, histtime, rc_projection_mode=mode, extra_tag=str(nPSbins) )

            if not plot_only:              
              TE.extract()
              TE.accumulate_data(histtime)
              TE.calculate_displacement_statistics()
              TE.write_histograms()
              TE.write_displacement_statistics()
              TE.write_running_averages()
              TE.write_basin_occupancies()
            
            p = plot()
                       
            #bg_colour = "dark-blue"
            bg_colour = "black"
            #bg_colour = "white"

            #splot_palette = "set palette; "   # The standard gnuplot black-to-yellow palette; use with black background colour
            #splot_palette = "set palette defined (0 \""+bg_colour+"\", 0.1 \"blue\", 1.0 \"red\"); "                         # use with white or black background
            #splot_palette = "set palette defined (0 \""+bg_colour+"\", 1 \"yellow\", 2 \"green\", 3 \"blue\", 4 \"red\"); "  # use with white or black background
            #splot_palette = "set palette defined (0 \""+bg_colour+"\", 1 \"light-blue\", 2 \"green\", 3 \"yellow\", 4 \"red\", 5 \"dark-red\"); "
            splot_palette = "set palette defined (0 \""+bg_colour+"\", 1 \"blue\", 2 \"dark-violet\", 3 \"magenta\", 4 \"orange\", 5 \"white\"); "
            
            splot_bg = "set object 1 rectangle from graph 0, graph 0 to graph 1, graph 1 behind fc rgbcolor '"+bg_colour+"' fs noborder; set tics out; "

            p.command = "splot"
            p.name = "hist_rc"
            p.definitions="set pm3d map; unset key; " + splot_palette + splot_bg + "set cbrange [-4:1]; set cblabel \"log(p)\";"
            p.rows = "($1/3600.0):2:(log($3*3600.0)/log(10))"
            p.xlabel = "t [h]"
            p.ylabel = "{/Symbol l}"
            p.xrange = "[0:*]"
            p.yrange = "[0:1]"
            p.xtics  = ""
            p.ytics  = "0, 0.2"
            if float(histtime) / 3600.0 < 1000:
                p.xformat = format_standard
            else:
                p.xformat = format_scientific
            TE.create_gnuplot_file(p)
            # unweighted version
            p.name = "hist_rc_unw"
            TE.create_gnuplot_file(p)
            # raw counts (unnormalized)
            p.name = "hist_rc_cnt"
            p.definitions="set pm3d map; unset key; " + splot_palette + splot_bg + "set cbrange [0:6]; set cblabel \"log(p)\";"
            TE.create_gnuplot_file(p)

            p.name = "hist_ps"
            p.definitions="set size square; set pm3d map; unset key; " + splot_palette + splot_bg + "set cbrange [-2:3]; set cblabel \"log(p)\";"
            p.rows = "1:2:(log($3)/log(10))"

            if mode=='Hb-Kni':
                p.xlabel = "[Hb]/[all]"
                p.ylabel = "[Kni]/[all]"
            elif mode=='Kr-Gt':
                p.xlabel = "[Kr]/[all]"
                p.ylabel = "[Gt]/[all]"
            elif mode=='differences':
                p.xlabel = "0.5 + 0.5 * ([Hb]-[Kni])/[all]"
                p.ylabel = "0.5 + 0.5 * ([Kr]-[Gt])/[all]"
            else: # 'default' or 'promoters'
                p.xlabel = "max([Hb],[Kni])/[all]"
                p.ylabel = "max([Kr],[Gt])/[all]"            
            p.xrange = "[0.2:0.7]"
            p.yrange = "[0.1:0.6]"
            #p.xtics  = "0, "+str(psmax/2)
            #p.ytics  = "0, "+str(psmax/2)
            #p.xformat = format_integer
            p.xtics  = "0.2, 0.1"
            p.ytics  = "0.1, 0.1"
            p.xformat = format_standard  # set to format_standard2 for float with 2 digits
            TE.create_gnuplot_file(p)
            # density snapshot histograms
            Ns = 4
            for s in range(1,Ns): # make sure Ns corresponds to no. of snapshots + 1 at least            
                p.name = "hist_ps_snapsh_"+str(s)
                TE.create_gnuplot_file(p)            
            # unweighted version
            p.name = "hist_ps_unw"
            TE.create_gnuplot_file(p)
            # raw counts (unnormalized)
            p.name = "hist_ps_cnt"
            p.definitions="set size square; set pm3d map; unset key; " + splot_palette + splot_bg + "set cbrange [0:4]; set cblabel \"log(p)\";"
            TE.create_gnuplot_file(p)

            q = plot()
            q.command = "plot"
            q.name = "basin_occupancies"
            q.definitions=""
            q.rows = "($2/3600.0):(1.0*$X/$3)"
            q.xlabel = "t [h]"
            q.ylabel = "basin wt / total wt"
            q.xrange = "[0:*]"
            q.yrange = "[0:*]"
            q.xtics  = ""
            q.ytics  = "0, 0.2"
            if float(histtime) / 3600.0 < 1000:
                q.xformat = format_standard
            else:
                q.xformat = format_scientific
            TE.create_gnuplot_file(q)

            q.name = "running_avg"
            q.definitions="set ytics nomirror; set y2tics nomirror; set y2label \'skewness\';"
            q.rows = "($2/3600.0):X"
            q.xlabel = "t [h]"
            q.ylabel = "average {/Symbol l}"
            q.xrange = "[0:*]"
            q.yrange = "[0:1]"
            TE.create_gnuplot_file(q)        


    ## Autocorrelations
    if str(taskname) == "correlations" :

        AC = AutocorrelationCalculator(filename, 100, dt_log) # std. (50, 300)

        if not plot_only:
          AC.iterate()
          AC.write()

        p = plot()
        plot.name = "corr_func"
        plot.definitions="unset key;"
        plot.rows = "($2/3600.0):($3/$4)"
        plot.xlabel = "{/Symbol t} [h]"
        plot.ylabel = "unnormalized C_{/Symbol l}({/Symbol t})"
        plot.xrange = "[0:*]"
        plot.yrange = "[*:*]"
        plot.xtics  = ""
        plot.ytics  = ""
        plot.xformat = format_standard        
        #AC.create_gnuplot_file(p)
        
        plot.rows = "($2/3600.0):( (($3/$4)-$8*$8)/$9 )"
        plot.ylabel = "normalized C_{/Symbol l}({/Symbol t})"
        AC.create_gnuplot_file(p)

        plot.name = "corr_func_local_mean"
        plot.rows = "($2/3600.0):( (($3/$4)-$6*$6)/$7 )"
        AC.create_gnuplot_file(p)
