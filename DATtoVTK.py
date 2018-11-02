# JUPITER .dat to .vtk file conversion class
#
# Evert Nasedkin, October 2018
# evertn@student.ethz.ch
#
# This is the class that converts the binary .dat and descriptor files
# output from JUPITER hydrodynamic simulations into .vtk format.
#
# The class can be run from the Convert.py script
#
# Process:
# Set up of  science factors and directories based on the user input
# Read in mesh vertices from descriptor file, convert to cartesian coordinates
# Calculate cell centers
# Read in science data from .dat file
# Setup VTK File format structure
#
# VTK File Format specs: https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
#
# The .dat data file contains a 1D array of doubles written to a binary file
# The descriptor file is structured as follows:
# 8 lines of simulation information, used by JUPITER, followed by
# 3 lines of grid data. In order: azimuthal angle, radius and polar angle
#
# The first and last two grid points in each dimension are overlap with the
# previous mesh resolution, or are just ghost points. They are deleted.
#
# To reconstruct the grid, the coordinates are iterated in column order,
# (ie azimuthal angle is iterated the fastest, polar angle the slowest)
# It is converted to Cartesian coordinates to be written to the VTK file.
#
# Notes:
# On importing the coordinate grid it is changed from left handed to write handed (azimuthal and polar angles)*-1
# Warning: Spherical coordinates are ordered (phi,r,theta). Just go with it. (phi is azimuthal, theta is polar)
#
#
# MAJOR FIXME:
# Need to fix AMR issues:
# Currently, cells at each mesh level are allowed to overlap
# In paraview, this leads to weird issues when viewing surface plots - trying to plot data from different cells to the same location.
# Need to figure out either:
# a) Write out only non-overlapping cells, choosing the finest grid to write
# b) Read in data as AMR, plotting only from finest grid.

import os,sys
import numpy as np
import astropy.units as u
import string
import Dialog
from multiprocessing import Pool

try:
    from pyvtk import *
except ImportError:
    print "Please install pyvtk. (pip install pyvtk)"

# Turns out we need pyvtk and vtk.
try:
    import vtk as vtk
    from vtk import vtkUnstructuredGrid, vtkStructuredGrid
    from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtkIdTypeArray
    from vtk.util.numpy_support import numpy_to_vtk
except:
    print "Please install vtk package (pip install vtk)"

class DATtoVTK:
    'Convert JUPITER .DAT files to binary VTK files for Paraview'
    def __init__(self):
        # Simulation Information
        self.simNumber = -1 # Which sim output number are we looking at
        self.nLevel = -1 # How many mesh levels are there?
        self.feature = 'notafeat'
        self.nLevelCoords = []
        self.nLevelh = []
        # density, temperature,energy, erad, opacity, potential, stheat, tau, taucell, velocity
        self.featlist = ['gasdensity', 'gastemperature','gasenergy',
                         'gaserad', 'gasopacity', 'gaspotential',
                         'gasstheat', 'gastau', 'gastaucell',
                         'gasvelocity','dustdensity', 'dusttemperature',
                         'dustenergy','dusterad', 'dustopacity',
                         'dustpotential','duststheat', 'dusttau',
                         'dusttaucell','dustvelocity']

        # Filepath information
        self.dataDir = 'notadir' # Where is the data from?
        self.dataOutPath = 'notapath' # Where do we output data
        self.inFilename = 'notaname' # Filename constructed from feature, outNumber and nLevels
        self.outFilename = 'notaname'
        self.descriptorName = 'notadesc'
        self.BASEPATH = os.getcwd() + '/'

        # Science information
        self.radius = -1.
        self.rcgs = -1.
        self.mass = -1.
        self.mcgs = -1.
        # Science Constants
        self.TEMP = -1.
        self.DENS = -1.
        self.PERIOD = -1.
        self.VEL = -1.

        # Grid Information
        self.sphere =[]# np.empty(0) # 3D Spherical array - edges (read in and filter)
        self.mesh = np.zeros((0,0,3),dtype=np.float64)   # 3D cartesian array - edges
        self.unfiltered = self.mesh = np.zeros((0,0,3),dtype=np.float64) # Unfiltered spherical coords
        self.ncell = 0 # How many cells are there (filtered)
        self.mlen = [0] # How many coordinates in the previous mesh level (unfiltered)
        self.mlenUF = [0]
        self.coordlist = np.zeros(0) # Which coordinates were filtered out
        self.cellist = [] # Which cells were filtered out
        self.mins = [] # Min boundaries of a each mesh level
        self.maxs = [] # Max boundaries of a each mesh level


    # ------------------------------------------------------------------------------------------------

    # Basic user Set functions
    def SetOutNumber(self, n):
        self.outNumber = n

    def SetLevel(self,level):
        self.nLevel = level

    def SetFeature(self, feat):
        if feat in self.featlist:
            self.feature = feat
            return 1
        else:
            print("Not an recognised input, may not be implemented. Continuing...")
            return 1
    def SetBasePath(path):
        self.BASEPATH = path
    def SetInDir(path):
        self.dataDir = path
    def SetOutDir(path):
        self.dataOutPath = path
    def SetInFilename(fname):
        self.inFilename = fname
    def SetRadius(self, rad):
        self.radius = rad*u.AU
        self.rcgs = self.radius.to(u.cm)
    def SetMass(self,mass):
        self.mass = mass*u.M_sun
        self.mcgs = self.mass.to(u.g)
        if(self.rcgs.value>0.):
            self.TEMP = ((self.rcgs.value)/(np.sqrt((self.rcgs.value)**3 / 6.67259e-8 / (self.mcgs.value))))**2 / 8.314e7
            self.DENS = self.mcgs.value/(self.rcgs.value)**3
            self.PERIOD = 2*np.pi*np.sqrt((self.rcgs.value)**3 / (6.67259e-8 * self.mcgs.value))
            self.VEL = self.rcgs.value/(self.PERIOD/2*np.pi)

    # ----------------------------------------------------------------------------------------

    # Directory and file Setup
    def SetupDirs( self ):
        # Error checking
        if self.outNumber < 0:
            print("Please set the simulation number")
            return
        if self.feature is 'notafeat':
            print("Please input a feature (e.g. velocity)")
            return


        # Create directory paths
        self.dataDir = self.BASEPATH + "output" + str(self.outNumber).zfill(5) + "/"
        self.dataOutPath =  self.BASEPATH + "VTK" + str(self.outNumber).zfill(5) +"/"
        self.descriptorName = "Descriptor" + str(self.outNumber) + ".dat"
        self.outFilename =  self.feature + str(self.outNumber) +"_" + str(self.nLevel) + ".vtk"

        # Check existance
        print("Data directory is: " + self.dataDir)
        print("If this is incorrect, please set a base path to the data directory")
        if not os.path.isdir(self.dataDir):
            print("ERROR: data directory " + self.dataDir + " does not exist!")
            return

        # Make out dir
        if not os.path.isdir(self.dataOutPath):
            os.mkdir(self.dataOutPath)


    def SetupNames( self, level ):
        if self.outNumber < 0:
            print("Please set the simulation number")
            return
        if self.feature is 'notafeat':
            print("Please input a feature (e.g. gasvelocity)")
            return
        if self.nLevel < 0:
            print("Please input the number of mesh refinement levels")
            return

        self.inFilename = self.feature + str(self.outNumber) +"_"+ str(level) + "_" + str(level) + ".dat"


    # ---------------------------------------------------------------------------------------------
    # Important part starts here

    #
    # This function wraps the binary .dat file reader for a given feature,
    # and output for the field in a VTK format. This is the only user facing function.
    # --------------------------------------------------------------------------------------------
    def ConvertFiles( self , binary = True ):
        data = []
        for i in range(self.nLevel):
            self.SetupNames(i)
            #feat = np.fromfile(self.dataDir + self.inFilename, dtype = 'double')

            if "velocity" in self.feature:
                data2 = np.zeros(0,dtype = np.float64)
                data2 = np.append(data2,feat.astype(np.float64))
                print self.VEL
                print data2
                data2 = np.reshape(data2,(3,-1))
                data3 = np.column_stack((data2[1], data2[0], data2[1]))
                data3[:,0] = data2[2]*self.rcgs*np.sin(data2[0])*np.cos(data2[1])
                data3[:,1] = data2[2]*self.rcgs*np.sin(data2[0])*np.sin(data2[1])
                data3[:,2] = data2[2]*self.rcgs*np.cos(data2[0])
                #print data3
                data3 = data3.transpose() # Velocity ordering is weird.
                print data3
                data3 = data3*self.VEL#.value
                print data3
                if i > 0:
                    data = np.concatenate((data,data3), axis = 0)
                else:
                    data = data3

            else:
                # Read in binary doubles into a 1D array
                feat = np.fromfile(self.dataDir + self.inFilename, dtype = 'double')
                data.extend(feat)

        # Compute all of the indices of cell vertices
        # Have to do this here to find out which cells
        # are overlapping, and should not be included
        self.ComputeStructuredIndices()
        inds = self.ComputeIndices()
        # Delete overlapping data points
        print len(data)
        data = self.FilterData(data,self.cellist)
        # Convert to CGS units
        if("density" in self.feature):
            data = [x*self.DENS for x in data]#.value
        if("temperature" in self.feature):
            data = [x*self.TEMP for x in data]#.value
        print str(len(data)) + " data points for " + self.feature
        self.WriteToVTK(data, inds, binary)


    # -------------------------------------------------------
    # GetCoordinates
    # This function reads in the descriptor file
    # and builds the VTK coordinate arrays
    #
    # It then Filters out the overlapping coordinates,
    # and transforms the spherical coordinates to
    # cartesian.
    # -------------------------------------------------------
    def GetCoordinates(self):

        phi = []
        r = []
        th = []
        coords = []
        for i in range(self.nLevel):
            dsc = open(self.dataDir + self.descriptorName)
            for j, line in enumerate(dsc):
                # Invert theta, phi so that coordinate system is right handed
                # (needed for cell orientation)
                if j == 8 + (i*11):
                    cur = [np.float64(x) for x in line.split()]
                    cur.pop(0) # First and last two points are 'ghost points'
                    cur.pop(0)
                    cur.pop()
                    cur.pop()
                    phi.append(cur)
                    cur = []
                if j == 9 + (i*11):
                    cur = [np.float64(x) for x in line.split()]
                    cur.pop(0)
                    cur.pop(0)
                    cur.pop()
                    cur.pop()
                    r.append(cur)
                    cur = []
                if j == 10 + (i*11):
                    cur = [np.float64(x) for x in line.split()]
                    cur.pop(0)
                    cur.pop(0)
                    cur.pop()
                    cur.pop()
                    th.append(cur)
                    cur = []
            self.nLevelCoords.append([len(phi[i]),len(r[i]),len(th[i])])
            dsc.close()
        self.unfiltered = self.BuildGrid(phi,r,th)
        self.sphere,self.coordlist = self.FilterCoords(phi,r,th)
        self.mesh = self.SphereToCart(self.sphere)
        print(str(self.sphere.shape[0]) + " Vertices in filtered grid.")


    # -------------------------------------------------------------------------
    # FilterCoords
    #
    # This function takes in a lists of AMR axis, ie
    #     x1s = [[phi level 1], [phi level 2], ... , [phi level n]]
    #
    # It then compares each level to the next finer level, to check if the outer
    # refinement level should be included or not
    #
    # This prevents issues in Paraview trying to render multiple
    # overlapping AMR levels.
    # ------------------------------------------------------------------------
    def FilterCoords(self, x1s, x2s, x3s):
        # Should have read in the same number of AMR levels for each axis
        try:
            assert(len(x1s) == len(x2s) and len(x1s) == len(x3s))
        except:
            print "Something went wrong reading in the coordinates. Please try again."
            sys.exit(1)
        coordlist = []
        newcoords = []
        newlev = []
        tcount = 0
        fcount = 0
        # Check if each element in a level is in the range of the next level
        for l in range(len(x1s)-1):
            nmin = [np.min(x1s[l+1]),np.min(x2s[l+1]),np.min(x3s[l+1])]
            nmax = [np.max(x1s[l+1]),np.max(x2s[l+1]),np.max(x3s[l+1])]
            self.mins.append(nmin)
            self.maxs.append(nmax)
            curlev,c = self.BuildOneLevel(x1s[l],x2s[l],x3s[l])
            mlen = 0
            for i in range(len(curlev)):
                # If yes, note the coordinate number (for data filtering)
                if (self.InRange(nmin,nmax,curlev[i])):
                    coordlist.append(tcount + i)
                # If not, that coordinate is not included.
                else:
                    newlev.append(curlev[i])
                    mlen+=1
            tcount += c
            fcount += mlen
            self.mlen.append(fcount) # This line is important - fcount if using filtered mesh, tcount else
            self.mlenUF.append(tcount)
            newcoords.extend(newlev)
            newlev = []
        curlev,c = self.BuildOneLevel(x1s[-1],x2s[-1],x3s[-1])
        #tcount += c
        fcount += c
        for i in range(len(curlev)):
            newcoords.append(curlev[i])
        # Return the coordinate array, and the coordinates to delete in reverse order (delete from end)
        return np.array(newcoords),np.sort(np.array(coordlist))

    # -------------------------------------------------------
    # In Range
    # This function is used in FilterCoords
    # It takes in two vectors with the minimum and maximum x,y,z values
    # These provide the boundaries to check if the coordinate provided
    # is in that range or not
    # -------------------------------------------------------
    def InRange(self, minvec, maxvec, coord):
        x = y = z = False
        if(coord[0] > minvec[0] and coord[0] < maxvec[0]):
            x = True
        if(coord[1] > minvec[1] and coord[1] < maxvec[1]):
            y = True
        if(coord[2] > minvec[2] and coord[2] < maxvec[2]):
            z = True
        return x and y and z

    def InPlane(self, minpair, maxpair, pair):
        # Assume thrid component is false
        # i.e., we're not in the forbidden region, bu
        # a line traced along an axis will intersect it.
        x = y = False
        if(pair[0] > minpair[0] and pair[0] < maxpair[0]):
            x = True
        if(pair[1] > minpair[1] and pair[1] < maxpair[1]):
            y = True
        return x and y
    # -------------------------------------------------------
    # BuildOneLevel
    # This function takes in 3 axis and builds an array of
    # 3 vectors, in column ordering (azimuthal, radial, polar)
    # -------------------------------------------------------
    def BuildOneLevel(self, x1, x2, x3):
        coords = []
        c = 0
        for az in x3:
            for ay in x2:
                for ax in x1:
                    coords.append([ax,ay,az])
                    c += 1
        return np.array(coords),c

    # -------------------------------------------------------
    # BuildGrid
    # Takes in multiple axis and builds each level, appending
    # into a single coordinate array.
    # -------------------------------------------------------
    def BuildGrid(self, x1s, x2s, x3s):
        grid = []
        tg = []
        lcount = []
        try:
            assert(len(x1s) == len(x2s) and len(x1s) == len(x3s))
        except:
            print "Something went wrong reading in the coordinates. Please try again."
            sys.exit(1)
        for i in range(len(x1s)):
            self.nLevelh.append([abs(x1s[i][1]-x1s[i][0]),abs(x2s[i][1]-x2s[i][0]),abs(x3s[i][1]-x3s[i][0])])
            tg,c = self.BuildOneLevel(x1s[i], x2s[i], x3s[i])
            grid.extend(tg)
            lcount.append(c)
        return np.array(grid)

    # -----------------------------------------------------
    # SphereToCart
    # Convert a list of spherical coordinates to cartisian
    # -----------------------------------------------------
    def SphereToCart(self, data):
        newcoords = np.zeros(data.shape)
        newcoords[:,0] = data[:,1]*self.rcgs*np.sin(data[:,2])*np.cos(data[:,0])
        newcoords[:,1] = data[:,1]*self.rcgs*np.sin(data[:,2])*np.sin(data[:,0])
        newcoords[:,2] = data[:,1]*self.rcgs*np.cos(data[:,2])
        return newcoords


    # --------------------------------------------------
    # WriteToVTK
    # Checks data formatting, prepares output file
    # --------------------------------------------------
    def WriteToVTK(self, data, inds, binary = True, append = False):
        # Data quality checking
        try:
            assert(len(data) == self.ncell)
        except ValueError:
            print "Error: Number of data points does not match number of mesh elements"
            return
        # File existance check
        if not os.path.isfile(self.dataOutPath + self.outFilename):
            self.VTKFormatOut(data,inds,  binary)
        else:
            userstr = self.ask("Output file already exists, do you want to overwrite or append data?",
                               ("Overwrite","Append","Quit"))
            if userstr is 'Overwrite':
                self.VTKFormatOut(data,inds,binary, False)
                return

            if userstr is 'Append': # not sure if append will work yet
                self.VTKFormatOut(data,inds, binary, True)
                return
            else:
                print "Nothing written to file!"
                return
            
    # -----------------------------------------
    # Compute Planes
    # Compute whether a point lies along a line
    # through the hole in the mesh along a given axis
    # ------------------------------------------
    def ComputePlanes(self,mins,maxs,aline,axis = 0):
        x = y = z = 0
        b0 = b3 = b4 = b7 = False
        if( axis == 0):
            y = 1
            z = 2
        elif(axis == 1):
            x = 1
            y = 0
            z = 2
        elif(axis == 2):
            x = 2
            y = 0
            z = 1
        else:
            print "Not a valid axis, returning False"
            return b0,b3,b4,b7

        b0 = self.InPlane([mins[y],mins[z]],
                          [maxs[y],maxs[z]],
                          [self.sphere[aline[0]][y],self.sphere[aline[0]][z]])
        b3 = self.InPlane([mins[y],mins[z]],
                          [maxs[y],maxs[z]],
                          [self.sphere[aline[3]][y],self.sphere[aline[3]][z]])
        b4 = self.InPlane([mins[y],mins[z]],
                          [maxs[y],maxs[z]],
                          [self.sphere[aline[4]][y],self.sphere[aline[4]][z]])
        b7 = self.InPlane([mins[y],mins[z]],
                          [maxs[y],maxs[z]],
                          [self.sphere[aline[7]][y],self.sphere[aline[7]][z]])
        return b0,b3,b4,b7

    # -----------------------------------------------
    # ComputeStructuredCell
    # Compute the indices of the vertices of a given
    # cell in the mesh, assuming a completed,
    # unfiltered grid
    #
    # Implements basic stride counting
    # -----------------------------------------------
    def ComputeStructuredCell(self,n,ix,iy,iz,filtered = True):
        if filtered:
            k1 = self.mlen[n]
        else:
            k1 = self.mlenUF[n]
        id0 = k1 + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*iy     + ix
        id3 = k1 + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*(iy+1) + ix
        id4 = k1 + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*iy     + ix
        id7 = k1 + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*(iy+1) + ix
        if(n == 0): # Only the base mesh level returns to the original location
            id1 = k1 + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*iy     + ((ix+1) % (self.nLevelCoords[n][0]-1))
            id2 = k1 + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*(iy+1) + ((ix+1) % (self.nLevelCoords[n][0]-1))
            id5 = k1 + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*iy     + ((ix+1) % (self.nLevelCoords[n][0]-1))
            id6 = k1 + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*(iy+1) + ((ix+1) % (self.nLevelCoords[n][0]-1))
        else:
            id1 = k1 + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*iy     + (ix+1)
            id2 = k1 + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*(iy+1) + (ix+1)
            id5 = k1 + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*iy     + (ix+1)
            id6 = k1 + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*(iy+1) + (ix+1)
        return id0,id1,id2,id3,id4,id5,id6,id7

    # -----------------------------------------------
    # ComputeCell
    # Compute the indices of the vertices of a given
    # cell in the mesh.
    # k1-8 allow for holes in the mesh when counting
    # indices
    # -----------------------------------------------
    def ComputeCell(self,n,ix,i3x,i4x,i7x,iy,iz,k1,k2,k3,k4,k5,k6,k7,k8):
        # ----------------------------------------------------------------------------------------
        # Write out the indices of the cell interface mesh that define a
        # hexahedron (VTK cell type #12)
        #
        # The indexing of a hexahedron is as follows
        #
        #                  7________6
        #                 /|      / |
        #                / |     /  |
        #               4_------5   |         z  th ^   ^ r y
        #               |  3____|___2               |  /
        #               | /     |  /                | /
        #               |/      | /                 |/
        #               0-------1                   0----->phi x
        #
        #
        # k1 - Number of points in previous planes up to current iz plane (not inclusive)
        # k2 - Remaining verticies of the current plane, plus the number of vertices up to index 4 in the above Plane
        # k3 - Number of completed, full length x axes in current plane
        # k4 - Length of the current x axis if left of the hole, length of the next x axis if right of the hole
        #    - Left and right are determined if strictly to the left of the leftmost border of the hole
        #    - gives number of points from i0 to i3
        # k5 - Length of a short x axis
        # k6 - k3, but for the top plane. Number of completed, full length x axes
        # k7 - k4, but for the top plane. The number of points from i4 to i7
        # k8[0-3] - Length of each x axis
        #         - only relevent for outer mesh level
        #         - used to wrap around to initial coord
        #
        # ---------------------------------------------------------------------------------------
        #print k1,k2,k3,k4

        id0 = self.mlen[n] + k1      + k3*(self.nLevelCoords[n][0]) + (iy-k3)*k5      + ix
        id1 = self.mlen[n] + k1      + k3*(self.nLevelCoords[n][0]) + (iy-k3)*k5      + ix +1
        id2 = self.mlen[n] + k1      + k3*(self.nLevelCoords[n][0]) + (iy-k3)*k5 + k4 + i3x+1
        id3 = self.mlen[n] + k1      + k3*(self.nLevelCoords[n][0]) + (iy-k3)*k5 + k4 + i3x
        id4 = self.mlen[n] + k1 + k2 + k6*(self.nLevelCoords[n][0]) + (iy-k6)*k5      + i4x
        id5 = self.mlen[n] + k1 + k2 + k6*(self.nLevelCoords[n][0]) + (iy-k6)*k5      + i4x+1
        id6 = self.mlen[n] + k1 + k2 + k6*(self.nLevelCoords[n][0]) + (iy-k6)*k5 + k7 + i7x+1
        id7 = self.mlen[n] + k1 + k2 + k6*(self.nLevelCoords[n][0]) + (iy-k6)*k5 + k7 + i7x
        if(n == 0): # Only the base mesh level returns to the original location
            if((ix+1)  % (k8[0]-1)) == 0:
                id1 = self.mlen[n] + k1      + k3*(self.nLevelCoords[n][0]) + (iy-k3)*k5
            if((i3x+1)  % (k8[1]-1)) == 0:
                id2 = self.mlen[n] + k1      + k3*(self.nLevelCoords[n][0]) + (iy-k3)*k5 + k4
            if((i4x+1)  % (k8[2]-1)) == 0:
                id5 = self.mlen[n] + k1 + k2 + k6*(self.nLevelCoords[n][0]) + (iy-k6)*k5
            if((i7x+1)  % (k8[3]-1)) == 0:
                id6 = self.mlen[n] + k1 + k2 + k6*(self.nLevelCoords[n][0]) + (iy-k6)*k5 + k7
        return id0,id1,id2,id3,id4,id5,id6,id7

    
    # ---------------------------------------------------
    # InitialCell
    # Computes the indices of the first cell of each mesh
    # refinement level.
    #
    # Effectively reimplements ComputeCell, but we haven't
    # computed all of the necessary quantities yet
    # ----------------------------------------------------
    def InitialCell(self,n,dtx,dty):
        id0 = self.mlen[n]
        id1 = self.mlen[n] + 1
        id2 = self.mlen[n] + self.nLevelCoords[n][0] + 1
        id3 = self.mlen[n] + self.nLevelCoords[n][0]
        id4 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1] 
        id5 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1] + 1
        id6 = id4 + self.nLevelCoords[n][0]
        id7 = id4 + self.nLevelCoords[n][0] + 1

        if (n<(self.nLevel-1)):
            if self.InPlane([self.mins[n][1],self.mins[n][2]],
                            [self.maxs[n][1],self.maxs[n][2]],
                            [self.sphere[id0][1],self.sphere[id0][2]]):
                id3 = self.mlen[n] + self.nLevelCoords[n][0] - int(dtx/self.nLevelh[n][0]-1)
                id2 = self.mlen[n] + self.nLevelCoords[n][0] - int(dtx/self.nLevelh[n][0]-1) + 1
            if self.sphere[id0][2] > (self.mins[n][2]-1.0e-5):
                id4 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1] - int((dtx/self.nLevelh[n][0] -1)*(dty/self.nLevelh[n][1] -1))
                id5 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1] - int((dtx/self.nLevelh[n][0] -1)*(dty/self.nLevelh[n][1] -1)) + 1
            if self.InPlane([self.mins[n][1],self.mins[n][2]],
                            [self.maxs[n][1],self.maxs[n][2]],
                            [self.sphere[id4][1],self.sphere[id4][2]]):
                id6 = id4 + self.nLevelCoords[n][0] - int(dtx/self.nLevelh[n][0]-1)
                id7 = id4 + self.nLevelCoords[n][0] - int(dtx/self.nLevelh[n][0]-1) + 1
                
        return id0,id1,id2,id3,id4,id5,id6,id7


    
    # ----------------------------------------------------
    # Increment Axes
    # For each of the four x (az) axis used in a cell,
    # we need to know whether to increment or not
    #
    # The main loop in ComputeIndices only loops over ix
    # so for a given ix, we need to figure out whether or
    # not to:
    #   - increment an axis by one
    #   - remain constant, because we're on a hole boundary
    #   - skip to the other end of the hole
    # -----------------------------------------------------
    def IncrementAxes(self,n,ix3,ix4,ix7,line,short):
        # ixj: the count along the x (az) axis for the jth cell index
        # bis[]: bool list of whether the index lies in line that
        #        passes through the hole
        # short: number of indices within the hole
        b0x=b3x=b4x=b7x = False
        b0y=b3y=b4y=b7y = False
        b0z=b3z=b4z=b7z = False
        if n<(self.nLevel -1):
            b0x,b3x,b4x,b7x = self.ComputePlanes(self.mins[n],self.maxs[n],line,0)
            b0y,b3y,b4y,b7y = self.ComputePlanes(self.mins[n],self.maxs[n],line,1)
            b0z,b3z,b4z,b7z = self.ComputePlanes(self.mins[n],self.maxs[n],line,2)
            # Remenant of when I passed the bools as args
        b0s = [b0x,b0y,b0z]
        b3s = [b3x,b3y,b3z]
        b4s = [b4x,b4y,b4z]
        b7s = [b7x,b7y,b7z]
        i3 = i4 = i7 = 0
        if (not b0s[0]): # 'Bottom' edge of hole
            if  b0s[1] and b3s[0]:
                i3 = ix3
            else:
                i3 = ix3+1
            if  b0s[2] and b4s[0]:
                i4 = ix4
            else:
                i4 = ix4+1
            if  b3s[2] and b7s[0]:
                i7 = ix7
            else:
                i7 = ix7+1
            return i3,i4,i7
        if b0s[0]: # 'Top' Edge of hole
            if (not b3s[0]) and self.InPlane([self.mins[n][0],self.mins[n][2]],
                                             [self.maxs[n][0],self.maxs[n][2]],
                                             [self.sphere[line[2]][0],self.sphere[line[2]][2]]):
                i3 = ix3 + short
            else:
                i3 = ix3 + 1
            if (not b4s[0]) and self.InPlane([self.mins[n][0],self.mins[n][1]],
                                             [self.maxs[n][0],self.maxs[n][1]],
                                             [self.sphere[line[5]][0],self.sphere[line[5]][1]]):
                i4 = ix4 + short
            else:
                i4 = ix4 + 1
            if (not b7s[0]) and self.InPlane([self.mins[n][0],self.mins[n][1]],
                                             [self.maxs[n][0],self.maxs[n][1]],
                                             [self.sphere[line[5]][0],self.sphere[line[5]][1]]):
                i7 = ix7 + short
            else:
                i7 = ix7 + 1
            return i3,i4,i7

    # ----------------------------------------
    # SkipCell
    # Do we use this cell or skip to the next one?
    #
    # Checks this by comparing the distance between
    # two consecutive indices and ensureing it is
    # less than or equal to the cell size
    def SkipCell(self,aline,n,ix,xlen):
        skip = False
        eps = 1.5*self.nLevelh[n][0] # just make sure that it's actually bigger than the step size
        if not(n==0 and ((ix+1)%(xlen-1) == 0)):
            #print eps, np.absolute(self.sphere[aline[6]][0] - self.sphere[aline[0]][0]),np.absolute(self.sphere[aline[2]][0] - self.sphere[aline[0]][0]),self.sphere[aline[6]][2],self.sphere[aline[0]][2]
            if (np.absolute(self.sphere[aline[1]][0] - self.sphere[aline[0]][0]) > eps):
                skip = True
            if (np.absolute(self.sphere[aline[2]][0] - self.sphere[aline[3]][0]) > eps):
                skip = True
            if (np.absolute(self.sphere[aline[5]][0] - self.sphere[aline[4]][0]) > eps):
                skip = True
            if (np.absolute(self.sphere[aline[6]][0] - self.sphere[aline[7]][0]) > eps):
                skip = True
            if (np.absolute(self.sphere[aline[6]][0] - self.sphere[aline[0]][0]) > eps):
                skip = True
        return skip

    # ---------------------------------------------------
    # ComputeIndices
    # For each cell in the mesh, this function computes
    # the index of the coordinate for each vertex
    # This is then formatted into a list to be output
    # to the VTK file.
    # ----------------------------------------------------
    def ComputeIndices(self):
        ##
        ## File Level Variables
        ##
        ls = []
        nn = xycount = nprev = 0
        id0=id1=id2=id3=id4=id5=id6=id7=0
        for n in range(self.nLevel):
            # ----------------------------------------------------------------
            ##
            ## Mesh refinement level Variables
            ##
            
            # Step Sizes
            hx = self.nLevelh[n][0]
            hy = self.nLevelh[n][1]
            hz = self.nLevelh[n][2]

            # Length of short axes
            # Unfortunately these have different meanings...
            shortx = int(self.nLevelCoords[n][0])  # Length of a short x axis
            shorty = 0#int(self.nLevelCoords[n][1])# Number of short x axis 
            
            dtx = dty = dtz = 0
            # Width of hole
            if (n < self.nLevel-1):
                dtx = np.absolute(self.maxs[n][0] - self.mins[n][0])
                dty = np.absolute(self.maxs[n][1] - self.mins[n][1])
                dtz = np.absolute(self.maxs[n][2] - self.mins[n][2])
                shortx = int((self.nLevelCoords[n][0])-int(round(dtx/hx -1))) #
                shorty = int(round(dty/hy -1))
            nfullx = self.nLevelCoords[n][1] - shorty    

            # Intialise cell
            id0,id1,id2,id3,id4,id5,id6,id7 = self.InitialCell(n,dtx,dty) #Cell vertex indices
            line = np.array([id0,id1,id2,id3,id4,id5,id6,id7])
            b0x=b3x=b4x=b7x = False
            b0y=b3y=b4y=b7y = False
            b0z=b3z=b4z=b7z = False
            if(n < self.nLevel-1):
                b0x,b3x,b4x,b7x = self.ComputePlanes(self.mins[n],self.maxs[n],line,0)
                b0y,b3y,b4y,b7y = self.ComputePlanes(self.mins[n],self.maxs[n],line,1)
                b0z,b3z,b4z,b7z = self.ComputePlanes(self.mins[n],self.maxs[n],line,2)
                
            # Counting
            nprev = xycount = 0     
            # ----------------------------------------------------------------
            ##############################################################
            #                                                            #
            #         Calculate the index of each of the vertices        #
            #                                                            #
            ##############################################################
            for iz in range(self.nLevelCoords[n][2]-1):
                # ------------------------------------------------------------
                ##
                ## Plane Level Variables
                ##
                nprev += xycount
                # How many points are in the current plane?
                xycount = self.nLevelCoords[n][0]*self.nLevelCoords[n][1]
                # print self.mins[n]
                # print self.sphere[id0]
                # print self.sphere[id0+2]
                # print self.sphere[id7+2]
                # print self.maxs[n]
                if n<self.nLevel-1:
                    if self.sphere[id0+2][2]> (self.mins[n][2]-1.0e-5) and self.sphere[id0+2][2]< (self.maxs[n][2]+1.0e-5):
                        xycount = (nfullx*self.nLevelCoords[n][0]) + (shorty*shortx)
                # Counters for completed lines
                nf = nf4 = 0
                for iy in range(self.nLevelCoords[n][1] - 1):
                    # Counters for each of the 4 x-axes
                    ix3 = ix4 = ix7 = 0
                    for ix in range(self.nLevelCoords[n][0]-1):
                        if n < (self.nLevel-1): 
                            nn+=1

                            # These are the values sent to the ComputeCell function
                            # The meaning of each is documented in the ComputeCell Function
                            # Yes I know I could just pass things to the function, but this keeps
                            #   me organized.
                            k1 = k2 = k3 = k4 = k5 = k6 = k7 = 0
                            k8 = [0,0,0,0]
                            #//////////////////////////////////
                            k1 = nprev #Total up to this iz plane
                            #if not b0x:
                            # k2 - Number of points from id0 to id4
                            # Decent chance of an off by one somewhere here FIXME
                            # - How many full remaining x axis on this plane
                            # - How many short axis remaining on this plane
                            # - How many full x axes up to this iy in the above plane
                            # - How many short axes up to this iy in the above plane
                            k2 = ((nfullx - nf)*self.nLevelCoords[n][0]) +\
                                 ((shorty - (iy - nf))*shortx)+\
                                 (nf4 * self.nLevelCoords[n][0]) +\
                                 (iy-nf4)*shortx
                                
                            # ///////////////////////////////////   
                            # Number of points from id0 to id3
                            if (self.sphere[id0+1][0] < self.mins[n][0]): # Check if we're left or right of the hole
                                if b0x:
                                    k4 = shortx # What's the length of the current x axis?
                                else:
                                    k4 = self.nLevelCoords[n][0]
                            else:
                                if b3x:
                                    k4 = shortx # what's the length of the next x axis?
                                else:
                                    k4 = self.nLevelCoords[n][0]
                                    
                            # ///////////////////////////////////
                            # Number of points from id4 to id7
                            if (self.sphere[id4+1][0] < self.mins[n][0]):
                                if b4x:
                                    k7 = shortx
                                else:
                                    k7 = self.nLevelCoords[n][0]
                            else:
                                if b7x:
                                    k7 = shortx
                                else:
                                    k7 = self.nLevelCoords[n][0]
                                    
                            # ////////////////////////////////////
                            k3 = nf     # Number of completed, full x axes up to and not including this iy
                            k5 = shortx # Length of a short x axis
                            k6 = nf4    # Number of completed, full x axes in the plane above

                            # ////////////////////////////////////
                            # Lengths of each x axis
                            # Yes I know this is probably redundant, but it's more understandable in the ComputeCell Function
                            if (n==0):
                                k8 = [0,0,0,0]
                                if b0x:
                                    k8[0] = shortx
                                else:
                                    k8[0] = self.nLevelCoords[n][0]
                                if b3x:
                                    k8[1] = shortx
                                else:
                                    k8[1] = self.nLevelCoords[n][0]
                                if b4x:
                                    k8[2] = shortx
                                else:
                                    k8[2] = self.nLevelCoords[n][0]
                                if b7x:
                                    k8[3] = shortx
                                else:
                                    k8[3] = self.nLevelCoords[n][0]
 
                            # //////////////////////////////////////
                            # Actually compute the indices of each vertex
                            id0,id1,id2,id3,id4,id5,id6,id7 = self.ComputeCell(n,ix,ix3,ix4,ix7,iy,iz,
                                                                               k1,k2,k3,k4,k5,k6,k7,k8)
                            line = np.array([id0,id1,id2,id3,id4,id5,id6,id7])


                            # Check positions relative to the hole
                            b0x,b3x,b4x,b7x = self.ComputePlanes(self.mins[n],self.maxs[n],line,0)
                            b0y,b3y,b4y,b7y = self.ComputePlanes(self.mins[n],self.maxs[n],line,1)
                            b0z,b3z,b4z,b7z = self.ComputePlanes(self.mins[n],self.maxs[n],line,2)
                            b1y = self.InPlane([self.mins[n][0],self.mins[n][2]],
                                               [self.maxs[n][0],self.maxs[n][2]],
                                               [self.sphere[line[1]][0],self.sphere[line[1]][2]])
                            # Check if we're skipping this one
                            if self.SkipCell(line,n,ix,k8[0]):
                                ix3,ix4,ix7 = self.IncrementAxes(n,ix3,ix4,ix7,
                                                                 line,
                                                                 self.nLevelCoords[n][0] - shortx - 1)
                                continue
                            # Append
                            ls.append(line)
                            self.ncell+=1

                            # ///////////////////////////////////////
                            # Count number of completed axes
                            if (ix+1)%(self.nLevelCoords[n][0]-1)==0:
                                nf+=1
                            if (ix4+1)%(self.nLevelCoords[n][0]-1)==0:
                                nf4+=1

                            # Check end of line
                            if b0x and ((ix+1)%shortx==0):
                                nn += self.nLevelCoords[n][0] - shortx -1
                                break
                            # Increment counters
                            ix3,ix4,ix7 = self.IncrementAxes(n,ix3,ix4,ix7,
                                                             line,
                                                             self.nLevelCoords[n][0] - shortx - 1)
                            # We made it!
                        else:
                            # There are simpler ways to do this,
                            # but wanted to ensure the full way works properly
                            nn += 1
                            k1 = nprev
                            k2 = ((nfullx - nf)*self.nLevelCoords[n][0]) +\
                                 ((shorty - (iy - nf))*shortx) +\
                                 (nf4 * self.nLevelCoords[n][0]) +\
                                 (iy-nf4)*shortx
                            k3 = nf
                            k4 = self.nLevelCoords[n][0]
                            k5 = shortx
                            k6 = nf4
                            k7 = self.nLevelCoords[n][0]
                            k8 = [self.nLevelCoords[n][0],self.nLevelCoords[n][0],self.nLevelCoords[n][0],self.nLevelCoords[n][0]]
                            #id0,id1,id2,id3,id4,id5,id6,id7 = self.ComputeStructuredCell(n,ix,iy,iz)
                            id0,id1,id2,id3,id4,id5,id6,id7 = self.ComputeCell(n,ix,ix3,ix4,ix7,iy,iz,
                                                                               k1,k2,k3,k4,k5,k6,k7,k8)
                            line = np.array([id0,id1,id2,id3,id4,id5,id6,id7])
                            if self.SkipCell(line,n,ix,k8[0]):
                               continue
                            ix3,ix4,ix7 = self.IncrementAxes(n,ix3,ix4,ix7,
                                                             line,
                                                             self.nLevelCoords[n][0] - shortx - 1)
                            if (ix+1)%(self.nLevelCoords[n][0] -1) == 0:
                                nf += 1
                                nf4 += 1
                            ls.append(line)
                            self.ncell+=1
        # Output to terminal, return the list of vertex indices
        print str(self.ncell) + " Cells in the mesh."
        print nn
        #print (self.nLevelCoords[1][0]-1)*(self.nLevelCoords[1][1]-1)*(self.nLevelCoords[1][2]-1),(self.nLevelCoords[0][0]-1)*(self.nLevelCoords[0][1]-1)*(self.nLevelCoords[0][2]-1)
        return ls


    # ------------------------------------------
    # ComputeStructuredIndices
    #
    # Reimplements ComputeIndices with simplifying
    # assumptions, so that we can easily compute
    # the indices for the unfiltered grid
    # ------------------------------------------
    def ComputeStructuredIndices(self):
        nn = 0
        c = 0
        ls = []
        for n in range(self.nLevel):
            totalus = nextus = 0 #k1,k2 for unfiltered grid
            for iz in range(self.nLevelCoords[n][2] -1):
                totalus+=nextus
                thisus = nextus
                nextus = int((self.nLevelCoords[n][0])*(self.nLevelCoords[n][1]))
                for iy in range(self.nLevelCoords[n][1] -1):
                    for ix in range(self.nLevelCoords[n][0] -1):
                            id0,id1,id2,id3,id4,id5,id6,id7 = self.ComputeStructuredCell(n,ix,iy,iz,False)
                            # Look at the unfiltered grid to build a list of cells to remove from the input data
                            line = np.array([id0,id1,id2,id3,id4,id5,id6,id7])
                            #ls.append(line)
                            cell = np.array([(self.unfiltered[id0][0] + self.unfiltered[id6][0])/2.,
                                             (self.unfiltered[id0][1] + self.unfiltered[id6][1])/2.,
                                             (self.unfiltered[id0][2] + self.unfiltered[id6][2])/2.])
                            if n < (self.nLevel-1):
                                if self.InRange(self.mins[n],self.maxs[n],cell):
                                    # Add to list of data points to remove
                                    c+=1
                                    self.cellist.append(nn)
                            nn+=1
        print str(c) + " Cells to remove"
        #self.ncell = nn - c
        return 


    #
    #  Format data into VTK structure
    #
    def VTKFormatOut(self, data, inds, binary = True, append = False):
        print "Writing to file..."
        # Write out based on RadMC3D WriteVTK()
        # and makes use of pyvtk
        #
        ncoords = self.mesh.shape[0]
        ncells = self.ncell
        if not append:
            if binary:
                if "velocity" not in self.feature:
                    vtk = VtkData(UnstructuredGrid(self.mesh,hexahedron = inds,),
                                  CellData(Scalars(data,self.feature)),
                                  name ="JUPITER Sim"+str(self.outNumber)+" "+self.feature+" field")
                else:
                    vtk = VtkData(UnstructuredGrid(self.mesh,hexahedron = inds,),
                                  CellData(Vectors(data,self.feature)),
                                  name ="JUPITER Sim"+str(self.outNumber)+" "+self.feature+" field")
                vtk.tofile(self.dataOutPath + self.outFilename, 'binary')
            else:
                if "velocity" not in self.feature:
                    vtk = VtkData(UnstructuredGrid(self.mesh,hexahedron = inds,),
                                  CellData(Scalars(data,self.feature)),
                                  name ="JUPITER Sim"+str(self.outNumber)+" "+self.feature+" field")

                else:
                    vtk = VtkData(UnstructuredGrid(self.mesh,hexahedron = inds,),
                                  CellData(Vectors(data,self.feature)),
                                  name ="JUPITER Sim"+str(self.outNumber)+" "+self.feature+" field")

                vtk.tofile(self.dataOutPath + self.outFilename)

        else:
            # -----------------------------------------------
            # Ascii out, append to file
            # -----------------------------------------------
            outfile = open(self.dataOutPath + self.outFilename, 'a+')
            if not binary:
                if "velocity" in self.feature:
                    print "Writing Vector Data"
                    outfile.write('\n%s %d\n'%('POINT_DATA', self.ncell))
                    outfile.write('%s\n'%('VECTORS ' + self.feature +' double'))
                    for i in range(self.ncell):
                        outfile.write('%.9e %.9e %.9e\n'%(data[i][0], data[i][1], data[0][2]))

                #
                # Write out the CELL CENTERED scalars
                #
                if "velocity" not in self.feature:
                    print "Writing scaler " + self.feature + " data..."
                    outfile.write('\n%s %d\n'%('CELL_DATA', self.ncell))
                    outfile.write('%s\n'%('SCALARS ' + self.feature + ' double'))
                    outfile.write('%s\n'%'LOOKUP_TABLE default')
                    for i in range(self.ncell):
                        self.printProgressBar(i,self.ncell)
                        outfile.write('%.9e\n'%data[i])

                outfile.close()
                return
            else:
                # ---------------------------------------------------
                # Binary out, append to file
                # ---------------------------------------------------
                print "This probably won't work. Recommend just overwriting or using ascii files."
                outfile = open(self.dataOutPath + self.outFilename, 'a+')
                if "velocity" in self.feature:
                    print "Writing Vector Data"
                    outfile.write('\n%s %d\n'%('CELL_DATA', self.ncell))
                    outfile.write('%s\n'%('VECTORS ' + self.feature +' double'))
                    outfile.close()
                    outfile = open(self.dataOutPath + self.outFilename, 'ab+')
                    for i in range(self.ncell):
                        data[i].tofile(outfile)

                #
                # Write out the CELL CENTERED scalars
                #
                if "velocity" not in self.feature:
                    print "Writing scaler " + self.feature + " data..."
                    outfile.close()
                    outfile = open(self.dataOutPath + self.outFilename, 'a+')
                    outfile.write('\n%s %d\n'%('CELL_DATA', self.ncell))
                    outfile.write('%s\n'%('SCALARS ' + self.feature + ' double'))
                    outfile.write('%s\n'%'LOOKUP_TABLE default')
                    outfile.close()
                    outfile = open(self.dataOutPath + self.outFilename, 'ab+')
                    for i in range(ncells):
                        self.printProgressBar(int(i),self.ncell)
                        data[i].tofile(outfile)

                # Close and exit
                outfile.close()
                return

    # --------------------------------------------------------------------------------------------------
    # Prettier outputs

    #
    # So I can use a nice dialog box
    #
    def ask(title, text, strings=('Yes', 'No'), bitmap='questhead', default=0):
        d = Dialog.Dialog(
            title=title, text=text, bitmap=bitmap, default=default, strings=strings)
        s = strings[d.num]
        d.destroy()
        return s

    #
    # This was nice before I used pyvtk...
    #
    def printProgressBar(self, iteration, total, prefix = '', suffix = '', decimals = 0, length = 67, fill = u'\u2588'):
        """
        Call in a loop to create terminal progress bar
        @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        """
        str_format = "{0:." + str(decimals) + "f}"
        percents = str_format.format(100 * (iteration / float(total-1)))
        filled_length = int(round(length * iteration / float(total-1)))
        bar = fill * filled_length + '-' * (length - filled_length)
        sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix))
        # Print New Line on Complete
        if iteration == total-1:
            print('\n')

    def ListLen(self, alist):
        t = 0
        for l in alist:
                t += l.size
        return

    # So this was an attempt to delete the cells within the
    # bounds of the next mesh layer. It's at least O(nlogn),
    # probably O(n^2) and is too slow to be useful.
    def FilterCells(self,inds,dels):
        print "Filtering Cells"
        delcount = 0
        curcount = 0
        out = []
        while curcount < len(inds):
            if(delcount<len(dels)):
                if(curcount == dels[delcount]):
                    inds.pop(curcount)
                    delcount +=1
                    continue
            for ind in range(len(inds[curcount])):
                d = self.find_nearest(dels,inds[curcount][ind])
                inds[curcount][ind] = inds[curcount][ind] - d
            curcount +=1
        return inds
    def FilterData(self,data,dels):
        print "Filtering Data"
        data = np.delete(data,dels)
        return data
    def find_nearest(self,array,value):
        idx = np.searchsorted(array, value, side="left")
        if idx > 0 and (idx == len(array) or (array[idx] > value)):
            return idx-1
        else:
            return idx
