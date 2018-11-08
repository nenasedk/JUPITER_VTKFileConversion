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
        inds = self.ComputeIndices()
        data = np.array([]) 
        for i in range(self.nLevel):
            self.SetupNames(i)
            feat = np.fromfile(self.dataDir + self.inFilename, dtype = 'double')  
            if "velocity" in self.feature:
                data2 = np.zeros(0,dtype = np.float64)
                data2 = np.append(data2,feat.astype(np.float64))
                data2 = np.reshape(data2,(3,-1))
                data3 = np.column_stack((data2[1], data2[0], data2[1]))
                data3 = data3*self.VEL#.value
                if i > 0:
                    data = np.concatenate((data,data3), axis = 0)
                else:
                    data = data3

            else:
                # Read in binary doubles into a 1D array
                #feat = np.fromfile(self.dataDir + self.inFilename, dtype = 'double')
                data = np.append(data,feat)
        # Compute all of the indices of cell vertices
        # Have to do this here to find out which cells
        # are overlapping, and should not be included
        # Delete overlapping data points       
        data = self.FilterData(data,self.cellist)
        if "velocity" in self.feature:
            data = self.SphereToCartDT(data,inds)
        # Convert to CGS units 
        print data.shape
        if("density" in self.feature):
            data = np.array([x*self.DENS for x in data])#.value
        if("temperature" in self.feature):
            data = np.array([x*self.TEMP for x in data])#.value
        print str(data.shape[0]) + " data points for " + self.feature
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
        #self.sphere,self.coordlist = self.FilterCoords(phi,r,th)
        self.mesh = self.SphereToCart(self.unfiltered)
        print(str(self.unfiltered.shape[0]) + " Vertices in filtered grid.")


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
            self.mins.append([np.min(x1s[i]),np.min(x2s[i]),np.min(x3s[i])])
            self.maxs.append([np.max(x1s[i]),np.max(x2s[i]),np.max(x3s[i])])
            tg,c = self.BuildOneLevel(x1s[i], x2s[i], x3s[i])
            grid.extend(tg)
            self.mlenUF.append(self.mlenUF[i] + tg.shape[0])
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

    def SphereToCartDT(self, data, inds):
        newcoords = np.zeros(data.shape)
        i = 0
        for ind in inds:
            phi = (self.unfiltered[ind[0]][0] + self.unfiltered[ind[6]][0])/2.0
            rad = (self.unfiltered[ind[0]][1] + self.unfiltered[ind[6]][1])/2.0 #* self.rcgs
            tht = (self.unfiltered[ind[0]][2] + self.unfiltered[ind[6]][2])/2.0

            xdot = np.cos(tht)*np.cos(phi)*data[i][2] +\
                   rad*np.cos(tht)*np.sin(phi)*data[i][1] +\
                   rad*np.sin(tht)*np.cos(phi)*data[i][0]
            ydot = np.cos(tht)*np.sin(phi)*data[i][1] -\
                   rad*np.cos(tht)*np.cos(phi)*data[i][1] +\
                   rad*np.sin(tht)*np.cos(phi)*data[i][0]
            zdot = np.sin(tht)*data[i][2] -\
                   rad*np.cos(tht)*data[i][0]

            newcoords[i][0] = xdot
            newcoords[i][1] = ydot
            newcoords[i][2] = zdot
            i+=1
        return newcoords
    # --------------------------------------------------
    # WriteToVTK
    # Checks data formatting, prepares output file
    # --------------------------------------------------
    def WriteToVTK(self, data, inds, binary = True, append = False):
        # Data quality checking
        try:
            assert(data.shape[0] == self.ncell)
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
            

    # -----------------------------------------------
    # ComputeStructuredCell
    # Compute the indices of the vertices of a given
    # cell in the mesh, assuming a completed,
    # unfiltered grid
    #
    # Implements basic stride counting
    # -----------------------------------------------
    def ComputeStructuredCell(self,n,ix,iy,iz,filtered = True):
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
        nn = 0
        id0=id1=id2=id3=id4=id5=id6=id7=0
        for n in range(self.nLevel):
            ##############################################################
            #                                                            #
            #         Calculate the index of each of the vertices        #
            #                                                            #
            ##############################################################
            for iz in range(self.nLevelCoords[n][2]-1):
                for iy in range(self.nLevelCoords[n][1] - 1):
                    for ix in range(self.nLevelCoords[n][0] - 1):
                        if n < (self.nLevel-1): 
                            id0,id1,id2,id3,id4,id5,id6,id7 = self.ComputeStructuredCell(n,ix,iy,iz,False)
                            line = np.array([id0,id1,id2,id3,id4,id5,id6,id7])
                            skip = False
                            for ind in line:
                                if self.InRange(self.mins[n+1],self.maxs[n+1],self.unfiltered[ind]):
                                    skip = True
                            if skip:
                                self.cellist.append(nn)
                                nn+=1
                                continue
                            ls.append(line)
                            self.ncell+=1
                            nn+=1

                            # Check end of line
                            #if b0x and ((ix+1)%shortx==0):
                            #    nn += self.nLevelCoords[n][0] - shortx -1
                            #    break
                            # We made it!                      
                        else:
                            nn += 1
                            id0,id1,id2,id3,id4,id5,id6,id7 = self.ComputeStructuredCell(n,ix,iy,iz,False)  
                            line = np.array([id0,id1,id2,id3,id4,id5,id6,id7])
                            ls.append(line)
                            self.ncell+=1
                #sys.exit(1)
                
        # Output to terminal, return the list of vertex indices
        print str(self.ncell) + " Cells in the mesh."
        return ls

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

    def FilterData(self,data,dels):
        print "Filtering Data"
        data = np.delete(data,dels,axis = 0)
        return data
