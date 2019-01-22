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


import os,sys
import numpy as np
import astropy.units as u
import astropy.constants as c
import string
import Dialog

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
    """
    Convert JUPITER .DAT files to binary VTK files for Paraview
    
    Reads in a Descriptor.dat file and a hydrodynamic field file (.dat)
    Outputs a vtk formatted file using an unstructured grid
    """
    def __init__(self):
        # Simulation Information
        self.simNumber = -1 # Which sim output number are we looking at
        self.nLevel = -1 # How many mesh levels are there?
        self.feature = 'notafeat' # what hydro field are we looking at
        self.nLevelCoords = [] # How many coordinates along each axis at this mesh refinement level?
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
        self.rcgs = -1.*u.g
        self.mass = -1.
        self.mcgs = -1.*u.g
        # Science Constants
        self.TEMP = -1.
        self.DENS = -1.
        self.PERIOD = -1.
        self.VEL = -1.
        self.units = "CGS"
        # Grid Information
        self.mesh = np.zeros((0,0,3),dtype=np.float64)   # 3D cartesian array - edges
        self.unfiltered = self.mesh = np.zeros((0,0,3),dtype=np.float64) # Unfiltered spherical coords
        self.ncell = 0 # How many cells are there
        self.mlenUF = [0] # Number of points in unfiltered in each mesh refinement level
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
            print("Not a recognised input, may not be implemented. Continuing...")
            self.feature = feat
            return 1
    def SetBasePath(self,path):
        self.BASEPATH = path
    def SetInDir(self,path):
        self.dataDir = path
    def SetOutDir(self,path):
        self.dataOutPath = path
    def SetInFilename(self,fname):
        self.inFilename = fname

    def SetUnits(self, units):
        ulist = ["CGS","cgs","AU","Au","au"]
        if units in ulist:
            self.units = units
        else:
            print "Please select a valid unit system."
        return
    # Science Set Functions
    def SetRadius(self, rad):
        """
        Orbital radius of companion
        """
        self.radius = rad*u.AU
        self.rcgs = self.radius.to(u.cm)
        if(self.mass>0. and self.TEMP <0):
            self.SetConstants(self.units)

    def SetMass(self,mass):
        """
        Mass of star in solar masses
        """
        self.mass = mass*u.M_sun
        self.mcgs = self.mass.to(u.g)
        if(self.radius > 0. and self.TEMP < 0):
            self.SetConstants(self.units)

    def SetConstants(self,units):
        # 6.67e-8 - cgs G
        # 8.314e7 - cgs gas constant
        if units == "AU" or units == "Au" or units == "au":
            self.TEMP = ((self.rcgs.value)/(np.sqrt((self.rcgs.value)**3 / 6.67259e-8 / (self.mcgs.value))))**2 / 8.314e7
            self.DENS = self.mcgs.value/(self.rcgs.value)**3
            self.PERIOD = 2*np.pi*np.sqrt((self.rcgs.value)**3 / (6.67259e-8 * self.mcgs.value))
            self.VEL = self.rcgs.value/(self.PERIOD/2*np.pi)
            self.rcgs = self.radius
            self.mcgs = self.mass
        else:    
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
        
        try:
            assert(os.path.isfile(self.dataDir + self.inFilename))
        except AssertionError:
            print self.dataDir + self.inFilename + " Does not exist. Please enter a valid filename or filepath."
            sys.exit(1)


    # ---------------------------------------------------------------------------------------------
    # Important part starts here

    #
    # This function wraps the binary .dat file reader for a given feature,
    # and output for the field in a VTK format. This is the only user facing function.
    # --------------------------------------------------------------------------------------------
    def ConvertFiles( self , binary = True , planet_centered = False):
        self.GetCoordinates()
        inds = self.ComputeIndices()
        data = np.array([]) 
        for i in range(self.nLevel):
            self.SetupNames(i)
            feat = np.fromfile(self.dataDir + self.inFilename, dtype = 'double')  
            if "velocity" in self.feature:
                data2 = np.array([],dtype = np.float64)
                data2 = np.append(data2,feat.astype(np.float64))
                data2 = np.reshape(data2,(3,-1))
                data3 = np.column_stack((data2[0], data2[1], data2[2]))
                data3 = data3
                if i>0:
                    data = np.concatenate((data,data3), axis = 0)
                else:
                    data = data3
            else:
                # Read in binary doubles into a 1D array
                data = np.append(data,feat)
                
        # Delete overlapping data points       
        data = self.FilterData(data,self.cellist)
        if "velocity" in self.feature:
            if planet_centered:
                data = self.PlanetCenteredVel(data,inds)
                data = data*self.VEL
            else:
                data = self.StarCenteredVel(data,inds)
                data = data*self.VEL
        # Convert to CGS units
        
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
        newcoords[:,0] = data[:,1]*self.rcgs.value*np.sin(data[:,2])*np.cos(data[:,0])
        newcoords[:,1] = data[:,1]*self.rcgs.value*np.sin(data[:,2])*np.sin(data[:,0])
        newcoords[:,2] = data[:,1]*self.rcgs.value*np.cos(data[:,2])
        return newcoords   
    def CartToSphere(self, data):
        newcoords = np.zeros(data.shape)
        newcoords[:,0] = np.arctan2(data[:,1],data[:,0])
        newcoords[:,1] = np.sqrt(data[:,0]**2 + data[:,1]**2 + data[:,2]**2)
        newcoords[:,2] = np.arccos(data[:,2],newcoords[:,1])
        return newcoords

    
    def PlanetCenteredVel(self, data, inds):
        """ 
        PlanetCenteredVelocities

        data: the 3D velocity data read in from .dat file

        inds: list of index of each vertex of each mesh cell

        PlanetCenteredVelocities converts the spherical velocities
        read in from file to a cartesian coordinate system.
        The orbital velocity of the planet is already subtracted from 
        the data read in from file.

        Note that the component ordering when reading in from file may change
        - this is a user adjustable parameter in JUPITER
        """
        
        newcoords = np.zeros(data.shape)
        i = 0
        for ind in inds: 
            x = (self.mesh[ind[0]][0] + self.mesh[ind[0]+1][0])/2.0 # x
            y = (self.mesh[ind[0]][1] + self.mesh[ind[3]][1])/2.0   # y
            z = (self.mesh[ind[0]][2] + self.mesh[ind[4]][2])/2.0   # z

            # Shift to planet centered reference frame
            x = x - self.rcgs.value
            sph = self.CartToSphere(np.array([[x,y,z]])/self.rcgs.value)

            # So I don't have to change any of the rest of the math
            phi = sph[0][0]
            rad = sph[0][1]
            rad = rad*u.AU.to(u.cm)
            tht = sph[0][2]
            
            # 0 = az (phi), 1 = rad, 2 = pol (tht)  Use this one, velocity ordering is(n't) weird
            # 0 = pol(tht), 1 = az (phi) 2 = rad

            # Using the coordinate conversion from RADMC3D - note no sin/cos theta dependance for x,y vel
            ydot = np.sin(phi)*np.sin(tht)*data[i][1]+\
                   rad*np.cos(phi)*data[i][0] +\
                   rad*np.cos(tht)*data[i][2]
            xdot = np.cos(phi)*np.sin(tht)*data[i][1] -\
                   rad*np.sin(phi)*data[i][0] +\
                   rad*np.cos(tht)*data[i][2]
            zdot = np.cos(tht)*data[i][1] -\
                   rad*np.sin(tht)*data[i][2]
            newcoords[i][0] = xdot
            newcoords[i][1] = ydot
            newcoords[i][2] = zdot

            i+=1
        return newcoords
    def StarCenteredVel(self, data, inds):
        '''
        Assuming circular keplerian orbital velocity for the planet
        '''
        newcoords = np.zeros(data.shape)
        i = 0
        for ind in inds:
            phi = (self.unfiltered[ind[0]][0] + self.unfiltered[ind[0]+1][0])/2.0 # Azimuth
            rad = (self.unfiltered[ind[0]][1] + self.unfiltered[ind[3]][1])/2.0   # Radial
            rad = rad*u.AU.to(u.cm)
            tht = (self.unfiltered[ind[0]][2] + self.unfiltered[ind[4]][2])/2.0   # Polar

            # 0 = az (phi), 1 = rad, 2 = pol (tht)  Use this one, velocity ordering is(n't) weird
            # 0 = pol(tht), 1 = az (phi) 2 = rad
            ydot = np.sin(phi)*np.sin(tht)*data[i][1]+\
                   rad*np.cos(phi)*data[i][0] +\
                   rad*np.cos(tht)*data[i][2]
            xdot = np.cos(phi)*np.sin(tht)*data[i][1] -\
                   rad*np.sin(phi)*data[i][0] +\
                   rad*np.cos(tht)*data[i][2]
            zdot = np.cos(tht)*data[i][1] -\
                   rad*np.sin(tht)*data[i][2]
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
                            # We made it!                      
                        else:
                            nn += 1
                            id0,id1,id2,id3,id4,id5,id6,id7 = self.ComputeStructuredCell(n,ix,iy,iz,False)  
                            line = np.array([id0,id1,id2,id3,id4,id5,id6,id7])
                            ls.append(line)
                            self.ncell+=1  
        # Output to terminal, return the list of vertex indices
        print str(self.ncell) + " Cells in the mesh."
        return ls


    def VTKFormatOut(self, data, inds, binary = True, append = False):
        print "Writing to file..."
        """
        Write out based on RadMC3D WriteVTK() function
        and makes use of vtk package

        data: the scaler or vector hydrodynamic data to be written to file
              so far velocity is the only vector quantity allowed

        inds: list of the indices of the vertices of each mesh cell
        
        binary: should the output file be binary (T) or ascii (F)

        append: should the data be appended to the file or overwrite existing?
        """
        ncoords = self.mesh.shape[0]
        ncells = self.ncell
        if not append:
            if binary:
                if "velocity" not in self.feature:
                    vtk = VtkData(UnstructuredGrid(self.mesh,hexahedron = inds,),
                                  CellData(Scalars(data,self.feature)),
                                  name ="JUPITER Sim "+str(self.outNumber)+" "+self.feature+" field")
                else:
                    vtk = VtkData(UnstructuredGrid(self.mesh,hexahedron = inds,),
                                  CellData(Vectors(data,self.feature)),
                                  name ="JUPITER Sim "+str(self.outNumber)+" "+self.feature+" field")
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
                        outfile.write('%.9e\n'%data[i])

                outfile.close()
                return
            else:
                # ---------------------------------------------------
                # Binary out, append to file
                #
                # FIXME - this probably doesn't work yet.
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

    def FilterData(self,data,dels):
        print "Filtering Data"
        if "velocity" in self.feature:
            data = np.delete(data,dels,axis = 0)
        else:
            data =  np.delete(data,dels)
        return data
    
