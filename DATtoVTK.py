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


import os,sys
import numpy as np
import astropy.units as u
import string
import Dialog
try:
    from pyvtk import *
except ImportError:
    print "Please install pyvtk. (pip install pyvtk)"
import struct

# Code Unit Constants, defined after mass and radius are given

class DATtoVTK:
    'Convert JUPITER .DAT files to binary VTK files for Paraview'
    def __init__(self):
        # Simulation Information
        self.simNumber = -1 # Which sim output number are we looking at
        self.nLevel = -1 # How many mesh levels are there?
        self.feature = 'notafeat'
        self.nLevelCoords = []
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
        self.sphere = np.zeros((0,0,0),dtype=np.float64) # 3D Spherical array - edges (read in)
        self.mesh = np.zeros((0,0,0),dtype=np.float64)   # 3D cartesian array - edges
        self.cent = np.zeros((0,0,0),dtype=np.float64)   # 3D cartesian array of cell centers
        self.nx = 0
        self.ny = 0
        self.nz = 0
        self.ncell = 0
        self.mlen = [0]

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
            self.TEMP = ((self.rcgs)/(np.sqrt((self.rcgs)**3 / 6.67259e-8 / (self.mcgs))))**2 / 8.314e7
            self.DENS = self.mcgs/(self.rcgs)**3
            self.PERIOD = 2*np.pi*np.sqrt((self.rcgs)**3 / (6.67259e-8 * self.mcgs))
            self.VEL = self.rcgs/(self.PERIOD/2*np.pi)
            
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
    #
    def ConvertFiles( self , binary = True ):
        data = np.zeros(0,dtype = np.float)
        for i in range(self.nLevel):
            self.SetupNames(i)
            feat = np.fromfile(self.dataDir + self.inFilename, dtype = 'double')
            
            if "velocity" in self.feature:
                data2 = np.zeros(0,dtype = np.float)
                data2 = np.append(data2,feat.astype(float))
                data2 = data2*self.VEL.value
                data2 = np.reshape(data2,(3,-1))
                data3 = np.zeros((3,len(data2[0])))
                data3[:,0] = data2[:,1]
                data3[:,1] = data2[:,2]
                data3[:,2] = data2[:,0]
                data3 = data3.transpose() # Velocity ordering is weird.
                if i > 0:
                    data = np.concatenate((data,data3), axis = 0)
                else:
                    data = data3
            else:
                data = np.append(data,feat.astype(float))
        # Convert to CGS units
        if("density" in self.feature):
            data = data*self.DENS.value
        if("temperature" in self.feature):
            data = data*self.TEMP.value

        print str(len(data)) + " data points for " + self.feature
        self.WriteToVTK(data, binary)
        
    #
    # This function reads in the descriptor file and builds the VTK coordinate arrays
    #
    def GetCoordinates(self):
        
        phi = []
        r = []
        th = []
        coords = []
        for i in range(0,self.nLevel):
            dsc = open(self.dataDir + self.descriptorName)
            for j, line in enumerate(dsc):
                if j == 8 + (i*11):
                    cur = [float(x) for x in line.split()]
                    cur.pop(0) # First and last two points are 'ghost points'
                    cur.pop(0)
                    cur.pop()
                    cur.pop() 
                    phi.extend(cur)
                    cur = []
                if j == 9 + (i*11):    
                    cur = [float(x) for x in line.split()]
                    cur.pop(0)
                    cur.pop(0)
                    cur.pop()
                    cur.pop()
                    r.extend(cur)
                    cur = []
                if j == 10 + (i*11):    
                    cur = [float(x) for x in line.split()]
                    cur.pop(0)
                    cur.pop(0)
                    cur.pop()
                    cur.pop()
                    th.extend(cur)
                    cur = []
            
            for ath in th:
                for ar in r:
                    for aphi in phi:
                        coords.append([-1.*aphi,ar,-1.*ath]) # Invert theta, phi so that coordinate system is right handed (needed for cell orientation)
            dsc.close()
            self.nLevelCoords.append([len(phi),len(r),len(th)])
            self.mlen.append(len(coords)) 
            phi = []
            r = []
            th = []

        # Convert from spherical to cartesian coordinates
        coords = np.array(coords)
        newcoords = np.hstack((coords,np.zeros(coords.shape)))
        newcoords[:,3] = coords[:,1]*self.rcgs*np.sin(coords[:,2])*np.cos(coords[:,0])
        newcoords[:,4] = coords[:,1]*self.rcgs*np.sin(coords[:,2])*np.sin(coords[:,0])
        newcoords[:,5] = coords[:,1]*self.rcgs*np.cos(coords[:,2])
        self.sphere = coords
        self.mesh = newcoords[:,[3,4,5]]
        for n in range(self.nLevel):
            self.nx += self.nLevelCoords[n][0]
            self.ny += self.nLevelCoords[n][1]
            self.nz += self.nLevelCoords[n][2]
            self.ncell += ((self.nLevelCoords[n][0]-1)*(self.nLevelCoords[n][1]-1)*(self.nLevelCoords[n][2]-1))
        #self.nx = len(self.x)
        print(str(len(self.mesh)) + " Vertices in grid.")
        print str(self.ncell) + " Cells in grid."

        #
        # Checks data formatting, prepares output file
        #
    def WriteToVTK(self, data, binary = True, append = False):
        # Data quality checking
        try:
            assert(len(data) == self.ncell) 
        except ValueError:
            print "Error: Number of data points does not match number of mesh elements"
            return
        # File existance check
        if not os.path.isfile(self.dataOutPath + self.outFilename):
            self.VTKFormatOut(data, binary)
        else:
            userstr = self.ask("Output file already exists, do you want to overwrite or append data?",("Overwrite","Append","Quit"))
            if userstr is 'Overwrite':
                self.VTKFormatOut(data,binary, False)
                return
                   
            if userstr is 'Append': # not sure if append will work yet
                self.VTKFormatOut(data, binary, True)
                return
            else:
                print "Nothing written to file!"
                return

    #                                  
    #  Format data into VTK structure  
    #                                        
    def VTKFormatOut(self, data, binary = True, append = False):
        # Write out based on RadMC3D WriteVTK()
        # and makes use of pyvtk
        #
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
        # ---------------------------------------------------------------------------------------
        ncoords = self.mesh.shape[0]
        ncells = self.ncell
        ls = []
        
        if not append:
            if binary:
                nn = nt= 0
                for n in range(self.nLevel): # mesh refinement levels are written out sequentially
                    for iz in range(self.nLevelCoords[n][2]-1):
                        for iy in range(self.nLevelCoords[n][1]-1):
                            for ix in range(self.nLevelCoords[n][0]-1):   # Calculate the index of each of the vertices of a given cell, and write to file.
                                nn +=1
                                id0 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*iy     + ix
                                if(n == 0):
                                    id1 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*iy     + ((ix+1) % (self.nLevelCoords[n][0]-1))
                                    id2 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*(iy+1) + ((ix+1) % (self.nLevelCoords[n][0]-1))
                                    id5 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*iy     + ((ix+1) % (self.nLevelCoords[n][0]-1))
                                    id6 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*(iy+1) + ((ix+1) % (self.nLevelCoords[n][0]-1))
                                else:
                                    id1 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*iy     + (ix+1)
                                    id2 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*(iy+1) + (ix+1)
                                    id5 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*iy     + (ix+1)
                                    id6 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*(iy+1)  + (ix+1)
                                id3 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*(iy+1) + ix 
                                id4 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*iy     + ix
                                id7 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*(iy+1) + ix
                                
                                line = np.array([id0,id1,id2,id3,id4,id5,id6,id7])
                                ls.append(line)

                    nt  = nn
                                
                if "velocity" not in self.feature:
                    vtk = VtkData(UnstructuredGrid(self.mesh,hexahedron = ls,),CellData(Scalars(data,self.feature)),name ="JUPITER Sim" + str(self.outNumber) + " " + self.feature + " field")
                else:
                    vtk = VtkData(UnstructuredGrid(self.mesh,hexahedron = ls,),CellData(Vectors(data,self.feature)),name="JUPITER Sim" + str(self.outNumber) + " " + self.feature + " field")
                vtk.tofile(self.dataOutPath + self.outFilename, 'binary')
            else:
                nn = nt= 0
                for n in range(self.nLevel): # mesh refinement levels are written out sequentially
                    for iz in range(self.nLevelCoords[n][2]-1):
                        for iy in range(self.nLevelCoords[n][1]-1):
                            for ix in range(self.nLevelCoords[n][0]-1):   # Calculate the index of each of the vertices of a given cell, and write to file.
                                nn +=1
                                id0 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*iy     + ix
                                if(n == 0):
                                    id1 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*iy     + ((ix+1) % (self.nLevelCoords[n][0]-1))
                                    id2 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*(iy+1) + ((ix+1) % (self.nLevelCoords[n][0]-1))
                                    id5 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*iy     + ((ix+1) % (self.nLevelCoords[n][0]-1))
                                    id6 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*(iy+1) + ((ix+1) % (self.nLevelCoords[n][0]-1))
                                else:
                                    id1 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*iy     + (ix+1)
                                    id2 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*(iy+1) + (ix+1)
                                    id5 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*iy     + (ix+1)
                                    id6 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*(iy+1)  + (ix+1)
                                id3 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*(iy+1) + ix 
                                id4 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*iy     + ix
                                id7 = self.mlen[n] + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*(iy+1) + ix
                                
                                line = np.array([id0,id1,id2,id3,id4,id5,id6,id7])
                                ls.append(line)

                    nt  = nn
                                
                if "velocity" not in self.feature:
                    vtk = VtkData(UnstructuredGrid(self.mesh,hexahedron = ls,),CellData(Scalars(data,self.feature)),name ="JUPITER Sim" + str(self.outNumber) + " " + self.feature + " field")
                else:
                    vtk = VtkData(UnstructuredGrid(self.mesh,hexahedron = ls,),CellData(Vectors(data,self.feature)),name="JUPITER Sim" + str(self.outNumber) + " " + self.feature + " field")
                vtk.tofile(self.dataOutPath + self.outFilename)
          
        else:
            # -----------------------------------------------
            # Ascii out, append to file
            # -----------------------------------------------
            outfile = open(self.dataOutPath + self.outFilename, 'a+')
            if not binary:
                if "velocity" in self.feature:
                    print "Writing Vector Data"
                    outfile.write('\n%s %d\n'%('POINT_DATA', ncoords))
                    outfile.write('%s\n'%('VECTORS ' + self.feature +' float'))
                    for i in range(0,ncoords):
                        outfile.write('%.9e %.9e %.9e\n'%(data[i][0], data[i][1], data[0][2]))

                # 
                # Write out the CELL CENTERED scalars
                #
                if "velocity" not in self.feature:
                    print "Writing scaler " + self.feature + " data..."
                    outfile.write('\n%s %d\n'%('CELL_DATA', ((self.nx-1)*(self.ny-1)*(self.nz-1))))
                    outfile.write('%s\n'%('SCALARS ' + self.feature + ' float'))
                    outfile.write('%s\n'%'LOOKUP_TABLE default')                   
                    for i in range(ncells):
                        self.printProgressBar(i,ncells)
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
                    outfile.write('\n%s %d\n'%('POINT_DATA', ncoords))
                    outfile.write('%s\n'%('VECTORS ' + self.feature +' float'))
                    outfile.close()
                    outfile = open(self.dataOutPath + self.outFilename, 'ab+')
                    for i in range(ncoords):
                        data[i].tofile(outfile)

                # 
                # Write out the CELL CENTERED scalars
                #
                if "velocity" not in self.feature:
                    print "Writing scaler " + self.feature + " data..."
                    outfile.close()
                    outfile = open(self.dataOutPath + self.outFilename, 'a+')
                    outfile.write('\n%s %d\n'%('CELL_DATA', ncells))
                    outfile.write('%s\n'%('SCALARS ' + self.feature + ' float'))
                    outfile.write('%s\n'%'LOOKUP_TABLE default')
                    outfile.close()
                    outfile = open(self.dataOutPath + self.outFilename, 'ab+')
                    for i in range(ncells):
                        self.printProgressBar(int(i),ncells)
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
        d.quit()
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

