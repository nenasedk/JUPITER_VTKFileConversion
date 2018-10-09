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
# VTK File Format specs:

import os,sys
import numpy as np
import astropy.units as u
import vtk
import string
import pyvtk
import Dialog

# Code Unit Constants, defined after mass and radius are given
TEMP = -1.
DENS = -1.
PERIOD = -1.
VEL = -1.
BASEPATH = os.getcwd() + '/'
class DATtoVTK:
    'Convert JUPITER .DAT files to binary VTK files for Paraview'
    def __init__(self):
        # Simulation Information
        self.simNumber = -1 # Which sim number are we looking at
        self.nLevel = -1 # How many mesh levels are there?
        self.feature = 'notafeat'
        self.nLevelCoords = []
        # density, energy, erad, opacity, potential, stheat, tau, taucell, velocity

        # Filepath information
        self.dataDir = 'notadir' # Where is the data from?
        self.dataOutPath = 'notapath' # Where do we output data
        self.inFilename = 'notaname' # Filename constructed from feature, simNumber and nLevels
        self.outFilename = 'notaname'
        self.descriptorName = 'notadesc'

        # Science information
        self.radius = -1.
        self.rcgs = -1.
        self.mass = -1.
        self.mcgs = -1.

        # Grid Information
        self.sphere = np.zeros((0,0,0),dtype=np.float64) # 3D Spherical array - edges (read in)
        self.mesh = np.zeros((0,0,0),dtype=np.float64)   # 3D cartesian array - edges
        self.cent = np.zeros((0,0,0),dtype=np.float64)   # 3D cartesian array of cell centers
        
        self.nx = 0
        self.ny = 0
        self.nz = 0
        self.ncell = 0

    # ------------------------------------------------------------------------------------------------
        
    # Basic user Set functions
    def SetSimNumber(self, n):
        self.simNumber = n

    def SetLevel(self,level):
        self.nLevel = level

    def SetFeature(self, feat):
        if feat in ['density', 'energy', 'erad', 'opacity', 'potential', 'stheat', 'tau', 'taucell', 'velocity']:
            self.feature = feat
            return 1
        else:
            print("Not an allowed input, please enter one of density, energy, erad, opacity, potential, stheat, tau, taucell or velocity")
            return -1
    def SetBasePath(path):
        BASEPATH = path

    def SetRadius(self, rad):
        self.radius = rad*u.AU
        self.rcgs = self.radius.to(u.cm)
    def SetMass(self,mass):
        self.mass = mass*u.M_sun
        self.mcgs = self.mass.to(u.g)
        if(self.radius>0):
            TEMP = ((self.rcgs)/(np.sqrt((self.rcgs)**3 / 6.67259e-8 / (self.mcgs))))**2 / 8.314e7
            DENS = self.mcgs/(self.rcgs)**3
            PERIOD = 2*np.pi*np.sqrt((self.rcgs)**3 / (6.67259e-8 * self.mcgs))
            VEL = self.rcgs/(PERIOD/2*np.pi)
            
    # ----------------------------------------------------------------------------------------
            
    # Directory and file Setup
    def SetupDirs( self ):
        # Error checking
        if self.simNumber < 0:
            print("Please set the simulation number")
            return
        if self.feature is 'notafeat':
            print("Please input a feature (e.g. velocity)")
            return
        

        # Create directory paths
        self.dataDir = BASEPATH + "output" + str(self.simNumber).zfill(5) + "/"
        self.dataOutPath =  BASEPATH + "VTK" + str(self.simNumber).zfill(5) +"/"
        self.descriptorName = "Descriptor" + str(self.simNumber) + ".dat"
        self.outFilename = "gas" + self.feature + str(self.simNumber) +"_" + str(self.nLevel) + ".vtk"

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
        if self.simNumber < 0:
            print("Please set the simulation number")
            return
        if self.feature is 'notafeat':
            print("Please input a feature (e.g. velocity)")
            return
        if self.nLevel < 0:
            print("Please input the number of mesh refinement levels")
            return

        self.inFilename = "gas" + self.feature + str(self.simNumber) +"_"+ str(level) + "_" + str(level) + ".dat"


    # ---------------------------------------------------------------------------------------------
    # Important part starts here

    #    
    # This function wraps the binary .dat file reader for a given feature,
    # and output for the field in a VTK format. This is the only user facing function.
    #
    def ConvertFiles( self , binary = True ):
        data = np.zeros(0,dtype = np.float)
        for i in range(self.nLevel):
            feat = np.zeros(0)
            self.SetupNames(i)
            if( self.feature is not 'velocity'): # velocity is a vector, so is treated differently, (and is vertex centered?)
                feat = np.zeros((self.nLevelCoords[i][0],self.nLevelCoords[i][1],self.nLevelCoords[i][2]))
            else:
                feat = np.zeros((self.nLevelCoords[i][0],self.nLevelCoords[i][1],self.nLevelCoords[i][2],3))
        

            feat = np.fromfile(self.dataDir + self.inFilename, dtype = 'double')
            data = np.append(data,feat.astype(float))

        # Convert to CGS units
        if(self.feature is 'density'):
            feat = feat0*DENS
        if(self.feature is 'temperature'):
            feat = feat0*DENS
        if(self.feature is 'velocity'): 
            feat = feat0*VEL
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
        if self.feature is not 'velocity':
            try:
                assert(len(data) == self.ncell) 
            except ValueError:
                print "Error: Number of data points does not match number of mesh elements"
                return
        else:
            try:
                assert(len(data) == len(self.mesh))
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
                   
            if userstr is 'Append': # not sure if append will work yet
                self.VTKFormatOut(data, binary, True)
            else:
                print "Nothing written to file!"
                return


    #                                  
    #  Format data into VTK structure  
    #                                        
    def VTKFormatOut(self, data, binary = True, append = False):
        
        ncoords = self.mesh.shape[0]
        ncells = self.ncell
        # VTK Header
        if not append:
            # -------------------------------------------
            # Ascii out
            # -------------------------------------------
            outfile = open(self.dataOutPath + self.outFilename, 'w+')
            outfile.write("# vtk DataFile Version 2.0\n")
            outfile.write("JUPITER Sim" + str(self.simNumber) + " " + self.feature + " field\n")
            if binary:
                outfile.write("BINARY\n")
            else:
                outfile.write("ASCII\n")
                outfile.write("\nDATASET UNSTRUCTURED_GRID\n")
            # Write out based on RadMC3D WriteVTK()
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
            ''' CELL COORDINATE OUT'''
            if not binary:
                outfile.write('%s\n'%('POINTS '+str(ncoords)+' float'))
                print "Writing points to file..."
                for i in range(ncoords):
                    self.printProgressBar(int(i),ncoords)
                    outfile.write('%.9e %9e %9e\n'%(self.mesh[i][0], self.mesh[i][1], self.mesh[i][2]))

                outfile.write('\n%s %d %d\n'%('CELLS ', ncells, ncells*9))
                print "Writing coordinates to file..."

                ntx = 0
                nty = 0
                ntz = 0
                npz = npy = npx = 0
                nn = nt= 0
                for n in range(self.nLevel): # mesh refinement levels are written out sequentially
                    ntx +=self.nLevelCoords[n][0]
                    nty +=self.nLevelCoords[n][1]
                    ntz +=self.nLevelCoords[n][2]
                    #print "tx: " + str(ntx) + ", cx: " + str(self.nLevelCoords[n][0]) +  ", px: " + str(npx)
                    #print "ty: " + str(nty) + ", cy: " + str(self.nLevelCoords[n][1]) +  ", py: " + str(npy)
                    #print "tz: " + str(ntz) + ", cz: " + str(self.nLevelCoords[n][2]) +  ", pz: " + str(npz)
                    for iz in range(self.nLevelCoords[n][2]-1):
                        for iy in range(self.nLevelCoords[n][1]-1):
                            for ix in range(self.nLevelCoords[n][0]-1):   # Calculate the index of each of the vertices of a given cell, and write to file.
                                nn +=1
                                id0 = nt + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*iy     + ix
                                id1 = nt + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*iy     + ((ix+1) % (self.nLevelCoords[n][0]-1))
                                id2 = nt + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*(iy+1) + ((ix+1) % (self.nLevelCoords[n][0]-1)) 
                                id3 = nt + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*(iy+1) + ix 
                                id4 = nt + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*iy     + ix
                                id5 = nt + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*iy     + ((ix+1) % (self.nLevelCoords[n][0]-1))
                                id6 = nt + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*(iy+1) + ((ix+1) % (self.nLevelCoords[n][0]-1))
                                id7 = nt + self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*(iy+1) + ix
                                
                                line = np.array([8,id0,id1,id2,id3,id4,id5,id6,id7])
                                line.tofile(outfile, sep=' ', format='%d')
                                outfile.write('\n')
                                self.printProgressBar(int(nn),ncells)
                    npx = ntx
                    npy = nty
                    npz = ntz
                    nt = nn
                print "\n" + str(nt) + " Cells (*8 for vertices) written"
                #
                # Now write out the type of each cell (#12)
                #
                outfile.write('\n%s %d\n'%('CELL_TYPES',ncells))
                for ix in range(ncells):     
                    outfile.write('%d\n'%12)

                # 
                # Now write out the CORNER CENTERED velocities, if existing (Vector outs)
                #
                if self.feature is 'velocity':
                    print "Writing Vector Data"
                    outfile.write('\n%s %d\n'%('POINT_DATA', ncoords))
                    outfile.write('%s\n'%('VECTORS ' + self.feature +' float'))
                    for i in range(ncoords):
                        outfile.write('%.9e %.9e %.9e\n'%(data[i][0], data[i][1], data[0][2]))

                # 
                # Write out the CELL CENTERED scalars
                #
                if self.feature is not 'velocity':
                    print "Writing scaler " + self.feature + " data..."
                    outfile.write('\n%s %d\n'%('CELL_DATA', ncells))
                    outfile.write('%s\n'%('SCALARS ' + self.feature + ' float'))
                    outfile.write('%s\n'%'LOOKUP_TABLE default')                   
                    for i in range(ncells):
                        self.printProgressBar(int(i),ncells)
                        outfile.write('%.9e\n'%data[i])

                outfile.close()
                return
            else:
                # -------------------------------------------------
                # Binary out (overwrite or new)
                # -------------------------------------------------
                outfile.write('%s\n'%('POINTS '+str(ncoords)+' float'))
                print "Writing points to file..."
                outfile.close()
                outfile = open(self.dataOutPath + self.outFilename, 'ab+')
                for i in range(ncoords):
                    self.printProgressBar(int(i),ncoords)
                    self.mesh[i].byteswap(True).tofile(outfile)
                outfile.close()
                outfile = open(self.dataOutPath + self.outFilename, 'a+')
                outfile.write('\n%s %d %d\n'%('CELLS ', ncells, ncells*9))
                print "Writing coordinates to file..."
                outfile.close()
                outfile = open(self.dataOutPath + self.outFilename, 'ab+')

                ntx = 0
                nty = 0
                ntz = 0
                npz = npy = npx = 0 
                for n in range(self.nLevel):
                    ntx +=self.nLevelCoords[n][0]
                    nty +=self.nLevelCoords[n][1]
                    ntz +=self.nLevelCoords[n][2]
                    for iz in range(npz, ntz-1):
                        self.printProgressBar(int(iz),self.nz-self.nLevel)
                        for iy in range(npy, nty-1):
                            for ix in range(npx, ntx-1):         
                                id0 = self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*iy     + ix
                                id1 = self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*iy     + ((ix+1) % (self.nLevelCoords[n][0]-1))
                                id2 = self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*(iy+1) + ((ix+1) % (self.nLevelCoords[n][0]-1)) 
                                id3 = self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*iz     + self.nLevelCoords[n][0]*(iy+1) + ix 
                                id4 = self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*iy     + ix
                                id5 = self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*iy     + ((ix+1) % (self.nLevelCoords[n][0]-1))
                                id6 = self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*(iy+1) + ((ix+1) % (self.nLevelCoords[n][0]-1))
                                id7 = self.nLevelCoords[n][0]*self.nLevelCoords[n][1]*(iz+1) + self.nLevelCoords[n][0]*(iy+1) + ix
                                
                                line = np.array([8,id0,id1,id2,id3,id4,id5,id6,id7])
                                line.byteswap(True).tofile(outfile)#, sep=' ', format='%d')
                                self.printProgressBar(int(nn),ncells)
                                #outfile.write('\n')
                    npx = ntx
                    npy = nty
                    npz = ntz
                #
                # Now write out the type of each cell (#12)
                #
                outfile.close()
                outfile = open(self.dataOutPath + self.outFilename, 'a+')
                outfile.write('\n%s %d\n'%('CELL_TYPES',ncells))
                outfile.close()
                outfile = open(self.dataOutPath + self.outFilename, 'ab+')
                for ix in range(ncells):
                    np.array([12]).byteswap(True).tofile(outfile)

                # 
                # Now write out the CORNER CENTERED velocities, if existing (Vector outs)
                #
                outfile.close()
                outfile = open(self.dataOutPath + self.outFilename, 'a+')
                if self.feature is 'velocity':
                    print "Writing Vector Data"
                    outfile.write('\n%s %d\n'%('POINT_DATA', ncoords))
                    outfile.write('%s\n'%('VECTORS ' + self.feature +' float'))
                    outfile.close()
                    outfile = open(self.dataOutPath + self.outFilename, 'ab+')
                    for i in range(ncoords):
                        data[i].byteswap(True).tofile(outfile)

                # 
                # Write out the CELL CENTERED scalars
                #
                if self.feature is not 'velocity':
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

                outfile.close()
                return
        else:
            # -----------------------------------------------
            # Ascii out, append to file
            # -----------------------------------------------
            outfile = open(self.dataOutPath + self.outFilename, 'a+')
            if not binary:
                if self.feature is 'velocity':
                    print "Writing Vector Data"
                    outfile.write('\n%s %d\n'%('POINT_DATA', ncoords))
                    outfile.write('%s\n'%('VECTORS ' + self.feature +' float'))
                    for i in range(0,ncoords):
                        outfile.write('%.9e %.9e %.9e\n'%(data[i][0], data[i][1], data[0][2]))

                # 
                # Write out the CELL CENTERED scalars
                #
                if self.feature is not 'velocity':
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
                outfile = open(self.dataOutPath + self.outFilename, 'a+')
                if self.feature is 'velocity':
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
                if self.feature is not 'velocity':
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
