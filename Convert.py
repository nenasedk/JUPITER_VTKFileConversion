# JUPITER .dat to .vtk file conversion script
#
# Evert Nasedkin, October 9, 2018
# evertn@student.ethz.ch
#
#  This script takes the user input from the terminal and
# passes it to the DATtoVTK Class to produce the converted files
#
# Additional inputs can be:
#     - dv.SetBasePath(dir), which sets the path where data directories are stored (defaults to cwd)
#     - Can set binary = True/False 
#
# Fixme: add feature so that all gas properties can be written to one vtk file

import DATtoVTK 
__CONT__ = True
BINARY = True # Set to false to write ASCII files
PATH = ''
while __CONT__:
    dv = DATtoVTK.DATtoVTK()
    dv.setBasePath(PATH) 
    while True:
        userstr = raw_input("Simulation Output number: ")
        try :
            userin = int(userstr)     
        except ValueError:
            print("Simulation number must be an integer.")
        else :
            break
    
    dv.SetOutNumber(userin)


    while True:
        userstr = raw_input("\nNumber of mesh levels (including base level): ")
        try :
            userin = int(userstr)
        except ValueError:
            print("Number of mesh levels must be an integer.")
        else :
            break
    dv.SetLevel(userin)
        

    err = -1
    while(err<0):
        userstr = raw_input("\nEnter a hydro field (e.g. gasdensity or gasvelocity): ")
        err = dv.SetFeature(userstr)
    
    while True:
        userstr = raw_input("\nOrbital radius of planet in AU: ")
        try :
            userin = float(userstr)     
        except ValueError:
            print("Radius must be a valid float.")
        else :
            break      
    dv.SetRadius(userin)

    while True:
        userstr = raw_input("\nMass of star in solar units: ")
        try :
            userin = float(userstr)     
        except ValueError:
            print("Radius must be a valid float.")
        else :
            break
    dv.SetMass(userin)
    print('\n')
    # Setup Directories
    dv.SetupDirs()

    # Extract mesh and data, output to VTK binary file
    dv.ConvertFiles(BINARY) # Can take binary = True/False (defaults to True)
    cont1 = True
    while cont1:
        userstr = raw_input("\nDo you want to convert another file? (y/n)")
        if userstr == 'y':
            __CONT__ = True
            cont1 = False
            print "\n\n\n\n\n"
        elif userstr == 'n':
            __CONT__ = False
            cont1 = False
        else:
            print "Please enter y or n... "

print "Exiting..."
                
                
            
