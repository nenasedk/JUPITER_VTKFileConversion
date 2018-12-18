# This is a script to automate the conversion of dat to vtk files
#
# default usage:
# python script_convert.py [first] [grid_level] [radius] [mass] -f [field list]
#
# with optional arguments:
# python script_convert.py [first] [grid level] [radius] [mass] -l [last] [-v] -b [binary/ascii] -d [dir] -f [field list]
#
# The list of fields must always be the final argument. 
import DATtoVTK
import argparse
import sys
import string

parser = argparse.ArgumentParser(description='dat to vtk file conversion')
parser.add_argument('o',action = 'store',nargs = 1,type = int,
                    metavar = 'first_out',
                    help= 'first simulation output to convert')
parser.add_argument('g',action = 'store',nargs = 1,type = int,
                    metavar = 'grid_level',
                    help= 'number of grid refinement levels')
parser.add_argument('-l', "--last",action ='store',nargs = 1,type = int,
                    metavar = 'last_out', required = False,
                    help= 'last simulation output to convert')
parser.add_argument('r',action = 'store',nargs = 1,type = float,
                    metavar = 'radius', help= 'planetary orbital radius')
parser.add_argument('m',action = 'store',nargs = 1,type = float,
                    metavar = 'mass', help= 'stellar mass')
parser.add_argument('-d',"--directory",action = 'store', nargs = 1,
		    required = False,	
                    help= "Directory to output folders")
parser.add_argument('-b',action = 'store', nargs = 1, metavar = 'binary', required = False,
                    help= '(b)inary or (a)scii')
parser.add_argument('-v', action = 'store_true',required = False,
                    help= 'Planet centred vels if included')
parser.add_argument('-f','--fields',action = 'append', nargs = argparse.REMAINDER,
                    help= 'list of hydrodynamic fields to convert',
	   	    required = True)

PVEL = False
args = parser.parse_args()

if args.last is None:
    args.last = [args.o[0]]
if args.b is None:
    args.b = ['b']
if args.v is not None:
    PVEL = True
print "Converting outputs from " + str(args.o[0]) + " to " + str(args.last[0])
print "Using fields " + str(args.fields[0])
print args.directory[0]
for i in range(args.o[0],args.last[0]+1):
    for index,f in enumerate(args.fields[0]):
        dv = DATtoVTK.DATtoVTK()
        dv.SetOutNumber(i)
        dv.SetLevel(args.g[0])
        dv.SetRadius(args.r[0])
        dv.SetMass(args.m[0])
        dv.SetFeature(f)
        if args.directory is None:
            dv.SetupDirs()
        else:
            dv.SetBasePath(args.directory[0])
            dv.SetupDirs()
        if args.b[0] == 'ascii' or args.b[0] == 'a':
            dv.ConvertFiles(False,PVEL)
        else:
            dv.ConvertFiles(True,PVEL)    

