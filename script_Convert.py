# This is a script to automate the conversion of dat to vtk files
#
# default usage:
# python script_convert.py [first] [grid_level] [radius] [mass] [field list]
#
# with optional arguments:
# python script_convert.py [first] [grid level] [radius] [mass] -l [last] -b [binary/ascii] -d [dir] [field list]
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
parser.add_argument('-l',action ='store',nargs = 1,type = int,
                    metavar = 'last_out',
                    help= 'last simulation output to convert')
parser.add_argument('r',action = 'store',nargs = 1,type = int,
                    metavar = 'radius', help= 'planetary orbital radius')
parser.add_argument('m',action = 'store',nargs = 1,type = int,
                    metavar = 'mass', help= 'stellar mass')
parser.add_argument('-d',action = 'store', nargs = 1,
                    help= 'directory to output folders')
parser.add_argument('-b',action = 'store', nargs = 1, metavar = 'binary',
                    help= '(b)inary or (a)scii')
parser.add_argument('fields',action = 'append', nargs = argparse.REMAINDER,
                    help= 'list of hydrodynamic fields to convert')

args = parser.parse_args()

if args.l is None:
    args.l = [args.o[0] + 1]
if args.b is None:
    args.b = ['b']

print "Converting outputs from " + str(args.o[0]) + " to " + str(args.l[0])
print "Using fields " + str(args.fields[0])
for i in range(args.o[0],args.l[0]):
    for index,f in enumerate(args.fields[0]):
        dv = DATtoVTK.DATtoVTK()
        dv.SetOutNumber(i)
        dv.SetLevel(args.g[0])
        dv.SetRadius(args.r[0])
        dv.SetMass(args.m[0])
        dv.SetFeature(f)
        if args.d is None:
            dv.SetupDirs()
        else:
            dv.SetBasePath(args.d[0])
            dv.SetupDirs()
        if args.b[0] == 'ascii' or args.b[0] == 'a':
            dv.ConvertFiles(False)
        else:
            dv.ConvertFiles(True)    

