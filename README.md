# JUPITER .dat to .vtk file conversion class

`v2.1`

Evert Nasedkin, Noveber 26, 2018

This is the class that converts the binary .dat and descriptor files
output from JUPITER hydrodynamic simulations into .vtk format.

## REQUIREMENTS:
* python 2.7 or greater
* numpy
* vtk `pip install vtk`
* astropy (for unit converions)

## Environment
A pipenv environment has been setup for this project. After pipenv has been installed with pip, the program can be run in the specified environemnt using ```pipenv shell'''. This is primarily for use with the matplotlib style plots of the VTK output files.

## USAGE:
The class can be run from the Convert.py script 

```python Convert.py```

Alternatively, the script_Convert program can be used from the command line. 
The following lists the required arguments. 
Note that the field list argument MUST be the final argument in the command.

```python script_Convert.py [first] [grid_level] [radius] [mass] [-f [field list]]```

or with full options:

```python script_Convert.py [first] [-l [last]] [grid level] [radius] [mass] [-v] [-b [binary/ascii]] [-d [dir]] [-u [units]] [-f [field list]]```

The directory should contain a folder labelled outputXXXXX, where XXXXX is the simulation output number padded to 5 digits. If multiple files are being converted, the -d argument should contain all of the output folders.
Units can be "CGS" or "AU"
For a full description of the arguments, type ```python script_Convert.py --help```.
An example, for output 280, with 5 mesh levels, at 5.2 AU around a solar mass star, located in the current directory, with planet centered velocities:

```python script_Convert.py 280 5 5.2 1.0 -v -d ./ -f gasdensity gasvelocity```

See script_Convert.py for more details.

## Process:
Set up of  science factors and directories based on the user input
Read in mesh vertices from descriptor file, convert to cartesian coordinates
Calculate cell centers
Read in science data from .dat file
Setup VTK File format structure

[VTK File Format specs](https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf)

The .dat data file contains a 1D array of doubles written to a binary file
The descriptor file is structured as follows:
8 lines of simulation information, used by JUPITER, followed by
3 lines of grid data. In order: azimuthal angle, radius and polar angle

The first and last two grid points in each dimension are overlap with the
previous mesh resolution, or are just ghost points. They are deleted.

To reconstruct the grid, the coordinates are iterated in column order,
(ie azimuthal angle is iterated the fastest, polar angle the slowest)
It is converted to Cartesian coordinates to be written to the VTK file.

On importing the coordinate grid it is changed from left handed to right 
handed (azimuthal and polar angles)*-1

## Matplotlib visualisation
The VTKtoNumpyPlotting.py module contains functions that can be used to read in VTK data and output a matplotlib style figure.
This module may be used as a template for numpy plotting, but note that the user will be required to modify the script so as to filter the data as necessary. 
It is recommended only to use the matplotlib visualisations for 2D contour plots, as vector and streamline plots do not produce useful results.

## Known Issues

## Changelog
v1.1
- Added filtering of vertices and cells to remove overlapping regions (not working properly)
- Refactored code to be more modular (ie SphereToCart was added) and remove duplicate code
- Improved scriptability by adding extra user facing functions, with a template python script

v1.2
- Changed calculation of vertex indices - count points in plane rather than multiplying axes
- Lots of minor +-1 changes to ensure index/coordinate computations are done correctly

v1.3
- Completely changed index calculation
  - documented in ComputeIndices() and ComputeCell()
  - individually count distinct axes, should deal with edge cases
- Added functions to check cell skipping, initialisation
- Separated Filtered and Unfiltered counting
- So far untested, likely off by one errors, but overall computation should be closer than prev versions.

v2.0
- Changed index counting again
- Paraview *can* deal with unused grid points, but not with unused cells
- Therefore, filtering the grid is unneccsary, we can just remove cells with points in the ROI
- This works, and is much easier and faster.

v2.1
- Velocity calculations correct
- Added star centered velocity function
- Started filling in extra documentation

v2.2
- Cleaned up script_Convert.py inputs and documentation
- Added VTKtoNumpyPlotting.py for matplotlib outputs
- Fixed default velocities settings in DAT_to_VTK
- Added pipenv environment
- Updates to report and documentation
