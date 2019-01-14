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

Alternatively

```python script_convert.py [first] [grid_level] [radius] [mass] -f [field list]```

or with full options:

```python script_convert.py [first] [grid level] [radius] [mass] -l [last] [-v] -b [binary/ascii] -d [dir] -f [field list]```

See script_convert.py for more details.

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
