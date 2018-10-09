# JUPITER .dat to .vtk file conversion class

'v1.0'

Evert Nasedkin, October 10 2018

This is the class that converts the binary .dat and descriptor files
output from JUPITER hydrodynamic simulations into .vtk format.

## REQUIREMENTS:
* python 2.7 or greater
* numpy
* pyvtk `pip install pyvtk`
* astropy (for unit converions)

## USAGE:
The class can be run from the Convert.py script 
```python Convert.py```

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

