"""
This module reads in a VTK data file and plots it using matplotlib.
The VTK data must be filtered into the correct format using the
vtk filter prior to being plotted with matplotlib. An example script 
is included.

Usage:
The user MUST use the pipenv set up in the local directory. The data
will not be plotted otherwise. Don't ask me why, it's magic.

pipenv shell
python VTKtoNumpyPlotting.py

The input variables are included at the top of the file.
The user may also need to modify the name of the hydro field
when calling GetVTKScalarOutput, and will also have to 
modify the plotting function inputs to ensure proper labelling and 
scaling.

author: Evert Nasedkin evertn@student.ethz.ch
"""
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm,colors
import os,sys

# VTK File and path
VTK_FILE = "gasdensity260_5.vtk"
VTK_DIR = "/home/evert/Documents/SemesterProject/output/opa_5jup_50AU/VTK00260/"
# Constants
AU = 1.496e+13 #cm


def GetVTKScalarOutput(filename,filepath,scalar_name):
    """
    GetVTKOutput
    :filename: name of VTK File
    :filepath: directory of VTK file
    :scalar_name: Hydrodynamic field of interest, typically specified in filename

    This function creates a VTK reader for an unstructured grid,
    and returns the output object of the reader for the specified
    hydrodynamic field.
    """
    if os.path.isfile(filepath + filename):
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(VTK_DIR + VTK_FILE)
        reader.SetScalarsName(scalar_name)
        reader.Update()
    else:
        print "File does not exist."
        sys.exit(1)
    return  reader

def GetVTKVectorOutput(filename,filepath,vector_name):
    """
    GetVTKOutput
    :filename: name of VTK File
    :filepath: directory of VTK file
    :vector_name: Hydrodynamic field of interest, typically specified in filename

    This function creates a VTK reader for an unstructured grid,
    and returns the output object of the reader for the specified
    vectorial hydrodynamic field.
    """
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(VTK_DIR + VTK_FILE)
    reader.SetVectorsName(vector_name)
    reader.Update()
    return  reader.GetOutput()

def GetDataAsNumpy(vtk_mapper, array_num = 0):
    """
    GetDataAsNumpy
    :vtk_mapper: VTKmapper object containing the filtered data
    :array_num: array number of interest if multiple fields stored in one file
    
    Converts the VTK formatted data into a numpy array.
    """
    return vtk_to_numpy(vtk_mapper.GetOutput().GetPointData().GetArray(array_num))

def GetNumpyMesh(points):
    """
    GetNumpyMesh
    :points: VTK object containing all point positions
    
    Returns the x,y, and z axes of the mesh
    """
    n_points = points.GetNumberOfPoints()
    x = np.zeros(n_points)
    y = np.zeros(n_points)
    z = np.zeros(n_points)
    for i in range(n_points):
        pt = points.GetPoint(i)
        x[i] = pt[0]
        y[i] = pt[1]
        z[i] = pt[2]
    return x,y,z

def PlotData(x,y,data,levels = 255,
             title = "", xlabel = "", ylabel = "", cbar_label = "",
             cmap = cm.plasma, log = True, outdir = "", dpi = 400,
             figsize = (15,15)):
    fig,ax = plt.subplots(figsize = figsize)
    if log:
        # Plot a triangulated, filled contour.
        res =( np.ceil(np.log10(data.max())+1) - np.floor(np.log10(data.min())-1) )/levels
        lev_exp = np.arange(np.floor(np.log10(data.min())-1),
                            np.ceil(np.log10(data.max())+1),
                            res)
        levs = np.power(10, lev_exp)
        tcf = ax.tricontourf(x,y,data,levs, norm=colors.LogNorm(), cmap = cmap)
    else:
        tcf = ax.tricontourf(x,y,data,levels,cmap = cmap)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    cbar = plt.colorbar(tcf,ax = ax)
    cbar.set_label("Log of Gas Density [$g/cm^{3}$]")
    plt.savefig(outdir + title + ".png", dpi = dpi)
    plt.show()
    if len(outdir)>0 and not outdir.endswith("/") :
        outdir += "/"



# Read in the data
reader = GetVTKScalarOutput(VTK_FILE,VTK_DIR,"gasdensity")
"""
VTK FILTERING

This part of the script filters the data using the vtk/paraview filters
The data must be filtered to the correct shape before plotting with numpy

The example script here creates a slice at the mideplane of the the disk,
and plots the data on the xy plane.
"""
# Select a Slice
#create a plane to cut,here it cuts in the XZ direction (xz normal=(1,0,0);XY =(0,0,1),YZ =(0,1,0)
plane=vtk.vtkPlane()
plane.SetOrigin(0,0,10)
plane.SetNormal(0,0,1)

# Create cut filter
cutter=vtk.vtkCutter()
cutter.SetCutFunction(plane)
cutter.SetInputConnection(reader.GetOutputPort())
cutter.Update()

# Map filter to data
cutMapper = vtk.vtkCellDataToPointData()
cutMapper.AddInputData(cutter.GetOutput())
cutMapper.Update()
"""
END OF SCRIPTING
"""
# Get Data and mesh coordinates
x,y,z = GetNumpyMesh(cutter.GetOutput().GetPoints())
data = GetDataAsNumpy(cutMapper)
# Convert Units (optional, assumes VTK File was stored in CGS units
x = x/AU
y = y/AU
z = z/AU
# Plot the data
# Plot a triangulated, filled contour
'''
fig, ax = plt.subplots(figsize = (16,20))
res =( np.ceil(np.log10(data.max())+1) - np.floor(np.log10(data.min())-1) )/255
lev_exp = np.arange(np.floor(np.log10(data.min())-1),
                    np.ceil(np.log10(data.max())+1),
                    res)
levs = np.power(10, lev_exp)
tcf = ax.tricontourf(x,y,data,levs, norm=colors.LogNorm())
plt.show()
'''
PlotData(x,y,data,
         levels = 255,
         title = "Gas Density of circumstellar disk around a 5 Jupiter mass planet at 50 AU",
         xlabel = "x [AU]",
         ylabel = "y [AU]",
         cbar_label = "Log of Gas Density [$g/cm^{3}$]",
         outdir = VTK_DIR + "Images/",
         log = True,
         dpi = 150
         figsize = (8,8))
