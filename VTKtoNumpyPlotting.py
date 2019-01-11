import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.colors as colors
from scipy.interpolate import griddata 
import sys

AU = 1.496e+13
VTK_DIR = "/home/evert/Documents/SemesterProject/output/opa_1jup_50AU/VTK00198/"
VTK_FILE = "gasdensity198_5.vtk"

# Read the source file.
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(VTK_DIR + VTK_FILE)
reader.SetScalarsName("gasdensity")
reader.Update() # Needed because of GetScalarRange
output = reader.GetOutput()

scalar_range = output.GetScalarRange()

# Select a Slice
#create a plane to cut,here it cuts in the XZ direction (xz normal=(1,0,0);XY =(0,0,1),YZ =(0,1,0)
plane=vtk.vtkPlane()
plane.SetOrigin(0,0,1)
plane.SetNormal(0,0,1)

#create cutter
cutter=vtk.vtkCutter()
cutter.SetCutFunction(plane)
cutter.SetInputConnection(reader.GetOutputPort())
cutter.Update()

# Data Array
data_vtk_array = cutter.GetOutput() #output.GetCells().GetData()#.GetScalars()

mapper = vtk.vtkCellDataToPointData()
mapper.AddInputData(data_vtk_array)
mapper.Update()
mapper_data = mapper.GetOutput()


# Mesh Coordinates
points = data_vtk_array.GetPoints()
n_points = points.GetNumberOfPoints()
# Node coordinates with data
#nodes_numpy_array = vtk_to_numpy(data_vtk_array.GetPoints())
x = np.zeros(n_points)
y = np.zeros(n_points)
z = np.zeros(n_points)
for i in range(n_points):
    pt = points.GetPoint(i)
    x[i] = pt[0]
    y[i] = pt[1]
    #z[i] = pt[2]
    
x = x/AU
y = y/AU
data_numpy_array = vtk_to_numpy(mapper_data.GetPointData().GetArray(0))
npts = len(x)
data = np.abs(data_numpy_array)
# Contours!
xmin,xmax = min(x),max(x)
ymin,ymax = min(y),max(y)

# Grid
data_grid = griddata((x,y), data, (x,y), method = 'nearest')
lvls = np.logspace(data.min(),data.max(),16)

plt.tricontourf(x,y,np.log10(data),255,)#, lvls)#cmap = cm.plasma, norm = colors.LogNorm())
plt.title("Log of Gas Density of circumstellar disk around a 1 Jupiter mass planet at 1AU")
plt.xlabel("x [AU]")
plt.ylabel("y [AU]")

cbar = plt.colorbar()
cbar.set_label("Log of Gas Density [$g/cm^{3}$]")
plt.show()
