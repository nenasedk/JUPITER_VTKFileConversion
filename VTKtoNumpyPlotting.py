import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata 
import sys


VTK_DIR = "/home/evert/Documents/SemesterProject/output/opa_1jup_50AU/VTK00198/"
VTK_FILE = "gasdensity198_5.vtk"

# Read the source file.
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(VTK_DIR + VTK_FILE)
reader.SetScalarsName("gasdensity")
reader.Update() # Needed because of GetScalarRange
output = reader.GetOutput()
#cellData = output.GetCellData()
#drange = cellData.GetRange()
scalar_range = output.GetScalarRange()

#cellCentersFilter = vtk.vtkCellCenters()
#cellCentersFilter.SetInput(output)
#cell

# Mesh Coordinates

nodes_vtk_array = output.GetPoints().GetData()

# Data Array
data_vtk_array = output.GetCells().GetData()#.GetScalars()
mapper = vtk.vtkCellDataToPointData()
mapper.AddInputData(output)
mapper.Update()
mapper_data = mapper.GetOutput()

# Node coordinates with data
nodes_numpy_array = vtk_to_numpy(nodes_vtk_array)
x,y,z = nodes_numpy_array[:,0], nodes_numpy_array[:,1], nodes_numpy_array[:,2]

data_numpy_array = vtk_to_numpy(mapper_data.GetPointData().GetArray(0))
data = data_numpy_array

# Contours!
npts = 100
xmin,xmax = min(x),max(x)
ymin,ymax = min(y),max(y)

# Grid
xi = np.linspace(xmin, xmax, npts)
yi = np.linspace(ymin, ymax, npts)
data_grid = griddata((x,y), data, (xi[None,:], yi[:,None]), method = 'cubic')

# Contour Plot
CS = plt.contour(xi,yi,data_grid, linewidths = 3, cmap = cm.jet)
plt.clabel(CS, inline=1,inline_spacing= 3, fontsize=12, colors='k', use_clabeltext=1)

plt.colorbar()
plt.show()
