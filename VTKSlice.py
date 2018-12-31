import vtk
from vtk.util.numpy_supper import vtk_to_numpy
import numpy as np

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

renderer = vtk.vtkRenderer()

# --- mappers, actors, render, etc. ---
# mapper and actor to view the cone
Mapper = vtk.vtkDataSetMapper()
Mapper.SetInputConnection(reader.GetOutputPort())
Mapper.SetScalarRange(scalar_range)
Mapper.GetLookupTable().SetScaleToLog10()

Actor = vtk.vtkActor()
Actor.SetMapper(Mapper)
Actor.GetProperty().SetOpacity(0.0)


#create a plane to cut,here it cuts in the XZ direction (xz normal=(1,0,0);XY =(0,0,1),YZ =(0,1,0)
plane=vtk.vtkPlane()
plane.SetOrigin(0,0,0)
plane.SetNormal(0,1,0)

#create cutter
cutter=vtk.vtkCutter()
cutter.SetCutFunction(plane)
cutter.SetInputConnection(reader.GetOutputPort())
cutter.Update()
cutterMapper=vtk.vtkDataSetMapper()
cutterMapper.SetScalarRange(scalar_range)
cutterMapper.GetLookupTable().SetScaleToLog10()
cutterMapper.SetInputConnection( cutter.GetOutputPort())

#create plane actor
planeActor=vtk.vtkActor()
planeActor.GetProperty().SetColor(1.0,1,0)
planeActor.GetProperty().SetLineWidth(2)
planeActor.SetMapper(cutterMapper)

# Create a Color Bar
scalarBar = vtk.vtkScalarBarActor()
scalarBar.SetLookupTable(cutterMapper.GetLookupTable())
scalarBar.SetTitle("Density")
#scalarBar.BoldOff()
#scalarBar.ShadowOff()
scalarBar.GetProperty().SetColor(0,0,0)

scalarBar.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
scalarBar.GetPositionCoordinate().SetValue(0.1,0.01)
scalarBar.SetOrientationToHorizontal()
scalarBar.SetWidth(0.8)
scalarBar.SetHeight(0.17)

# AXES
cubeAxesActor = vtk.vtkCubeAxesActor()
cubeAxesActor.SetBounds(reader.GetOutput().GetBounds())
cubeAxesActor.SetCamera(renderer.GetActiveCamera())
cubeAxesActor.GetTitleTextProperty(0).SetColor(1.0, 0.0, 0.0)
cubeAxesActor.GetLabelTextProperty(0).SetColor(1.0, 0.0, 0.0)

cubeAxesActor.GetTitleTextProperty(1).SetColor(0.0, 1.0, 0.0)
cubeAxesActor.GetLabelTextProperty(1).SetColor(0.0, 1.0, 0.0)

cubeAxesActor.GetTitleTextProperty(2).SetColor(0.0, 0.0, 1.0)
cubeAxesActor.GetLabelTextProperty(2).SetColor(0.0, 0.0, 1.0)

cubeAxesActor.DrawXGridlinesOn()
cubeAxesActor.DrawYGridlinesOn()
cubeAxesActor.DrawZGridlinesOn()
if vtk.VTK_MAJOR_VERSION > 5:
    cubeAxesActor.SetGridLineLocation(vtk.VTK_GRID_LINES_FURTHEST)

cubeAxesActor.XAxisMinorTickVisibilityOff()
cubeAxesActor.YAxisMinorTickVisibilityOff()
cubeAxesActor.ZAxisMinorTickVisibilityOff()

# A renderer and render window
renderer.SetBackground(0.2, 0.75, 0.85)

# add the actors
renderer.AddActor(Actor)
renderer.AddActor(planeActor)
renderer.AddActor(scalarBar)
renderer.AddActor(cubeAxesActor)

renderer.GetActiveCamera().Azimuth(30)
renderer.GetActiveCamera().Elevation(30)

renderer.ResetCamera()

renwin = vtk.vtkRenderWindow()
renwin.AddRenderer(renderer)

# An interactor
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renwin)

# Start
interactor.Initialize()
interactor.Start()
