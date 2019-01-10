import vtk
from vtk.util.misc import vtkGetDataRoot

VTK_DATA_ROOT = vtkGetDataRoot()
 
def MakeLUTFromCTF(tableSize,scalar_range):
    '''
    Use a color transfer Function to generate the colors in the lookup table.
    See: http://www.vtk.org/doc/nightly/html/classvtkColorTransferFunction.html
    :param: tableSize - The table size
    :return: The lookup table.
    '''
    ctf = vtk.vtkColorTransferFunction()
    ctf.SetColorSpaceToDiverging()
    # Green to tan.
    ctf.AddRGBPoint(scalar_range[0], 0.085, 0.532, 0.201)
    ctf.AddRGBPoint((scalar_range[0]+scalar_range[1])/2.0, 0.865, 0.865, 0.865)
    ctf.AddRGBPoint(scalar_range[1], 0.677, 0.492, 0.093)
 
    lut = vtk.vtkLogLookupTable()
    lut.SetNumberOfTableValues(tableSize)
    lut.Build()
 
    for i in range(0,tableSize):
        rgb = list(ctf.GetColor(float(i)/tableSize))+[1]
        lut.SetTableValue(i,rgb)
 
    return lut


# The source file
file_name = "VTK00262/gasdensity262_5.vtk"
#file_name = "VTK00262/uGrid.vtk"

# Read the source file.
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(file_name)
reader.Update()  # Needed because of GetScalarRange
reader.GetOutput().GetPointData().SetActiveScalars("Density")
scalar_range = reader.GetOutput().GetScalarRange()
lut = MakeLUTFromCTF(10,scalar_range)

isoContour = vtk.vtkContourFilter()
isoContour.SetInputConnection( reader.GetOutputPort() )
isoContour.SetValue(0,scalar_range[1])

rMapper = vtk.vtkDataSetMapper()
rMapper.SetInputConnection(reader.GetOutputPort())
rMapper.SetScalarRange(scalar_range)
rMapper.ScalarVisibilityOn()
rMapper.Update()

#create a plane to cut,here it cuts in the XZ direction (xz normal=(1,0,0);XY =(0,0,1),YZ =(0,1,0)
plane=vtk.vtkPlane()
plane.SetOrigin(0,0,0)
plane.SetNormal(-1,0,0)

#create cutter
cutter=vtk.vtkCutter()
cutter.SetCutFunction(plane)
cutter.SetInputConnection(reader.GetOutputPort())
cutter.Update()

cutterMapper=vtk.vtkPolyDataMapper()
cutterMapper.SetInputConnection( cutter.GetOutputPort() )
cutterMapper.SetLookupTable(lut)
cutterMapper.SetScalarRange(scalar_range)

#create plane actor
cutActor=vtk.vtkActor()
#planeActor.GetProperty().SetColor(1.0,1.0,0)
#planeActor.GetProperty().SetLineWidth(2)
cutActor.SetMapper(cutterMapper)

#create  actor
dActor=vtk.vtkActor()
#dActor.GetProperty().SetColor(0.5,1,0.5)
#dActor.GetProperty().SetOpacity(0.5)
dActor.SetMapper( rMapper )
dActor.GetProperty().SetRepresentationToSurface()


sagittal = vtk.vtkMatrix4x4()
sagittal.DeepCopy((0, 0,-1, 0,
                   1, 0, 0, 0,
                   0,-1, 0, 0,
                   0, 0, 0, 1))

# Extract a slice in the desired orientation
reslice = vtk.vtkImageReslice()
reslice.SetInputConnection(reader.GetOutputPort())
reslice.SetOutputDimensionality(2)
reslice.SetResliceAxes(sagittal)
reslice.SetInterpolationModeToLinear()

# Map the image through the lookup table
color = vtk.vtkImageMapToColors()
color.SetLookupTable(lut)
color.SetInputConnection(reslice.GetOutputPort())

# Display the image
colorActor = vtk.vtkImageActor()
colorActor.SetInput(color.GetOutput())


#create renderers and add actors of plane and data
ren = vtk.vtkRenderer()
ren.AddActor(colorActor)

#Add renderer to renderwindow and render
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(600, 600)
interactorStyle = vtk.vtkInteractorStyleImage()
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
ren.SetBackground(0.76,0.87,0.95)
renWin.Render()

# Create callbacks for slicing the image
actions = {}
actions["Slicing"] = 0

def ButtonCallback(obj, event):
    if event == "LeftButtonPressEvent":
        actions["Slicing"] = 1
    else:
        actions["Slicing"] = 0

def MouseMoveCallback(obj, event):
    (lastX, lastY) = interactor.GetLastEventPosition()
    (mouseX, mouseY) = interactor.GetEventPosition()
    if actions["Slicing"] == 1:
        deltaY = mouseY - lastY
        reslice.GetOutput().UpdateInformation()
        sliceSpacing = reslice.GetOutput().GetSpacing()[2]
        matrix = reslice.GetResliceAxes()
        # move the center point that we are slicing through
        center = matrix.MultiplyPoint((0, 0, sliceSpacing*deltaY, 1))
        matrix.SetElement(0, 3, center[0])
        matrix.SetElement(1, 3, center[1])
        matrix.SetElement(2, 3, center[2])
        window.Render()
    else:
        interactorStyle.OnMouseMove()
        

interactorStyle.AddObserver("MouseMoveEvent", MouseMoveCallback)
interactorStyle.AddObserver("LeftButtonPressEvent", ButtonCallback)
interactorStyle.AddObserver("LeftButtonReleaseEvent", ButtonCallback)

# Start interaction
iren.Start()
'''
output = reader.GetOutput()
scalar_range = output.GetScalarRange()
lut.SetTableRange(scalar_range)
lut.Build()
# Create the mapper that corresponds the objects of the vtk.vtk file
# into graphics elements
mapper = vtk.vtkDataSetMapper()
mapper.SetInputData(output)
mapper.SetInputConnection(reader.GetOutputPort())
mapper.SetScalarModeToUseCellData()
mapper.SetScalarRange(scalar_range)
mapper.SetLookupTable(lut)
# Create the Actor
actor = vtk.vtkActor()
actor.SetMapper(mapper)

# Create the Renderer
renderer = vtk.vtkRenderer()
renderer.AddActor(actor)
renderer.SetBackground(1, 1, 1)  # Set background to white

# Create the RendererWindow
renderer_window = vtk.vtkRenderWindow()
renderer_window.AddRenderer(renderer)

# Create the RendererWindowInteractor and display the vtk_file
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renderer_window)
interactor.Initialize()
interactor.Start()
'''
