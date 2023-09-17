import random
import numpy as np

from vtkmodules.vtkRenderingAnnotation import (
    vtkAnnotatedCubeActor,
    vtkAxesActor
)
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkPolyDataMapper,
    vtkPropAssembly,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkRenderer
)
import vtk
import xml.etree.ElementTree as ET

# Erstellen ein VTK-Array aus Ihrem NumPy-Array
# vtkDataArray = vtk.vtkFloatArray()
# vtkDataArray.SetName("MyData")
# vtkDataArray.SetNumberOfComponents(3)
# for i in range(gridToPolyData.GetNumberOfPoints()):
#     vtkDataArray.InsertNextTuple([1.0, 0.0,  0.0])
# pd = gridToPolyData.GetPointData()
# pd.AddArray(vtkDataArray)

def make_axes_actor(scale, xyz_labels):
    """
    :param scale: Sets the scale and direction of the axes.
    :param xyz_labels: Labels for the axes.
    :return: The axes actor.
    """
    axes = vtkAxesActor()
    axes.SetScale(scale)
    axes.SetShaftTypeToCylinder()
    axes.SetCylinderRadius(1.5 * axes.GetCylinderRadius())
    axes.SetConeRadius(2.0 * axes.GetConeRadius())
    axes.SetSphereRadius(1.5 * axes.GetSphereRadius())
    axes.AxisLabelsOff()
    return axes


def make_annotated_cube_actor(cube_labels, colors):
    """
    :param cube_labels: The labels for the cube faces.
    :param colors: Used to determine the cube color.
    :return: The annotated cube actor.
    """
    # A cube with labeled faces.
    cube = vtkAnnotatedCubeActor()
    cube.SetXPlusFaceText(cube_labels[0])
    cube.SetXMinusFaceText(cube_labels[1])
    cube.SetYPlusFaceText(cube_labels[2])
    cube.SetYMinusFaceText(cube_labels[3])
    cube.SetZPlusFaceText(cube_labels[4])
    cube.SetZMinusFaceText(cube_labels[5])
    cube.SetFaceTextScale(0.5)
    cube.GetCubeProperty().SetColor(colors.GetColor3d('Gainsboro'))

    cube.GetTextEdgesProperty().SetColor(colors.GetColor3d('LightSlateGray'))

    # Change the vector text colors.
    cube.GetXPlusFaceProperty().SetColor(colors.GetColor3d('Tomato'))
    cube.GetXMinusFaceProperty().SetColor(colors.GetColor3d('Tomato'))
    cube.GetYPlusFaceProperty().SetColor(colors.GetColor3d('DeepSkyBlue'))
    cube.GetYMinusFaceProperty().SetColor(colors.GetColor3d('DeepSkyBlue'))
    cube.GetZPlusFaceProperty().SetColor(colors.GetColor3d('SeaGreen'))
    cube.GetZMinusFaceProperty().SetColor(colors.GetColor3d('SeaGreen'))
    return cube


def make_cube_actor(label_selector, colors):
    """
    :param label_selector: The selector used to define labels for the axes and cube.
    :param colors: Used to set the colors of the cube faces.
    :return: The combined axes and annotated cube prop.
    """
    if label_selector == 'sal':
        # xyz_labels = ['S', 'A', 'L']
        xyz_labels = ['+X', '+Y', '+Z']
        cube_labels = ['S', 'I', 'A', 'P', 'L', 'R']
        scale = [1.5, 1.5, 1.5]
    elif label_selector == 'rsp':
        # xyz_labels = ['R', 'S', 'P']
        xyz_labels = ['+X', '+Y', '+Z']
        cube_labels = ['R', 'L', 'S', 'I', 'P', 'A']
        scale = [1.5, 1.5, 1.5]
    elif label_selector == 'lsa':
        # xyz_labels = ['L', 'S', 'A']
        xyz_labels = ['+X', '+Y', '+Z']
        cube_labels = ['L', 'R', 'S', 'I', 'A', 'P']
        scale = [1.5, 1.5, 1.5]
    else:
        xyz_labels = ['+X', '+Y', '+Z']
        cube_labels = ['+X', '-X', '+Y', '-Y', '+Z', '-Z']
        scale = [1.5, 1.5, 1.5]

    # We are combining a vtkAxesActor and a vtkAnnotatedCubeActor
    # into a vtkPropAssembly
    cube = make_annotated_cube_actor(cube_labels, colors)
    axes = make_axes_actor(scale, xyz_labels)

    # Combine orientation markers into one with an assembly.
    assembly = vtkPropAssembly()
    assembly.AddPart(axes)
    assembly.AddPart(cube)
    return assembly


# Funktion zum Erstellen eines vtkColorTransferFunction aus XML
def create_color_transfer_function_from_xml(xml_filename):
    colorTransferFunction = vtk.vtkColorTransferFunction()
    
    tree = ET.parse(xml_filename)
    root = tree.getroot()
    
    for point in root.iter("Point"):
        x = float(point.get("x"))
        r = float(point.get("r"))
        g = float(point.get("g"))
        b = float(point.get("b"))
        colorTransferFunction.AddRGBPoint(x, r, g, b)
    
    color_steps = 8
    # Erstellen einer VTK-Lookup-Tabelle (LookupTable) basierend auf der Farbtransferfunktion
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(color_steps)  # Anzahl der Werte in der LUT (z.B. 256 für 8-Bit Farbtiefe)
    lut.Build()

    # Anwenden der Farbtransferfunktion auf die Lookup-Tabelle
    for i in range(color_steps):
        x = i / (color_steps - 1.0)  # Skalieren auf den Bereich [0, 1]
        color = colorTransferFunction.GetColor(x)
        lut.SetTableValue(i, color[0], color[1], color[2], 1.0)  # Alpha auf 1.0 (undurchsichtig) festlegen

    colorTransferFunction.SetVectorModeToMagnitude()
    lut.SetVectorModeToMagnitude()
    return colorTransferFunction, lut


def iglToVtkPolydata(sf, sv):
    # Erstellen Sie ein VTK-PolyData-Objekt für die Dreiecke
    triangle_polydata = vtk.vtkPolyData()

    # Erstellen Sie Punkte und Zellen für die Dreiecke
    points = vtk.vtkPoints()
    for v in sv:
        points.InsertNextPoint(v)

    cells = vtk.vtkCellArray()
    for T in sf:
        triangle = vtk.vtkTriangle()
        triangle.GetPointIds().SetId(0, T[0])
        triangle.GetPointIds().SetId(1, T[1])
        triangle.GetPointIds().SetId(2, T[2])
        cells.InsertNextCell(triangle)

    # Fügen Sie die Punkte und Zellen zur PolyData hinzu
    triangle_polydata.SetPoints(points)
    triangle_polydata.SetPolys(cells)

    return triangle_polydata


def addRandomCellData(triangle_polydata):
    tgroupping = vtk.vtkFloatArray()
    tgroupping.SetNumberOfComponents(1)
    tgroupping.SetName("group")

    n = triangle_polydata.GetNumberOfCells()
    for i in range(n):
        tgroupping.InsertNextTuple1(random.random())
    
    triangle_polydata.GetCellData().SetScalars(tgroupping)

    return triangle_polydata


def addScalarCellData(triangle_polydata, cell_data):
    tgroupping = vtk.vtkFloatArray()
    tgroupping.SetNumberOfComponents(1)
    tgroupping.SetName("group")

    scale = len(np.unique(cell_data))
    n = triangle_polydata.GetNumberOfPoints()

    print("polys     : " + str(triangle_polydata.GetNumberOfPolys()))
    print("points    : " + str(triangle_polydata.GetNumberOfPoints()))
    print("cells     : " + str(triangle_polydata.GetNumberOfCells()))
    print("verts     : " + str(triangle_polydata.GetNumberOfVerts()))
    print("celldata  : " + str(len(cell_data)))

    for i in range(n):
        tgroupping.InsertNextTuple1(float(cell_data[i])/float(scale))
    
    triangle_polydata.GetPointData().SetScalars(tgroupping)

    return triangle_polydata
