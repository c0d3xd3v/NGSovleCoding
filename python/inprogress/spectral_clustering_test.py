#!/usr/bin/env python

import numpy as np
import vtkmodules.all as vtk
import igl

# noinspection PyUnresolvedReferences
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonCore import vtkIdTypeArray
from vtkmodules.vtkCommonDataModel import (
    vtkSelection,
    vtkSelectionNode,
    vtkUnstructuredGrid
)
from vtkmodules.vtkFiltersCore import vtkTriangleFilter
from vtkmodules.vtkFiltersExtraction import vtkExtractSelection
from vtkmodules.vtkFiltersSources import vtkPlaneSource
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkCellPicker,
    vtkDataSetMapper,
    vtkPolyDataMapper,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkRenderer
)


# Catch mouse events
class MouseInteractorStyle(vtkInteractorStyleTrackballCamera):
    def __init__(self, data):
        self.AddObserver('LeftButtonPressEvent', self.left_button_press_event)
        self.data = data
        self.selected_mapper = vtkDataSetMapper()
        self.selected_actor = vtkActor()

    def left_button_press_event(self, obj, event):
        colors = vtkNamedColors()

        # Get the location of the click (in window coordinates)
        pos = self.GetInteractor().GetEventPosition()

        picker = vtkCellPicker()
        picker.SetTolerance(0.0005)

        # Pick from this location.
        picker.Pick(pos[0], pos[1], 0, self.GetDefaultRenderer())

        world_position = picker.GetPickPosition()
        print(f'Cell id is: {picker.GetCellId()}')

        if picker.GetCellId() != -1:
            print(f'Pick position is: ({world_position[0]:.6g}, {world_position[1]:.6g}, {world_position[2]:.6g})')

            ids = vtkIdTypeArray()
            ids.SetNumberOfComponents(1)
            ids.InsertNextValue(picker.GetCellId())

            selection_node = vtkSelectionNode()
            selection_node.SetFieldType(vtkSelectionNode.CELL)
            selection_node.SetContentType(vtkSelectionNode.INDICES)
            selection_node.SetSelectionList(ids)

            selection = vtkSelection()
            selection.AddNode(selection_node)

            extract_selection = vtkExtractSelection()
            extract_selection.SetInputData(0, self.data)
            extract_selection.SetInputData(1, selection)
            extract_selection.Update()

            # In selection
            selected = vtkUnstructuredGrid()
            selected.ShallowCopy(extract_selection.GetOutput())

            print(f'Number of points in the selection: {selected.GetNumberOfPoints()}')
            print(f'Number of cells in the selection : {selected.GetNumberOfCells()}')

            self.selected_mapper.SetInputData(selected)
            self.selected_actor.SetMapper(self.selected_mapper)
            self.selected_actor.GetProperty().EdgeVisibilityOn()
            self.selected_actor.GetProperty().SetColor(colors.GetColor3d('Tomato'))

            self.selected_actor.GetProperty().SetLineWidth(3)

            self.GetInteractor().GetRenderWindow().GetRenderers().GetFirstRenderer().AddActor(self.selected_actor)

        # Forward events
        self.OnLeftButtonDown()


def main(argv):
    colors = vtkNamedColors()

    plane_source = vtkPlaneSource()
    plane_source.Update()

    sv, sf = igl.read_triangle_mesh("data/tuningfork.stl")

    # Erstellen Sie ein VTK-PolyData-Objekt für die Dreiecke
    triangle_polydata = vtk.vtkPolyData()

    # Erstellen Sie Punkte und Zellen für die Dreiecke
    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()

    for i, triangle in enumerate(sf):

        # Extrahieren Sie die Punkte des Dreiecks aus dem STL-Mesh
        point1 = sv[triangle[0]]
        point2 = sv[triangle[1]]
        point3 = sv[triangle[2]]

        # Fügen Sie die Punkte zur Punktwolke hinzu
        point_id1 = points.InsertNextPoint(point1)
        point_id2 = points.InsertNextPoint(point2)
        point_id3 = points.InsertNextPoint(point3)

        # Erstellen Sie ein Dreieck und fügen Sie es zur Zellenliste hinzu
        triangle = vtk.vtkTriangle()
        triangle.GetPointIds().SetId(0, point_id1)
        triangle.GetPointIds().SetId(1, point_id2)
        triangle.GetPointIds().SetId(2, point_id3)
        cells.InsertNextCell(triangle)

    # Fügen Sie die Punkte und Zellen zur PolyData hinzu
    triangle_polydata.SetPoints(points)
    triangle_polydata.SetPolys(cells)

    tgroupping = vtk.vtkFloatArray()
    tgroupping.SetNumberOfComponents(1)
    tgroupping.SetName("group")

    for i, label in enumerate(sf):
        tgroupping.InsertNextTuple1(i)
    
    triangle_polydata.GetCellData().SetScalars(tgroupping)

    triangle_filter = vtkTriangleFilter()
    triangle_filter.SetInputData(triangle_polydata)
    triangle_filter.Update()

    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(triangle_filter.GetOutputPort())
    mapper.ScalarVisibilityOn()
    mapper.SelectColorArray('group')
    mapper.SetScalarModeToUsePointFieldData()
    mapper.SetColorModeToMapScalars()
    mapper.InterpolateScalarsBeforeMappingOn()

    actor = vtkActor()
    #actor.GetProperty().SetColor(colors.GetColor3d('SeaGreen'))
    actor.SetMapper(mapper)

    renderer = vtkRenderer()
    ren_win = vtkRenderWindow()
    ren_win.AddRenderer(renderer)
    ren_win.SetWindowName('CellPicking')
    iren = vtkRenderWindowInteractor()
    iren.SetRenderWindow(ren_win)

    renderer.AddActor(actor)
    # renderer.ResetCamera()
    renderer.SetBackground(colors.GetColor3d('PaleTurquoise'))

    # Add the custom style.
    style = MouseInteractorStyle(triangle_filter.GetOutput())
    style.SetDefaultRenderer(renderer)
    iren.SetInteractorStyle(style)

    ren_win.Render()
    iren.Initialize()
    iren.Start()


if __name__ == '__main__':
    import sys

    main(sys.argv)