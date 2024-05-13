import vtkmodules.all as vtk
from Visualization.vtkhelper import *
import numpy as np
import math
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk

class VtkFEMActor():
    def __init__(self):
        self.dataset = None
        self.center = (0., 0., 0.)
        self.clipPlane = vtk.vtkPlane()
        self.gridToPolyData = vtk.vtkGeometryFilter()
        self.normals = vtk.vtkPolyDataNormals()
        self.colorTransferFunction, self.lut = create_color_transfer_function_from_xml("Visualization/color_theme.xml")
        self.colorTransferFunction.SetVectorModeToMagnitude()
        self.mapper = vtk.vtkOpenGLPolyDataMapper()
        self.actor = vtk.vtkActor()
        self.function_name = "eigenmode0"
        self.texture_image_reader = None
        self.clipper = None
        self.texture = None
        self.tcoords = None
        self.prgfilter = None
        self.pdata = None

    def setDataset(self, dataset):
        self.dataset = dataset
        self.center = dataset.GetCenter()
        xnorm = [1.0, 0.0, 0.0]
        self.clipPlane = vtk.vtkPlane()
        self.clipPlane.SetOrigin(self.center)
        self.clipPlane.SetNormal(xnorm)
        self.clipper = self.setupClipping(dataset, self.center)
        self.gridToPolyData = self.toPolyData(self.clipper.GetOutput())
        self.normals = self.setupNormals(self.gridToPolyData.GetOutputPort())
        self.pdata = self.normals.GetOutput()

        self.prgfilter = vtk.vtkProgrammableFilter()
        self.prgfilter.SetInputData(self.pdata)
        self.time = 0
        def animateFilter():
            input_data = self.prgfilter.GetInputDataObject(0, 0)
            output_data = self.prgfilter.GetOutputDataObject(0)

            points = vtk_to_numpy(input_data.GetPoints().GetData())
            real = vtk_to_numpy(input_data.GetPointData().GetArray(self.function_name))
            #imag = vtk_to_numpy(input_data.GetPointData().GetArray("eigenmode1"))
            self.time = self.time + 1
            t = self.time / 10.0
            self.time = self.time % 10
            amp = 0.1
            result = points + real*np.cos(t*2.0*math.pi)*amp
            #result = real*np.cos(t*2.0*math.pi) + imag*np.sin(t*2.0*math.pi)
            result_array_vtk = numpy_to_vtk(result, deep=True)
            output_data.GetPoints().SetData(result_array_vtk)

        self.prgfilter.SetExecuteMethod(animateFilter)
        self.prgfilter.Update()

        self.mapper = self.setupMapper(self.prgfilter.GetOutput())
        #self.mapper = self.setupMapper(self.pdata)
        self.actor = self.setupActor(self.mapper)

        self.disableClipping()

    def enableClipping(self):
        self.gridToPolyData.SetInputData(self.clipper.GetOutput())
        self.gridToPolyData.Update()
        self.normals.Update()

    def disableClipping(self):
        self.gridToPolyData.SetInputData(self.dataset)
        self.gridToPolyData.Update()
        self.normals.Update()

    def enableEdges(self):
        p = self.actor.GetProperty() # vtk.vtkProperty()

        colors = vtk.vtkNamedColors()
        edgeColor = colors.GetColor3d("Black")
        p.SetEdgeColor(edgeColor)
        #self.mapper.SetResolveCoincidentTopologyToPolygonOffset()
        self.mapper.SetResolveCoincidentTopologyToPolygonOffset()
        self.mapper.SetResolveCoincidentTopologyLineOffsetParameters(100.0, 10.0)
        p.EdgeVisibilityOn()
        #p.SetEdgeOpacity(1.0)
        #self.mapper.SetColorModeToDefault()

    def disableEdges(self):
        self.actor.GetProperty().EdgeVisibilityOff()

    def setupClipping(self, dataset, center):
        clipper = vtk.vtkClipDataSet()
        clipper.SetInputData(dataset)
        clipper.SetClipFunction(self.clipPlane)
        clipper.SetValue(0.0)
        clipper.GenerateClippedOutputOff()
        clipper.GenerateClipScalarsOff()
        clipper.Update()

        return clipper

    def toPolyData(self, datasetport):
        gridToPolyData = vtk.vtkGeometryFilter()
        gridToPolyData.SetInputData(datasetport)
        gridToPolyData.Update()

        return gridToPolyData

    def setupNormals(self, datasetport):
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputConnection(datasetport)
        normals.AutoOrientNormalsOn()
        normals.ConsistencyOn()
        normals.ComputePointNormalsOn()
        normals.SplittingOn()
        #normals.ComputeCellNormalsOn()
        normals.Update()

        return normals

    def select_function(self, name):
        self.function_name = name
        self.mapper.SelectColorArray(self.function_name)
        real = vtk_to_numpy(self.dataset.GetPointData().GetArray(self.function_name))
        print("range : ", np.min(real), np.max(real))
        self.mapper.SetScalarRange([0., np.max(real)])

    def setupMapper(self, dataset):
        mapper = vtk.vtkOpenGLPolyDataMapper()
        mapper.SetInputData(self.pdata)
        mapper.SetScalarModeToDefault()
        #self.pdata = self.setupMatCapsTextureCoords(self.pdata)
        mapper.SelectColorArray(self.function_name)
        mapper.ScalarVisibilityOn()
        mapper.SetScalarModeToUsePointFieldData()
        mapper.InterpolateScalarsBeforeMappingOn()
        mapper.SetScalarRange([0., 1.])
        mapper.SetLookupTable(self.lut)

        return mapper

    def setupActor(self, mapper):
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        prop: vtk.vtkProperty = actor.GetProperty()
        prop.SetInterpolationToFlat()
        prop.SetMetallic(1.0)
        prop.SetRoughness(0.2)
        prop.SetDiffuse(0.9)
        prop.SetBaseIOR(0.5)
        prop.SetCoatStrength(0.2)
        prop.SetCoatIOR(0.5)
        white = vtk.vtkNamedColors().GetColor3d('White')
        prop.SetColor(white)
        bgcolor = vtk.vtkNamedColors().HTMLColorToRGB("#363737")
        bgcolor = [bgcolor[0]/255., bgcolor[1]/255., bgcolor[2]/255.]
        prop.SetEdgeColor(white)
        prop.SetCoatRoughness(0.0)
        prop.SetCoatColor(vtk.vtkNamedColors().GetColor3d('White'))
        return actor

    def GetCenter(self):
        return self.center

    def setTimer(self, time):
        self.normals.Modified()
        self.normals.Update()
        self.prgfilter.Modified()
        self.prgfilter.Update()

        self.mapper.Modified()
        self.mapper.Update()
