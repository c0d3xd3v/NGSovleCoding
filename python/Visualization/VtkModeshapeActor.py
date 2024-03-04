import vtkmodules.all as vtk
from Visualization.vtkhelper import *
import numpy as np
import math
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk

class VtkModeShapeActor():
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
        self.function_name = "eigenmode8"

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

        self.prgfilter = vtk.vtkProgrammableFilter()
        self.prgfilter.SetInputData(self.normals.GetOutput())
        self.time = 0
        def animateFilter():
            input_data = self.prgfilter.GetInputDataObject(0, 0)
            output_data = self.prgfilter.GetOutputDataObject(0)

            points = vtk_to_numpy(input_data.GetPoints().GetData())
            real = vtk_to_numpy(input_data.GetPointData().GetArray(self.function_name))
            #imag = vtk_to_numpy(input_data.GetPointData().GetArray("eigenmode1"))

            self.time = self.time + 1
            t = self.time / 100.0
            self.time = self.time % 100

            amp = 0.5
            result = points + (real*np.cos(t*2.0*math.pi))*amp;
            #result = real*np.cos(t*2.0*math.pi) + imag*np.sin(t*2.0*math.pi)
            result_array_vtk = numpy_to_vtk(result, deep=True)
            output_data.GetPoints().SetData(result_array_vtk)

        self.prgfilter.SetExecuteMethod(animateFilter)
        self.prgfilter.Update()

        #self.mapper = self.setupMapper(self.prgfilter.GetOutput())
        self.mapper = self.setupMapper(self.normals.GetOutput())
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
        self.actor.GetProperty().EdgeVisibilityOn()

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
        normals.ComputePointNormalsOff()
        normals.Update()

        return normals

    def select_function(self, name):
        self.function_name = name
        self.mapper.SelectColorArray(self.function_name)
        self.mapper.MapDataArrayToVertexAttribute("real", self.function_name, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS)
        real = vtk_to_numpy(self.dataset.GetPointData().GetArray(self.function_name))
        print(np.min(real), np.max(real))
        self.mapper.SetScalarRange([0., np.max(real)])

    def setupMapper(self, dataset):
        mapper = vtk.vtkOpenGLPolyDataMapper()
        mapper.SetInputData(dataset)
        mapper.SelectColorArray(self.function_name)
        mapper.MapDataArrayToVertexAttribute("real", self.function_name, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS)
        #mapper.MapDataArrayToVertexAttribute("imag", "eigenmode1", vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS)
        #mapper.ScalarVisibilityOn()
        mapper.SetScalarModeToUsePointFieldData()
        mapper.InterpolateScalarsBeforeMappingOn()
        mapper.SetScalarRange([0., 1.])
        mapper.SetLookupTable(self.lut)

        return mapper

    def setupActor(self, mapper):
        actor = vtk.vtkActor()
        sp = actor.GetShaderProperty()
        uniforms = sp.GetVertexCustomUniforms()
        uniforms.SetUniformf("time_value", 0.0)

        actor.GetShaderProperty().AddVertexShaderReplacement(
            "//VTK::System::Dec\n",
            True,
            "//VTK::System::Dec\n"
            "in vec3 real;\n"
            "in vec3 imag;\n",
            False
        )

        actor.GetShaderProperty().AddVertexShaderReplacement(
            "//VTK::ValuePass::Impl",
            True, # before the standard replacements
            "//VTK::ValuePass::Impl\n" # we still want the default
            "float pi = 3.14159;"
            "float scale = 15.;"
            "vec3 d = real*sin(2.*pi*time_value);"
            "vec3 r = vertexMC.xyz;\n"
            "vertexVCVSOutput = MCVCMatrix * vec4(r + 5.0*d, 1.0);"
            "gl_Position = MCDCMatrix * vec4(r + 5.0*d, 1.0);\n",
            False # only do it once
        )
        actor.GetProperty().EdgeVisibilityOn()
        actor.SetMapper(mapper)
        return actor

    def GetCenter(self):
        return self.center

    def setTimer(self, time):
        #sp = self.actor.GetShaderProperty()
        #uniforms = sp.GetVertexCustomUniforms()
        #uniforms.SetUniformf("time_value", 0.)
        self.prgfilter.Modified()
        self.prgfilter.Update()

