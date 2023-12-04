import vtkmodules.all as vtk
from Visualization.vtkhelper import *


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
        self.mapper = self.setupMapper(self.normals.GetOutput())
        self.actor = self.setupActor(self.mapper)

    def enableClipping(self):
        self.gridToPolyData.SetInputData(self.clipper.GetOutput())
        self.gridToPolyData.Update()
        self.normals.Update()

    def disableClipping(self):
        self.gridToPolyData.SetInputData(self.dataset)
        self.gridToPolyData.Update()
        self.normals.Update()

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
        normals.Update()

        return normals

    def setupMapper(self, dataset):
        mapper = vtk.vtkOpenGLPolyDataMapper()
        mapper.SetInputData(dataset)
        mapper.ScalarVisibilityOn()
        #mapper.SelectColorArray('eigenmode0')
        #mapper.SetScalarModeToUsePointFieldData()
        #mapper.SetColorModeToMapScalars()
        #mapper.InterpolateScalarsBeforeMappingOn()
        mapper.MapDataArrayToVertexAttribute("real", "eigenmode0", vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, 3)
        mapper.MapDataArrayToVertexAttribute("imag", "eigenmode1", vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, 3)
        #sr = 0.1
        #mapper.SetScalarRange([0.0, sr])
        #mapper.SetLookupTable(self.lut)

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
            "vec3 displacement = real*cos(6.28*time_value) + imag*sin(6.28*time_value);"
            "vec3 r = vertexMC.xyz + displacement*3.;\n"
            "vertexVCVSOutput = MCVCMatrix * vec4(r, 1.0);\n"
            "gl_Position = MCDCMatrix * vec4(r, 1.0);\n",
            False # only do it once
        )
        actor.GetProperty().EdgeVisibilityOff()
        actor.SetMapper(mapper)
        return actor

    def GetCenter(self):
        return self.center

    def setTimer(self, time):
        sp = self.actor.GetShaderProperty()
        uniforms = sp.GetVertexCustomUniforms()
        uniforms.SetUniformf("time_value", time)

