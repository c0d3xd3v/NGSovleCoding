
import numpy as np
import math
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from Visualization.vtkhelper import *


class VtkPressureFieldActor():
    def __init__(self):
        self.dataset = None
        self.timestepper = None
        self.prgfilter = None
        self.center = (0., 0., 0.)
        self.clipPlane = vtk.vtkPlane()
        self.gridToPolyData = vtk.vtkGeometryFilter()
        self.colorTransferFunction, self.lut = create_color_transfer_function_from_xml("Visualization/color_theme.xml")
        self.colorTransferFunction.SetVectorModeToMagnitude()
        self.mapper = vtk.vtkOpenGLPolyDataMapper()
        self.actor = vtk.vtkActor()


    def setDataset(self, dataset):
        self.dataset = dataset
        self.center = dataset.GetCenter()
        xnorm = [0.0, 1.0, 0.0]
    
        self.clipPlane = vtk.vtkPlane()
        self.clipPlane.SetOrigin(self.center)
        self.clipPlane.SetNormal(xnorm)
        self.clipper = self.setupClipping(self.dataset, self.center)
        self.gridToPolyData = self.toPolyData(self.clipper.GetOutput())

        self.prgfilter = vtk.vtkProgrammableFilter()
        self.prgfilter.SetInputData(self.gridToPolyData.GetOutput())
        self.time = 0
        def animateFilter():
            input_data = self.prgfilter.GetInputDataObject(0, 0)
            output_data = self.prgfilter.GetOutputDataObject(0)

            real = vtk_to_numpy(input_data.GetPointData().GetArray("pressurefield0Real"))
            imag = vtk_to_numpy(input_data.GetPointData().GetArray("pressurefield0Imag"))
 
            self.time = self.time + 1
            t = self.time / 100.0
            self.time = self.time % 100

            result = real*np.cos(t*2.0*math.pi) + imag*np.sin(t*2.0*math.pi)
            result_array_vtk = numpy_to_vtk(result, deep=True)
            result_array_vtk.SetName("Result")
            output_data.GetPointData().AddArray(result_array_vtk)
        self.prgfilter.SetExecuteMethod(animateFilter)
        self.prgfilter.Update()
        
        self.mapper = self.setupMapper(self.prgfilter.GetOutput())
        self.actor = self.setupActor(self.mapper)

    def enableClipping(self):
        self.gridToPolyData.SetInputData(self.clipper.GetOutput())
        self.gridToPolyData.Update()

    def disableClipping(self):
        self.gridToPolyData.SetInputData(self.dataset)
        self.gridToPolyData.Update()

    def toPolyData(self, datasetport):
        gridToPolyData = vtk.vtkGeometryFilter()
        gridToPolyData.SetInputData(datasetport)
        gridToPolyData.Update()

        return gridToPolyData

    def setupClipping(self, dataset, center):
        clipper = vtk.vtkClipDataSet()
        clipper.SetInputData(dataset)
        clipper.SetClipFunction(self.clipPlane)
        clipper.SetValue(0.0)
        clipper.GenerateClippedOutputOn()
        clipper.GenerateClipScalarsOn()
        clipper.Update()

        return clipper

    def setupMapper(self, dataset):
        mapper = vtk.vtkOpenGLPolyDataMapper()
        mapper.SetInputData(dataset)
        mapper.ScalarVisibilityOn()
        mapper.SelectColorArray('Result')
        mapper.SetScalarModeToUsePointFieldData()
        mapper.SetColorModeToMapScalars()
        mapper.InterpolateScalarsBeforeMappingOn()
        sr = 0.00001
        mapper.SetScalarRange([-sr, sr])
        mapper.SetLookupTable(self.lut)

        return mapper

    def setupActor(self, mapper):
        actor = vtk.vtkActor()
        actor.GetProperty().EdgeVisibilityOff()
        actor.SetMapper(mapper)
        return actor

    def setTimer(self, time):
        self.prgfilter.Modified()
        self.prgfilter.Update()

