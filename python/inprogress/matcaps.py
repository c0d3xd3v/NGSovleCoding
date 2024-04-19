import sys
import vtk

# Load the STL mesh
reader = vtk.vtkSTLReader()
reader.SetFileName(sys.argv[1])
reader.Update()

polydata = reader.GetOutput()

normals = vtk.vtkPPolyDataNormals()
normals.SetInputData(polydata)
normals.ComputePointNormalsOn()
normals.Update()

polydata = normals.GetOutput()

tcoords = vtk.vtkFloatArray()
tcoords.SetNumberOfComponents(2)

pointdata = polydata.GetPointData()
normals = polydata.GetPointData().GetNormals()
'''
perlinNoise = vtk.vtkPerlinNoise()
perlinNoise.SetFrequency(200, 150.25, 150.5)
perlinNoise.SetPhase(0, 10, 0)
perlinNoise.SetAmplitude(0.01)
'''
bounds = polydata.GetBounds()
print(bounds)
for i in range(polydata.GetNumberOfPoints()):
    p = polydata.GetPoint(i)
    n = normals.GetTuple3(i)
    nx = p[0]*0.005
    ny = n[1]*0.005
    #pn1 = perlinNoise.EvaluateFunction(p[0], p[1], p[2])
    #pn2 = perlinNoise.EvaluateFunction(p[2], p[1], p[0])
    tcoords.InsertNextTuple2(nx, ny)
    #print(pn1, pn2)

polydata.GetPointData().SetTCoords(tcoords)
'''
smoother = vtk.vtkAttributeSmoothingFilter()
smoother.SetInputData(polydata)
smoother.SetRelaxationFactor(0.1)
smoother.SetNumberOfIterations(4)
smoother.Update()
'''
# Create a mapper
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(polydata)

# Create an actor
actor = vtk.vtkActor()
actor.SetMapper(mapper)
#actor.GetProperty().SetInterpolationToPBR()

texture_image_reader = vtk.vtkPNGReader()
texture_image_reader.SetFileName("/home/kai/chrome2.png")  # Replace "texture.jpg" with the path to your texture image
texture_image_reader.Update()

# Create a texture object
texture = vtk.vtkTexture()
texture.SetInputConnection(texture_image_reader.GetOutputPort())  # Assuming you have already loaded the texture image

# Set texture properties
texture.InterpolateOff()

actor.SetTexture(texture)
actor.GetShaderProperty().AddVertexShaderReplacement(
            "//VTK::Normal::Dec",
            True,
            "//VTK::Normal::Dec\n",
            False
        )
actor.GetShaderProperty().AddVertexShaderReplacement(
            "//VTK::ValuePass::Impl",
            True,
            "//VTK::ValuePass::Impl\n"
            "vec3 normalDir = normalize(normalMatrix*normalMC);\n"
            "tcoordVCVSOutput =  0.5*(1. + normalDir.xy);\n",
            False
        )
actor.GetShaderProperty().AddFragmentShaderReplacement(
              "//VTK::Normal::Dec",
              True,
              "//VTK::Normal::Dec\n",
              False
            )
actor.GetShaderProperty().AddFragmentShaderReplacement(
              "//VTK::Light::Impl",
              True,
              "//VTK::Light::Impl\n"
              "fragOutput0 = vec4(fragOutput0.rgb, fragOutput0.a);\n",
              False
          )

actor.GetProperty().SetAmbient(1.0)
actor.GetProperty().SetDiffuse(0.0)
# Render the scene
renderer = vtk.vtkRenderer()

bgcolor = vtk.vtkNamedColors().HTMLColorToRGB("#363737")
bgcolor = [bgcolor[0]/255., bgcolor[1]/255., bgcolor[2]/255.]

renderer.SetBackground(bgcolor)
renderer.AddActor(actor)

render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)

render_window_interactor = vtk.vtkRenderWindowInteractor()
render_window_interactor.SetRenderWindow(render_window)

render_window.Render()
render_window_interactor.Start()
