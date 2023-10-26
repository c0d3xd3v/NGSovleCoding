#include<string>
#include<stdio.h>
#include<stdlib.h>
#include <iostream>

#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <QQuickWindow>

#include <vtk/vtkRenderWindow.h>
#include <vtk/vtkLight.h>
#include <vtk/vtkOpenGLRenderer.h>
#include <vtk/QQuickVTKRenderItem.h>
#include <vtk/vtkRenderWindowInteractor.h>
#include <vtk/vtkOrientationMarkerWidget.h>
#include <vtk/vtkCamera.h>
#include <vtk/vtkSequencePass.h>
#include <vtk/vtkCameraPass.h>
#include <vtk/vtkShadowMapBakerPass.h>
#include <vtk/vtkShadowMapPass.h>
#include <vtk/vtkRenderPassCollection.h>

#include "vtkorientationindicator.h"
#include "instantrendering.h"
#include "meshrendercontroller.h"
#include "cellpicking.h"
#include "colormaphelper.h"

int renderTriangleMesh(Eigen::MatrixXf &nodes, Eigen::MatrixXi &tris)
{
    QQuickVTKRenderWindow::setupGraphicsBackend();
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
    int argc = 0;
    char **argv = 0;
    QGuiApplication app(argc, argv);
    QQmlApplicationEngine engine;
    engine.addImportPath("/usr/lib/qml/");
    engine.load(QUrl("qrc:InstantRendering.qml"));

    QObject* topLevel = engine.rootObjects().value(0);
    QQuickWindow* window = qobject_cast<QQuickWindow*>(topLevel);
    window->show();

    QQuickVTKRenderItem* qquickvtkItem = topLevel->findChild<QQuickVTKRenderItem*>("ConeView");
    vtkRenderWindow *renderWindow = qquickvtkItem->renderWindow()->renderWindow();
    vtkRenderWindowInteractor *iRen = renderWindow->GetInteractor();
    renderWindow->SetMultiSamples(8);
/*
    vtkNew<vtkNamedColors> light_colors;
    // Color temp. 5400k.
    light_colors->SetColor("HighNoonSun", 1.0, 1.0, .9843, 1.0);
    // Color temp. 2850k.
    light_colors->SetColor("100W Tungsten", 1.0, .8392, .6667, 1.0);

    vtkNew<vtkLight> light1;
    light1->SetFocalPoint(0, 0, 0);
    light1->SetPosition(0, 1, 0.2);
    light1->SetColor(light_colors->GetColor3d("HighNoonSun").GetData());
    light1->SetIntensity(0.2);
    qquickvtkItem->renderer()->AddLight(light1);

    vtkNew<vtkLight> light2;
    light2->SetFocalPoint(0, 0, 0);
    light2->SetPosition(1.0, 1.0, 1.0);
    light2->SetColor(light_colors->GetColor3d("100W Tungsten").GetData());
    light2->SetIntensity(0.2);
    qquickvtkItem->renderer()->AddLight(light2);

    vtkNew<vtkShadowMapPass> shadows;
    vtkNew<vtkSequencePass> seq;
    vtkNew<vtkRenderPassCollection> passes;
    passes->AddItem(shadows->GetShadowMapBakerPass());
    passes->AddItem(shadows);
    seq->SetPasses(passes);

    vtkNew<vtkCameraPass> cameraP;
    cameraP->SetDelegatePass(seq);

    vtkOpenGLRenderer* glr = dynamic_cast<vtkOpenGLRenderer*>(qquickvtkItem->renderer());
    glr->SetPass(cameraP);
*/
    auto labels{"sal"};
    vtkNew<vtkNamedColors> colors;
    colors->SetColor("ParaViewBkg",
                     std::array<unsigned char, 4>{82, 87, 110, 255}.data());
    auto axes = MakeCubeActor(labels, colors);
    vtkNew<vtkOrientationMarkerWidget> om;
    om->SetOrientationMarker(axes);
    om->SetViewport(0, 0, 0.2, 0.2);
    om->SetInteractor(iRen);
    om->EnabledOn();

    //ColorMap cm = loadColorMapFromXML("color_theme.xml");
    ColorMap cm = loadIdColorMap(0, 6000);
    MeshRenderController rc(nodes, tris);
    rc.setColormap(cm);

    // Set the custom stype to use for interaction.
    vtkNew<MouseInteractorStyle> style;
    style->SetDefaultRenderer(qquickvtkItem->renderer());
    iRen->SetInteractorStyle(style);

    qquickvtkItem->renderer()->AddActor(rc.getActor());
    qquickvtkItem->renderer()->ResetCamera();
    qquickvtkItem->renderer()->SetBackground(0.5, 0.5, 0.7);
    qquickvtkItem->renderer()->SetBackground2(0.7, 0.7, 0.7);
    qquickvtkItem->renderer()->SetGradientBackground(true);
    qquickvtkItem->update();

    int a = app.exec();
    return a;
}
