#include<string>
#include<stdio.h>
#include<stdlib.h>
#include <iostream>

#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <QQuickWindow>

#include <vtk/vtkRenderWindow.h>
#include <vtk/QQuickVTKRenderItem.h>
#include <vtk/vtkRenderWindowInteractor.h>
#include <vtk/vtkOrientationMarkerWidget.h>

#include "vtkorientationindicator.h"
#include "instantrendering.h"
#include "rendercontroller.h"
#include "cellpicking.h"
#include "colormaphelper.h"

int renderTriangleMesh(Eigen::MatrixXd &nodes, Eigen::MatrixXi &tris)
{
    QQuickVTKRenderWindow::setupGraphicsBackend();
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
    int argc = 0;
    char **argv = 0;
    QGuiApplication app(argc, argv);
    QQmlApplicationEngine engine;
    engine.addImportPath("/usr/lib/qml/");
    engine.load(QUrl("qrc:main.qml"));

    QObject* topLevel = engine.rootObjects().value(0);
    QQuickWindow* window = qobject_cast<QQuickWindow*>(topLevel);
    window->show();

    QQuickVTKRenderItem* qquickvtkItem = topLevel->findChild<QQuickVTKRenderItem*>("ConeView");
    vtkRenderWindow *renderWindow = qquickvtkItem->renderWindow()->renderWindow();
    vtkRenderWindowInteractor *iRen = renderWindow->GetInteractor();
    renderWindow->SetMultiSamples(16);

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
    RenderController rc(nodes, tris);
    rc.setColormap(cm);

    // Set the custom stype to use for interaction.
    vtkNew<MouseInteractorStyle> style;
    style->SetDefaultRenderer(qquickvtkItem->renderer());
    style->rc = &rc;
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

