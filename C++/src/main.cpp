#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <QQuickWindow>

#include <vtk/vtkRenderWindow.h>
#include <vtk/vtkRenderWindowInteractor.h>
#include <vtk/vtkOrientationMarkerWidget.h>

#include <vtk/QQuickVTKRenderItem.h>

#include <Ui/vtkorientationindicator.h>

int main(int argc, char *argv[])
{
    QQuickVTKRenderWindow::setupGraphicsBackend();
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
    QGuiApplication app(argc, argv);

    QQmlApplicationEngine engine;
    engine.addImportPath("/usr/lib/qml/");
    engine.load(QUrl("qrc:main.qml"));

    QObject* topLevel = engine.rootObjects().value(0);
    QQuickWindow* window = qobject_cast<QQuickWindow*>(topLevel);
    window->show();

    QQuickVTKRenderItem* qquickvtkItem = topLevel->findChild<QQuickVTKRenderItem*>("ConeView");
    vtkRenderWindowInteractor *iRen = qquickvtkItem->renderWindow()->renderWindow()->GetInteractor();

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

    qquickvtkItem->renderer()->ResetCamera();
    qquickvtkItem->renderer()->SetBackground(0.5, 0.5, 0.7);
    qquickvtkItem->renderer()->SetBackground2(0.7, 0.7, 0.7);
    qquickvtkItem->renderer()->SetGradientBackground(true);
    qquickvtkItem->update();

    if (engine.rootObjects().isEmpty())
        return -1;

    return app.exec();
}
