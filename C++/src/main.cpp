#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <QQuickWindow>
#include <QQmlContext>
//#include <QQuickStyle>

#include <vtkRenderWindow.h>
#include <QQuickVTKRenderItem.h>
#include <vtkRenderWindowInteractor.h>

#include "Ui/vtkqtcontroller.h"
#include "Ui/qmlpanemeshinterface.h"
#include "Ui/meshingcontroller.h"

static void registerQmlTypes() {
    qRegisterMetaType<QmlPaneMeshInterface*>("QmlPaneMeshInterface*");
}
Q_COREAPP_STARTUP_FUNCTION(registerQmlTypes)

int main (int argc, char ** argv)
{
    std::cout << "EIGEN_MAJOR_VERSION : " << EIGEN_MAJOR_VERSION << std::endl;
    std::cout << "EIGEN_MINOR_VERSION : " << EIGEN_MINOR_VERSION << std::endl;
    QQuickVTKRenderWindow::setupGraphicsBackend();
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
    //QQuickStyle::setStyle("Fusion");

    QGuiApplication app(argc, argv);

    VtkQtController vtkqtcontroller;
    QQmlApplicationEngine engine;
    engine.addImportPath("/usr/lib/qml/");

    engine.rootContext()->setContextProperty("vtkqtcontroller", &vtkqtcontroller);
    engine.load(QUrl("qrc:/main.qml"));

    QObject* topLevel = engine.rootObjects().value(0);
    QQuickWindow* window = qobject_cast<QQuickWindow*>(topLevel);
    window->show();

    QQuickVTKRenderItem* qquickvtkItem = topLevel->findChild<QQuickVTKRenderItem*>("ConeView");
    vtkRenderWindow *renderWindow = qquickvtkItem->renderWindow()->renderWindow();
    vtkRenderWindowInteractor *iRen = renderWindow->GetInteractor();
    vtkRenderer *renderer = qquickvtkItem->renderer();

    vtkqtcontroller.init(qquickvtkItem, renderer, renderWindow, iRen);

    qquickvtkItem->update();

    return app.exec();
}
