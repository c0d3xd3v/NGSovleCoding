#ifndef VTKQTCONTROLLER_H
#define VTKQTCONTROLLER_H

#include <QObject>

#include <vtk/vtkRenderer.h>
#include <vtk/vtkLight.h>
#include <vtk/vtkRenderWindow.h>
#include <vtk/vtkRenderWindowInteractor.h>
#include <vtk/vtkOrientationMarkerWidget.h>
#include <vtk/QQuickVTKRenderItem.h>

class VtkQtController : public QObject
{
    Q_OBJECT
private:
    QQuickVTKRenderItem* qquickvtkItem;
    vtkRenderer* renderer;
    vtkRenderWindow *renderWindow;
    vtkRenderWindowInteractor *interactor;
    vtkSmartPointer<vtkOrientationMarkerWidget> om;
public:
    explicit VtkQtController(QObject *parent = nullptr);
    void init(QQuickVTKRenderItem* qquickvtkItem,
              vtkRenderer* renderer,
              vtkRenderWindow *renderWindow,
              vtkRenderWindowInteractor *interactor);

public slots:
    void resize(int width, int height);

signals:

};

#endif // VTKQTCONTROLLER_H
