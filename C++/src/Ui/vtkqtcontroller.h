#ifndef VTKQTCONTROLLER_H
#define VTKQTCONTROLLER_H

#include <map>
#include <tuple>

#include <QObject>

#include <vtk/vtkRenderer.h>
#include <vtk/vtkLight.h>
#include <vtk/vtkRenderWindow.h>
#include <vtk/vtkRenderWindowInteractor.h>
#include <vtk/vtkOrientationMarkerWidget.h>
#include <vtk/QQuickVTKRenderItem.h>

#include "Ui/meshrendercontroller.h"
#include "Ui/cellpicking.h"

class VtkQtController : public QObject
{
    Q_OBJECT
private:
    vtkSmartPointer<MouseInteractorStyle> selectionStyle;
    vtkInteractorObserver* defaultInteractionStyle;
    QQuickVTKRenderItem* qquickvtkItem;
    vtkRenderer* renderer;
    vtkRenderWindow *renderWindow;
    vtkRenderWindowInteractor *interactor;
    vtkSmartPointer<vtkOrientationMarkerWidget> om;
    std::map<std::string, MeshRenderController *> meshRenderController;

public:
    explicit VtkQtController(QObject *parent = nullptr);
    void init(QQuickVTKRenderItem* qquickvtkItem,
              vtkRenderer* renderer,
              vtkRenderWindow *renderWindow,
              vtkRenderWindowInteractor *interactor);

public slots:
    void resize(int width, int height);
    QString loadFile(QString filepath);
    void removeObject(QString hashString);

signals:

};

#endif // VTKQTCONTROLLER_H
