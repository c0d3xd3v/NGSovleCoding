#ifndef VTKQTCONTROLLER_H
#define VTKQTCONTROLLER_H

#include <map>
#include <tuple>

#include <QObject>

#include <vtkRenderer.h>
#include <vtkLight.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOrientationMarkerWidget.h>
#include <QQuickVTKRenderItem.h>

#include "Ui/meshingcontroller.h"
#include "Ui/qmlpanemeshinterface.h"
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
    //std::map<std::string, MeshRenderController *> meshRenderController;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    explicit VtkQtController(QObject *parent = nullptr);
    void init(QQuickVTKRenderItem* qquickvtkItem,
              vtkRenderer* renderer,
              vtkRenderWindow *renderWindow,
              vtkRenderWindowInteractor *interactor);

public slots:
    void resize(int width, int height);
    QmlPaneMeshInterface *loadFile(QString filepath);
    void removeObject(QString hashString, QmlPaneMeshInterface* pane);

signals:

};

#endif // VTKQTCONTROLLER_H
