#ifndef QMLPANEMESHINTERFACE_H
#define QMLPANEMESHINTERFACE_H

#include <QObject>
#include <QThread>

class vtkRenderer;
class MeshingController;
class MeshingWorker;

class QmlPaneMeshInterface : public QObject
{
    Q_OBJECT
private:
    vtkRenderer *renderer;
    MeshingController* controller;
    MeshingWorker *worker;
public slots:
    void handleResults();
signals:
    void meshingFinished();
public:
    QmlPaneMeshInterface(vtkRenderer *renderer);
    ~QmlPaneMeshInterface();
    QmlPaneMeshInterface(QmlPaneMeshInterface &controller);
    void set(MeshingController* controller);

    Q_INVOKABLE void doMeshing(double stop_energy, double rel_edge_length, double rel_eps);
    Q_INVOKABLE QString getHashString();
    Q_INVOKABLE void setRepresentation(QString mode);
    Q_INVOKABLE QString getRepresentation();
    Q_INVOKABLE void setVisualization(QString mode);
    Q_INVOKABLE void cleanupRendering();
};

#endif // QMLPANEMESHINTERFACE_H
