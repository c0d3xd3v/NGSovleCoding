#ifndef QMLPANEMESHINTERFACE_H
#define QMLPANEMESHINTERFACE_H

#include <QObject>
#include <QThread>

class MeshRenderController;
class MeshingWorker;

class QmlPaneMeshInterface : public QObject
{
    Q_OBJECT
private:
    MeshRenderController* controller;
    MeshingWorker *worker;
public slots:
    void handleResults();
public:
    QmlPaneMeshInterface();
    QmlPaneMeshInterface(QmlPaneMeshInterface &controller);
    void set(MeshRenderController* controller);

    Q_INVOKABLE void doMeshing(double stop_energy, double rel_edge_length, double rel_eps);
    Q_INVOKABLE QString getHashString();
    Q_INVOKABLE void setVisualization(QString mode);
};

#endif // QMLPANEMESHINTERFACE_H
