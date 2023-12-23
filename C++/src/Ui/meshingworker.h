#ifndef MESHINGWORKER_H
#define MESHINGWORKER_H

#include <QThread>

class MeshingController;

class MeshingWorker : public QThread
{
private:
    Q_OBJECT
    void run() override;
    double stop_energy;
    double rel_edge_length;
    double rel_eps;
    MeshingController* controller;
public slots:
signals:
    void resultReady();
public:
    MeshingWorker(MeshingController* controller, double stop_energy, double rel_edge_length, double rel_eps);
};

#endif // MESHINGWORKER_H
