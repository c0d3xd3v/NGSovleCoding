#include <QDebug>

#include "meshingworker.h"

#include "Ui/meshrendercontroller.h"

MeshingWorker::MeshingWorker(MeshRenderController *controller, double stop_energy, double rel_edge_length, double rel_eps) :
    controller(controller),
    stop_energy(stop_energy), rel_edge_length(rel_edge_length), rel_eps(rel_eps)
{
}

void MeshingWorker::run()
{
    controller->domeshing(stop_energy, rel_edge_length, rel_eps);
    emit resultReady();
}
