#include <QDebug>

#include "Ui/meshrendercontroller.h"
#include "Ui/meshingworker.h"

#include "qmlpanemeshinterface.h"

QmlPaneMeshInterface::QmlPaneMeshInterface()
    : controller(nullptr)
{
}

QmlPaneMeshInterface::QmlPaneMeshInterface(QmlPaneMeshInterface &controller)
{
    this->controller = controller.controller;
}

void QmlPaneMeshInterface::set(MeshRenderController *controller)
{
    this->controller = controller;
}

void QmlPaneMeshInterface::doMeshing(double stop_energy, double rel_edge_length, double rel_eps)
{
    if(controller != nullptr)
    {
        MeshingWorker *workerThread = new MeshingWorker(controller, stop_energy, rel_edge_length, rel_eps);
        connect(workerThread, &MeshingWorker::resultReady, this, &QmlPaneMeshInterface::handleResults);
        connect(workerThread, &MeshingWorker::finished, workerThread, &QObject::deleteLater);
        workerThread->start();
    }
}

void QmlPaneMeshInterface::handleResults()
{
    qDebug() << "meshing ready";
}

QString QmlPaneMeshInterface::getHashString()
{
    std::stringstream ss;
    if(controller != nullptr)
        ss << controller->getActor();
    return QString(ss.str().c_str());
}

void QmlPaneMeshInterface::setVisualization(QString mode)
{
    if(mode == "Surface")
    {
        controller->renderSurface();
    }else if(mode == "Surface edges")
    {
        controller->renderSurfaceWithEdges();
    }
}
