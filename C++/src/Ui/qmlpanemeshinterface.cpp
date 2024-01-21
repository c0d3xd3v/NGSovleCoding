#include <QDebug>

#include <vtk/vtkRenderer.h>

#include "Ui/meshingcontroller.h"
#include "Ui/meshingworker.h"

#include "qmlpanemeshinterface.h"

QmlPaneMeshInterface::QmlPaneMeshInterface(vtkRenderer* renderer)
    : controller(nullptr), renderer(renderer)
{
}

QmlPaneMeshInterface::~QmlPaneMeshInterface()
{
}

void QmlPaneMeshInterface::cleanupRendering()
{
    controller->setCurrentMeshToResult();
    renderer->RemoveActor(controller->getActor());

    controller->setCurrentMeshToSource();
    renderer->RemoveActor(controller->getActor());

    setParent(nullptr);
    deleteLater();
}

QmlPaneMeshInterface::QmlPaneMeshInterface(QmlPaneMeshInterface &controller)
{
    this->controller = controller.controller;
}

void QmlPaneMeshInterface::set(MeshingController *controller)
{
    this->controller = controller;
}

void QmlPaneMeshInterface::doMeshing(double stop_energy, double rel_edge_length, double rel_eps)
{
    if(controller != nullptr)
    {
        controller->setCurrentMeshToResult();
        if(controller->getActor() != nullptr)
            renderer->RemoveActor(controller->getActor());
        controller->setCurrentMeshToSource();
        controller->showCurrentMesh();
        MeshingWorker *workerThread = new MeshingWorker(controller, stop_energy, rel_edge_length, rel_eps);
        connect(workerThread, &MeshingWorker::resultReady, this, &QmlPaneMeshInterface::handleResults);
        connect(workerThread, &MeshingWorker::finished, workerThread, &QObject::deleteLater);
        workerThread->start();
    }
}

void QmlPaneMeshInterface::handleResults()
{
    qDebug() << "meshing ready";
    controller->setCurrentMeshToSource();
    controller->hideCurrentMesh();
    controller->setCurrentMeshToResult();
    renderer->AddActor(controller->getActor());

    emit meshingFinished();
}

QString QmlPaneMeshInterface::getHashString()
{
    std::stringstream ss;
    if(controller != nullptr)
        ss << controller->getActor();
    return QString(ss.str().c_str());
}

void QmlPaneMeshInterface::setRepresentation(QString mode)
{
    if(controller != nullptr)
        if(mode == "Surface")
        {
            controller->renderSurface();
        }else if(mode == "Surface edges")
        {
            controller->renderSurfaceWithEdges();
        }
}

QString QmlPaneMeshInterface::getRepresentation()
{
    if(controller != nullptr)
        if(controller->getRenderSurfaceWidthEdges())
            return "Surface edges";
    return "Surface";
}

void QmlPaneMeshInterface::setVisualization(QString mode)
{
    if(controller != nullptr)
        if(mode == "Source Mesh")
        {
            controller->setCurrentMeshToResult();
            controller->hideCurrentMesh();
            controller->setCurrentMeshToSource();
            controller->showCurrentMesh();
        }
        else if(mode == "Result Mesh")
        {
            controller->setCurrentMeshToSource();
            controller->hideCurrentMesh();
            controller->setCurrentMeshToResult();
            controller->showCurrentMesh();
        }
}
