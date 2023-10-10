#include <QDebug>
#include <QQuickWindow>
#include <igl/read_triangle_mesh.h>

#include "Ui/vtkorientationindicator.h"

#include "vtkqtcontroller.h"

VtkQtController::VtkQtController(QObject *parent)
    : QObject{parent}
{
    renderer = nullptr;
    renderWindow = nullptr;
    interactor = nullptr;
    om = nullptr;
    qquickvtkItem = nullptr;
    //meshRenderController = nullptr;
    selectionStyle = nullptr;
}

void VtkQtController::init(QQuickVTKRenderItem* qquickvtkItem, vtkRenderer *renderer, vtkRenderWindow *renderWindow, vtkRenderWindowInteractor *interactor)
{
    if(this->renderer == nullptr &&
       this->renderWindow  == nullptr &&
       this->interactor == nullptr &&
       this->qquickvtkItem == nullptr)
    {
        this->qquickvtkItem = qquickvtkItem;
        this->renderer = renderer;
        this->renderWindow = renderWindow;
        this->interactor = interactor;

        defaultInteractionStyle = this->interactor->GetInteractorStyle();

        auto labels{"sal"};
        vtkNew<vtkNamedColors> colors;
        colors->SetColor("ParaViewBkg", std::array<unsigned char, 4>{82, 87, 110, 255}.data());
        vtkSmartPointer<vtkPropAssembly> axes = MakeCubeActor(labels, colors);
        om = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
        om->SetOrientationMarker(axes);
        om->SetViewport(0, 0, 0.2, 0.2);
        om->SetInteractor(interactor);
        om->EnabledOn();
/*
        vtkNew<vtkLight> light1;
        light1->SetFocalPoint(0, 0, 0);
        light1->SetPosition(0, 1, 0.2);
        light1->SetColor(colors->GetColor3d("HighNoonSun").GetData());
        light1->SetIntensity(0.3);
        renderer->AddLight(light1);

        vtkNew<vtkLight> light2;
        light2->SetFocalPoint(0, 0, 0);
        light2->SetPosition(1.0, 1.0, 1.0);
        light2->SetColor(colors->GetColor3d("100W Tungsten").GetData());
        light2->SetIntensity(0.8);
        renderer->AddLight(light2);
*/
        renderWindow->SetMultiSamples(16);
        renderer->ResetCamera();
        renderer->SetBackground(0.5, 0.5, 0.7);
        renderer->SetBackground2(0.7, 0.7, 0.7);
        renderer->SetGradientBackground(true);

        selectionStyle = vtkSmartPointer<MouseInteractorStyle>::New();
        selectionStyle->SetDefaultRenderer(renderer);
        interactor->SetInteractorStyle(selectionStyle);
    }
}

void VtkQtController::resize(int width, int height)
{
    //qDebug() << width << ", " << height;
    if(om != nullptr)
    {
        om->SetViewport(0.0, 0.0, 100.0/float(width), 100.0/float(height));
        om->Modified();
        renderer->Modified();
        renderWindow->Modified();
        interactor->Modified();
        qquickvtkItem->update();
        if(qquickvtkItem->window())
        qquickvtkItem->window()->update();
    }
}

QString VtkQtController::loadFile(QString filepath)
{
    QString hashString;
    const QUrl url(filepath);
    if (url.isLocalFile())
    {
        std::string path = url.path().toStdString();

        const std::string oldLocale = std::setlocale(LC_NUMERIC, nullptr);
        std::setlocale(LC_NUMERIC, "C");

        Eigen::MatrixXf V;
        Eigen::MatrixXi F;
        if(igl::read_triangle_mesh(path, V, F))
        {
            /*
            if(meshRenderController != nullptr)
            {
                renderer->RemoveActor(meshRenderController->getActor());
                delete meshRenderController;
            }
            */
            MeshRenderController* mc = new MeshRenderController(V, F);

            std::stringstream ss;
            ss << mc->getActor();

            typedef std::pair<std::string, MeshRenderController*> pair;
            meshRenderController.insert(pair(ss.str(), mc));

            hashString = ss.str().c_str();

            renderer->AddActor(mc->getActor());
            renderer->ResetCamera();
        }

        std::setlocale(LC_NUMERIC, oldLocale.c_str());
    }

    return hashString;
}

void VtkQtController::removeObject(QString hashString)
{
    typedef std::map<std::string, MeshRenderController*>::iterator Iterator;
    Iterator itr = meshRenderController.find(hashString.toStdString());
    if(itr != meshRenderController.end())
    {
        qDebug() << "Delete " << itr->first.c_str();
        renderer->RemoveActor(itr->second->getActor());
        delete itr->second;
        meshRenderController.erase(itr);
    }
}
