#ifndef VTKSCENEEVENTPROCESSOR_H
#define VTKSCENEEVENTPROCESSOR_H

#include <QObject>

class VTKSceneEventProcessor : public QObject
{
    Q_OBJECT
public:
    explicit VTKSceneEventProcessor(QObject *parent = nullptr);

signals:

};

#endif // VTKSCENEEVENTPROCESSOR_H
