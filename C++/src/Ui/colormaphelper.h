#ifndef COLORMAPHELPER_H
#define COLORMAPHELPER_H

#include <string>

#include <vtk/vtkColorTransferFunction.h>
#include <vtk/vtkLookupTable.h>

typedef struct {
    vtkSmartPointer<vtkLookupTable> lut;
    vtkSmartPointer<vtkColorTransferFunction> ctf;
}ColorMap;

ColorMap loadColorMapFromXML(std::string path);
ColorMap loadIdColorMap(int minId, int maxId);

#endif // COLORMAPHELPER_H
