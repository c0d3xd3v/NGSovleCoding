#include <iostream>
#include <random>
#include "tinyxml2/tinyxml2.h"

#include "colormaphelper.h"

ColorMap loadColorMapFromXML(std::string path)
{
    tinyxml2::XMLDocument doc;
    if (doc.LoadFile(path.c_str()) != tinyxml2::XML_SUCCESS) {
        std::cerr << "Fehler beim Laden der XML-Datei." << std::endl;
    }

    tinyxml2::XMLElement* colorMapElement = doc.FirstChildElement("ColorMap");
    if (!colorMapElement) {
        std::cerr << "Fehler: ColorMap-Element nicht gefunden." << std::endl;
    }

    const char* name = colorMapElement->Attribute("name");
    const char* space = colorMapElement->Attribute("space");
    std::cout << "Name: " << name << std::endl;
    std::cout << "Space: " << space << std::endl;

    tinyxml2::XMLElement* pointElement = colorMapElement->FirstChildElement("Point");
    vtkNew<vtkColorTransferFunction> colorTransferFunction;

    // we need this fix to convert strings to double from the xml file.
    // https://stackoverflow.com/a/54145018
    // std::stod is defined in terms of std::strtod, which is inherited from
    // the C standard library. The C function strtod works in terms of the
    // C locale, accessible via the setlocale function from the <locale.h>
    // header.                                                                                                                                                In C++, the C locale is still accessible via std::setlocale function in the <clocale> header, and it does influence both std::strtod and std::stod.
    // Qt's QApplication uses std::setlocale to set the user-chosen locale.
    // Thus whenever you use a C-locale-dependent function in a
    // GUI Qt application, you'll have locale-dependent radix point.
    // To force a particular locale for numbers, you can use std::setlocale.
    const std::string oldLocale = std::setlocale(LC_NUMERIC, nullptr);
    std::setlocale(LC_NUMERIC, "C");

    while (pointElement) {
        const char* xStr = pointElement->Attribute("x");
        //const char* oStr = pointElement->Attribute("o");
        const char* rStr = pointElement->Attribute("r");
        const char* gStr = pointElement->Attribute("g");
        const char* bStr = pointElement->Attribute("b");

        double x = std::atof(xStr);
        //double o = std::atof(oStr);
        double r = std::atof(rStr);
        double g = std::atof(gStr);
        double b = std::atof(bStr);

        /*
        std::cout << "load colormap point - x: " << x
                  << ", o: " << o
                  << ", r: " << r
                  << ", g: " << g
                  << ", b: " << b
                  << std::endl;
        */

        colorTransferFunction->AddRGBPoint(x, r, g, b);

        pointElement = pointElement->NextSiblingElement("Point");
    }

    std::setlocale(LC_NUMERIC, oldLocale.c_str());

    int color_steps = 8;
    vtkNew<vtkLookupTable> lut;
    lut->SetNumberOfTableValues(color_steps);
    lut->Build();

    for(int i = 0; i < color_steps; i++)
    {
        float x = float(i) / float(color_steps - 1.0);
        double rgb[3];
        colorTransferFunction->GetColor(x, rgb);
        lut->SetTableValue(i, rgb[0], rgb[1], rgb[2]);
    }

    colorTransferFunction->SetVectorModeToMagnitude();
    lut->SetVectorModeToMagnitude();

    return ColorMap{lut, colorTransferFunction};
}

ColorMap loadIdColorMap(int minId, int maxId)
{
    vtkNew<vtkColorTransferFunction> colorTransferFunction;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int i = minId; i < maxId; i++) {
        double x = float(i);
        double o = 1.0;
        double r = dis(gen);//std::atof(rStr);
        double g = dis(gen);//std::atof(gStr);
        double b = dis(gen);//std::atof(bStr);

        /*
        std::cout << "load colormap point - x: " << x
                  << ", o: " << o
                  << ", r: " << r
                  << ", g: " << g
                  << ", b: " << b
                  << std::endl;
        */

        colorTransferFunction->AddRGBPoint(x, r, g, b);

    }

    int color_steps = maxId;
    vtkNew<vtkLookupTable> lut;
    //lut->SetNumberOfTableValues(color_steps);
    lut->Build();
    /*
    for(int i = 0; i < color_steps; i++)
    {
        float x = float(i) / float(color_steps - 1.0);
        double rgb[3];
        colorTransferFunction->GetColor(x, rgb);
        lut->SetTableValue(i, rgb[0], rgb[1], rgb[2]);
    }
    */
    colorTransferFunction->SetVectorModeToMagnitude();
    lut->SetVectorModeToMagnitude();

    return ColorMap{lut, colorTransferFunction};
}
