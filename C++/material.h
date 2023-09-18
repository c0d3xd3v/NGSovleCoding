#ifndef MATERIAL_H
#define MATERIAL_H

#include <cmath>

class Material
{
    public:
        typedef enum {
            LINEAR,
            END_LIST
        } MATERIAL_TYPE;
        double E;
        double nu;
        double mu;
        double lam;
        double rho;
        MATERIAL_TYPE type;

        Material(MATERIAL_TYPE type, double E=0, double nu=0, double rho=0);
};

extern Material aluminium;
extern Material steel;

#endif // MATERIAL_H
