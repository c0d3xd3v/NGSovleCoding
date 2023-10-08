#include "material.h"

Material aluminium = Material(Material::LINEAR, 69.0  * pow(10.0, 6.0), 0.33, 2700.0 * pow(10.0, -9.0));
Material steel     = Material(Material::LINEAR, 220.0 * pow(10.0, 6.0), 0.28, 7700.0 * pow(10.0, -9.0));

Material::Material(MATERIAL_TYPE type, double E, double nu, double rho)
    : E(E)
    , nu(nu)
    , type(type)
    , mu(0)
    , lam(0)
    , rho(0)
{
    if(fabs(E)>10.0e-15  &&
            fabs(nu)>10.0e-15 &&
            fabs(nu+1.0)>10.0e-15 &&
            fabs(nu-1.0/2.0)>10.0e-15 &&
            fabs(rho)>10.0e-15)
    {
        mu        = E / 2.0 / ( 1.0 + nu);
        lam       = E * nu  / ((1.0 + nu) * (1.0 - 2.0 * nu));
        this->rho = rho;
    }
}
