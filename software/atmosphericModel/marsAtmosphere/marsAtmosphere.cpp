#include "math.h"
void simpleMarsAtmosphere(const double height,const double speed,double *T_atm, double *p_atm,double *rho, double *M, double *g){
    double a;
  // ==================%
  //  Atmosphere model %
  // ==================%

  if (height <= 7000) {
    *T_atm = -23.4 - 0.00222 * height;
    //  Temperature
    
    *p_atm = 699.0 * std::exp(-0.00009*height);
    //  Pressure
  } else {
    *T_atm = -31.0 - 0.000998 * height;
    //  Temperature
    *p_atm = 699.0 * std::exp(-0.00009*height);
    //  Pressure
  }
  
  *rho = *p_atm / (.1921 * (273.1 + *T_atm));
  //  Density
  
  // Gamma = 1.2, R = 188
  *M = speed / std::sqrt(225.6 * *T_atm);
  //  Mach number

  a = 3,389.5 / (height / 1000.0 + 3,389.5);
  *g = 3.721 * (a * a);
  //  Accleration due to gravity at height h
}