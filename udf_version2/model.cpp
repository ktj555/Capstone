#include "model.h"
#include "udf.h"

constant figure;

figure.sigma = 0.072;
figure.T_ref = 303.15;
figure.c_sf = 0.006;
figure.D = 0.06;
figure.thickness = 0.008;
figure.init_dt = 0.1;
figure.reservoir_temp = 303.15;
figure.reservoir_enthalpy = 303.15 * 4182;

real Porosity(cell_t c, Thread* t){
    real cen[ND_ND];
    C_CENTROID(cen,c,t);
    return 0.3
}
real Permeability(cell_t c, Thread* t){
    return pow(Particle_Diameter(c,t),2) * pow(Porosity(c,t),3) / (150 * pow(1-Porosity(c,t),2));
}
real Particle_Diameter(cell_t c, Thread* t){
    real cen[ND_ND];
    C_CENTROID(cen,c,t);
    return 100e-6;
}
real Porosity(face_t c, Thread* t){
    real cen[ND_ND];
    C_CENTROID(cen,c,t);
    return 0.3
}
real Particle_Diameter_face(face_t f, Thread* t){
    real cen[ND_ND];
    F_CENTROID(cen,f,t);
    return 100e-6;
}