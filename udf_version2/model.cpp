#include "model.h"
#include "udf.h"

constant model;

model.sigma = 0.072;
model.T_ref = 303.15;
model.c_sf = 0.006;
model.D = 0.06;
model.thickness = 0.008;
model.init_dt = 0.1;
model.reservoir_temp = 303.15;
model.reservoir_enthalpy = 303.15 * 4182;

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
real Particle_Diameter(face_t f, Thread* t){
    real cen[ND_ND];
    F_CENTROID(cen,f,t);
    return 100e-6;
}