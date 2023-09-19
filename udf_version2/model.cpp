#include "model.h"
#include "udf.h"

extern struct constant models;

real Porosity(cell_t c, Thread* t){
    real cen[ND_ND];
    C_CENTROID(cen,c,t);
    return 0.3;
}
real Permeability(cell_t c, Thread* t){
    return pow(Particle_Diameter(c,t),2) * pow(Porosity(c,t),3) / (150 * pow(1-Porosity(c,t),2));
}
real Particle_Diameter(cell_t c, Thread* t){
    real cen[ND_ND];
    C_CENTROID(cen,c,t);
    return 100e-6;
}
real Porosity_face(face_t c, Thread* t){
    real cen[ND_ND];
    C_CENTROID(cen,c,t);
    return 0.3;
}
real Particle_Diameter_face(face_t f, Thread* t){
    real cen[ND_ND];
    F_CENTROID(cen,f,t);
    return 100e-6;
}