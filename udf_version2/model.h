#ifndef MODEL_H
#define MODEL_H

#include "udf.h"

typedef struct constant{
    real sigma;
    real T_ref;
    real c_sf;
    real D;
    real thickness;
    real init_dt;
    real reservoir_temp;
    real reservoir_enthalpy;
}constant;

enum uds{
    enthalpy,
    temp_s
};

enum udm{
    noting_else
};

enum state_of_cell{
    liquid,
    vapor,
    mixture
};

real Porosity(cell_t c, Thread* t);
real Permeability(cell_t c, Thread* t);
real Particle_Diameter(cell_t c, Thread* t);
real Porosity(face_t f,Thread* t);
real Particle_Diameter_face(face_t f,Thread* t);

extern constant figure;

#endif