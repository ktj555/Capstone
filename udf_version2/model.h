#ifndef MODEL_H
#define MODEL_H

#include "udf.h"

struct constant{
    real sigma = 0.072;
    real T_ref = 303.15;
    real c_sf = 0.006;
    real D = 0.06;
    real thickness = 0.008;
    real init_dt = 0.1;
    real reservoir_temp = 303.15;
    real reservoir_enthalpy = 303.15 * 4182;
};

enum uds{
    enthalpy,
    temp_s
};

enum udm{

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
real Particle_Diameter(face_t f,Thread* t);

constant model;

#endif