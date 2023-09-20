#ifndef MODEL_H
#define MODEL_H

#include "udf.h"

struct constant{
    real sigma = 0.072;
    real T_ref = 303.15;
    real c_sf = 0.006;
    real D = 0.06;
    real thickness = 0.008;
    real init_dt = 10;
    real reservoir_temp = 303.15;
    real reservoir_enthalpy = 0;
};

enum uds{
    enthalpy,
    temp_s
};

enum udm{
    data
};

enum state_of_cell{
    liquid,
    vapor,
    mixture
};

real Porosity(cell_t c, Thread* t);
real Permeability(cell_t c, Thread* t);
real Particle_Diameter(cell_t c, Thread* t);
real Porosity_face(face_t f,Thread* t);
real Particle_Diameter_face(face_t f,Thread* t);


#endif