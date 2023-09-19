#ifndef INFO_H
#define INFO_H

#include "udf.h"

int state(cell_t c, Thread* t);


real T_f(cell_t c, Thread* t);
real S_(cell_t c, Thread* t);
real dT_dH(cell_t c,Thread* t);
real dS_dH(cell_t c,Thread* t);

real T_sat(cell_t c, Thread* t);
real H_fg(cell_t c, Thread* t);


real H_sat_l(cell_t c, Thread* t);
real H_sat_v(cell_t c, Thread* t);


real J_function(cell_t c,Thread* t);
real J_function_derivative(cell_t c,Thread* t);
real beta_current(cell_t c,Thread* t);


// past
int state_past(cell_t c,Thread* t);
real T_f_past(cell_t c,Thread* t);
real S_past(cell_t c,Thread* t);
real T_sat_past(cell_t c, Thread* t);
real H_fg_past(cell_t c,Thread* t);
real H_sat_l_past(cell_t c, Thread* t);
real H_sat_v_past(cell_t c,Thread* t);
real beta_past(cell_t c,Thread* t);

// face
int state(face_t f, Thread* t);
real T_f(face_t f,Thread* t);
real S_(face_t f,Thread* t);
real T_sat(face_t f,Thread* t);
real H_fg(face_t f,Thread* t);
real H_sat_l(face_t f,Thread* t);
real H_sat_v(face_t f,Thread* t);
real beta_current(face_t f,Thread* t);

#endif