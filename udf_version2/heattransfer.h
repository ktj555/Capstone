#ifndef HEATTRANSFER_H
#define HEATTRANSFER_H

#include "udf.h"

real alpha_sf(cell_t c, Thread* t);

// Dimensionless Number
real Re_l(cell_t c, Thread* t);
real Re_v(cell_t c, Thread* t);
real Pr_l(cell_t c, Thread* t);
real Pr_v(cell_t c, Thread* t);

real Nu_l(cell_t c, Thread* t);
real Nu_v(cell_t c, Thread* t);
real h_l(cell_t c, Thread* t);
real h_v(cell_t c, Thread* t);

real q_l(cell_t c, Thread* t);
real q_v(cell_t c, Thread* t);
real q_boil(cell_t c, Thread* t);

// derivative of T_f
real dRel_dT(cell_t c,Thread* t);
real dRev_dT(cell_t c,Thread* t);
real dPrl_dT(cell_t c,Thread* t);
real dPrv_dT(cell_t c,Thread* t);
real dNul_dT(cell_t c,Thread* t);
real dNuv_dT(cell_t c,Thread* t);
real dhl_dT(cell_t c,Thread* t);
real dhv_dT(cell_t c,Thread* t);
real dql_dT(cell_t c,Thread* t);
real dqv_dT(cell_t c,Thread* t);

real dqboil_dS(cell_t c,Thread* t);

// derivative of T_s
real dql_dTs(cell_t c, Thread* t);
real dqv_dTs(cell_t c, Thread* t);
real dqboil_dTs(cell_t c, Thread* t);

// BC
real inlet_enthalpy_l(face_t f,Thread* t);
real inlet_enthalpy_v(face_t f,Thread* t);
real inlet_enthalpy_m(face_t f,Thread* t);

// face
real Re_l_face(face_t f, Thread* t);
real Re_v_face(face_t f, Thread* t);
real Pr_l_face(face_t f, Thread* t);
real Pr_v_face(face_t f, Thread* t);
real Nu_l_face(face_t f, Thread* t);
real Nu_v_face(face_t f, Thread* t);
real h_l_face(face_t f, Thread* t);
real h_v_face(face_t f, Thread* t);

#endif