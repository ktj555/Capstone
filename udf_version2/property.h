#ifndef PROPERTY_H
#define PROPERTY_H

#include "udf.h"

// heat coefficient
real Conductivity_l(cell_t c, Thread* t);
real Conductivity_v(cell_t c, Thread* t);
real Conductivity_s(cell_t c, Thread* t);
real Specific_Heat_l(cell_t c, Thread* t);
real Specific_Heat_v(cell_t c, Thread* t);
real Specific_Heat_s(cell_t c, Thread* t);

real Rho_s(cell_t c, Thread* t);

// flow coefficient
real Rho_l(cell_t c, Thread* t);
real Rho_v(cell_t c,Thread* t);
real Viscosity_l(cell_t c, Thread* t);
real Viscosity_v(cell_t c, Thread* t);
real Kinematic_Viscosity_l(cell_t c, Thread* t);
real Kinematic_Viscosity_v(cell_t c, Thread* t);



// mixture
real K_rl(cell_t c, Thread* t);
real K_rv(cell_t c, Thread* t);
real L(cell_t c, Thread* t);
real Rho_m(cell_t c, Thread* t);
real Viscosity_m(cell_t c, Thread* t);
real Kinematic_Viscosity_m(cell_t c,Thread* t);
real Conductivity_m(cell_t c,Thread* t);



// derivative of T_f
real dkl_dT(cell_t c,Thread* t);
real dkv_dT(cell_t c,Thread* t);
real dcpl_dT(cell_t c,Thread* t);
real dcpv_dT(cell_t c,Thread* t);
real drhol_dT(cell_t c,Thread* t);
real drhov_dT(cell_t c,Thread* t);
real dmul_dT(cell_t c,Thread* t);
real dmuv_dT(cell_t c,Thread* t);
real dnul_dT(cell_t c,Thread* t);
real dnuv_dT(cell_t c,Thread* t);

// derivative of S
real dkrl_dS(cell_t c,Thread* t);
real dkrv_dS(cell_t c,Thread* t);
real dl_dS(cell_t c,Thread* t);
real drho_dS(cell_t c,Thread* t);
real dmu_dS(cell_t c,Thread* t);
real dnu_dS(cell_t c,Thread* t);

// past
real Specific_Heat_l_past(cell_t c,Thread* t);
real Specific_Heat_v_past(cell_t c,Thread* t);
real Specific_Heat_s_past(cell_t c,Thread* t);
real Rho_l_past(cell_t c, Thread* t);
real Rho_v_past(cell_t c, Thread* t);
real Rho_s_past(cell_t c,Thread* t);
real Viscosity_l_past(cell_t c,Thread* t);
real Viscosity_v_past(cell_t c,Thread* t);
real Kinematic_Viscosity_l_past(cell_t c,Thread* t);
real Kinematic_Viscosity_v_past(cell_t c,Thread* t);
real K_rl_past(cell_t c, Thread* t);
real K_rv_past(cell_t c, Thread* t);
real L_past(cell_t c, Thread* t);
real Rho_m_past(cell_t c,Thread* t);
real Viscosity_m_past(cell_t c, Thread* t);
real Kinematic_Viscosity_m_past(cell_t c,Thread* t);

// face
real Specific_Heat_l_face(face_t f,Thread* t);
real Specific_Heat_v_face(face_t f,Thread* t);
real Specific_Heat_s_face(face_t f,Thread* t);
real Specific_Heat_m_face(face_t f,Thread* t);
real Conductivity_l_face(face_t f, Thread* t);
real Conductivity_v_face(face_t f, Thread* t);
real Conductivity_s_face(face_t f, Thread* t);
real Conductivity_m_face(face_t f, Thread*t);
real Rho_l_face(face_t f, Thread* t);
real Rho_v_face(face_t f, Thread* t);
real Rho_m_face(face_t f, Thread* t);
real Viscosity_l_face(face_t f,Thread* t);
real Viscosity_v_face(face_t f,Thread* t);
real Kinematic_Viscosity_l_face(face_t f,Thread* t);
real Kinematic_Viscosity_v_face(face_t f,Thread* t);
real K_rl_face(face_t f,Thread* t);
real L_face(face_t f,Thread* t);

#endif