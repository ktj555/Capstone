#include "property.h"
#include "info.h"
#include "udf.h"

// heat coefficient
real Conductivity_l(cell_t c, Thread* t){
    return 0.68;
}
real Conductivity_v(cell_t c, Thread* t){
    real a0, a1, T;
    a0 = -21.99433e-3;
    a1 = 118.42e-6;
    T = T_f(c,t);
    return a0 + a1 * T;
}
real Conductivity_s(cell_t c, Thread* t){
    return 21.7;
}
real Specific_Heat_l(cell_t c, Thread* t){
    return 4182;
}
real Specific_Heat_v(cell_t c, Thread* t){
    return 2030;
}
real Specific_Heat_s(cell_t c,Thread* t){
    return 625;
}



real Rho_s(cell_t c, Thread* t){
    return 8.4e3;
}

// flow coefficient
real Rho_l(cell_t c, Thread* t){
    return 998;
}
real Rho_v(cell_t c,Thread* t){
    real R = 8.314, M_V = 0.018;
    real p;
    if(NNULLP(THREAD_STORAGE(t,SV_P))){
        p = (C_P(c,t) + RP_Get_Real("operating-pressure")) / 1e6;
    }
    else{
        p = RP_Get_Real("operating-pressure") / 1e6;
    }
    real T = T_f(c,t);
    return p / (T * R / M_V);
}
real Viscosity_l(cell_t c, Thread* t){
    real T;
    T = T_f(c,t);
    return 24.141e-3 * pow(10,247.8 / (T - 140));
}
real Viscosity_v(cell_t c, Thread* t){
    real T;
    T = T_f(c,t);
    return -2.77567e-3 + 40.35e-6 * T;
}
real Kinematic_Viscosity_l(cell_t c, Thread* t){
    return Viscosity_l(c,t) / Rho_l(c,t);
}
real Kinematic_Viscosity_v(cell_t c, Thread* t){
    return Viscosity_v(c,t) / Rho_v(c,t);
}

// mixture
real K_rl(cell_t c, Thread* t){
    return pow(S_(c,t),3);
}
real K_rv(cell_t c, Thread* t){
    return pow(1-S_(c,t),3);
}
real L(cell_t c, Thread* t){
    return K_rl(c,t) * Kinematic_Viscosity_m(c,t) / Kinematic_Viscosity_l(c,t);
}
real Rho_m(cell_t c, Thread* t){
    return S_(c,t) * Rho_l(c,t) + (1 - S_(c,t)) * Rho_v(c,t);
}
real Viscosity_m(cell_t c, Thread* t){
    return Rho_m(c,t) * Kinematic_Viscosity_m(c,t);
}
real Kinematic_Viscosity_m(cell_t c, Thread* t){
    return 1 / (K_rl(c,t) / Kinematic_Viscosity_l(c,t) + K_rv(c,t) / Kinematic_Viscosity_v(c,t));
}
real Conductivity_m(cell_t c,Thread* t){
    return S_(c,t) * Conductivity_l(c,t) + (1-S_(c,t)) * Conductivity_v(c,t);
}

// derivative of T_f
real dkl_dT(cell_t c,Thread* t){
    return 0;
}
real dkv_dT(cell_t c,Thread* t){
    real a1 = 118.42e-6;
    return a1;
}
real dcpl_dT(cell_t c,Thread* t){
    return 0;
}
real dcpv_dT(cell_t c,Thread* t){
    return 0;
}
real drhol_dT(cell_t c,Thread* t){
    return 0;
}
real drhov_dT(cell_t c,Thread* t){
    return - Rho_v(c,t) / T_f(c,t);
}
real dmul_dT(cell_t c,Thread* t){
    real T;
    T = T_f(c,t);
    return -Viscosity_l(c,t) * 247.8 / pow(T - 140,2) * log(10);
}
real dmuv_dT(cell_t c,Thread* t){
    real a1 = 40.35e-6;
    return a1;
}
real dnul_dT(cell_t c,Thread* t){
    return (dmul_dT(c,t) * Rho_l(c,t) - Viscosity_l(c,t) * drhol_dT(c,t)) / pow(Rho_l(c,t),2);
}
real dnuv_dT(cell_t c,Thread* t){
    return (dmuv_dT(c,t) * Rho_v(c,t) - Viscosity_v(c,t) * drhov_dT(c,t)) / pow(Rho_v(c,t),2);
}

// derivative of S
real dkrl_dS(cell_t c,Thread* t){
    return 3 * pow(S_(c,t),2);
}
real dkrv_dS(cell_t c,Thread* t){
    return -3 * pow(1-S_(c,t),2);
}
real dl_dS(cell_t c,Thread* t){
    return (dkrl_dS(c,t) * Kinematic_Viscosity_m(c,t) + K_rl(c,t) * dnu_dS(c,t)) / Kinematic_Viscosity_l(c,t);
}
real drho_dS(cell_t c,Thread* t){
    return Rho_l(c,t) - Rho_v(c,t);
}
real dmu_dS(cell_t c,Thread* t){
    return drho_dS(c,t) * Kinematic_Viscosity_m(c,t) + Rho_m(c,t) * dnu_dS(c,t);
}
real dnu_dS(cell_t c,Thread* t){
    return -pow(Kinematic_Viscosity_m(c,t),2) * (dkrl_dS(c,t) / Kinematic_Viscosity_l(c,t) + dkrv_dS(c,t) / Kinematic_Viscosity_v(c,t));
}

// past
real Specific_Heat_l_past(cell_t c,Thread* t){
    return 4182;
}
real Specific_Heat_v_past(cell_t c,Thread* t){
    return 2030;
}
real Specific_Heat_s_past(cell_t c,Thread* t){
    return 625;
}
real Rho_l_past(cell_t c,Thread* t){
    return 998;
}
real Rho_s_past(cell_t c,Thread* t){
    return 8.4e3;
}
real Rho_v_past(cell_t c,Thread* t){
    real R = 8.314, M_V = 0.018;
    real p;
    if(NNULLP(THREAD_STORAGE(t,SV_P))){
        p = (C_P_M1(c,t) + RP_Get_Real("operating-pressure")) / 1e6;
    }
    else{
        p = RP_Get_Real("operating-pressure") / 1e6;
    }
    real T = T_f_past(c,t);
    return p / (T * R / M_V);
}
real Viscosity_l_past(cell_t c, Thread* t){
    real T;
    T = T_f_past(c,t);
    return 24.141e-3 * pow(10,247.8 / (T - 140));
}
real Viscosity_v_past(cell_t c, Thread* t){
    real T;
    T = T_f_past(c,t);
    return -2.77567e-3 + 40.35e-6 * T;
}
real Kinematic_Viscosity_l_past(cell_t c, Thread* t){
    return Viscosity_l_past(c,t) / Rho_l_past(c,t);
}
real Kinematic_Viscosity_v_past(cell_t c, Thread* t){
    return Viscosity_v_past(c,t) / Rho_v_past(c,t);
}
real K_rl_past(cell_t c, Thread* t){
    return pow(S_past(c,t),3);
}
real K_rv_past(cell_t c, Thread* t){
    return pow(1-S_past(c,t),3);
}
real L_past(cell_t c, Thread* t){
    return K_rl_past(c,t) * Kinematic_Viscosity_m_past(c,t) / Kinematic_Viscosity_l_past(c,t);
}
real Rho_m_past(cell_t c,Thread* t){
    return S_past(c,t) * Rho_l_past(c,t) + (1 - S_past(c,t)) * Rho_v_past(c,t);
}
real Viscosity_m_past(cell_t c, Thread* t){
    return Rho_m_past(c,t) * Kinematic_Viscosity_m_past(c,t);
}
real Kinematic_Viscosity_m_past(cell_t c, Thread* t){
    return 1 / (K_rl_past(c,t) / Kinematic_Viscosity_l_past(c,t) + K_rv_past(c,t) / Kinematic_Viscosity_v_past(c,t));
}


// face
real Conductivity_l_face(face_t f, Thread* t){
    return 0.68;
}
real Conductivity_v_face(face_t f, Thread* t){
    real a0, a1, T;
    a0 = -21.99433e-3;
    a1 = 118.42e-6;
    T = T_f_face(f,t);
    return a0 + a1 * T;
}
real Conductivity_s_face(face_t f, Thread* t){
    return 21.7;
}
real Conductivity_m_face(face_t f,Thread* t){
    return S_face(f,t) * Conductivity_l_face(f,t) + (1-S_face(f,t)) * Conductivity_v_face(f,t);
}
real Specific_Heat_l_face(face_t f, Thread* t){
    return 4182;
}
real Specific_Heat_v_face(face_t f, Thread* t){
    return 2030;
}
real Specific_Heat_s_face(face_t f,Thread* t){
    return 625;
}
real Specific_Heat_m_face(face_t f,Thread* t){
    return (S_face(f,t) * Rho_l_face(f,t) * Specific_Heat_l_face(f,t) + (1-S_face(f,t)) * Rho_v_face(f,t) * Specific_Heat_v_face(f,t)) / Rho_m_face(f,t);
}
real Rho_l_face(face_t f, Thread* t){
    return 998;
}
real Rho_v_face(face_t f,Thread* t){
    real R = 8.314, M_V = 0.018;
    real p;
    if(NNULLP(THREAD_STORAGE(t,SV_P))){
        p = (F_P(f,t) + RP_Get_Real("operating-pressure")) / 1e6;
    }
    else{
        p = RP_Get_Real("operating-pressure") / 1e6;
    }
    real T = T_f_face(f,t);
    return p / (T * R / M_V);
}
real Rho_m_face(face_t f, Thread* t){
    return S_face(f,t) * Rho_l_face(f,t) + (1 - S_face(f,t)) * Rho_v_face(f,t);
}
real Viscosity_l_face(face_t f, Thread* t){
    real T;
    T = T_f_face(f,t);
    return 24.141e-3 * pow(10,247.8 / (T - 140));
}
real Viscosity_v_face(face_t f, Thread* t){
    real T;
    T = T_f_face(f,t);
    return -2.77567e-3 + 40.35e-6 * T;
}
real Kinematic_Viscosity_l_face(face_t f, Thread* t){
    return Viscosity_l_face(f,t) / Rho_l_face(f,t);
}
real Kinematic_Viscosity_v_face(face_t f, Thread* t){
    return Viscosity_v_face(f,t) / Rho_v_face(f,t);
}
real Kinematic_Viscosity_m_face(face_t f, Thread* t){
    return 1 / (K_rl_face(f,t) / Kinematic_Viscosity_l_face(f,t) + K_rv_face(f,t) / Kinematic_Viscosity_v_face(f,t));
}
real K_rl_face(face_t f, Thread* t){
    return pow(S_face(f,t),3);
}
real K_rv_face(face_t f, Thread* t){
    return pow(1-S_face(f,t),3);
}
real L_face(face_t f, Thread* t){
    return K_rl_face(f,t) * Kinematic_Viscosity_m_face(f,t) / Kinematic_Viscosity_l_face(f,t);
}