#include "info.h"
#include "property.h"
#include "model.h"
#include "udf.h"

int state(cell_t c, Thread* t){
    real h_l_sat = H_sat_l(c,t);
    real h_v_sat = H_sat_v(c,t);
    real h;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(uds::enthalpy)))){
        h = C_UDSI(c,t,uds::enthalpy);
    }
    else{
        h = 0;
    }
    
    if(h<=h_l_sat) return state_of_cell::liquid;
    else if(h>=h_v_sat) return state_of_cell::vapor;
    else return state_of_cell::mixture;
}



real T_f(cell_t c, Thread* t){
    real h;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(uds::enthalpy)))){
        h = C_UDSI(c,t,uds::enthalpy);
    }
    else{
        h = 0;
    }
    switch(state(c,t)){
    case state_of_cell::liquid:
        return model.T_ref + h / Specific_Heat_l(c,t);
    case state_of_cell::vapor:
        return T_sat(c,t) + (h - H_v_sat(c,t)) / Specific_Heat_v(c,t);
    case state_of_cell::mixture:
        return T_sat(c,t);
    }
}
real S_(cell_t c,Thread* t){
    real h;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(uds::enthalpy)))){
        h = C_UDSI(c,t,uds::uds::enthalpy);
    }
    else{
        h = 0;
    }
    switch(state(c,t)){
    case state_of_cell::liquid:
        return 1;
    case state_of_cell::vapor:
        return 0;
    case state_of_cell::mixture:
        real nu_r, h_r;
        nu_r = Kinematic_Viscosity_v(c,t) / Kinematic_Viscosity_l(c,t);
        h_r = H_fg(c,t) / (H_sat_v(c,t) - h);
        return 1 / (pow(nu_r * (h_r - 1), 0.33) + 1);
    }
}
real dT_dH(cell_t c, Thread* t){
    switch(state(c,t)){
    case state_of_cell::liquid:
        return 1/(Specific_Heat_l(c,t) + dcpl_dT(c,t) * T_f(c,t));
    case state_of_cell::vapor:
        return 1/(Specific_Heat_v(c,t) + dcpv_dT(c,t) * T_f(c,t));
    case state_of_cell::mixture:
        return dS_dH(c,t);
    }
}
real dS_dH(cell_t c,Thread* t){
    return -1/(H_fg(c,t) * dl_dS(c,t))
}

real T_sat(cell_t c, Thread* t){
    real P;
    real a0,a1,a2,a3;
    if(NNULLP(THREAD_STORAGE(t,SV_P))){
        P = (C_P(c,t); + RP_Get_Real("operating-pressure")) / 1e6;
    }
    else{
        P = RP_Get_Real("operating-pressure") / 1e6;
    }
    a0 = 429.69474687605754;
	a1 = 11.983108526790891;
	a2 = -0.2940193901122584;
	a3 = 25.197563642164386;
    return a0 + a1 * p + a2 * pow(p, 2) + a3 * log(p);
}
real H_fg(cell_t c, Thread* t){
    real T;
    T = T_sat(c,t);
    if (T > 650) { return 0; }
	real a0 = -2511.1589762597114;
	real a1 = 1.790037469059945;
	real a2 = 0.0028969464729942194;
	real a3 = 222.03437607975306;
	return (a0 + a1 * T + a2 * pow(T, 2) + a3 * sqrt(650 - T))*1000;
}


real H_sat_l(cell_t c, Thread* t){
    return Specific_Heat_l(c,t) * (T_sat(c,t) - model.T_ref);
}
real H_sat_v(cell_t c, Thread* t){
    return H_sat_l(c,t) + H_fg(c,t);
}


real J_function(cell_t c, Thread* t){
    return 1.417*(1-S_(c,t)) - 2.120*pow(1-S_(c,t),2) + 1.263*pow(1-S_(c,t),3);
}
real J_function_derivative(cell_t c,Thread* t){
    return - (1.417 - 2*2.120 * (1-S_(c,t)) + 3 * 1.263 * pow(1-S_(c,t),2));
}
real beta_current(cell_t c,Thread* t){
    switch(state(c,t)){
    case state_of_cell::liquid:
        return 1;
    case state_of_cell::vapor:
        return 1;
    case state_of_cell::mixture:
        return (Rho_v(c,t)/ Rho_m(c,t) * (1-S_(c,t)) * H_fg(c,t) + Specific_Heat_l(c,t)*T_sat(c,t)) / ((1-L(c,t))*H_fg(c,t)+Specific_Heat_l(c,t)*T_sat(c,t))
    }
}


// past
int state_past(cell_t c,Thread* t){
    real h_l_sat_past = H_sat_l_past(c,t);
    real h_v_sat_past = H_sat_v_past(c,t);
    real h;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(uds::enthalpy)))){
        h = C_UDSI_M1(c,t,uds::enthalpy);
    }
    else{
        h = 0;
    }

    if(h<=h_l_sat) return state_of_cell::liquid;
    else if(h>=h_v_sat) return state_of_cell::vapor;
    else return state_of_cell::mixture;
}
real T_f_past(cell_t c, Thread* t){
    real h;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(uds::enthalpy)))){
        h = C_UDSI_M1(c,t,uds::enthalpy);
    }
    else{
        h = 0;
    }
    switch(state_past(c,t)){
    case state_of_cell::liquid:
        return model.T_ref + h / Specific_Heat_l_past(c,t);
    case state_of_cell::vapor:
        return T_sat_past(c,t) + (h - H_v_sat_past(c,t)) / Specific_Heat_v_past(c,t);
    case state_of_cell::mixture:
        return T_sat(c,t);
    }
}
real S_past(cell_t c, Thread* t){
    real h;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(uds::enthalpy)))){
        h = C_UDSI_M1(c,t,uds::enthalpy);
    }
    else{
        h = 0;
    }
    switch(state_past(c,t)){
    case state_of_cell::liquid:
        return 1;
    case state_of_cell::vapor:
        return 0;
    case state_of_cell::mixture:
        real nu_r, h_r;
        nu_r = Kinematic_Viscosity_v_past(c,t) / Kinematic_Viscosity_l_past(c,t);
        h_r = H_fg_past(c,t) / (H_sat_v_past(c,t) - h);
        return 1 / (pow(nu_r * (h_r - 1), 0.33) + 1);
    }
}
real T_sat_past(cell_t c, Thread* t){
    real P;
    real a0,a1,a2,a3;
    if(NNULLP(THREAD_STORAGE(t,SV_P))){
        P = (C_P_M1(c,t); + RP_Get_Real("operating-pressure")) / 1e6;
    }
    else{
        P = RP_Get_Real("operating-pressure") / 1e6;
    }
    a0 = 429.69474687605754;
	a1 = 11.983108526790891;
	a2 = -0.2940193901122584;
	a3 = 25.197563642164386;
    return a0 + a1 * p + a2 * pow(p, 2) + a3 * log(p);
}
real H_fg_past(cell_t c, Thread* t){
    real T;
    T = T_sat_past(c,t);
    if (T > 650) { return 0; }
	real a0 = -2511.1589762597114;
	real a1 = 1.790037469059945;
	real a2 = 0.0028969464729942194;
	real a3 = 222.03437607975306;
	return (a0 + a1 * T + a2 * pow(T, 2) + a3 * sqrt(650 - T))*1000;
}
real H_sat_l_past(cell_t c, Thread* t){
    return Specific_Heat_l_past(c,t) * T_sat_past(c,t);
}
real H_sat_v_past(cell_t c, Thread* t){
    return H_sat_l_past(c,t) + H_fg_past(c,t);
}
real beta_past(cell_t c, Thread* t){
    switch(state_past(c,t)){
    case state_of_cell::liquid:
        return 1;
    case state_of_cell::vapor:
        return 1;
    case state_of_cell::mixture:
        return (Rho_v_past(c,t) / Rho_m_past(c,t) * (1-S_past(c,t)) * H_fg_past(c,t) + Specific_Heat_l_past(c,t) * T_sat_past(c,t)) / ((1-L_past(c,t))*H_fg_past(c,t) + Specific_Heat_l_past(c,t) * T_sat_past(c,t))
    }
}

// face
int state(face_t f, Thread* t){
    real h_l_sat = H_sat_l(f,t);
    real h_v_sat = H_sat_v(f,t);
    real h;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(uds::enthalpy)))){
        h = F_UDSI(f,t,uds::uds::enthalpy);
    }
    else{
        h = 0;
    }
    
    if(h<=h_l_sat) return state_of_cell::liquid;
    else if(h>=h_v_sat) return state_of_cell::vapor;
    else return state_of_cell::mixture;
}
real T_f(face_t f, Thread* t){
    real h;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(uds::enthalpy)))){
        h = F_UDSI(f,t,uds::enthalpy);
    }
    else{
        h = 0;
    }
    switch(state(f,t)){
    case state_of_cell::liquid:
        return model.T_ref + h / Specific_Heat_l(f,t);
    case state_of_cell::vapor:
        return T_sat(f,t) + (h - H_v_sat(f,t)) / Specific_Heat_v(f,t);
    case state_of_cell::mixture:
        return T_sat(f,t);
    }
}
real S_(face_t f,Thread* t){
    real h;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(uds::enthalpy)))){
        h = F_UDSI(f,t,uds::enthalpy);
    }
    else{
        h = 0;
    }
    switch(state(f,t)){
    case state_of_cell::liquid:
        return 1;
    case state_of_cell::vapor:
        return 0;
    case state_of_cell::mixture:
        real nu_r, h_r;
        nu_r = Kinematic_Viscosity_v(f,t) / Kinematic_Viscosity_l(f,t);
        h_r = H_fg(f,t) / (H_sat_v(f,t) - h);
        return 1 / (pow(nu_r * (h_r - 1), 0.33) + 1);
    }
}
real T_sat(face_t f, Thread* t){
    real P;
    real a0,a1,a2,a3;
    if(NNULLP(THREAD_STORAGE(t,SV_P))){
        P = (F_P(f,t); + RP_Get_Real("operating-pressure")) / 1e6;
    }
    else{
        P = RP_Get_Real("operating-pressure") / 1e6;
    }
    a0 = 429.69474687605754;
	a1 = 11.983108526790891;
	a2 = -0.2940193901122584;
	a3 = 25.197563642164386;
    return a0 + a1 * p + a2 * pow(p, 2) + a3 * log(p);
}
real H_fg(face_t f, Thread* t){
    real T;
    T = T_sat(f,t);
    if (T > 650) { return 0; }
	real a0 = -2511.1589762597114;
	real a1 = 1.790037469059945;
	real a2 = 0.0028969464729942194;
	real a3 = 222.03437607975306;
	return (a0 + a1 * T + a2 * pow(T, 2) + a3 * sqrt(650 - T))*1000;
}
real H_sat_l(face_t f, Thread* t){
    return Specific_Heat_l(f,t) * T_sat(f,t);
}
real H_sat_v(cell_t c, Thread* t){
    return H_sat_l(f,t) + H_fg(f,t);
}
real beta_current(face_t f,Thread* t){
    switch(state(f,t)){
    case state_of_cell::liquid:
        return 1;
    case state_of_cell::vapor:
        return 1;
    case state_of_cell::mixture:
        return (Rho_v(f,t)/ Rho_m(f,t) * (1-S_(f,t)) * H_fg(f,t) + Specific_Heat_l(f,t)*T_sat(f,t)) / ((1-L(f,t))*H_fg(f,t)+Specific_Heat_l(f,t)*T_sat(f,t))
    }
}