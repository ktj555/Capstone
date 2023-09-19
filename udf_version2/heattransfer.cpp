#include "heattransfer.h"
#include "property.h"
#include "model.h"
#include "info.h"

real alpha_sf(cell_t c, Thread* t){
    return 6 * (1-Porosity(c,t)) / Particle_Diameter(c,t);
}

real Re_l(cell_t c,Thread* t){
    real NV_VEC(psi);
    real v;
    if(NNULLP(THREAD_STORAGE(t,SV_U))){
        NV_D(psi, =, C_U(c,t), C_V(c,t), C_W(c,t));
        v = NV_MAG(psi);
    }
    else{
        v = 0;
    }
    return v * Particle_Diameter(c,t) / Kinematic_Viscosity_l(c,t);
}
real Re_v(cell_t c, Thread* t){
    real NV_VEC(psi);
    real v;
    if(NNULLP(THREAD_STORAGE(t,SV_U))){
        NV_D(psi, =, C_U(c,t), C_V(c,t), C_W(c,t));
        v = NV_MAG(psi);
    }
    else{
        v = 0;
    }
    return v * Particle_Diameter(c,t) / Kinematic_Viscosity_v(c,t);
}
real Pr_l(cell_t c, Thread* t){
    return Viscosity_l(c,t) * Specific_Heat_l(c,t) / Conductivity_l(c,t);
}
real Pr_v(cell_t c, Thread* t){
    return Viscosity_v(c,t) * Specific_heat_v(c,t) / Conductivity_v(c,t);
}

real Nu_l(cell_t c, Thread* t){
    return 2 + 1.1 * pow(Pr_l(c,t), 0.33) * pow(Re_l(c,t), 0.6);
}
real Nu_v(cell_t c, Thread* t){
    return 2 + 1.1 * pow(Pr_v(c,t),0.33) * pow(Re_v(c,t), 0.6);
}
real h_l(cell_t c, Thread* t){
    return Nu_l(c,t) * Conductivity_l(c,t) / Particle_Diameter(c,t);
}
real h_v(cell_t c, Thread* t){
    return Nu_v(c,t) * Conductivity_v(c,t) / Particle_Diameter(c,t);
}

real q_l(cell_t c,Thread* t){
    real T_s;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(temp_s)))){
        T_s = C_UDSI(c,t,temp_s);
    }
    else{
        T_s = T_f(c,t);
    }
    return h_l(c,t) * alpha_sf(c,t) * (T_s - T_f(c,t));
}
real q_v(cell_t c,Thread* t){
    real T_s;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(temp_s)))){
        T_s = C_UDSI(c,t,temp_s);
    }
    else{
        T_s = T_f(c,t);
    }
    return h_v(c,t) * alpha_sf(c,t) * (T_s - T_f(c,t));
}
real q_boil(cell_t c,Thread* t){
    real T_s;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(temp_s)))){
        T_s = C_UDSI(c,t,temp_s);
    }
    else{
        T_s = T_f(c,t);
    }
    return alpha_sf(c,t) * Viscosity_m(c,t) * H_fg(c,t) * sqrt(9.81 * (Rho_l(c,t) - Rho_v(c,t)) / model.sigma) * pow(Specific_Heat_l(c,t) * (T_s - T_sat(c,t)) / (model.c_sf * H_fg(c,t) * Pr_l(c,t)), 3);
}


// derivative
real dRel_dT(cell_t c,Thread* t){
    real NV_VEC(psi);
    real v;
    if(NNULLP(THREAD_STORAGE(t,SV_U))){
        NV_D(psi, =, C_U(c,t), C_V(c,t), C_W(c,t));
        v = NV_MAG(psi);
    }
    else{
        v = 0;
    }
    return -v * Particle_Diameter(c,t) / pow(Kinematic_Viscosity_l(c,t),2) * dnul_dT(c,t);
}
real dRev_dT(cell_t c,Thread* t){
    real NV_VEC(psi);
    real v;
    if(NNULLP(THREAD_STORAGE(t,SV_U))){
        NV_D(psi, =, C_U(c,t), C_V(c,t), C_W(c,t));
        v = NV_MAG(psi);
    }
    else{
        v = 0;
    }
    reteurn -v * Particle_Diameter(c,t) / pow(Kinematic_Viscosity_v(c,t),2) * dnuv_dT(c,t);
}
real dPrl_dT(cell_t c,Thread* t){
    return (dmul_dT(c,t) * Specific_Heat_l(c,t) + Viscosity_l(c,t) * dcpl_dT(c,t) - Pr_l(c,t) * dkl_dT(c,t)) / Conductivity_l(c,t);
}
real dPrv_dT(cell_t c,Thread* t){
    return (dmuv_dT(c,t) * Specific_Heat_v(c,t) + Viscosity_v(c,t) * dcpv_dT(c,t) - Pr_v(c,t) * dkv_dT(c,t)) / Conductivity_v(c,t);
}
real dNul_dT(cell_t c,Thread* t){
    return 1.1 * (0.33 * pow(Pr_l(c,t),-0.67) * dPrl_dT(c,t) * pow(Re_l(c,t),0.6) + 0.6 * pow(Pr_l(c,t),0.33) * pow(Re_l(c,t),-0.4) * dRel_dT(c,t));
}
real dNuv_dT(cell_t c,Thread* t){
    return 1.1 * (0.33 * pow(Pr_v(c,t),-0.67) * dPrv_dT(c,t) * pow(Re_v(c,t),0.6) + 0.6 * pow(Pr_v(c,t),0.33) * pow(Re_v(c,t),-0.4) * dRev_dT(c,t));
}
real dhl_dT(cell_t c,Thread* t){
    return (dNul_dT(c,t) * Conductivity_l(c,t) + Nu_l(c,t) * dkl_dT(c,t)) / Particle_Diameter(c,t);
}
real dhv_dT(cell_t c,Thread* t){
    return (dNuv_dT(c,t) * Conductivity_v(c,t) + Nu_v(c,t) * dkv_dT(c,t)) / Particle_Diameter(c,t);
}
real dql_dT(cell_t c,Thread* t){
    real T_s;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(temp_s)))){
        T_s = C_UDSI(c,t,temp_s);
    }
    else{
        T_s = T_f(c,t);
    }
    return dhl_dT(c,t) * alpha_sf(c,t) * (T_s - T_f(c,t)) - h_l(c,t) * alpha_sf(c,t);
}
real dqv_dT(cell_t c,Thread* t){
    real T_s;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(temp_s)))){
        T_s = C_UDSI(c,t,temp_s);
    }
    else{
        T_s = T_f(c,t);
    }
    return dhv_dT(c,t) * alpha_sf(c,t) * (T_s - T_f(C,t)) - h_v(c,t) * alpha_sf(c,t);
}

real dqboil_dS(cell_t c,Thread* t){
    real T_s;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(temp_s)))){
        T_s = C_UDSI(c,t,temp_s);
    }
    else{
        T_s = T_f(c,t);
    }
    return alpha_sf(c,t) * dnu_dS(c,t) * H_fg(c,t) * sqrt(9.81 * (Rho_l(c,t) - Rho_v(c,t)) / model.sigma) * pow(Specific_Heat_l(c,t) * (T_s - T_sat(c,t)) / (model.c_sf * H_fg(c,t) * Pr_l(c,t)) , 3)
}

// BC
real inlet_enthalpy_l(face_t f,Thread* t){
	real q_in, mass_in, m_flux;

	q_in = RP_Get_Real("myudf/heat");
	mass_in = RP_Get_Real("myudf/mass");
	m_flux = mass_in / (pow(model.D,2) / 4);

    return (model.reservoir_enthalpy + h_inlet(f,t) / m_flux * Solid_temp(f,t) - Conductivity_m(f,t) / m_flux / Specific_Heat_l(f,t) * Gradient_Enthalpy(f,t)) / (1 + h_inlet(f,t) / m_flux / Specific_Heat_l(f,t));
}
real inlet_enthalpy_v(face_t f,Thread* t){
    real q_in, mass_in, m_flux;

	q_in = RP_Get_Real("myudf/heat");
	mass_in = RP_Get_Real("myudf/mass");
	m_flux = mass_in / (pow(model.D,2) / 4);

    return (model.reservoir_enthalpy + h_inlet(f,t) / m_flux * (Solid_temp(f,t) - T_sat(f,t) + H_sat_v(f,t) / Specific_Heat_v(f,t)) - Conductivity_m(f,t) * Porosity(f,t) / Specific_Heat_v(f,t) * Gradient_Enthalpy(f,t)) / (1+h_inlet(f,t) / m_flux / Specific_Heat_v(f,t));
}
real inlet_enthalpy_m(face_t f,Thread* t){
    real q_in, mass_in, m_flux;

	q_in = RP_Get_Real("myudf/heat");
	mass_in = RP_Get_Real("myudf/mass");
	m_flux = mass_in / (pow(model.D,2) / 4);

    return (model.reservoir_enthalpy + h_inlet(f,t) / m_flux * (Solid_temp(f,t) - T_sat(f,t)) - Conductivity_m(f,t) * Porosity(f,t) / m_flux * dTsat_dP(f,t) * Gradient_Pressure(f,t)) / beta_current(f,t);
}

// face
real Re_l(face_t f,Thread* t){
    real NV_VEC(psi);
    real v;
    if(NNULLP(THREAD_STORAGE(t,SV_U))){
        NV_D(psi, =, F_U(f,t), F_V(f,t), C_W(f,t));
        v = NV_MAG(psi);
    }
    else{
        v = 0;
    }
    return v * Particle_Diameter(f,t) / Kinematic_Viscosity_l(f,t);
}
real Re_v(face_t f, Thread* t){
    real NV_VEC(psi);
    real v;
    if(NNULLP(THREAD_STORAGE(t,SV_U))){
        NV_D(psi, =, F_U(f,t), F_V(f,t), C_W(f,t));
        v = NV_MAG(psi);
    }
    else{
        v = 0;
    }
    return v * Particle_Diameter(f,t) / Kinematic_Viscosity_v(f,t);
}
real Pr_l(face_t f, Thread* t){
    return Viscosity_l(f,t) * Specific_Heat_l(f,t) / Conductivity_l(f,t);
}
real Pr_v(face_t f, Thread* t){
    return Viscosity_v(f,t) * Specific_heat_v(f,t) / Conductivity_v(f,t);
}
real Nu_l(face_t f, Thread* t){
    return 2 + 1.1 * pow(Pr_l(f,t), 0.33) * pow(Re_l(f,t), 0.6);
}
real Nu_v(face_t f, Thread* t){
    return 2 + 1.1 * pow(Pr_v(f,t),0.33) * pow(Re_v(f,t), 0.6);
}
real h_l(face_t f, Thread* t){
    return Nu_l(f,t) * Conductivity_l(f,t) / Particle_Diameter(f,t);
}
real h_v(face_t f, Thread* t){
    return Nu_v(f,t) * Conductivity_v(f,t) / Particle_Diameter(f,t);
}