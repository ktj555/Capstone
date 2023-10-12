#include "heattransfer.h"
#include "property.h"
#include "model.h"
#include "info.h"

extern constant models;

real alpha_sf(cell_t c, Thread* t){
    return 6 * (1-Porosity(c,t)) / Particle_Diameter(c,t);
}

real Re_l(cell_t c,Thread* t){
    real NV_VEC(psi);
    real v;
    real mass_in, m_flux;

    mass_in = RP_Get_Real("myudf/mass");
	m_flux = mass_in / (pow(models.D,2) * M_PI / 4);

    if(NNULLP(THREAD_STORAGE(t,SV_U))){
        NV_D(psi, =, C_U(c,t), C_V(c,t), C_W(c,t));
        v = NV_MAG(psi);
    }
    else{
        v = m_flux / Rho_l(c,t);
    }
    return v * Particle_Diameter(c,t) / Kinematic_Viscosity_l(c,t);
}
real Re_v(cell_t c, Thread* t){
    real NV_VEC(psi);
    real v;
    real mass_in, m_flux;

    mass_in = RP_Get_Real("myudf/mass");
	m_flux = mass_in / (pow(models.D,2) * M_PI / 4);

    if(NNULLP(THREAD_STORAGE(t,SV_U))){
        NV_D(psi, =, C_U(c,t), C_V(c,t), C_W(c,t));
        v = NV_MAG(psi);
    }
    else{
        v = m_flux / Rho_v(c,t);
    }
    return v * Particle_Diameter(c,t) / Kinematic_Viscosity_v(c,t);
}
real Pr_l(cell_t c, Thread* t){
    return Viscosity_l(c,t) * Specific_Heat_l(c,t) / Conductivity_l(c,t);
}
real Pr_v(cell_t c, Thread* t){
    return Viscosity_v(c,t) * Specific_Heat_v(c,t) / Conductivity_v(c,t);
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
    return alpha_sf(c,t) * Viscosity_m(c,t) * H_fg(c,t) * sqrt(9.81 * (Rho_l(c,t) - Rho_v(c,t)) / models.sigma) * pow(Specific_Heat_l(c,t) * (T_s - T_sat(c,t)) / (models.c_sf * H_fg(c,t) * Pr_l(c,t)), 3);
}


// derivative
real dRel_dT(cell_t c,Thread* t){
    real NV_VEC(psi);
    real v;
    real mass_in, m_flux;
    mass_in = RP_Get_Real("myudf/mass");
	m_flux = mass_in / (pow(models.D,2) * M_PI / 4);
    if(NNULLP(THREAD_STORAGE(t,SV_U))){
        NV_D(psi, =, C_U(c,t), C_V(c,t), C_W(c,t));
        v = NV_MAG(psi);
    }
    else{
        v = m_flux / Rho_l(c,t);
    }
    return -v * Particle_Diameter(c,t) / pow(Kinematic_Viscosity_l(c,t),2) * dnul_dT(c,t);
}
real dRev_dT(cell_t c,Thread* t){
    real NV_VEC(psi);
    real v;
    real mass_in, m_flux;
    mass_in = RP_Get_Real("myudf/mass");
	m_flux = mass_in / (pow(models.D,2) * M_PI / 4);
    if(NNULLP(THREAD_STORAGE(t,SV_U))){
        NV_D(psi, =, C_U(c,t), C_V(c,t), C_W(c,t));
        v = NV_MAG(psi);
    }
    else{
        v = m_flux / Rho_v(c,t);
    }
    return -v * Particle_Diameter(c,t) / pow(Kinematic_Viscosity_v(c,t),2) * dnuv_dT(c,t);
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
    return dhv_dT(c,t) * alpha_sf(c,t) * (T_s - T_f(c,t)) - h_v(c,t) * alpha_sf(c,t);
}

real dqboil_dS(cell_t c,Thread* t){
    real T_s;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(temp_s)))){
        T_s = C_UDSI(c,t,temp_s);
    }
    else{
        T_s = T_f(c,t);
    }
    return alpha_sf(c,t) * dnu_dS(c,t) * H_fg(c,t) * sqrt(9.81 * (Rho_l(c,t) - Rho_v(c,t)) / models.sigma) * pow(Specific_Heat_l(c,t) * (T_s - T_sat(c,t)) / (models.c_sf * H_fg(c,t) * Pr_l(c,t)) , 3);
}

// derivative of T_s
real dql_dTs(cell_t c, Thread* t){
    return h_l(c,t) * alpha_sf(c,t);
}
real dqv_dTs(cell_t c, Thread* t){
    return h_v(c,t) * alpha_sf(c,t);
}
real dqboil_dTs(cell_t c, Thread* t){
    real T_s;
    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(temp_s)))){
        T_s = C_UDSI(c,t,temp_s);
    }
    else{
        T_s = T_f(c,t);
    }
    return q_boil(c,t) * 3 / (T_s - T_sat(c,t));
}

// BC
real inlet_enthalpy_l(face_t f,Thread* t){
	real q_in, mass_in, m_flux, T_s, dhdx, h_c;
    real dhdx_lim;
    cell_t c;
    Thread* tc;

	q_in = RP_Get_Real("myudf/heat");
	mass_in = RP_Get_Real("myudf/mass");
	m_flux = mass_in / (pow(models.D,2) * M_PI / 4);

    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(temp_s)))){
        T_s = F_UDSI(f,t,temp_s);
    }
    else{
        T_s = T_f_face(f,t);
    }
    c = F_C0(f,t);
    tc = THREAD_T0(t);
    if (NNULLP(T_STORAGE_R_NV(tc, SV_UDSI_G(enthalpy)))) {
        dhdx = C_UDSI_G(c, tc, enthalpy)[0];
    }
    else {
        dhdx = 0.0;
    }

    h_c = Nu_l_face(f,t) * Conductivity_l_face(f,t) / models.D;

    dhdx_lim = (Specific_Heat_l_face(f,t) * (T_s - models.T_ref) - models.reservoir_enthalpy) / (Porosity_face(f,t) * Conductivity_l_face(f,t) / Specific_Heat_l_face(f,t) / m_flux);
    if(dhdx > dhdx_lim) dhdx = dhdx_lim / 2;

    return (models.reservoir_enthalpy + h_c / m_flux * (T_s - models.T_ref) + Porosity_face(f,t) * Conductivity_l_face(f,t) / m_flux / Specific_Heat_l_face(f,t) * dhdx) / (1 + h_c / m_flux / Specific_Heat_l_face(f,t));
}
real inlet_enthalpy_v(face_t f,Thread* t){
    real q_in, mass_in, m_flux, T_s, dhdx, h_c;
    real dhdx_lim;
    cell_t c;
    Thread* tc;

	q_in = RP_Get_Real("myudf/heat");
	mass_in = RP_Get_Real("myudf/mass");
	m_flux = mass_in / (pow(models.D,2) * M_PI / 4);

    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(temp_s)))){
        T_s = F_UDSI(f,t,temp_s);
    }
    else{
        T_s = T_f_face(f,t);
    }
    c = F_C0(f,t);
    tc = THREAD_T0(t);
    if (NNULLP(T_STORAGE_R_NV(tc, SV_UDSI_G(enthalpy)))) {
        dhdx = C_UDSI_G(c, tc, enthalpy)[0];
    }
    else {
        dhdx = 0.0;
    }

    h_c = Nu_v_face(f,t) * Conductivity_v_face(f,t) / models.D;

    dhdx_lim = m_flux * (Specific_Heat_v_face(f,t) * (T_s - T_sat_face(f,t)) + H_sat_v_face(f,t) - models.reservoir_enthalpy) / (Porosity_face(f,t) * Conductivity_v_face(f,t) / Specific_Heat_v_face(f,t));
    if(dhdx > dhdx_lim) dhdx = dhdx_lim / 2;

    return (models.reservoir_enthalpy + h_c / m_flux * (T_s - T_sat_face(f,t) + H_sat_v_face(f,t) / Specific_Heat_v_face(f,t)) + Conductivity_m_face(f,t) * Porosity_face(f,t) / Specific_Heat_v_face(f,t) * dhdx) / (1+h_c / m_flux / Specific_Heat_v_face(f,t));
}
real inlet_enthalpy_m(face_t f,Thread* t){
    real q_in, mass_in, m_flux, h_inlet, T_s, dPdx, h_cl, h_cv;
    cell_t c;
    Thread* tc;

	q_in = RP_Get_Real("myudf/heat");
	mass_in = RP_Get_Real("myudf/mass");
	m_flux = mass_in / (pow(models.D,2) * M_PI / 4);

    h_cl = Nu_l_face(f,t) * Conductivity_l_face(f,t) / models.D;
    h_cv = Nu_v_face(f,t) * Conductivity_v_face(f,t) / models.D;
    h_inlet = S_face(f,t) * h_cl + (1-S_face(f,t)) * h_cv;

    if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(temp_s)))){
        T_s = F_UDSI(f,t,temp_s);
    }
    else{
        T_s = T_f_face(f,t);
    }
    c = F_C0(f,t);
    tc = THREAD_T0(t);
    if(NNULLP(THREAD_STORAGE(tc,SV_P_G))){
        dPdx = C_P_G(c,tc)[0];
    }
    else{
        dPdx = 0;
    }

    return (models.reservoir_enthalpy + h_inlet / m_flux * (T_s - T_sat_face(f,t)) + Conductivity_m_face(f,t) * Porosity_face(f,t) / m_flux * dTsat_dP_face(f,t) * dPdx) / beta_current_face(f,t);
}

// face
real Re_l_face(face_t f,Thread* t){
    real NV_VEC(psi);
    real v;
    real mass_in, m_flux;
    mass_in = RP_Get_Real("myudf/mass");
	m_flux = mass_in / (pow(models.D,2) * M_PI / 4);
    if(NNULLP(THREAD_STORAGE(t,SV_U))){
        NV_D(psi, =, F_U(f,t), F_V(f,t), F_W(f,t));
        v = NV_MAG(psi);
    }
    else{
        v = m_flux / Rho_l_face(f,t);
    }
    return v * Particle_Diameter_face(f,t) / Kinematic_Viscosity_l_face(f,t);
}
real Re_v_face(face_t f, Thread* t){
    real NV_VEC(psi);
    real v;
    real mass_in, m_flux;
    mass_in = RP_Get_Real("myudf/mass");
	m_flux = mass_in / (pow(models.D,2) * M_PI / 4);
    if(NNULLP(THREAD_STORAGE(t,SV_U))){
        NV_D(psi, =, F_U(f,t), F_V(f,t), F_W(f,t));
        v = NV_MAG(psi);
    }
    else{
        v = m_flux / Rho_v_face(f,t);
    }
    return v * Particle_Diameter_face(f,t) / Kinematic_Viscosity_v_face(f,t);
}
real Pr_l_face(face_t f, Thread* t){
    return Viscosity_l_face(f,t) * Specific_Heat_l_face(f,t) / Conductivity_l_face(f,t);
}
real Pr_v_face(face_t f, Thread* t){
    return Viscosity_v_face(f,t) * Specific_Heat_v_face(f,t) / Conductivity_v_face(f,t);
}
real Nu_l_face(face_t f, Thread* t){
    return 2 + 1.1 * pow(Pr_l_face(f,t), 0.33) * pow(Re_l_face(f,t), 0.6);
}
real Nu_v_face(face_t f, Thread* t){
    return 2 + 1.1 * pow(Pr_v_face(f,t),0.33) * pow(Re_v_face(f,t), 0.6);
}
real h_l_face(face_t f, Thread* t){
    return Nu_l_face(f,t) * Conductivity_l_face(f,t) / Particle_Diameter_face(f,t);
}
real h_v_face(face_t f, Thread* t){
    return Nu_v_face(f,t) * Conductivity_v_face(f,t) / Particle_Diameter_face(f,t);
}