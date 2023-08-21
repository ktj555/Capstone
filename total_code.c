#include "udf.h"

// User defined values
#define MY_SIGMA 0.072
#define MY_C_SF 0.006
#define INLET_TEMP_F 303.15
#define D_P 100e-6
#define DIAMETER 0.06
#define THICKNESS 0.008

// Gravity
#define G_X -9.81
#define G_Y 0

// User Define Scalar Index
#define MODIFIED_ENTHALPY 0
#define TEMP_S 1

// User Define Memory Index
#define MY_POROSITY 0
#define PERMEABILITY 1
#define SATURATION 2
#define TEMP_F 3
#define SOURCE_F 4
#define SOURCE_S 5

static real MASS_IN=0.0;
static real Q_IN=0.0;

// fluid dynamics properties for each phases
real RHO_L(real T, real P);
real dRHO_LdT(real T, real P);
real RHO_V(real T, real P);
real dRHO_VdT(real T, real P);
real RHO_S(real T);

real MU_L(real T, real P);
real dMU_LdT(real T, real P);
real MU_V(real T, real P);
real dMU_VdT(real T, real P);
real NU_L(real T, real P);
real dNU_LdT(real T, real P);
real NU_V(real T, real P);
real dNU_VdT(real T, real P);

// thermodynamics properties for each phases
real CP_L(real T, real P);
real dCP_LdT(real T, real P);
real CP_V(real T, real P);
real dCP_VdT(real T, real P);
real CP_S(real T);
real K_L(real T, real P);
real dK_LdT(real T, real P);
real K_V(real T, real P);
real dK_VdT(real T, real P);
real K_S(real T);
real T_SAT(real P);
real dT_SATdP(real P);
real H_L_SAT(real P);
real H_FG(real P);
real H_V_SAT(real P);


real H_to_T(real h, real P);
real H_to_S(real h, real P);


real K_RL(real S);
real dK_RLdS(real S);
real K_RV(real S);
real dK_RLdS(real S);

// mixture properties
real RHO(real H, real P);
real dRHOdS(real H, real P);
real NU(real H, real P);
real dNUdS(real H, real P);
real MU(real H, real P);
real dMUdS(real H, real P);
real CP_F(real H, real P);
real dCP_FdS(real H, real P);
real K_F(real H, real P);
real dK_FdS(real H, real P);

real LAMBDA_L(real H, real P);
real dLAMBDA_LdS(real H, real P);
real LAMBDA_V(real H, real P);
real dLAMBDA_VdS(real H, real P);

real MY_BETA(real H, real P);

real Re_L(real T, real P, real v);
real dRe_LdT(real T, real P, real v);
real Re_V(real T, real P, real v);
real dRe_VdT(real T, real P, real v);
real Pr_L(real T, real P);
real dPr_LdT(real T, real P);
real Pr_V(real T, real P);
real dPr_VdT(real T, real P);
real Nu_L(real T, real P, real v);
real dNu_LdT(real T, real P, real v);
real Nu_V(real T, real P, real v);
real dNu_VdT(real T, real P, real v);
real h_sl(real T, real P, real v);
real dh_sldT(real T, real P, real v);
real h_sv(real T, real P, real v);
real dh_svdT(real T, real P, real v);
real q_l(real T, real P, real v, real alpha, real T_s);
real dq_ldT(real T, real P, real v, real alpha, real T_s);
real dq_ldT_s(real T, real P, real v, real alpha, real T_s);
real q_v(real T, real P, real v, real alpha, real T_s);
real dq_vdT(real T, real P, real v, real alpha, real T_s);
real dq_vdT_s(real T, real P, real v, real alpha, real T_s);
real q_boil(real H, real P, real alpha, real T_s);
real dq_boildS(real H, real P, real alpha, real T_s);
real dq_boildT_s(real H, real P, real alpha, real T_s);
real qsf(real H, real P, real v, real alpha, real T_s);
real dqsfdH(real H, real P, real v, real alpha, real T_s);
real dqsfdT_s(real H, real P, real v, real alpha, real T_s);

real h_inject(real H, real P, real v);
void check_list_message();

// Set UDS & UDM names for GUI
DEFINE_EXECUTE_ON_LOADING(set_load, libname)
{
	Set_User_Scalar_Name(MODIFIED_ENTHALPY, "Modified_enthalpy");
	Set_User_Scalar_Name(TEMP_S, "solid_temperature");

	Set_User_Memory_Name(MY_POROSITY, "local_porosity");
	Set_User_Memory_Name(PERMEABILITY, "local_permeability");
	Set_User_Memory_Name(TEMP_F, "fluid_temperature");
	Set_User_Memory_Name(SATURATION, "Saturation");
	Set_User_Memory_Name(SOURCE_F, "Source_in_fiuld");
	Set_User_Memory_Name(SOURCE_S, "Source_in_solid");
}

DEFINE_EXECUTE_FROM_GUI(mass_and_heat,libname,mode){
    MASS_IN = RP_Get_Real("myudf/mass");
    Message("Set mass flow rate : %g [kg/s]\n",MASS_IN);
    Q_IN = RP_Get_Real("myudf/heat");
    Message("Set heat flux : %g [W/m^2]\n",Q_IN);
}

DEFINE_INIT(initialize_saturation_and_temperature, d)
{
#if !RP_HOST
	cell_t c;
	Thread* t;

	real q_in, mass_in;
	real h_in, h_out, A;
	real x;
	real cen[ND_ND];

	q_in = RP_Get_Real("myudf/heat");
	mass_in = RP_Get_Real("myudf/mass");
	A = M_PI * pow(DIAMETER, 2) / 4;
	h_in = CP_L(INLET_TEMP_F, 0) * INLET_TEMP_F;
	h_out = h_in + q_in/(mass_in/A);


	thread_loop_c(t, d) {
		begin_c_loop_all(c, t)
		{
			C_CENTROID(cen,c,t);
			x = cen[0];
			C_UDSI(c, t, MODIFIED_ENTHALPY) = h_in + (h_out - h_in) / THICKNESS * (x + THICKNESS * 0.5);
			C_UDSI(c, t, TEMP_S) = H_to_T(h_in + (h_out - h_in) / THICKNESS * (x + THICKNESS * 0.5), 0);

			C_UDMI(c, t, MY_POROSITY) = 0.3;
			C_UDMI(c, t, PERMEABILITY) = 3.67e-12;

			C_UDMI(c, t, SATURATION) = 0;
			C_UDMI(c, t, TEMP_F) = 0;
			C_UDMI(c, t, SOURCE_F) = 0;
			C_UDMI(c, t, SOURCE_S) = 0;
		}
		end_c_loop_all(c, t)
	}
#endif
}

DEFINE_ADJUST(adjust_variables, d) 
{
#if !RP_HOST
	cell_t c;
	Thread* t;

	real P;
	real h_m;

	thread_loop_c(t, d) {
		begin_c_loop(c, t)
		{
			h_m = C_UDSI(c, t, MODIFIED_ENTHALPY);
			P = C_P(c, t);
			C_UDMI(c, t, TEMP_F) = H_to_T(h_m, P);
			C_UDMI(c, t, SATURATION) = H_to_S(h_m, P);
		}
		end_c_loop(c, t)
	}

#endif
}

DEFINE_PROFILE(inlet_enthalpy, t, i)
{
	Thread* tc;
	face_t f;
	cell_t c;

	real h_m, h_in, P, v, T_s, T;
	real m_flux, dhdx, dPdx;
	real h_c;
	real e, k_f_eff;
	real h_m_s0, h_m_s1;

	m_flux = MASS_IN / (M_PI * pow(DIAMETER, 2) / 4);

	begin_f_loop(f, t)
	{
		c = F_C0(f, t);
		tc = THREAD_T0(t);
		
		if(NNULLP(THREAD_STORAGE(t,SV_P))){
			P = F_P(f, t);
		}
		else{
			P = 0;
		}
		if(NNULLP(THREAD_STORAGE(tc,SV_P_G))){
			dPdx = C_P_G(c,tc)[0];
		}
		else{
			dPdx = 0;
		}
		if(NNULLP(THREAD_STORAGE(t,SV_U))){
			v = F_U(f,t);
		}
		else{
			v = m_flux / 998.0;
		}
		if (NNULLP(THREAD_STORAGE(t, SV_UDS_I(TEMP_S)))) {
			T_s = F_UDSI(f, t, TEMP_S);
		}
		else {
			T_s = INLET_TEMP_F;
		}
		if (NNULLP(THREAD_STORAGE(t, SV_UDM_I))) {
			e = F_UDMI(f, t, MY_POROSITY);
		}
		else {
			e = 1.0;
		}
		if (NNULLP(THREAD_STORAGE(t, SV_UDS_I(MODIFIED_ENTHALPY)))) {
			h_m = F_UDSI(f, t, MODIFIED_ENTHALPY);
		}
		else {
			h_m = h_in;
		}
		if (NNULLP(T_STORAGE_R_NV(tc, SV_UDSI_G(MODIFIED_ENTHALPY)))) {
			dhdx = C_UDSI_G(c, tc, MODIFIED_ENTHALPY)[0];
		}
		else {
			dhdx = 0.0;
		}
		h_in = CP_L(INLET_TEMP_F, P) * INLET_TEMP_F;
		T = H_to_T(h_m, P);
		k_f_eff = e * K_F(h_m, P);
		h_c = h_inject(h_m, P, v);
		h_m_s0 = H_V_SAT(P);
		h_m_s1 = H_L_SAT(P);

		if (h_m <= h_m_s1) {
			F_PROFILE(f, t, i) = (h_in + h_c / m_flux * T_s - k_f_eff / m_flux / CP_L(T, P) * dhdx) / (1 + h_c / m_flux / CP_L(T, P));
		}
		else if (h_m >= h_m_s0) {
			F_PROFILE(f, t, i) = (h_in + h_c / m_flux * (T_s - T_SAT(P) + H_V_SAT(P) / CP_V(T, P)) - k_f_eff / CP_V(T, P) * dhdx) / (1 + h_c / m_flux / CP_V(T, P));
		}
		else {
			F_PROFILE(f, t, i) = (h_in + h_c / m_flux * (T_s - T_SAT(P)) - k_f_eff / m_flux * dT_SATdP(P) * dPdx) / MY_BETA(h_m,P);
		}
	}
	end_f_loop(f, t)
}


DEFINE_PROFILE(inlet_temp_s_flux, t, i)
{
	face_t f;

	real h_m, P, v, T_s, T;
	real m_flux;
	real h_c;
	real e, k_s_eff;

	m_flux = MASS_IN / (M_PI * pow(DIAMETER, 2) / 4);

	begin_f_loop(f, t)
	{
		if(NNULLP(THREAD_STORAGE(t,SV_P))){
			P = F_P(f, t);
		}
		else{
			P = 0;
		}
		if(NNULLP(THREAD_STORAGE(t,SV_U))){
			v = F_U(f,t);
		}
		else{
			v = m_flux / 998.0;
		}
		if (NNULLP(THREAD_STORAGE(t, SV_UDS_I(TEMP_S)))) {
			T_s = F_UDSI(f, t, TEMP_S);
		}
		else {
			T_s = INLET_TEMP_F;
		}
		if (NNULLP(THREAD_STORAGE(t, SV_UDM_I))) {
			e = F_UDMI(f, t, MY_POROSITY);
		}
		else {
			e = 1.0;
		}
		if (NNULLP(THREAD_STORAGE(t, SV_UDS_I(MODIFIED_ENTHALPY)))) {
			h_m = F_UDSI(f, t, MODIFIED_ENTHALPY);
		}
		else {
			h_m = CP_L(INLET_TEMP_F, P) * INLET_TEMP_F;
		}
		T = H_to_T(h_m, P);
		h_c = h_inject(h_m, P, v);
		k_s_eff = (1 - e) * K_S(T_s);

		F_PROFILE(f, t, i) = -h_c * (T_s - T);
	}
	end_f_loop(f, t)
}

DEFINE_PROFILE(inlet_velocity, t, i)
{
	face_t f;

	real T, P, h_m, rho;
	real h_m_s0, h_m_s1;

	begin_f_loop(f, t)
	{
		if(NNULLP(THREAD_STORAGE(t,SV_P))){
			P = F_P(f,t);
		}
		else{
			P = 0;
			inlet_velocity_check_list[0] = 1;
		}
		if (NNULLP(THREAD_STORAGE(t, SV_UDS_I(MODIFIED_ENTHALPY)))) {
			h_m = F_UDSI(f,t,MODIFIED_ENTHALPY);
			T = H_to_T(h_m,P);
		}
		else {
			h_m = INLET_TEMP_F * CP_L(INLET_TEMP_F,P);
			T = INLET_TEMP_F;
		}
		h_m_s1 = H_L_SAT(P);
		h_m_s0 = H_V_SAT(P);

		if (h_m <= h_m_s1) {
			rho = RHO_L(T,P);
		}
		else if (h_m >= h_m_s0) {
			rho = RHO_V(T,P);
		}
		else {
			rho = RHO(h_m,P);
		}
		F_PROFILE(f, t, i) = MASS_IN / (M_PI * pow(DIAMETER, 2) / 4) / rho;
	}
	end_f_loop(f, t)
}

DEFINE_PROFILE(heat_flux, t, i)
{
	face_t f;

	begin_f_loop(f, t)
	{
		F_PROFILE(f, t, i) = Q_IN;
	}
	end_f_loop(f, t)
}

DEFINE_PROFILE(porosity, t, i) {
	cell_t c;
	begin_c_loop(c, t) {
		if (NNULLP(THREAD_STORAGE(t, SV_UDM_I))) {
			C_PROFILE(c, t, i) = C_UDMI(c, t, MY_POROSITY);
		}
		else {
			C_PROFILE(c, t, i) = 1.0;
		}
	}
	end_c_loop(c,t)
}

DEFINE_PROFILE(permeability_resistence, t, i) {
	cell_t c;
	begin_c_loop(c, t) {
		if (NNULLP(THREAD_STORAGE(t, SV_UDM_I))) {
			C_PROFILE(c, t, i) = 1 / C_UDMI(c, t, PERMEABILITY);
		}
		else {
			C_PROFILE(c, t, i) = 2.11e10;
		}
	}
	end_c_loop(c,t)
}


DEFINE_PROPERTY(Mixture_Rho, c, t) {
	real H, T, P;
	P = C_P(c, t);
	H = C_UDSI(c, t, MODIFIED_ENTHALPY);
	return RHO(H, P);
}
DEFINE_PROPERTY(Mixture_viscsity, c, t) {
	real H, P;
	P = C_P(c, t);
	H = C_UDSI(c, t, MODIFIED_ENTHALPY);
	return MU(H, P);
}

DEFINE_UDS_UNSTEADY(Unsteady_for_solid, c, t, i, apu, su)
{
	real physical_dt, vol, rho, rho_old, phi_old, epsilon, c_p, c_p_old;

	physical_dt = RP_Get_Real("physical-time-step");
	vol = C_VOLUME(c, t);
	epsilon = C_UDMI(c, t, MY_POROSITY);

	c_p = CP_S(C_UDSI(c, t, TEMP_S));
	c_p_old = CP_S(C_UDSI_M1(c, t, TEMP_S));
	rho = RHO_S(C_UDSI(c, t, TEMP_S));
	rho_old = RHO_S(C_UDSI_M1(c, t, TEMP_S));

	phi_old = C_STORAGE_R(c, t, SV_UDSI_M1(TEMP_S));

	*apu = -rho * vol * (1 - epsilon) * c_p / physical_dt;
	*su = rho_old * vol * (1 - epsilon) * c_p_old * phi_old / physical_dt;
}
DEFINE_UDS_UNSTEADY(Unsteady_for_fluid, c, t, i, apu, su)
{
	real physical_dt, vol, rho, rho_old, phi_old, epsilon, beta, beta_old;

	physical_dt = RP_Get_Real("physical-time-step");
	vol = C_VOLUME(c, t);
	epsilon = C_UDMI(c, t, MY_POROSITY);

	beta = MY_BETA(C_UDSI(c, t, MODIFIED_ENTHALPY), C_P(c, t));
	beta_old = MY_BETA(C_UDSI_M1(c, t, MODIFIED_ENTHALPY), C_P_M1(c, t));
	rho = C_R(c, t);
	rho_old = C_R_M1(c, t);

	phi_old = C_STORAGE_R(c, t, SV_UDSI_M1(MODIFIED_ENTHALPY));

	*apu = -epsilon * rho * vol * beta / physical_dt;
	*su = epsilon * rho_old * vol * beta_old * phi_old / physical_dt;
}

DEFINE_DIFFUSIVITY(Diffusivity_for_solid, c, t, i)
{
	real epsilon, k_s;
	epsilon = C_UDMI(c, t, MY_POROSITY);
	k_s = K_S(C_UDSI(c, t, TEMP_S));
	return k_s * (1 - epsilon);
}
DEFINE_DIFFUSIVITY(Diffusivity_for_fluid, c, t, i)
{
	real T_sat, D, h_m_s0, h_m_s1, e;
	real h_m, P, T, S;
	T_sat = T_SAT(C_P(c, t));
	h_m_s1 = H_L_SAT(C_P(c, t));
	h_m_s0 = H_V_SAT(C_P(c, t));
	P = C_P(c, t);
	h_m = C_UDSI(c, t, MODIFIED_ENTHALPY);
	T = H_to_T(h_m, P);
	S = H_to_S(h_m, P);
	e = C_UDMI(c, t, MY_POROSITY);

	if (h_m <= h_m_s1) {
		return e * K_F(h_m, P) / CP_L(T, P);
	}
	else if (h_m >= h_m_s0) {
		return e * K_F(h_m, P) / CP_V(T, P);
	}
	else {
		real nu, lambda, K, dlds, dnuds;
		K = C_UDMI(c, t, PERMEABILITY);
		nu = NU(h_m, P);
		lambda = LAMBDA_L(h_m, P);
		D = K / nu * lambda * (1 - lambda) * sqrt(e / K) * MY_SIGMA * (1.417 - 2 * 2.120 * (1 - S) + 3 * 1.263 * pow(1 - S, 2));
		dnuds = 3 * pow(nu, 2) * (pow(1 - S, 2) / NU_V(T, P) - pow(S, 2) / NU_L(T, P));
		dlds = (dnuds * pow(S, 3) + 3 * nu * pow(S, 2)) / NU_L(T, P);
		return D / dlds;
	}
}

DEFINE_SOURCE(Source_for_solid, c, t, dS, eqn)
{
	real NV_VEC(psi);
	real alpha_sf, T_s;
	real source;
	real h_m, h_m_s0, h_m_s1, P, v;

	P = C_P(c, t);
	h_m = C_UDSI(c, t, MODIFIED_ENTHALPY);
	T_s = C_UDSI(c, t, TEMP_S);
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);

	alpha_sf = 6 * (1 - C_UDMI(c, t, MY_POROSITY)) / D_P;
	NV_D(psi, =, C_U(c, t), C_V(c, t), C_W(c, t));
	v = NV_MAG(psi);

	source = -qsf(h_m, P, v, alpha_sf, T_s);
	dS[eqn] = -dqsfdT_s(h_m, P, v, alpha_sf, T_s);
	C_UDMI(c, t, SOURCE_S) = source;
	return source;
}
DEFINE_SOURCE(Source_for_fluid, c, t, dS, eqn)
{
	real NV_VEC(psi);
	real alpha_sf, T_s;
	real source;
	real h_m, h_m_s0, h_m_s1, P, v;

	P = C_P(c, t);
	h_m = C_UDSI(c, t, MODIFIED_ENTHALPY);
	T_s = C_UDSI(c, t, TEMP_S);
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);

	alpha_sf = 6 * (1 - C_UDMI(c, t, MY_POROSITY)) / D_P;
	NV_D(psi, =, C_U(c, t), C_V(c, t), C_W(c, t));
	v = NV_MAG(psi);

	source = qsf(h_m, P, v, alpha_sf, T_s);
	dS[eqn] = dqsfdH(h_m, P, v, alpha_sf, T_s);
	C_UDMI(c, t, SOURCE_F) = source;
	return source;
}


real RHO_L(real T, real P) {
	return 998.0;
}
real dRHO_LdT(real T, real P) {
	return 0;
}
real RHO_V(real T, real P) {
	real R = 8.314, M_V = 18;
	real p = P + RP_Get_Real("operating-pressure");
	return p / (T * R / M_V * 1000);
}
real dRHO_VdT(real T, real P) {
	return -RHO_V(T, P) / T;
}
real RHO_S(real T) {
	return 8.4e3;
}
real MU_L(real T, real P) {
	return 24.141e-6 * pow(10, 247.8 / (T - 140));
}
real dMU_LdT(real T, real P) {
	return -MU_L(T, P) * 247.8 / pow(T - 140, 2) * log(10);
}
real MU_V(real T, real P) {
	real a0 = -2.77567e-6;
	real a1 = 40.35e-9;
	return a0 + a1 * T;
}
real dMU_VdT(real T, real P) {
	real a1 = 40.35e-6;
	return a1;
}
real NU_L(real T, real P) {
	return MU_L(T, P) / RHO_L(T, P);
}
real dNU_LdT(real T, real P) {
	return (dMU_LdT(T, P) * RHO_L(T, P) - MU_L(T, P) * dRHO_LdT(T, P)) / pow(RHO_L(T, P), 2);
}
real NU_V(real T, real P) {
	return MU_V(T, P) / RHO_V(T, P);
}
real dNU_VdT(real T, real P) {
	return (dMU_VdT(T, P) * RHO_V(T, P) - MU_V(T, P) * dRHO_VdT(T, P)) / pow(RHO_V(T, P), 2);
}
real CP_L(real T, real P) {
	return 4182.0;
}
real dCP_LdT(real T, real P) {
	return 0;
}
real CP_V(real T, real P) {
	return 2030.0;
}
real dCP_VdT(real T, real P) {
	return 0;
}
real CP_S(real T) {
	return 625;
}
real K_L(real T, real P) {
	return 0.68;
}
real dK_LdT(real T, real P) {
	return 0;
}
real K_V(real T, real P) {
	real a0 = -21.994433e-3;
	real a1 = 118.42e-6;
	return a0 + a1 * T;
}
real dK_VdT(real T, real P) {
	real a1 = 118.42e-6;
}
real K_S(real T) {
	return 21.7;
}
real T_SAT(real P) {
	real p = (P + RP_Get_Real("operating-pressure")) / 1e6;
	real a0 = 429.69474687605754;
	real a1 = 11.983108526790891;
	real a2 = -0.2940193901122584;
	real a3 = 25.197563642164386;
	return a0 + a1 * p + a2 * pow(p, 2) + a3 * log(p);
}
real dT_SATdP(real P){
	real p = (P + RP_Get_Real("operating-pressure")) / 1e6;
	real a1 = 11.983108526790891;
	real a2 = -0.2940193901122584;
	real a3 = 25.197563642164386;
	return (a1 + 2 * a2 * p + a3 / p) / 1e6;
}
real H_L_SAT(real P) {
	return CP_L(T_SAT(P), P) * T_SAT(P);
}
real H_FG(real P) {
	real T;
	T = T_SAT(P);
	if (T > 650) { return 0; }
	real a0 = -2511.1589762597114;
	real a1 = 1.790037469059945;
	real a2 = 0.0028969464729942194;
	real a3 = 222.03437607975306;
	return a0 + a1 * T + a2 * pow(T, 2) + a3 * sqrt(650 - T);
}
real H_V_SAT(real P) {
	return H_FG(P) + H_L_SAT(P);
}
real H_to_T(real h, real P) {
	real h_m_s0, h_m_s1, T;
	h_m_s1 = H_L_SAT(P);
	h_m_s0 = H_V_SAT(P);
	T = T_SAT(P);
	// if cp_l is not constant, need modification
	if (h < h_m_s1) { return h / CP_L(T, P); }
	else if (h > h_m_s0) { return T + (h - h_m_s0) / CP_V(T, P); }
	else { return T; }
}
real H_to_S(real h, real P) {
	real h_m_s0, h_m_s1, T;
	h_m_s1 = H_L_SAT(P);
	h_m_s0 = H_V_SAT(P);
	T = T_SAT(P);
	if (h <= h_m_s1) { return 1; }
	else if (h >= h_m_s0) { return 0; }
	else {
		real nu_r, h_r;
		nu_r = NU_V(T, P) / NU_L(T, P);
		h_r = H_FG(P) / (h_m_s0 - h);
		return 1 / (pow(nu_r * (h_r - 1), 0.33) + 1);
	}
}
real K_RL(real S) {
	return pow(S, 3);
}
real dK_RLdS(real S) {
	return 3 * pow(S, 2);
}
real K_RV(real S) {
	return pow(1 - S, 3);
}
real dK_RVdS(real S) {
	return -3 * pow(1 - S, 2);
}
real RHO(real H, real P) {
	real T, S;
	real h_m_s0, h_m_s1;
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);
	T = H_to_T(H, P);
	S = H_to_S(H, P);
	if (H <= h_m_s1) {
		return RHO_L(T, P);
	}
	else if (H >= h_m_s0) {
		return RHO_V(T, P);
	}
	else {
		return S * RHO_L(T, P) + (1 - S) * RHO_V(T, P);
	}
}
real dRHOdS(real H, real P) {
	real T, S;
	real h_m_s0, h_m_s1;
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);
	T = H_to_T(H, P);
	S = H_to_S(H, P);
	if (H <= h_m_s1) {
		return 0;
	}
	else if (H >= h_m_s0) {
		return 0;
	}
	else {
		return RHO_L(T, P) - RHO_V(T, P);
	}
}
real NU(real H, real P) {
	real T, S;
	real h_m_s0, h_m_s1;
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);
	T = H_to_T(H, P);
	S = H_to_S(H, P);
	if (H <= h_m_s1) {
		return NU_L(T, P);
	}
	else if (H >= h_m_s0) {
		return NU_V(T, P);
	}
	else {
		return 1 / (K_RL(S) / NU_L(T, P) + K_RV(S) / NU_V(T, P));
	}
}
real dNUdS(real H, real P) {
	real T, S;
	real h_m_s0, h_m_s1;
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);
	T = H_to_T(H, P);
	S = H_to_S(H, P);
	if (H <= h_m_s1) {
		return 0;
	}
	else if (H >= h_m_s0) {
		return 0;
	}
	else {
		return -pow(NU(H, P), 2) * (dK_RLdS(S) / NU_L(T, P) + dK_RVdS(S) / NU_V(T, P));
	}
}
real MU(real H, real P) {
	return RHO(H, P) * NU(H, P);
}
real dMUdS(real H, real P) {
	return dRHOdS(H, P) * NU(H, P) + RHO(H, P) * dNUdS(H, P);
}
real CP_F(real H, real P) {
	real T, S;
	real h_m_s0, h_m_s1;
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);
	T = H_to_T(H, P);
	S = H_to_S(H, P);
	if (H <= h_m_s1) {
		return CP_L(T, P);
	}
	else if (H >= h_m_s0) {
		return CP_V(T, P);
	}
	else {
		return (S * RHO_L(T, P) * CP_L(T, P) + (1 - S) * RHO_V(T, P) * CP_V(T, P)) / RHO(H, P);
	}
}
real dCP_FdS(real H, real P) {
	real T, S;
	real h_m_s0, h_m_s1;
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);
	T = H_to_T(H, P);
	S = H_to_S(H, P);
	if (H <= h_m_s1) {
		return 0;
	}
	else if (H >= h_m_s0) {
		return 0;
	}
	else {
		return (RHO_L(T, P) * CP_L(T, P) - RHO_V(T, P) * CP_V(T, P) - CP_F(H, P) * dRHOdS(H, P)) / RHO(H, P);
	}
}
real K_F(real H, real P) {
	real T, S;
	real h_m_s0, h_m_s1;
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);
	T = H_to_T(H, P);
	S = H_to_S(H, P);
	if (H <= h_m_s1) {
		return K_L(T, P);
	}
	else if (H >= h_m_s0) {
		return K_V(T, P);
	}
	else {
		return S * K_L(T, P) + (1 - S) * K_V(T, P);
	}
}
real dK_FdS(real H, real P) {
	real T, S;
	real h_m_s0, h_m_s1;
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);
	T = H_to_T(H, P);
	S = H_to_S(H, P);
	if (H <= h_m_s1) {
		return 0;
	}
	else if (H >= h_m_s0) {
		return 0;
	}
	else {
		return K_L(T, P) - K_V(T, P);
	}
}
real LAMBDA_L(real H, real P) {
	real T;
	T = H_to_T(H,P);
	return NU(H, P) / NU_L(T, P) * K_RL(S);
}
real dLAMBDA_LdS(real H, real P) {
	real T, S;
	real h_m_s0, h_m_s1;
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);
	T = H_to_T(H, P);
	S = H_to_S(H, P);
	if (H <= h_m_s1) { return 0; }
	else if (H >= h_m_s0) { return 0; }
	else {
		return (dK_RLdS(S) * NU(H, P) + K_RL(S) * dNUdS(H, P)) / NU_L(T, P);
	}
}
real LAMBDA_V(real H, real P) {
	real T, S;
	real h_m_s0, h_m_s1;
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);
	T = H_to_T(H, P);
	S = H_to_S(H, P);
	if (H <= h_m_s1) { return 0; }
	else if (H >= h_m_s0) { return 1; }
	else {
		return NU(H, P) / NU_V(T, P) * K_RV(S);
	}
}
real dLAMBDA_VdS(real H, real P) {
	real T, S;
	real h_m_s0, h_m_s1;
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);
	T = H_to_T(H, P);
	S = H_to_S(H, P);
	if (H <= h_m_s1) { return 0; }
	else if (H >= h_m_s0) { return 0; }
	else {
		return (dK_RVdS(S) * NU(H, P) + K_RV(S) * dNUdS(H, P)) / NU_V(T, P);
	}
}
real MY_BETA(real H, real P) {
	real T_sat;
	real T, S;
	real h_m_s0, h_m_s1;
	T_sat = T_SAT(P);
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);
	T = H_to_T(H, P);
	S = H_to_S(H, P);
	if (H <= h_m_s1) {
		return 1;
	}
	else if (H >= h_m_s0) {
		return 1;
	}
	else {
		return (RHO_V(T, P) / RHO(H, P) * (1 - S) * H_FG(P) + CP_L(T_sat, P) * T_sat) / ((1 - LAMBDA_L(H, P)) * H_FG(P) + CP_L(T_sat, P) * T_sat);
	}
}
real Re_L(real T, real P, real v) {
	return v * D_P / NU_L(T, P);
}
real dRe_LdT(real T, real P, real v) {
	return -v * D_P / pow(NU_L(T, P), 2) * dNU_LdT(T, P);
}
real Re_V(real T, real P, real v) {
	return v * D_P / NU_V(T, P);
}
real dRe_VdT(real T, real P, real v) {
	return -v * D_P / pow(NU_V(T, P), 2) * dNU_VdT(T, P);
}
real Pr_L(real T, real P) {
	return MU_L(T, P) * CP_L(T, P) / K_L(T, P);
}
real dPr_LdT(real T, real P) {
	return (dMU_LdT(T, P) * CP_L(T, P) + MU_L(T, P) * dCP_LdT(T, P) - Pr_L(T, P) * dK_LdT(T, P)) / K_L(T, P);
}
real Pr_V(real T, real P) {
	return MU_V(T, P) * CP_V(T, P) / K_V(T, P);
}
real dPr_VdT(real T, real P) {
	return (dMU_VdT(T, P) * CP_V(T, P) + MU_V(T, P) * dCP_VdT(T, P) - Pr_V(T, P) * dK_VdT(T, P)) / K_V(T, P);
}
real Nu_L(real T, real P, real v) {
	return 2 + 1.1 * pow(Pr_L(T, P), 0.33) * pow(Re_L(T, P, v), 0.6);
}
real dNu_LdT(real T, real P, real v) {
	return 1.1 * (0.33 * pow(Pr_L(T, P), -0.67) * dPr_LdT(T, P) * pow(Re_L(T, P, v), 0.6) + 0.6 * pow(Pr_L(T, P), 0.33) * pow(Re_L(T, P, v), -0.4) * dRe_LdT(T, P, v));
}
real Nu_V(real T, real P, real v) {
	return 2 + 1.1 * pow(Pr_V(T, P), 0.33) * pow(Re_V(T, P, v), 0.6);
}
real dNu_VdT(real T, real P, real v) {
	return 1.1 * (0.33 * pow(Pr_V(T, P), -0.67) * dPr_VdT(T, P) * pow(Re_V(T, P, v), 0.6) + 0.6 * pow(Pr_V(T, P), 0.33) * pow(Re_V(T, P, v), -0.4) * dRe_VdT(T, P, v));
}
real h_sl(real T, real P, real v) {
	return Nu_L(T, P, v) * K_L(T, P) / D_P;
}
real dh_sldT(real T, real P, real v) {
	return (dNu_LdT(T, P, v) * K_L(T, P) + Nu_L(T, P, v) * dK_LdT(T, P)) / D_P;
}
real h_sv(real T, real P, real v) {
	return Nu_V(T, P, v) * K_V(T, P) / D_P;
}
real dh_svdT(real T, real P, real v) {
	return (dNu_VdT(T, P, v) * K_V(T, P) + Nu_V(T, P, v) * dK_VdT(T, P)) / D_P;
}
real q_l(real T, real P, real v, real alpha, real T_s) {
	return h_sl(T, P, v) * alpha * (T_s - T);
}
real dq_ldT(real T, real P, real v, real alpha, real T_s) {
	return dh_sldT(T, P, v) * alpha * (T_s - T) - h_sl(T, P, v) * alpha;
}
real dq_ldT_s(real T, real P, real v, real alpha, real T_s) {
	return h_sl(T, P, v) * alpha;
}
real q_v(real T, real P, real v, real alpha, real T_s) {
	return h_sv(T, P, v) * alpha * (T_s - T);
}
real dq_vdT(real T, real P, real v, real alpha, real T_s) {
	return dh_svdT(T, P, v) * alpha * (T_s - T) - h_sv(T, P, v) * alpha;
}
real dq_vdT_s(real T, real P, real v, real alpha, real T_s) {
	return h_sv(T, P, v) * alpha;
}
real q_boil(real H, real P, real alpha, real T_s) {
	real T, S;
	T = H_to_T(H, P);
	S = H_to_S(H, P);
	return S * alpha * MU(H, P) * H_FG(P) * sqrt(sqrt(G_X * G_X + G_Y * G_Y) * (RHO_L(T, P) - RHO_V(T, P)) / MY_SIGMA) * pow(CP_L(T, P) * (T_s - T_SAT(P)) / (MY_C_SF * H_FG(P) * Pr_L(T, P)), 3);
}
real dq_boildS(real H, real P, real alpha, real T_s) {
	return q_boil(H, S, alpha, T_s) / S + q_boil(H, S, alpha, T_s) / MU(H, P) * dMUdS(H, P);
}
real dq_boildT_s(real H, real P, real alpha, real T_s) {
	return q_boil(H, P, alpha, T_s) * 3 / (T_s - T_SAT(P));
}
real qsf(real H, real P, real v, real alpha, real T_s) {
	real T_sat;
	real T, S;
	real h_m_s0, h_m_s1;
	T_sat = T_SAT(P);
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);
	T = H_to_T(H, P);
	S = H_to_S(H, P);
	if (H <= h_m_s1) {
		return q_l(T, P, v, alpha, T_s);
	}
	else if (H >= h_m_s0) {
		return q_v(T, P, v, alpha, T_s);
	}
	else {
		return q_boil(H, P, alpha, T_s) + (1 - S) * h_sv(T, P, v) * alpha * (T_s - T_sat);
	}
}
real dqsfdH(real H, real P, real v, real alpha, real T_s) {
	real T_sat;
	real T, S;
	real h_m_s0, h_m_s1;
	T_sat = T_SAT(P);
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);
	T = H_to_T(H, P);
	S = H_to_S(H, P);
	if (H <= h_m_s1) {
		return dq_ldT(T, P, v, alpha, T_s) / (dCP_LdT(T, P) * T + CP_L(T, P));
	}
	else if (H >= h_m_s0) {
		return dq_vdT(T, P, v, alpha, T_s) / (dCP_VdT(T, P) * (T - T_sat) + CP_V(T, P));
	}
	else {
		return (dq_boildS(H, P, alpha, T_s) - h_sv(T, P, v) * alpha * (T_s - T_sat)) / (dLAMBDA_LdS(H, P) * H_FG(P));
	}
}
real dqsfdT_s(real H, real P, real v, real alpha, real T_s) {
	real T_sat;
	real T, S;
	real h_m_s0, h_m_s1;
	T_sat = T_SAT(P);
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);
	T = H_to_T(H, P);
	S = H_to_S(H, P);
	if (H <= h_m_s1) {
		return dq_ldT_s(T, P, v, alpha, T_s);
	}
	else if (H >= h_m_s0) {
		return dq_vdT_s(T, P, v, alpha, T_s);
	}
	else {
		return dq_boildT_s(H, P, alpha, T_s) - (1 - S) * h_sv(T, P, v) * alpha;
	}
}
real h_inject(real H, real P, real v) {
	real h_m_s0,h_m_s1;
	real T,S;
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);
	T = H_to_T(H,P);
	S = H_to_S(H,P);
	if(H<=h_m_s1){
		return h_sl(T,P,v);
	}
	else if(H >= h_m_s0){
		return h_sv(T,P,v);
	}
	else{
		return S * h_sl(T,P,v) + (1-S) * h_sv(T,P,v);
	}
}