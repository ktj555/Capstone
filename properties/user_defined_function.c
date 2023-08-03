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

void check_list_message(){

	int check = 0;

	if(inlet_enthalpy_check_list[0]==1){
		Message("\nAccess Error : Pressure | inlet | during inlet enthalpy\n");
		check = 1;
	}
	if(inlet_enthalpy_check_list[1]==1){
		Message("\nAccess Error : Velocity | inlet | during inlet enthalpy\n");
		check = 1;
	}
	if(inlet_enthalpy_check_list[2]==1){
		Message("\nAccess Error : Solid temperature | inlet | during inlet enthalpy\n");
		check = 1;
	}
	if(inlet_enthalpy_check_list[3]==1){
		Message("\nAccess Error : UDM | inlet | during inlet enthalpy\n");
		check = 1;
	}
	if(inlet_enthalpy_check_list[4]==1){
		Message("\nAccess Error : Modified Enthalpy | inlet | during inlet enthalpy\n");
			check = 1;
	}
	if(inlet_enthalpy_check_list[5]==1){
		Message("\nAccess Error : Enthalpy Gradient | inlet | during inlet enthalpy\n");
		check = 1;
	}
	if(inlet_temp_s_check_list[0]==1){
		Message("\nAccess Error : Pressure | inlet | during inlet flux solid\n");
		check = 1;
	}
	if(inlet_temp_s_check_list[1]==1){
		Message("\nAccess Error : Velocity | inlet | during inlet flux solid\n");
		check = 1;
	}
	if(inlet_temp_s_check_list[2]==1){
		Message("\nAccess Error : Solid temperature | inlet | during flux solid\n");
		check = 1;
	}
	if(inlet_temp_s_check_list[3]==1){
		Message("\nAccess Error : UDM | inlet | during flux solid\n");
		check = 1;
	}
	if(inlet_temp_s_check_list[4]==1){
		Message("\nAccess Error : Modified Enthalpy | inlet | during flux solid\n");
		check = 1;
	}
	if(inlet_velocity_check_list[0]==1){
		Message("\nAccess Error : Pressure | inlet | during inlet velocity\n");
		check = 1;
	}
	if(inlet_velocity_check_list[1]==1){
		Message("\nAccess Error : Modified Enthalpy | inlet | during inlet velocity\n");
		check = 1;
	}
	if(porosity_and_permeability_check_list[0]==1){
		Message("\nAccess Error : Porosity | cell zone | during cell zone condition\n");
		check = 1;
	}
	if(porosity_and_permeability_check_list[1]==1){
		Message("\nAccess Error : Resistence | cell zone | during cell zone condition\n");
		check = 1;
	}
	
	if(check == 0){
		Message("\nAll access \n");
	}

	for(int i = 0; i<6;i++){
		inlet_enthalpy_check_list[i]=0;
	}
	for(int i = 0; i<5;i++){
		inlet_temp_s_check_list[i]=0;
	}
	for(int i = 0; i<2;i++){
		inlet_velocity_check_list[i]=0;
	}
	for(int i = 0; i<2;i++){
		porosity_and_permeability_check_list[i]=0;
	}
}