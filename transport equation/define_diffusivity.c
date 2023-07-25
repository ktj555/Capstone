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