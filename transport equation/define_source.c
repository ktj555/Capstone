DEFINE_SOURCE(Source_for_solid, c, t, dS, eqn)
{
	real NV_VEC(psi);
	real s, alpha_sf, T_s;
	real source;
	real h_m, h_m_s0, h_m_s1, P, v;

	P = C_P(c, t);
	h_m = C_UDSI(c, t, MODIFIED_ENTHALPY);
	T_s = C_UDSI(c, t, TEMP_S);
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);

	s = H_to_S(h_m, P);

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
	real s, alpha_sf, T_s;
	real source;
	real h_m, h_m_s0, h_m_s1, P, v;

	P = C_P(c, t);
	h_m = C_UDSI(c, t, MODIFIED_ENTHALPY);
	T_s = C_UDSI(c, t, TEMP_S);
	h_m_s0 = H_V_SAT(P);
	h_m_s1 = H_L_SAT(P);

	s = H_to_S(h_m, P);

	alpha_sf = 6 * (1 - C_UDMI(c, t, MY_POROSITY)) / D_P;
	NV_D(psi, =, C_U(c, t), C_V(c, t), C_W(c, t));
	v = NV_MAG(psi);

	source = qsf(h_m, P, v, alpha_sf, T_s);
	dS[eqn] = dqsfdH(h_m, P, v, alpha_sf, T_s);
	C_UDMI(c, t, SOURCE_F) = source;
	return source;
}