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