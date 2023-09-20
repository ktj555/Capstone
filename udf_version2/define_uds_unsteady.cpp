#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

/*
DEFINE_SOURCE 매크로는 해당 셀에서의 unsteady term을 반환하는 것이 목표
Argument : macro name, cell thread, thread (Group of cells), index, implicit term, explicit term

{ macro name을 제외한 인자들은 fluent에서 자동으로 매칭}
*/

DEFINE_UDS_UNSTEADY(unsteady_for_fluid, c, t, i, apu, su){
    real dt, vol, rho, rho_old, phi_old, e, beta, beta_old;

	dt = RP_Get_Real("physical-time-step");
	vol = C_VOLUME(c, t);
	e = Porosity(c,t);

	beta = beta_current(c,t);
	beta_old = beta_past(c,t);

	rho = C_R(c, t);
	rho_old = C_R_M1(c, t);

	phi_old = C_STORAGE_R(c, t, SV_UDSI_M1(enthalpy));

	*apu = -e * rho * vol * beta / dt;
	*su = e * rho_old * vol * beta_old * phi_old / dt;
}

DEFINE_UDS_UNSTEADY(unsteady_for_solid, c, t, i, apu, su)
{
	real dt, vol, rho, rho_old, phi_old, e, c_p, c_p_old;

	dt = RP_Get_Real("physical-time-step");
	vol = C_VOLUME(c, t);
	e = Porosity(c,t);

	c_p = Specific_Heat_s(c,t);
	c_p_old = Specific_Heat_s_past(c,t);
	rho = Rho_s(c,t);
	rho_old = Rho_s_past(c,t);

	phi_old = C_STORAGE_R(c, t, SV_UDSI_M1(temp_s));

	*apu = -rho * vol * (1 - e) * c_p / dt;
	*su = rho_old * vol * (1 - e) * c_p_old * phi_old / dt;
}