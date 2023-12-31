#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

extern constant models;

real min(real x, real y){
	if(x>=y) return y;
	return x;
}

real max(real x, real y){
	if(x>=y) return x;
	return y;
}

DEFINE_PROFILE(inlet_enthalpy, t, i)
{
	face_t f;

	real q_in, mass_in, m_flux;
	mass_in = RP_Get_Real("myudf/mass");
	q_in = RP_Get_Real("myudf/heat");
	m_flux = mass_in / (pow(models.D,2) * M_PI / 4);

	begin_f_loop(f,t){
		switch(state(f,t)){
		case liquid:
			F_PROFILE(f,t,i) = max(min(H_sat_l_face(f,t),inlet_enthalpy_l(f,t)),models.reservoir_enthalpy);
			break;
		case vapor:
			F_PROFILE(f,t,i) = max(min(models.reservoir_enthalpy + q_in / m_flux,inlet_enthalpy_v(f,t)),H_sat_v_face(f,t));
			break;
		case mixture:
			F_PROFILE(f,t,i) = max(min(H_sat_v_face(f,t),inlet_enthalpy_m(f,t)),H_sat_l_face(f,t));
			break;
		}
	}
	end_f_loop(f,t)
}

DEFINE_PROFILE(inlet_temp_s_flux, t, i)
{
	face_t f;

	real T_s, T, h_c;

	begin_f_loop(f, t)
	{	
		if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(temp_s)))){
        	T_s = F_UDSI(f,t,temp_s);
		}
		else{
			T_s = T_f_face(f,t);
		}
		T = T_f_face(f,t);

		switch(state(f,t)){
		case liquid:
			h_c = Nu_l_face(f,t) * Conductivity_l_face(f,t) / models.D;
			break;
		case vapor:
			h_c = Nu_v_face(f,t) * Conductivity_v_face(f,t) / models.D;
			break;
		case mixture:
			h_c = S_face(f,t) * (Nu_l_face(f,t) * Conductivity_l_face(f,t) / models.D) + (1-S_face(f,t)) * (Nu_v_face(f,t) * Conductivity_v_face(f,t) / models.D);
			break;
		}

		F_PROFILE(f, t, i) = -h_c * (T_s - T);
	}
	end_f_loop(f, t)
}

DEFINE_PROFILE(inlet_velocity, t, i)
{
	face_t f;
	real mass_in, m_flux;

	mass_in = RP_Get_Real("myudf/mass");
	m_flux = mass_in / (pow(models.D,2) * M_PI / 4);

	begin_f_loop(f, t)
	{
		switch(state(f,t)){
		case liquid:
			F_PROFILE(f,t,i) = m_flux / Rho_l_face(f,t);
			break;
		case vapor:
			F_PROFILE(f,t,i) = m_flux / Rho_v_face(f,t);
			break;
		case mixture:
			F_PROFILE(f,t,i) = m_flux / Rho_m_face(f,t);
			break;
		}
	}
	end_f_loop(f, t)
}

DEFINE_PROFILE(heat_flux, t, i)
{
	face_t f;

	real q_in;
	q_in = RP_Get_Real("myudf/heat");

	begin_f_loop(f, t)
	{
		F_PROFILE(f, t, i) = q_in;
	}
	end_f_loop(f, t)
}

DEFINE_PROFILE(porosity, t, i) {
	cell_t c;
	begin_c_loop(c, t) {
		C_PROFILE(c, t, i) = Porosity(c,t);
	}
	end_c_loop(c,t)
}

DEFINE_PROFILE(permeability_resistence, t, i) {
	cell_t c;
	begin_c_loop(c, t) {
		C_PROFILE(c,t,i) = 1/ Permeability(c,t);
	}
	end_c_loop(c,t)
}