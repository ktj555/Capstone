#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

DEFINE_PROFILE(inlet_enthalpy, t, i)
{
	face_t f;

	begin_f_loop(f,t){
		switch(state(f,t)){
		case liquid:
			F_PROFILE(f,t,i) = inlet_enthalpy_l(f,t);
			break;
		case vapor:
			F_PROFILE(f,t,i) = inlet_enthalpy_v(f,t);
			break;
		case mixture:
			F_PROFILE(f,t,i) = inlet_enthalpy_m(f,t);
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
			T_s = T_f(f,t);
		}
		T = T_f(f,t);

		switch(state(f,t)){
		case liquid:
			h_c = h_l(f,t);
			break;
		case vapor:
			h_c = h_v(f,t);
			break;
		case mixture:
			h_c = S_(f,t) * h_l(f,t) + (1-S_(f,t)) * h_v(f,t);
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
	m_flux = mass_in / (pow(constant::D,2) / 4);

	begin_f_loop(f, t)
	{
		switch(state(f,t)){
		case liquid:
			F_PROFILE(f,t,i) = m_flux / Rho_l(f,t);
			break;
		case vapor:
			F_PROFILE(f,t,i) = m_flux / Rho_v(f,t);
			break;
		case mixture:
			F_PROFILE(f,t,i) = m_flux / Rho_m(f,t);
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
		C_PROFILE(c,t,i) = Permeability(c,t);
	}
	end_c_loop(c,t)
}