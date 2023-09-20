#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

extern constant models; // 계산에 필요한 여러 상수들(ex. 시편 두께 및 직경)을 모아둔 구조체
                        // model.h 파일에 정의됨

/*
DEFINE_PROFILE 매크로는 해당 face에서 지정된 값을 반환하는 것이 목표
Argument list : macro name, thread (Group of faces), index

{ thread, index는 fluent에서 자동으로 매칭시킴 }
*/

DEFINE_PROFILE(inlet_enthalpy, t, i)	// 입구에서 들어오는 엔탈피
{										// Specific value type BC (특정 값 반환)
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

DEFINE_PROFILE(inlet_temp_s_flux, t, i)		// 입구에서 부여되는 고체 열유속 (고체와 유체 사이에서 일어나는 대류 열전달)
{											// Specific flux type BC (확산계수 * 공간좌표에 대한 미분{= 확산항} 반환)
	face_t f;

	real T_s, T, h_c;

	begin_f_loop(f, t)
	{	
		if(NNULLP(THREAD_STORAGE(t, SV_UDS_I(temp_s)))){	// uds1(solid temperature)가 메모리에 할당되어서 접근이 가능한지 판단
        	T_s = F_UDSI(f,t,temp_s);						// 접근 가능한 경우 해당 고체 온도 반환
		}
		else{
			T_s = T_f(f,t);									// 접근이 되지 않을 경우 해당 유체 온도 반환
		}
		T = T_f(f,t);

		switch(state(f,t)){		// 입구단의 유체 상태에 따라 대류 계수 산출
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

DEFINE_PROFILE(inlet_velocity, t, i)	// 입구 유속
{
	face_t f;
	real mass_in, m_flux;

	mass_in = RP_Get_Real("myudf/mass");
	m_flux = mass_in / (pow(models.D,2) * M_PI / 4);

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

DEFINE_PROFILE(heat_flux, t, i)		// 출구에서 부여되는 열유속
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

DEFINE_PROFILE(porosity, t, i) {	// 시편 전체의 Porosity
	cell_t c;
	begin_c_loop(c, t) {
		C_PROFILE(c, t, i) = Porosity(c,t);
	}
	end_c_loop(c,t)
}

DEFINE_PROFILE(permeability_resistence, t, i) {		// 시편 전체의 Permeability의 역수
	cell_t c;
	begin_c_loop(c, t) {
		C_PROFILE(c,t,i) = 1/ Permeability(c,t);
	}
	end_c_loop(c,t)
}