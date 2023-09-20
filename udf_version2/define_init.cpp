#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

extern constant models;	// 계산에 필요한 여러 상수들(ex. 시편 두께 및 직경)을 모아둔 구조체
                        // model.h 파일에 정의됨

/*
DEFINE_INIT 매크로는 domain에서 어떻게 scalar 값들을 초기화 할지 결정
Argument list : macro name, domain pointer

{ domain pointer는 fluent에서 자동으로 매칭시킴 }
*/

DEFINE_INIT(init_function,d){
#if !RP_HOST
	cell_t c;
	Thread* t;

	real q_in, mass_in, m_flux;
	real h_in, h_out, A;
	real x;
	real cen[ND_ND];

	// fluent에서 heat flux와 mass flow rate을 불러와 h_in, h_out을 결정
	q_in = RP_Get_Real("myudf/heat");
	mass_in = RP_Get_Real("myudf/mass");
	A = M_PI * pow(models.D, 2) / 4;
	h_in = models.reservoir_enthalpy;
	m_flux = mass_in / A;
    h_out = h_in + q_in/m_flux;

	thread_loop_c(t, d) {
		begin_c_loop_all(c, t)
		{
			C_CENTROID(cen,c,t);	// cen변수에 cell의 중심좌료를 대입
			x = cen[0];
			C_UDSI(c, t, enthalpy) = h_in + (h_out - h_in) / models.thickness * (x + models.thickness * 0.5);	// 입구 h_in에서 출구 h_out까지 선형적인 분포 가정
			C_UDSI(c, t, temp_s) = T_f(c,t) + models.init_dt;		// 유체 온도보다 미리 설정한 온도만큼 높게 설정
		}
		end_c_loop_all(c, t)
	}
	// cell 뿐만 아니라 face에 대해서도 수행
	thread_loop_f(t,d){
		begin_f_loop(f,t){
			F_CENTRIOD(cen,c,t);
			x = cen[0];
			F_UDSI(c,t,enthalpy) = h_in + (h_out - h_in) / models.thickness * (x + models.thickness) * 0.5;
			F_UDSI(c,t,temp_s) = T_f_face(f,t) + models.init_dt;
		}
	}


#endif
}