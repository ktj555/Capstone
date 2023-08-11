// DEFINE_INIT 매크로는 Fluent에서 사용하는 변수들을 초기화하기 위해서 사용
// Fluent에서 수행하는 기본적인 initialization을 수행한 후 DEFINE_INIT 매크로 호출

// 해당 DEFINE_INIT 매크로는 UDS를 초기화하는 것을 목적으로 함
// UDM의 값 중 Porosity와 Permeabiity는 계산에 사용되어야 하는 값으로 설정
// UDM의 나머지 값들은 계산 과정 중 덮어씌워지는 값이므로 0으로 초기화만 진행

DEFINE_INIT(initialize_saturation_and_temperature, d) // parameter 1 : macro name | 2 : domain pointer
{
#if !RP_HOST	// Host cpu가 아닌 Node cpu 라면
	cell_t c;	// Cell 데이터 구조를 저장할 cell_t 타입 c 선언
	Thread* t;	// Cell Thread 데이터 구조를 저장할 Thread pointer 타입 t 선언

	real q_in, mass_in;
	real h_in, h_out, A;
	real x;
	real cen[ND_ND];

	q_in = RP_Get_Real("myudf/mass");
	mass_in = RP_Get_Real("myudf/heat");
	A = M_PI * pow(DIAMETER, 2) / 4;
	h_in = CP_L(INLET_TEMP_F, 0) * INLET_TEMP_F;
	h_out = h_in + q_in/(mass_in/A);


	thread_loop_c(t, d) {		// domain에 포함된 cell thread를 하나씩 불러와 t에 부여하고 루프를 진행
		begin_c_loop_all(c, t) 	// cell thread에 포함된 cell을 하나씩 불러와 c에 부여하고 루프를 진행
		{
			C_CENTROID(cen,c,t);
			x = cen[0];
			C_UDSI(c, t, MODIFIED_ENTHALPY) = h_in + h_out * x / THICKNESS;	// Enthalpy 초기값을 300K의 물의 Enthalpy로 초기화
			C_UDSI(c, t, TEMP_S) = H_to_T(h_in + h_out * x / THICKNESS, 0);			// 고체 온도를 300K으로 초기화
															// 이러한 온도 profile은 열이 가해지지 않았을 경우 외기온도(300K)으로 평형을 이루는 것을 가정했기 때문
			C_UDMI(c, t, MY_POROSITY) = 0.3;		// 모든 Porosity를 0.3으로 설정 
			C_UDMI(c, t, PERMEABILITY) = 3.67e-12;	// 모든 Permeability를 3.67e-12로 설정

			// 공간 할당을 위해서 사용
			C_UDMI(c, t, SATURATION) = 0;
			C_UDMI(c, t, TEMP_F) = 0;
			C_UDMI(c, t, SOURCE_F) = 0;
			C_UDMI(c, t, SOURCE_S) = 0;
		}
		end_c_loop_all(c, t)
	}

	check_list_message();
#endif
}