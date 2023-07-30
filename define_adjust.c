// DEFINE_ADJUST 매크로는 Fluent에서 사용하는 변수들을 계산마다 업데이트 해주기 위해서 사용
// UDS는 자동적으로 업데이트가 이루어지므로 DEFINE_ADJUST를 통해 업데이트 할 필요가 없음
// 다만 UDM는 자동적으로 업데이트 되지 않기 때문에 코드 작성 필요

// UDM에 저장되는 값은 Porosity와 Permeability를 제외하고는 계산에 사용되는 값이 아님
// UDM에 저장되는 값은 원하는 값들을 시각화하기 위해 메모리에 값을 저장하기 위해서 사용

// 따라서 DEFINE_ADJUST 매크로가 수행하는 일은 매 계산이 끝나고 시작하기 전,
// 현재의 UDS(Modified Enthalpy)를 기준으로 Saturation과 Fluid temperature를 저장하는 일을 수행

DEFINE_ADJUST(adjust_variables, d) // parameter 1 : macro name | 2 : domain pointer
{
#if !RP_HOST		// Host cpu가 아닌 Node cpu 라면
	cell_t c;		// Cell 데이터 구조를 저장할 cell_t 타입 c 선언
	Thread* t;		// Cell Thread 데이터 구조를 저장할 Thread pointer 타입 t 선언

	real P;			// 계산에 필요한 Cell의 압력을 저장하기 위한 변수
	real h_m;		// 계산에 필요한 Enthalpy를 저장하기 위한 변수

	thread_loop_c(t, d) {	// domain에 포함된 cell thread를 하나씩 불러와 t에 부여하고 루프를 진행
		begin_c_loop(c, t)	// cell thread에 포함된 cell을 하나씩 불러와 c에 부여하고 루프를 진행
		{
			h_m = C_UDSI(c, t, MODIFIED_ENTHALPY);		// Thread t에 포함된 c가 지정하는 Cell에 저장되어있는 UDS 중 MODIFIED_ENTALPY번째에 저장된 값 호출
			P = C_P(c, t);								// Thread t에 포함된 c가 지정하는 Cell에 저장되어있는 압력값 호출
			C_UDMI(c, t, TEMP_F) = H_to_T(h_m, P);		// Thread t에 포함된 c가 지정하는 Cell에 저장되어 있는 UDM 중 TEMP_F번째에 저장된 값을 설정
			C_UDMI(c, t, SATURATION) = H_to_S(h_m, P);	// Thread t에 포함된 c가 지정하는 Cell에 저장되어 있는 UDM 중 SATURATION번째에 저장된 값을 설정
		}
		end_c_loop(c, t)
	}

#endif
}