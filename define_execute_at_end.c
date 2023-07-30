// DEFINE_EXECUTE_AT_END 매크로는 Fluent에서 매 계산이 끝나고 난 후 할 작업들을 정의하기 위해 사용
// DEFINE_ADJUST나 DEFINE_INIT과는 다르게 인자로 domain pointer를 Ansys로부터 넘겨받지 못하므로
// Domain에 접근할 필요성이 있을 경우 직접 사용자가 이를 명시해줘야 함

// 해당 DEFINE_EXECUTE_AT_END 매크로는 매 계산 직후 고체를 기준으로 출입하는 열량을 계산하여 출력

DEFINE_EXECUTE_AT_END(calculate_energy_balance_vr_transient)	// paramter 1 : macro name
{
	// 출력하기 위해 계산해야하는 값들을 정의하고 0으로 초기화
	real energy_in_solid=0.0, energy_out_solid=0.0, energy_get_solid=0.0, energy_to_fluid=0.0, error=0.0;

	int BC_in = 5;	// inlet boundary face ID
	int BC_out = 6;	// outlet boundary face ID

#if !RP_HOST	// Host cpu가 아니라면 계산 수행
	Domain* d;	// Domain 타입의 데이터 구조를 저장하기 위한 domain pointer 타입 d 선언
	Thread* t;	// Cell Thread 데이터 구조를 저장할 Thread pointer 타입 t 선언
	Thread* ti;	// inlet boundary face의 Face Thread 데이터 구조를 저장할 Thread pointer 타입 t 선언
	Thread* to;	// outlet boundary face의 Face Thread 데이터 구조를 저장할 Thread pointer 타입 t 선언
	Thread* tc;	// face와 맞닿은 cell Thread 데이터 구조를 저장할 Thread pointer 타입 t 선언

	cell_t c;	// Cell 데이터 구조를 저장할 cell_t 타입 c 선언
	face_t f;	// Face 데이터 구조를 저장할 Face_t 타입 f 선언

	real P;		// 계산에 필요한 압력값을 저장하기 위한 변수
	real h_m;	// 계산에 필요한 Enthapy값을 저장하기 위한 변수

	real dt;	// 현재 해석에서 사용하고 있는 time step을 저장하기 위한 변수
	real sum_of_solid_U = 0.0, sum_of_solid_U_old = 0.0;	// 현재와 이전 time에서 계산된 고체의 내부에너지를 저장하기 위한 변수들

	real dtdx;	// 경계에서 온도 그레디언트를 저장하기 위한 변수
	real area;	// 경계면의 면적을 저장하기 위한 변수
	real A[ND_ND];	// 경계면의 면적 벡터를 저장하기 위한 배열, ND_ND는 2D 해석의 경우 2이고 3D 해석의 경우 3
#endif

#if !RP_HOST
	dt = RP_Get_Real("physical-time-step");		// RP_Get_Real은 지정된 변수를 real 타입으로 반환해주는 함수
												// 이를 통해 현재 해석에서 사용하는 time step을 불러와 dt에 저장
	
	d = Get_Domain(1);				// Get_domain함수는 지정된 인덱스의 domain pointer를 반환, domain의 인덱스는 mixture = 1을 가장 상위 계층으로 함
									// Multiphase flow에서 각 상의 domain은 2,3, ... 등으로 fluent gui에서 id 확인 가능
	ti = Lookup_Thread(d, BC_in);	// Lookup_Thread함수는 지정된 domain에서 해당 index에 face thread를 반환
	to = Lookup_Thread(d, BC_out);	// 각각 inlet의 face thread와 outlet의 face thread를 반환하여 ti, to 에 저장
#endif

#if !RP_HOST
	thread_loop_c(t, d) {	// domain에 포함된 cell thread를 하나씩 불러와 t에 부여하고 루프를 진행
		begin_c_loop(c, t)	// cell thread에 포함된 cell을 하나씩 불러와 c에 부여하고 루프를 진행
		{
			// C_UDSI 함수를 통해 지정된 UDS의 값을 불러오고 C_UDSI_M1 함수를 통해 지정된 UDS의 이전 time에서의 값을 불러옴
			// C_VOLUME 함수를 통해 지정된 Cell의 Volume을 불러오고 이를 사용하여 현재 시간 고체의 내부에너지와 이전 시간 고체의 내부에너지를 게산
			sum_of_solid_U += C_UDSI(c, t, TEMP_S) * C_VOLUME(c, t) * CP_S(C_UDSI(c, t, TEMP_S)) * (1 - C_UDMI(c,t,MY_POROSITY));
			sum_of_solid_U_old += C_UDSI_M1(c, t, TEMP_S) * C_VOLUME(c, t) * CP_S(C_UDSI_M1(c, t, TEMP_S)) * (1 - C_UDMI(c, t, MY_POROSITY));

			// fluid로 흘러들어가는 source term을 합산
			energy_to_fluid += C_UDMI(c, t, SOURCE_F) * C_VOLUME(c, t);
		}
		end_c_loop(c,t)
	}

	energy_get_solid = sum_of_solid_U - sum_of_solid_U_old;	// 현재와 이전 내부 에너지 차이가 고체가 얻은 에너지 [J]
	energy_get_solid /= dt;	// 이를 time step으로 나눔 [W]
#endif

#if !RP_HOST
	begin_f_loop(f, to)	// face thread에 포함된 face를 하나씩 불러와 f에 부여하고 루프를 진행
	{
		if(PRINCIPAL_FACE_P(f,to)){					// 해당 face가 principal face인지 판별
			c = F_C0(f, to);						// 기울기 연산은 face에 대해서 수행 불가, 오로지 cell에 대해서 가능하기 때문에 f에 접하고 있는 Cell의 정보를 불러와 c에 저장
			tc = THREAD_T0(to);						// c를 이용하여 값을 불러올 경우, c를 포함하고 있는 cell thread가 필요하기 때문에 to에 있는 face들과 접하는 cell들의 정보를 포함하는 cell thread를 불러와 tc에 저장
			F_AREA(A, f, to);						// f의 면적 벡터를 A에 저장
			area = NV_MAG(A);						// 면적 벡터 A의 크기는 해당 f의 넓이
			dtdx = C_UDSI_G(c, tc, TEMP_S)[0];		// 경계면에서의 온도 그레디언트를 구하기 위해 C_UDSI_G함수로 고체 온도 그레디언트 벡터를 불러오고 그 중 x축 방향의 기울기만 택함
			energy_in_solid += area * K_S(F_UDSI(f, ti, TEMP_S)) * (1 - C_UDMI(c, t, MY_POROSITY)) * dtdx;
		}	
	}
	end_f_loop(f, to)
#endif

#if !RP_HOST
	begin_f_loop(f, ti) 
	{
		if(PRINCIPAL_FACE_P(f,ti)){
			c = F_C0(f, ti);
			tc = THREAD_T0(ti);
			F_AREA(A, f, ti);
			area = NV_MAG(A);
			dtdx = C_UDSI_G(c, tc, TEMP_S)[0];
			energy_out_solid += area * K_S(F_UDSI(f, to, TEMP_S)) * dtdx * (1 - C_UDMI(c, t, MY_POROSITY));
		}

	}
	end_f_loop(f, ti)
#endif

#if !RP_HOST
	// 고체에 출입하는 에너지의 총량과 고체에 내부에너지 변화량을 다 더하면 0이 되어야 함
	error = energy_in_solid + energy_out_solid + energy_get_solid - energy_to_fluid;

	Message("\n\nCompute each node value\n");
	Message("\n- - - - - - - - - - - - - - - - - - - - - -\n");
	Message("outlet    : %.8g [W]\n",energy_in_solid);
	Message("inlet    : %.8g [W]\n", energy_out_solid);
	Message("exchange to fluid : %.8g [W]\n", energy_to_fluid);
	Message("get               : %.8g [W]\n", energy_get_solid);
	Message("Error                   : %.8g [W]\n",error);
	Message("\n- - - - - - - - - - - - - - - - - - - - - -\n\n\n");

	energy_in_solid = PRF_GRSUM1(energy_in_solid);		// 각각의 node cpu에서 계산된 값들을 합산하여 최종 값을 반환
	energy_out_solid = PRF_GRSUM1(energy_out_solid);
	energy_get_solid = PRF_GRSUM1(energy_get_solid);
	energy_to_fluid = PRF_GRSUM1(energy_to_fluid);
	error = PRF_GRSUM1(error);

	Message("\n\nAfter Summing\n");
	Message("\n- - - - - - - - - - - - - - - - - - - - - -\n");
	Message("outlet    : %.8g [W]\n",energy_in_solid);
	Message("inlet    : %.8g [W]\n", energy_out_solid);
	Message("exchange to fluid : %.8g [W]\n", energy_to_fluid);
	Message("get               : %.8g [W]\n", energy_get_solid);
	Message("Error                   : %.8g [W]\n",error);
	Message("\n- - - - - - - - - - - - - - - - - - - - - -\n\n\n");

#endif
	// node cpu에서 계산된 값을 host cpu로 반환
	node_to_host_real_5(energy_in_solid, energy_out_solid, energy_get_solid, energy_to_fluid, error);

#if RP_HOST	// host cpu라면 값을 출력
	Message("\n-------------------------------------------\n");
	Message("Energy balance at solid\n\n");
	Message("Total in from outlet    : %.8g [W]\n",energy_in_solid);
	Message("Total out from inlet    : %.8g [W]\n", energy_out_solid);
	Message("Total exchange to fluid : %.8g [W]\n", energy_to_fluid);
	Message("Total get               : %.8g [W]\n", energy_get_solid);
	Message("Error                   : %.8g [W]\n",error);
	Message("\n-------------------------------------------\n");
#endif

}

DEFINE_EXECUTE_AT_END(calculate_energy_balance_vr_steady) {

	real energy_in_solid=0.0, energy_out_solid=0.0, energy_to_fluid=0.0, error=0.0;

	int BC_in = 5;
	int BC_out = 6;

#if !RP_HOST
	Domain* d;
	Thread* t;
	Thread* ti;
	Thread* to;
	Thread* tc;

	cell_t c;
	face_t f;

	real P;
	real h_m;

	real dtdx;
	real area;
	real A[ND_ND];
#endif

#if !RP_HOST
	d = Get_Domain(1);
	ti = Lookup_Thread(d, BC_in);
	to = Lookup_Thread(d, BC_out);
#endif

#if !RP_HOST
	thread_loop_c(t, d) {
		begin_c_loop(c, t)
		{
			energy_to_fluid += C_UDMI(c, t, SOURCE_F) * C_VOLUME(c, t);
		}
		end_c_loop(c, t)
	}
#endif

#if !RP_HOST
	begin_f_loop(f, to)
	{
		if(PRINCIPAL_FACE_P(f,to)){
			c = F_C0(f, to);
			tc = THREAD_T0(to);
			F_AREA(A, f, to);
			area = NV_MAG(A);
			dtdx = C_UDSI_G(c, tc, TEMP_S)[0];
			energy_in_solid += area * K_S(F_UDSI(f, ti, TEMP_S)) * dtdx * (1 - C_UDMI(c, t, MY_POROSITY));
		}

	}
	end_f_loop(f, to)
#endif

#if !RP_HOST
	begin_f_loop(f, ti)
	{
		if(PRINCIPAL_FACE_P(f,ti)){
			c = F_C0(f, ti);
			tc = THREAD_T0(ti);
			F_AREA(A, f, ti);
			area = NV_MAG(A);
			dtdx = C_UDSI_G(c, tc, TEMP_S)[0];
			energy_out_solid += area * K_S(F_UDSI(f, to, TEMP_S)) * dtdx * (1 - C_UDMI(c, t, MY_POROSITY));
		}

	}
	end_f_loop(f, ti)
#endif

#if !RP_HOST
	error = energy_in_solid + energy_out_solid - energy_to_fluid;

	Message("\n\nCompute each node value\n");
	Message("\n- - - - - - - - - - - - - - - - - - - - - -\n");
	Message("outlet    : %.8g [W]\n",energy_in_solid);
	Message("inlet    : %.8g [W]\n", energy_out_solid);
	Message("exchange to fluid : %.8g [W]\n", energy_to_fluid);
	Message("Error                   : %.8g [W]\n",error);
	Message("\n- - - - - - - - - - - - - - - - - - - - - -\n\n\n");

	energy_in_solid = PRF_GRSUM1(energy_in_solid);
	energy_out_solid = PRF_GRSUM1(energy_out_solid);
	energy_to_fluid = PRF_GRSUM1(energy_to_fluid);
	error = PRF_GRSUM1(error);

	Message("\n\nAfter Summing\n");
	Message("\n- - - - - - - - - - - - - - - - - - - - - -\n");
	Message("outlet    : %.8g [W]\n",energy_in_solid);
	Message("inlet    : %.8g [W]\n", energy_out_solid);
	Message("exchange to fluid : %.8g [W]\n", energy_to_fluid);
	Message("Error                   : %.8g [W]\n",error);
	Message("\n- - - - - - - - - - - - - - - - - - - - - -\n\n\n");

#endif

	node_to_host_real_4(energy_in_solid, energy_out_solid, energy_to_fluid, error);

#if RP_HOST
	Message("\n-------------------------------------------\n");
	Message("Energy balance at solid\n\n");
	Message("Total in from outlet    : %.8g [W]\n", energy_in_solid);
	Message("Total out from inlet    : %.8g [W]\n", energy_out_solid);
	Message("Total exchange to fluid : %.8g [W]\n", energy_to_fluid);
	Message("Error                   : %.8g [W]\n", error);
	Message("\n-------------------------------------------\n");
#endif

}