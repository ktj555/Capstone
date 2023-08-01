// DEFINE_PROFILE 매크로는 Fluent에서 필요한 profile들을 정의하기 위해 사용
// 일반적으로 경계조건을 부여하기 위한 용도로 많이 사용됨
// DEFNE_PROFILE 매크로는 해당 데이터에 solver가 접근이 가능한 상황인지 체크하는 것이 중요
// 만약 sovler에서 해당 데이터에 접근하는 것이 허용되지 않는 상황이면 Fluent는 강제 종료됨

// 해당 DEFINE_PROFILE 매크로는 입구 Enthalpy를 정의하기 위해 사용

DEFINE_PROFILE(inlet_enthalpy, t, i)	// paramter 1 : macro name | 2 : thread pointer | 3 : index
{
	Thread* tc;	// Cell Thread 데이터 구조를 저장할 Thread pointer 타입 t 선언
	face_t f;	// Face 데이터 구조를 저장할 Face_t 타입 f 선언
	cell_t c;	// Cell 데이터 구조를 저장할 cell_t 타입 c 선언

	real h_m, h_in, P, v, T_s, T;
	real m_flux, dhdx, dPdx;
	real h_c;
	real e, k_f_eff;
	real h_m_s0, h_m_s1;

	m_flux = MASS_IN / (M_PI * pow(DIAMETER, 2) / 4);	// mass flux = mass flow rate / Area

	begin_f_loop(f, t)	// face thread에 포함된 face를 하나씩 불러와 f에 부여하고 루프를 진행
	{
		c = F_C0(f, t);		// 기울기 연산은 face에 대해서 수행 불가, 오로지 cell에 대해서 가능하기 때문에 f에 접하고 있는 Cell의 정보를 불러와 c에 저장
		tc = THREAD_T0(t);	// c를 이용하여 값을 불러올 경우, c를 포함하고 있는 cell thread가 필요하기 때문에 to에 있는 face들과 접하는 cell들의 정보를 포함하는 cell thread를 불러와 tc에 저장
		
		if(NNULLP(THREAD_STORAGE(t,SV_P))){
			P = F_P(f, t);		// F_P함수를 통해 해당 면에서 압력 호출
		}
		else{
			P = 0;
			inlet_enthalpy_check_list[0] = 1;
		}
		if(NNULLP(THREAD_STORAGE(tc,SV_P_G))){
			dPdx = C_P_G(c,tc)[0];
		}
		else{
			dPdx = 0;
		}
		if(NNULLP(THREAD_STORAGE(t,SV_U))){
			v = F_U(f,t);		// F_U함수를 통해 해당 면에서 x방향 속도 호출	
		}
		else{
			v = m_flux / 998.0;
			inlet_enthalpy_check_list[1] = 1;
		}
		if (NNULLP(THREAD_STORAGE(t, SV_UDS_I(TEMP_S)))) {	// face thread t에서 UDS의 주소를 반한화여 Null인지 아닌지 판단
			T_s = F_UDSI(f, t, TEMP_S);						// Null이 아니라면 메모리에 할당되어 접근 가능한 상태이므로 값을 호출
		}
		else {												// Null 이라면 메모리에 할당되지 않아 접근 불가능한 상태이므로 기본값 할당
			T_s = INLET_TEMP_F;
			inlet_enthalpy_check_list[2] = 1;
		}
		if (NNULLP(THREAD_STORAGE(t, SV_UDM_I))) {			// face thread t에서 UDM의 주소를 반환하여 Null인지 아닌지 판단
			e = F_UDMI(f, t, MY_POROSITY);					// Null이 아니라면 메모리에 할당되어 접근 가능한 상태이므로 값을 호출
		}
		else {												// Null 이라면 메모리에 할당되지 않아 접근 불가능한 상태이므로 기본값 할당
			e = 1.0;
			inlet_enthalpy_check_list[3] = 1;
		}
		if (NNULLP(THREAD_STORAGE(t, SV_UDS_I(MODIFIED_ENTHALPY)))) {
			h_m = F_UDSI(f, t, MODIFIED_ENTHALPY);
		}
		else {
			h_m = h_in;
			inlet_enthalpy_check_list[4] = 1;
		}
		if (NNULLP(T_STORAGE_R_NV(tc, SV_UDSI_G(MODIFIED_ENTHALPY)))) {
			dhdx = C_UDSI_G(c, tc, MODIFIED_ENTHALPY)[0];
		}
		else {
			dhdx = 0.0;
			inlet_enthalpy_check_list[5] = 1;
		}
		h_in = CP_L(INLET_TEMP_F, P) * INLET_TEMP_F;	// reservoir coolant Enthalpy
		// face에서 불러온 변수들을 바탕으로 계산에 필요한 값들을 저장
		T = H_to_T(h_m, P);
		k_f_eff = e * K_F(h_m, P);
		h_c = h_inject(h_m, P, v);
		h_m_s0 = H_V_SAT(P);
		h_m_s1 = H_L_SAT(P);

		// 각 상의 상태에 따라 반환할 값을 계산하여 F_PROFILE로 반환
		// UDS specified value 조건의 경우 해당 UDS가 가져야 할 값을 반환
		// UDS specified flux 조건의 경우 해당 UDS의 diffusive flux 값을 반환
		if (h_m <= h_m_s1) {
			F_PROFILE(f, t, i) = (h_in + h_c / m_flux * T_s - k_f_eff / m_flux / CP_L(T, P) * dhdx) / (1 + h_c / m_flux / CP_L(T, P));
		}
		else if (h_m >= h_m_s0) {
			F_PROFILE(f, t, i) = (h_in + h_c / m_flux * (T_s - T_SAT(P) + H_V_SAT(P) / CP_V(T, P)) - k_f_eff / CP_V(T, P) * dhdx) / (1 + h_c / m_flux / CP_V(T, P));
		}
		else {
			F_PROFILE(f, t, i) = (h_in + h_c / m_flux * (T_s - T_SAT(P)) - k_f_eff / m_flux * dT_SATdP(P) * dPdx) / MY_BETA(h_m,P);
		}
	}
	end_f_loop(f, t)
}

// 해당 DEFINE_PROFILE 매크로는 입구에서 유체와의 대류열전달에 의해 빠져나가는 고체 표면의 열유속을 부여하기 위해 사용

DEFINE_PROFILE(inlet_temp_s_flux, t, i)
{
	face_t f;

	real h_m, P, v, T_s, T;
	real m_flux;
	real h_c;
	real e, k_s_eff;

	m_flux = MASS_IN / (3.14 * pow(DIAMETER, 2) / 4);

	begin_f_loop(f, t)
	{
		if(NNULLP(THREAD_STORAGE(t,SV_P))){
			P = F_P(f, t);
		}
		else{
			P = 0;
			inlet_temp_s_check_list[0]=1;
		}
		if(NNULLP(THREAD_STORAGE(t,SV_U))){
			v = F_U(f,t);
		}
		else{
			v = m_flux / 998.0;
			inlet_temp_s_check_list[1]=1;
		}
		if (NNULLP(THREAD_STORAGE(t, SV_UDS_I(TEMP_S)))) {
			T_s = F_UDSI(f, t, TEMP_S);
		}
		else {
			T_s = INLET_TEMP_F;
			inlet_temp_s_check_list[2]=1;
		}
		if (NNULLP(THREAD_STORAGE(t, SV_UDM_I))) {
			e = F_UDMI(f, t, MY_POROSITY);
		}
		else {
			e = 1.0;
			inlet_temp_s_check_list[3]=1;
		}
		if (NNULLP(THREAD_STORAGE(t, SV_UDS_I(MODIFIED_ENTHALPY)))) {
			h_m = F_UDSI(f, t, MODIFIED_ENTHALPY);
		}
		else {
			h_m = CP_L(INLET_TEMP_F, P) * INLET_TEMP_F;
			inlet_temp_s_check_list[4]=1;
		}
		T = H_to_T(h_m, P);
		h_c = h_inject(h_m, P, v);
		k_s_eff = (1 - e) * K_S(T_s);

		F_PROFILE(f, t, i) = -h_c * (T_s - T);
	}
	end_f_loop(f, t)
}

// 해당 DEFINE_PROFILE 매크로는 입구 속력을 정의하기 위해 사용

DEFINE_PROFILE(inlet_velocity, t, i)
{
	face_t f;

	real T, P, h_m, rho;

	begin_f_loop(f, t)
	{
		if(NNULLP(THREAD_STORAGE(t,SV_P))){
			P = F_P(f,t);
		}
		else{
			P = 0;
			inlet_velocity_check_list[0] = 1;
		}
		if (NNULLP(THREAD_STORAGE(t, SV_UDS_I(MODIFIED_ENTHALPY)))) {
			h_m = F_UDSI(f,t,MODIFIED_ENTHALPY);
			T = H_to_T(h_m,P);
		}
		else {
			h_m = INLET_TEMP_F * CP_L(INLET_TEMP_F,P);
			T = INLET_TEMP_F;
			inlet_velocity_check_list[1] = 1;
		}

		rho = RHO(h_m,P);

		F_PROFILE(f, t, i) = MASS_IN / (M_PI * pow(DIAMETER, 2) / 4) / rho;
	}
	end_f_loop(f, t)
}

// 해당 DEFINE_PROFILE 매크로는 출구단에서 고체에 조사되는 열유속을 부여하기 위해 사용

DEFINE_PROFILE(heat_flux, t, i)
{
	face_t f;

	begin_f_loop(f, t)
	{
		F_PROFILE(f, t, i) = Q_IN;
	}
	end_f_loop(f, t)
}

// 해당 DEFINE_PROFILE 매크로는 porosity를 정의하기 위해 사용

DEFINE_PROFILE(porosity, t, i) {
	cell_t c;
	begin_c_loop(c, t) {
		if (NNULLP(THREAD_STORAGE(t, SV_UDM_I))) {
			C_PROFILE(c, t, i) = C_UDMI(c, t, MY_POROSITY);
		}
		else {
			C_PROFILE(c, t, i) = 1.0;
			porosity_and_permeability_check_list[0] = 1;
		}
	}
	end_c_loop(c,t)
}

// 해당 DEFINE_PROFILE 매크로는 Permeability를 정의하기 위해 사용

DEFINE_PROFILE(permeability_resistence, t, i) {
	cell_t c;
	begin_c_loop(c, t) {
		if (NNULLP(THREAD_STORAGE(t, SV_UDM_I))) {
			C_PROFILE(c, t, i) = 1 / C_UDMI(c, t, PERMEABILITY);
		}
		else {
			C_PROFILE(c, t, i) = 2.11e10;
			porosity_and_permeability_check_list[1] = 1;
		}
	}
	end_c_loop(c,t)
}