// DEFINE_PROPERTY 매크로는 Fluent에서 사용하는 물성을 정의하기 위해서 사용
// 해당 DEFINE_PROPERTY 매크로는 mixture의 밀도와 점성을 정의하기 위해서 사용

DEFINE_PROPERTY(Mixture_Rho, c, t) {
	real H, T, P;
	P = C_P(c, t);
	H = C_UDSI(c, t, MODIFIED_ENTHALPY);
	return RHO(H, P);
}
DEFINE_PROPERTY(Mixture_viscsity, c, t) {
	real H, P;
	P = C_P(c, t);
	H = C_UDSI(c, t, MODIFIED_ENTHALPY);
	return MU(H, P);
}