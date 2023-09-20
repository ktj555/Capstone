#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

extern constant models; // 계산에 필요한 여러 상수들(ex. 시편 두께 및 직경)을 모아둔 구조체
                        // model.h 파일에 정의됨

/*
DEFINE_DIFFUSIVITY 매크로는 해당 셀에서 확산 계수를 반환하는 것이 목표
Argument list : macro name, cell thread, thread (Group of cells), index

{ cell thread와 thread, index는 fluent에서 자동으로 매칭시킴 }
*/

DEFINE_DIFFUSIVITY(diffusivity_for_fluid, c, t, i){
    real e;
    real k_l, cp_l, k_v, cp_v;
    real K, nu, D, lambda, dj, dlambda;

    e = Porosity(c,t);

    switch(state(c,t)){                 // state 함수는 입력된 cell이 single phase인지, tow phase인지 판정하여 결과를 반환
    case liquid:
        k_l = Conductivity_l(c,t);
        cp_l = Specific_Heat_l(c,t);
        return e * k_l / cp_l;

    case vapor:
        k_v = Conductivity_v(c,t);
        cp_v = Specific_Heat_v(c,t);
        return e * k_v / cp_v;
        
    case mixture:
        K = Permeability(c,t);
        nu = Kinematic_Viscosity_m(c,t);
        lambda = L(c,t);
        dlambda = dl_dS(c,t);
        dj = -J_function_derivative(c,t);
        D = K / nu * lambda * (1-lambda) * sqrt(e / K) * models.sigma * dj;
        return D / dlambda ; 
    }
}

DEFINE_DIFFUSIVITY(diffusivity_for_solid, c, t, i){
    real e, k_s;

    e = Porosity(c,t);
    k_s = Conductivity_s(c,t);

    return (1-e)*k_s;       // k_s_eff
}