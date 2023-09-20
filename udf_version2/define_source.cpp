#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

/*
DEFINE_SOURCE 매크로는 해당 셀에서의 source term을 반환하는 것이 목표
Argument : macro name, cell thread, thread (Group of cells), derivative array of each source, equation index

{ macro name을 제외한 인자들은 fluent에서 자동으로 매칭}
*/

DEFINE_SOURCE(source_for_fluid,c,t,dS,eqn){
    switch(state(c,t)){         // 해당 셀의 유체 상태를 확인하여 각 상태별로 적절한 source term 값을 반환
    case liquid:
        dS[eqn] = dql_dT(c,t) * dT_dH(c,t);     // dS array에 eqn index에는 해당 source term을 해당 scalar로 미분한 값을 반환해 줘야함
        return q_l(c,t);
    case vapor:
        dS[eqn] = dqv_dT(c,t) * dT_dH(c,t);
        return q_v(c,t);
    case mixture:
        dS[eqn] = q_boil(c,t) + S_(c,t) * dqboil_dS(c,t) - q_v(c,t);
        return S_(c,t) * q_boil(c,t) + (1-S_(c,t)) * q_v(c,t);      // S만큼은 boiling, 1-S 만큼은 vapor에 의한 대류
    }
}

DEFINE_SOURCE(source_for_solid,c,t,dS,eqn){     // solid에 대한 source term은 fluid와 부호 반대
    switch(state(c,t)){
    case liquid:
        dS[eqn] = - dql_dT(c,t) * dT_dH(c,t);
        return -q_l(c,t);
    case vapor:
        dS[eqn] = -dqv_dT(c,t) * dT_dH(c,t);
        return -q_v(c,t);
    case mixture:
        dS[eqn] = -q_boil(c,t) - S_(c,t) * dqboil_dS(c,t) + q_v(c,t);
        return -S_(c,t) * q_boil(c,t) - (1-S_(c,t)) * q_v(c,t);
    }
}