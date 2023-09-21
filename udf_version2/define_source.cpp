#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

DEFINE_SOURCE(source_for_fluid,c,t,dS,eqn){
    switch(state(c,t)){
    case liquid:
        dS[eqn] = dql_dT(c,t) * dT_dH(c,t);
        return q_l(c,t);
    case vapor:
        dS[eqn] = dqv_dT(c,t) * dT_dH(c,t);
        return q_v(c,t);
    case mixture:
        dS[eqn] = (q_boil(c,t) + S_(c,t) * dqboil_dS(c,t) - q_v(c,t)) * dS_dH(c,t);
        return S_(c,t) * q_boil(c,t) + (1-S_(c,t)) * q_v(c,t);
    }
}

DEFINE_SOURCE(source_for_solid,c,t,dS,eqn){
    switch(state(c,t)){
    case liquid:
        dS[eqn] = - dql_dTs(c,t);
        return -q_l(c,t);
    case vapor:
        dS[eqn] = -dqv_dTs(c,t);
        return -q_v(c,t);
    case mixture:
        dS[eqn] = -S_(c,t) * dqboil_dTs(c,t) - (1-S) * dqv_dTs(c,t);
        return -S_(c,t) * q_boil(c,t) - (1-S_(c,t)) * q_v(c,t);
    }
}