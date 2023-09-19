#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

DEFINE_SOURCE(source_for_fluid,c,t,dS,eqn){
    switch(state(c,t)){
    case state_of_cell::liquid:
        dS[eqn] = dql_dT(c,t) * dT_dH(c,t);
        return q_l(c,t);

    case state_of_cell::vapor:
        dS[eqn] = dqv_dT(c,t) * dT_dH(c,t);
        return q_v(c,t);
        
    case state_of_cell::mixture:
        dS[eqn] = q_boil(c,t) + S(c,t) * dqboil_dS(c,t) - q_v(c,t);
        return S(c,t) * q_boil(c,t) + (1-S(c,t)) * q_v(c,t);
    }
}

DEFINE_SOURCE(source_for_solid,c,t,ds,eqn){
    switch(state(c,t)){
    case state_of_cell::liquid:
        dS[eqn] = - dql_dT(c,t) * dT_dH(c,t);
        return -q_l(c,t);
    case state_of_cell::vapor:
        dS[eqn] = -dqv_dT(c,t) * dT_dH(c,t);
        return -q_v(c,t);
    case state_of_cell::mixture:
        dS[eqn] = -q_boil(c,t) - S(c,t) * dqboil_dS(c,t) + q_v(c,t);
        return -S(c,t) * q_boil(c,t) - (1-S(c,t)) * q_v(c,t);
    }
}