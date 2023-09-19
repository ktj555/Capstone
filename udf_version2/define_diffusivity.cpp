#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

DEFINE_DIFFUSIVITY(diffusivity_for_fluid, c, t, i){
    real e;

    real k_l, cp_l, k_v, cp_v;
    real K, nu, D, lambda, dj, dlambda;

    e = Porosity(c,t);

    switch(state(c,t)){
    case state_of_cell::liquid:
        k_l = Conductivity_l(c,t);
        cp_l = Specific_Heat_l(c,t);
        return e * k_l / cp_l;

    case state_of_cell::vapor:
        k_v = Conductivity_v(c,t);
        cp_v = Specific_Heat_v(c,t);
        return e * k_v / cp_v;
        
    case state_of_cell::mixture:
        K = Permeability(c,t);
        nu = Kinematic_Viscosity_m(c,t);
        lambda = L(c,t);
        dlambda = dl_dS(c,t);
        dj = -J_function_derivative(c,t);
        D = K / nu * lambda * (1-lambda) * sqrt(e / K) * model.sigma * dj;
        return D / dlambda ; 
    }
}

DEFINE_DIFFUSIVITY(diffusivity_for_solid, c, t, i){
    real e, k_s;

    e = Porosity(c,t);
    k_s = Condcutivity_s(c,t);

    return (1-e)*k_s;
}