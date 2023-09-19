#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

DEFINE_DIFFUSIVITY(diffusivity_for_fluid, c, t, i){
    real e;

    e = Porosity(c,t);

    switch(state(c,t)){
    case state_of_cell::liquid:
        real k_l = Conductivity_l(c,t);
        real cp_l = Specific_Heat_l(c,t);
        return e * k_l / cp_l;

    case state_of_cell::vapor:
        real k_v = Conductivity_v(c,t);
        real cp_v = Specific_Heat_v(c,t);
        return e * k_v / cp_v;
        
    case state_of_cell::mixture:
        real K, nu, lambda, dj, dlambda;
        K = Permeability(c,t);
        nu = Kinematic_Viscosity(c,t);
        lambda = L(c,t);
        dlambda = dl_dS(c,t);
        dj = -J_function_derivative(c,t);
        D = K / nu * lambda * (1-lmabda) * sqrt(e / K) * model.sigma * dj
        return D / dlambda ; 
    }
}

DEFINE_DIFFUSIVITY(diffusivity_for_solid, c, t, i){
    real e, k_s;

    e = Porosity(c,t);
    k_s = Condcutivity_s(c,t);

    return (1-e)*k_s;
}