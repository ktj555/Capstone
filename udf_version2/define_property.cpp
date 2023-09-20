#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

DEFINE_PROPERTY(Mixture_Rho, c, t) {
	switch(state(c,t)){
    case liquid:
        return Rho_l(c,t);
    case vapor:
        return Rho_v(c,t);
    case mixture:
        return Rho_m(c,t);
    }
}
DEFINE_PROPERTY(Mixture_viscsity, c, t) {
	switch(state(c,t)){
    case liquid:
        return Viscosity_l(c,t);
    case vapor:
        return Viscosity_v(c,t);
    case mixture:
        return Viscosity_m(c,t);
    }
}