#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

/*
DEFINE_PROPERTY 매크로는 fluent에서 물성을 불러올 때 사용하는 매크로, 해당 셀에서의 물성을 반환하는 것이 목표
*/

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