#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

extern constant models;

DEFINE_ADJUST(adjust_variables, d)
{
#if !RP_HOST
	cell_t c;
	Thread* t;

	real P;
	real h_m;
	real T_f, T_s;

	thread_loop_c(t, d) {
		begin_c_loop(c, t)
		{
			C_UDMI(c,t,saturation) = S_(c,t);
            C_UDMI(c,t,fluid_temperature) = T_f(c,t);
		}
		end_c_loop(c, t)
	}

#endif
}