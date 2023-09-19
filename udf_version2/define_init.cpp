#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

DEFINE_INIT(init_function,d){
#if !RP_HOST
	cell_t c;
	Thread* t;

	real q_in, mass_in;
	real h_in, h_out, A;
	real x;
	real cen[ND_ND];

	q_in = RP_Get_Real("myudf/heat");
	mass_in = RP_Get_Real("myudf/mass");
	A = M_PI * pow(constant::D, 2) / 4;
    h_in = Specific_Heat_l(c,t) * constant::reservoir_temp;
    h_out = h_in + q_in/(mass_in/A);

	thread_loop_c(t, d) {
		begin_c_loop_all(c, t)
		{
			C_CENTROID(cen,c,t);
			x = cen[0];
			C_UDSI(c, t, uds::enthalpy) = h_in + (h_out - h_in) / constant::thickness * (x + constant::thickness * 0.5);
			C_UDSI(c, t, uds::temp_s) = T_f(c,t) + constant::init_dt;
		}
		end_c_loop_all(c, t)
	}
#endif
}