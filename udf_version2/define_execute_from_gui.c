#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

DEFINE_EXECUTE_FROM_GUI(mass_and_heat,libname,mode){
    MASS_IN = RP_Get_Real("myudf/mass");
    Q_IN = RP_Get_Real("myudf/heat");
#if RP_HOST
    Message("Set mass flow rate : %g [kg/s]\n",MASS_IN);
    Message("Set heat flux : %g [W/m^2]\n",Q_IN);
#endif
}