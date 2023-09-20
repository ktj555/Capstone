#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

/*
DEFINE_EXECUTE_FROM_GUI 매크로는 그래픽 창에서 사용자가 설정해둔 버튼에 기능을 부여함

해당 매크로는 GUI 창에서 질량유량 및 열유속을 설정하고 OK를 눌렀을 시 fluent가 정상적으로 적용되었는지 확인하기 위함
*/

DEFINE_EXECUTE_FROM_GUI(mass_and_heat,libname,mode){
    real MASS_IN, Q_IN;
    MASS_IN = RP_Get_Real("myudf/mass");
    Q_IN = RP_Get_Real("myudf/heat");
#if RP_HOST
    Message("Set mass flow rate : %g [kg/s]\n",MASS_IN);
    Message("Set heat flux : %g [W/m^2]\n",Q_IN);
#endif
}