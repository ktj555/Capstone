#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

/*
DEFINE_EXECTURE_ON_LOADING 매크로는 udf library가 fluent로 load될 때에 한 번만 실행됨

해당 매크로는 uds가 fluent GUI에서 uds0, uds1으로 표기되는 대신에 설정된 이름으로 표기되도록 함
*/

DEFINE_EXECUTE_ON_LOADING(set_load, libname)
{
	Set_User_Scalar_Name(enthalpy, "Modified_enthalpy");
	Set_User_Scalar_Name(temp_s, "solid_temperature");
}