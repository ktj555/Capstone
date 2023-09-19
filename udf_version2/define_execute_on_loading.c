#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

DEFINE_EXECUTE_ON_LOADING(set_load, libname)
{
	Set_User_Scalar_Name(enthalpy, "Modified_enthalpy");
	Set_User_Scalar_Name(temp_s, "solid_temperature");
}