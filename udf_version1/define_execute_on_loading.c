// DEFINE_EXECUTE_ON_LOADING 매크로는 Fluent에서 udf가 로드될 때 수행되는 작업을 정의하기 위해 사용
// 해당 DEFINE_EXECUTE_ON_LOADING 매크로는 GUI에서 UDM과 UDS의 이름을 설정

DEFINE_EXECUTE_ON_LOADING(set_load, libname)	// paramter 1 : macro name | 2 : library name loaded by Fluent
{
	// UDS 이름 설정
	Set_User_Scalar_Name(MODIFIED_ENTHALPY, "Modified_enthalpy");
	Set_User_Scalar_Name(TEMP_S, "solid_temperature");

	// UDM 이름 설정
	Set_User_Memory_Name(MY_POROSITY, "local_porosity");
	Set_User_Memory_Name(PERMEABILITY, "local_permeability");
	Set_User_Memory_Name(TEMP_F, "fluid_temperature");
	Set_User_Memory_Name(SATURATION, "Saturation");
	Set_User_Memory_Name(SOURCE_F, "Source_in_fiuld");
	Set_User_Memory_Name(SOURCE_S, "Source_in_solid");
}