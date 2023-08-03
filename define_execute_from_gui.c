DEFINE_EXECUTE_FROM_GUI(mass_and_heat,libname,mode){
    MASS_IN = RP_Get_Real("myudf/mass");
    Message("Set mass flow rate : %g [kg/s]",MASS_IN);
    Q_IN = RP_Get_Real("myudf/heat");
    Message("Set heat flux : %g [W/m^2]",Q_IN);
}