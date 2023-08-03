DEFINE_EXECUTE_FROM_GUI(mass_and_heat,libname,mode){
    MASS_IN = RP_Get_Real("myudf/mass");
    Q_IN = RP_Get_Real("myudf/heat");
    Message("%g,%g\n",MASS_IN,Q_IN);
}