#include "udf.h"
#include "property.h"
#include "model.h"
#include "heattransfer.h"
#include "info.h"

extern constant models;

DEFINE_ON_DEMAND(get_data){
    Domain* d;
    Thread* t;
    cell_t c;

    thread_loop_c(t,d){
        begin_c_loop(c,t){
            C_UDMI(c,t,data) = 0; // you can change
        }
        end_c_loop(c,t)
    }
}