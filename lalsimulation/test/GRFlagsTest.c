#include  <lal/LALSimInspiralWaveformGRTestFlags.h>
#include <lal/LALStatusMacros.h>

int main(int argc , char *argv[])
{
    /* Set DebugLevel to print all info, warnings, errors */
    lalDebugLevel = 7;

    (void) argc;
    (void) argv;
    LALSimGRTestParam *test=XLALSimCreateGRParam("first_param",10.0); ;
    printf("value:%lf\n",XLALSimGetGRParamValue(test, "first_param"));
    XLALSimAddGRParam(test,"second_param",20.0);
    XLALSimPrintGRParamStruct(stderr,test);
    printf("value:%lf\n",XLALSimGetGRParamValue(test, "second_param"));
    if( !(XLALSimSetGRParamValue(test, "first_param",1000.0) == XLAL_SUCCESS) )
        XLAL_ERROR(XLAL_EFUNC);
    printf("new value:%lf\n",XLALSimGetGRParamValue(test, "first_param"));
    printf("first_param:%d third_param:%d\n",XLALSimGRParamExists(test, "first_param"),XLALSimGRParamExists(test, "third_param"));
    XLALSimAddGRParam(test,"third_param",12.0);
    printf("first_param:%d third_param:%d\n",XLALSimGRParamExists(test, "first_param"),XLALSimGRParamExists(test, "third_param"));
    XLALSimAddGRParam(test,"third_param",12.0);
    XLALSimPrintGRParamStruct(stderr,test);
    XLALSimDestroyGRParam(test);
}
