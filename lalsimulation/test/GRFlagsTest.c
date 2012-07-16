#include  <lal/LALSimInspiralWaveformGRTestFlags.h>

int main(int argc , char *argv[])
{
    (void) argc;
    (void) argv;
    LALSimGRTestParam *test=XLALSimCreateGRParam("first_param",10.0); ;
    printf("value:%lf\n",XLALSimGetGRParamValue(test, "first_param"));
    XLALSimAddGRParam(test,"second_param",20.0);
    XLALSimPrintGRParamStruct(stderr,test);
    printf("value:%lf\n",XLALSimGetGRParamValue(test, "second_param"));
    XLALSimSetGRParamValue(test, "first_param",1000.0);
    printf("new value:%lf\n",XLALSimGetGRParamValue(test, "first_param"));
    printf("first_param:%d third_param:%d\n",XLALSimGRParamExists(test, "first_param"),XLALSimGRParamExists(test, "third_param"));
    XLALSimAddGRParam(test,"third_param",12.0);
    printf("first_param:%d third_param:%d\n",XLALSimGRParamExists(test, "first_param"),XLALSimGRParamExists(test, "third_param"));
    XLALSimAddGRParam(test,"third_param",12.0);
    XLALSimPrintGRParamStruct(stderr,test);
    XLALSimFreeGRParam(test);
}