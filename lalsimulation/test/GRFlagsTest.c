#include  <lal/LALSimInspiralWaveformGRTestFlags.h>
#include <lal/LALStatusMacros.h>

int main(int argc , char *argv[])
{
    /* Set lalDebugLevel to print all info, warnings, errors */
    lalDebugLevel = 7;

    (void) argc;
    (void) argv;
    int errnum;

    /* Create a new struct */
    printf("Creating a new struct with one parameter...\n");
    LALSimGRTestParam *test=XLALSimCreateGRParam("first_param",10.0); ;
    printf("value:%lf\n",XLALSimGetGRParamValue(test, "first_param"));

    /* Add a second parameter to the struct */
    printf("Adding a second parameter to the struct...\n");
    if( XLALSimAddGRParam(test,"second_param",20.0) != XLAL_SUCCESS)
        XLAL_ERROR(XLAL_EFUNC);
    printf("Printing the struct after adding second parameter...\n");
    if( XLALSimPrintGRParamStruct(stderr,test) != XLAL_SUCCESS)
        XLAL_ERROR(XLAL_EFUNC);
    printf("Get second parameter and print its value...\n");
    printf("value:%lf\n",XLALSimGetGRParamValue(test, "second_param"));

    /* Set first parameter to new value */
    printf("Setting first parameter to new value...\n");
    if( XLALSimSetGRParamValue(test, "first_param",1000.0) != XLAL_SUCCESS )
        XLAL_ERROR(XLAL_EFUNC);
    printf("new value:%lf\n",XLALSimGetGRParamValue(test, "first_param"));

    /* Check for existing and non-existent parameters */
    printf("Checking for existing first parameter and non-existent third parameter...\n");
    printf("first_param:%d second_param:%d third_param:%d\n",
            XLALSimGRParamExists(test, "first_param"),
            XLALSimGRParamExists(test, "second_param"),
            XLALSimGRParamExists(test, "third_param"));
    printf("Now add a third parameter...\n");
    if( XLALSimAddGRParam(test,"third_param",12.0) != XLAL_SUCCESS )
        XLAL_ERROR(XLAL_EFUNC);
    printf("first_param:%d second_param:%d third_param:%d\n",
            XLALSimGRParamExists(test, "first_param"),
            XLALSimGRParamExists(test, "second_param"),
            XLALSimGRParamExists(test, "third_param"));

    /* Print the params struct */
    printf("We print the struct as it appears now...\n");
    if( XLALSimPrintGRParamStruct(stderr,test) != XLAL_SUCCESS )
        XLAL_ERROR(XLAL_EFUNC);

    /* Try to add a parameter that already exists */
    printf("Trying to add a parameter that already exists...\n");
    XLAL_TRY( XLALSimAddGRParam(test,"third_param",12.0), errnum );
    printf("This throws the above message and error code %d\n", errnum);

    /* Destroy the params struct */
    XLALSimDestroyGRParam(test);

    return 0;
}
