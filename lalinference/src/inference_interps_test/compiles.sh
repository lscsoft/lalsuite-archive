gcc -g $(pkg-config --libs lal) LALInferenceInterps_unit_test.c -I./ -I/home/rory.smith/gstlocal/include --std=gnu99
#gcc -l $(pkg-config --libs gsl lal) cheby_interp.c -I./ -I/home/channa/gstopt/include --std=gnu99
