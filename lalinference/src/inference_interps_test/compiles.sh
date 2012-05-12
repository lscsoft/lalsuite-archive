gcc -g $(pkg-config --libs lal) cheby_interp.c -I./ -I/home/channa/gstopt/include --std=gnu99
#gcc -l $(pkg-config --libs gsl lal) cheby_interp.c -I./ -I/home/channa/gstopt/include --std=gnu99
