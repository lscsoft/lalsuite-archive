#ifndef __UTIL_H__
#define __UTIL_H__

/* locate argument number arg in the character line of length length */
void locate_arg(char *line, int length, int arg, int *arg_start, int *arg_stop);
float hann_response(float delta);
void fill_hann_filterN(float *coeffs, int length, int middle, float mismatch);
void fill_hann_filter7(float *coeffs, float mismatch);
void fill_diff_hann_filter7(float *coeffs, float mismatch);
void tabulate_hann_filter7(void);
void tabulated_fill_hann_filter7(float *coeffs, float mismatch);

#endif
