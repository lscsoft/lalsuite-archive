//bayespputils.h

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*
 * ============================================================================
 *
 *                                burnin.c
 *
 * ============================================================================
 */

#define EOF_MARKER 26    /* Decimal code of DOS end-of-file marker */
#define MAX_REC_LEN 2048 /* Maximum size of input buffer */
#define MAX_IN_FILES 128

#define SPIN 1
#define NOSPIN 0

#define SPIN_PAR_NUM 15 //+cycles,posterior,prior
#define NOSPIN_PAR_NUM 9  //+cycles, posterior, prior

double findMaxLogL(FILE * input, double	maxLogL);
void printBurnIn(FILE * input, FILE * output, double maxLogL, int chain, double * bayes, int * numpoints);

typedef struct tagBurnInInput{
	char** files;
	int nfiles;
	int spin;
	float deltaLogL;
	char* output_file_name;
}BurnInInput;

typedef struct tagBurnInOutput{
    double** pos;
    int nSamples;
    double bayesFactorFromHarmonicMean;
}BurnInOutput;

int BurnIn(BurnInInput*,BurnInOutput*);
