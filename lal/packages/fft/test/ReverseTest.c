#include <lal/LALStdlib.h>
#include <lal/LALNoiseModels.h>
#include <lal/Random.h>

NRCSID( REVERSETESTC, "$Id$" );

#define REVERSETESTC_ESUB   0
#define REVERSETESTC_MSGESUB   "Subroutine failed"

#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
     "        %s %s\n", (code),  __FILE__,                           \
     __LINE__, REVERSETESTC,                                         \
     statement ? statement : "", (msg) );                            \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( REVERSETESTC_ESUB,                                          \
   REVERSETESTC_MSGESUB,                                             \
   "Function call \"" #func "\" failed:" );                          \
  return REVERSETESTC_ESUB;                                          \
}                                                                    \
while (0)

#define N (4096) 

int main()
  {
  static LALStatus stat;
  FILE *fp = NULL;

  REAL4Vector     *time = NULL;

  COMPLEX8Vector  *freq  = NULL;
  COMPLEX8Vector  *freq2 = NULL;

  REAL4Vector  *freq_re = NULL;
  REAL4Vector  *freq_im = NULL;

  RandomParams *random = NULL;

  REAL4FFTPlan    *fwd = NULL;
  REAL4FFTPlan    *rev = NULL;

  INT4  tLen = N;
  INT4  fLen = N/2+1;

  INT4  j, k;


  /*
   * Allocate Memory
   */
  SUB( LALCCreateVector( &stat, &freq, fLen ), &stat );
  SUB( LALCCreateVector( &stat, &freq2, fLen ), &stat );
  SUB( LALSCreateVector( &stat, &freq_re, fLen ), &stat );
  SUB( LALSCreateVector( &stat, &freq_im, fLen ), &stat );
  SUB( LALSCreateVector( &stat, &time, tLen ), &stat );


  /*
   * Create Original Frequency Series
   */
  SUB( LALCreateRandomParams( &stat, &random, 0 ), &stat );
  SUB( LALNormalDeviates( &stat, freq_re, random ), &stat );  
  SUB( LALDestroyRandomParams( &stat, &random ), &stat );
  SUB( LALCreateRandomParams( &stat, &random, 1 ), &stat );
  SUB( LALNormalDeviates( &stat, freq_im, random ), &stat );  
  SUB( LALDestroyRandomParams( &stat, &random ), &stat );

  fp = fopen( "FDdata1.dat", "w" );
  for( k=0; k<fLen; k++ )
  {
    if( k == 0 || k == (fLen-1) )
    {    
      freq->data[k].re = 0.0;
      freq->data[k].im = 0.0;
    }
    else
    {
      freq->data[k].re = freq_re->data[k];
      freq->data[k].im = freq_im->data[k];
    }
    fprintf( fp, "%e %e\n",
             freq->data[k].re,
             freq->data[k].im );
  }
  fclose( fp );


  /*
   * FFT to time domain
   */
  SUB( LALCreateReverseRealFFTPlan( &stat, &rev, tLen, 0), &stat );
  SUB( LALReverseRealFFT( &stat, time, freq, rev ), &stat );

  fp = fopen( "TDdata.dat", "w" );
  for( j=0; j<tLen; ++j )
  {
    /* Normalise */
    time->data[j] /= tLen;
    fprintf( fp, "%e\n", time->data[j] );
  }

  SUB( LALDestroyRealFFTPlan( &stat, &rev ), &stat );


  /*
   * FFT back to frequency domain
   */
  SUB( LALCreateForwardRealFFTPlan( &stat, &fwd, tLen, 0), &stat );
  SUB( LALForwardRealFFT( &stat, freq2, time, fwd ), &stat );

  fp = fopen( "FDdata2.dat", "w" );
  for( k=0; k<fLen; k++ )
  {
    if( k == 0 || k == (fLen-1) )
    {    
      freq2->data[k].re = 0.0;
      freq2->data[k].im = 0.0;
    }
    fprintf( fp, "%e %e\n",
             freq2->data[k].re,
             freq2->data[k].im );
  }
  fclose( fp );

  SUB( LALDestroyRealFFTPlan( &stat, &fwd ), &stat );
  

  /*
   * CLEAN UP MEMORY
   */
  SUB( LALSDestroyVector( &stat, &time ), &stat );
  SUB( LALCDestroyVector( &stat, &freq ), &stat );
  SUB( LALCDestroyVector( &stat, &freq2 ), &stat );
  SUB( LALSDestroyVector( &stat, &freq_re ), &stat );
  SUB( LALSDestroyVector( &stat, &freq_im ), &stat );

  LALCheckMemoryLeaks(); 

  return 0;
}
