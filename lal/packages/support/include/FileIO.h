/* <lalVerbatim file="FileIOHV">
$Id$
</lalVerbatim> */

/* <lalLaTeX>
\section{Header \texttt{FileIO.h}}
\label{s:FileIO.h}

Provides standard LAL support IO functions.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
\end{verbatim}
\noindent Only use \texttt{FileIO.h} in test code that links to
the \texttt{lalsupport} library.

\vfill{\footnotesize\input{FileIOHV}}
\newpage\input{FileIOC}
</lalLaTeX> */ 

#ifndef _FILEIO_H
#define _FILEIO_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdarg.h>
#include <lal/LALRCSID.h>

NRCSID( FILEIOH, "$Id$" );

#ifndef LALFopen
#define LALFopen fopen
#endif

#ifndef LALFclose
#define LALFclose fclose
#endif

/* maximum string size to print with LAL Printf routines */
#define LAL_PRINTF_BUFSIZE 4096

FILE *
LALOpenDataFile( const char * );

typedef struct tagLALFILE { int compression; void *fp; } LALFILE;
int XLALFileIsCompressed( const char *path );
LALFILE * XLALFileOpenRead( const char *path );
LALFILE * XLALFileOpenWrite( const char *path, int compression );
LALFILE * XLALFileOpen( const char *path, const char *mode );
int XLALFileClose( LALFILE * file );
size_t XLALFileRead( void *ptr, size_t size, size_t nobj, LALFILE *file );
size_t XLALFileWrite( const void *ptr, size_t size, size_t nobj, LALFILE *file );
int XLALFileGetc( LALFILE *file );
int XLALFilePutc( int c, LALFILE *file );
char * XLALFileGets( char * s, int size, LALFILE *file );
int XLALFilePuts( const char * s, LALFILE *file );
int XLALFileVPrintf( LALFILE *file, const char *fmt, va_list ap );
int XLALFilePrintf( LALFILE *file, const char *fmt, ... );
int XLALFileFlush( LALFILE *file );
int XLALFileSeek( LALFILE *file, long offset, int whence );
long XLALFileTell( LALFILE *file );
void XLALFileRewind( LALFILE *file );
int XLALFileEOF( LALFILE *file );


#ifdef __cplusplus
}
#endif

#endif /* _FILEIO_H */
