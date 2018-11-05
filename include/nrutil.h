/* ---------------- */
/* --- nrutil.h --- */
/* ---------------- */

#ifndef __NRUTIL_H__
#define __NRUTIL_H__

//#pragma message("  include  nrutil.h")
#include "nrdef.h"

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif


#define NR_END 0
#define FREE_ARG char*

extern long nr_end;


/* ------------------------------- */
/* -- PGM and PNM binary format -- */
/* ------------------------------- */

uint8** LoadPGM_ui8matrix(char *filename, long *nrl, long *nrh, long *ncl, long *nch);
void   MLoadPGM_ui8matrix(char *filename, int nrl, int nrh, int ncl, int nch, uint8 **m);
void    SavePGM_ui8matrix(uint8 **m,       long  nrl, long  nrh, long  ncl, long  nch, char *filename);

rgb8 ** LoadPPM_rgb8matrix(char *filename, long *nrl, long *nrh, long *ncl, long *nch);
void    SavePPM_rgb8matrix(rgb8 **m,       long  nrl, long  nrh, long  ncl, long  nch, char *filename);

#ifdef __cplusplus
}
#endif

#endif // __NRUTL_H__
