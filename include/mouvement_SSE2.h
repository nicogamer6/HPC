#ifndef __MOUVEMENT_SSE_H__
#define __MOUVEMENT_SSE_H__

#include "nrutil.h"
#include "nrdef.h"
#include "vnrdef.h"
#include "vnrutil.h"

//#define NB_IMAGE 299
//#define SEUILFD 25


void routine_FrameDifference_SSE2(vuint8** It, vuint8** It_1, vuint8** Et, long nrl, long nrh, long ncl, long nch, int seuil);

void SigmaDelta_step0_SSE2 (vuint8** M, vuint8** V, vuint8** It,long nrl, long nrh, long ncl, long nch);
void SigmaDelta_step0_SSE2_OMP (vuint8** M, vuint8** V,  vuint8** It,long nrl, long nrh, long ncl, long nch);

void SigmaDelta_1step_SSE2 (vuint8** V,vuint8** Vtm1, vuint8** M, vuint8** Mtm1, vuint8** It, vuint8** Et,long nrl, long nrh, long ncl, long nch);
void SigmaDelta_1step_SSE2_OMP (vuint8** V,vuint8** Vtm1, vuint8** M, vuint8** Mtm1, vuint8** It, vuint8** Et,long nrl, long nrh, long ncl, long nch);

void routine_FrameDifference_SSE2_OMP(vuint8** It, vuint8** It_1, vuint8** Et, long nrl, long nrh, long ncl, long nch, int seuil);
void SigmaDelta_SSE2_AoSoA(uint8 **It, uint8 **It_1, uint8 **Et, uint8 **Vt, uint8 **Vt_1, uint8 **Mt, uint8 **Mt_1, uint8 **Ot,long nrl, long nrh, long ncl, long nch);


#endif // __#endif // __MOUVEMENT_SSE_H___H__
