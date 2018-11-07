#ifndef __MOUVEMENT_H__
#define __MOUVEMENT_H__

#include "nrutil.h"
#include "nrdef.h"


#define VMIN 1
#define VMAX 254
#define N 2

uint8** routine_FrameDifference(uint8 **in, uint8 **out,  long nrl, long nrh, long ncl, long nch, int seuil);
uint8** routine_SigmaDelta_step0(uint8 **V, uint8 **M, uint8 **I, long nrl, long nrh, long ncl, long nch);
uint8** routine_SigmaDelta_1step(uint8 **V, uint8 **Vtm1, uint8 **M, uint8 **Mtm1, uint8 **I, uint8 **Et, uint8 **Ot, long nrl, long nrh, long ncl, long nch);

int min(int a, int b);
int max(int a, int b);


#endif // __MOUVEMENT_H__

