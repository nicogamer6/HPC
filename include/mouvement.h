#ifndef __MOUVEMENT_H__
#define __MOUVEMENT_H__

#include "nrutil.h"
#include "nrdef.h"


//#define VMIN 1
//#define VMAX 254
//#define N 4

void routine_FrameDifference(uint8 **in1, uint8 **in2, uint8 **res, long nrl, long nrh, long ncl, long nch, int seuil);
void routine_SigmaDelta_step0(uint8 **V, uint8 **M, uint8 **I, long nrl, long nrh, long ncl, long nch);
void routine_SigmaDelta_step0OMP(uint8 **V, uint8 **M, uint8 **I, long nrl, long nrh, long ncl, long nch);
void routine_SigmaDelta_1step(uint8 **V, uint8 **Vtm1, uint8 **M, uint8 **Mtm1, uint8 **I, uint8 **Et, long nrl, long nrh, long ncl, long nch);
void routine_SigmaDelta_1step_opti(uint8 **V, uint8 **Vtm1, uint8 **M, uint8 **Mtm1, uint8 **I, uint8 **Et, long nrl, long nrh, long ncl, long nch);
void routine_SigmaDelta_soa(SoA Vm, uint8 **Vtm1, uint8 **M, uint8 **Mtm1, uint8 **Et, long nrl, long nrh, long ncl, long nch);
void routine_SigmaDelta_1step_SOA(SoA3 Vm, uint8 **V, uint8 **M, uint8 **Et, long nrl, long nrh, long ncl, long nch);

void routine_SigmaDelta_1stepOMP(uint8 **V, uint8 **Vtm1, uint8 **M, uint8 **Mtm1, uint8 **I, uint8 **Et, long nrl, long nrh, long ncl, long nch);


int min(int a, int b);
int max(int a, int b);


#endif // __MOUVEMENT_H__

