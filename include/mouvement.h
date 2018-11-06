#ifndef __MOUVEMENT_H__
#define __MOUVEMENT_H__

#include "nrutil.h"
#include "nrdef.h"


#define VMIN 1
#define VMAX 254

uint8 ** routine_FrameDifference(uint8 **in, uint8 **out,  long nrl, long nrh, long ncl, long nch, int seuil);
uint8** SigmaDelta_step0(uint8** V, uint8 ** M	, uint8** I, long nrl, long nrh, long ncl, long nch);




#endif // __MOUVEMENT_H__

