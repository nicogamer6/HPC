#ifndef __MOUVEMENT_H__
#define __MOUVEMENT_H__

#include "nrutil.h"
#include "nrdef.h"


uint8 ** routine_FrameDifference(uint8 **in, uint8 **out,  long nrl, long nrh, long ncl, long nch, int seuil);


#endif // __MOUVEMENT_H__

