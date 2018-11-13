#ifndef __MORPHO_SSE_H__
#define __MORPHO_SSE_H__

#include "nrutil.h"
#include "nrdef.h"


void erosion3SSE(vuint8 ** Et, vuint8 **EtE, long nrl, long nrh, long ncl, long nch);
void dilatation3SSE(vuint8 ** Et, vuint8 **EtD, long nrl, long nrh, long ncl, long nch);
void ouverture3SSE(vuint8 ** Et, vuint8 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture3SSE(vuint8 ** Et, vuint8 **Etout, long nrl, long nrh, long ncl, long nch);

void erosion5SSE(vuint8 ** Et, vuint8 **EtE, long nrl, long nrh, long ncl, long nch);
void dilatation5SSE(vuint8 ** Et, vuint8 **EtD, long nrl, long nrh, long ncl, long nch);
void ouverture5SSE(vuint8 ** Et, vuint8 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture5SSE(vuint8 ** Et, vuint8 **Etout, long nrl, long nrh, long ncl, long nch);



#endif // __MORPHO_SSE_H__

