#ifndef __MORPHO_H__
#define __MORPHO_H__

#include "nrutil.h"
#include "nrdef.h"


void erosion3(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch);
void erosion3_opti_lu_rr(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch);
void dilatation3(uint8 ** Et, uint8 **EtD, long nrl, long nrh, long ncl, long nch);
void dilatation3_opti_lu_rr(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch);
void ouverture3(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);
void ouverture3_opti(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture3(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture3_opti(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);

void erosion5(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch);
void dilatation5(uint8 ** Et, uint8 **EtD, long nrl, long nrh, long ncl, long nch);
void ouverture5(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture5(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);



#endif // __MORPHO_H__

