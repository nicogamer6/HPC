#ifndef __MORPHO_H__
#define __MORPHO_H__

#include "nrutil.h"
#include "nrdef.h"


void erosion3(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch);
void erosion3_opti_lu_rr(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch);
void erosion3_column(uint8 ** Et, uint8 **EtE, int j, long nrl, long nrh, long ncl, long nch);
void dilatation3(uint8 ** Et, uint8 **EtD, long nrl, long nrh, long ncl, long nch);
void dilatation3_opti_lu_rr(uint8 ** Et, uint8 **EtD, long nrl, long nrh, long ncl, long nch);
void dilatation3_column(uint8 ** Et, uint8 **EtD, int j, long nrl, long nrh, long ncl, long nch);

void erosion3_bin(ulong64 ** Et, ulong64 **EtE, long nrl, long nrh, long ncl, long nch);
void dilatation3_bin(ulong64 ** Et, ulong64 **EtD, long nrl, long nrh, long ncl, long nch);

void erosion3_line_bin(ulong64 ** Et, ulong64 **EtE, int i, long nrl, long nrh, long ncl, long nch);
void dilatation3_line_bin(ulong64 ** Et, ulong64 **EtD, int i, long nrl, long nrh, long ncl, long nch);

void ouverture3(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);
void ouverture3_opti(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);
void ouverture3_pipe(uint8 ** Et, uint8 ** tmp, uint8 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture3(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture3_opti(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture3_pipe(uint8 ** Et, uint8 ** tmp, uint8 **Etout, long nrl, long nrh, long ncl, long nch);

void ouverture3_bin(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture3_bin(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch);

void ouverture3_pipe_bin(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture3_pipe_bin(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch);


void erosion5(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch);
void erosion5_opti_lu_rr(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch);
void dilatation5(uint8 ** Et, uint8 **EtD, long nrl, long nrh, long ncl, long nch);
void dilatation5_opti_lu_rr(uint8 ** Et, uint8 **EtD, long nrl, long nrh, long ncl, long nch);

void ouverture5(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);
void ouverture5_opti(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture5(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture5_opti(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);


#endif // __MORPHO_H__

