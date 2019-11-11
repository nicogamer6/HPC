#ifndef __MORPHO_SSE_H__
#define __MORPHO_SSE_H__

#include "nrutil.h"
#include "nrdef.h"


void erosion3SSE(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch);
void dilatation3SSE(uint8 ** Et, uint8 **EtD, long nrl, long nrh, long ncl, long nch);
void ouverture3SSE(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture3SSE(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);

void erosion5SSE(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch);
void dilatation5SSE(uint8 ** Et, uint8 **EtD, long nrl, long nrh, long ncl, long nch);
void ouverture5SSE(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture5SSE(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch);


void erosion3SSE_bin(ulong64 ** Et, ulong64 **EtE, long nrl, long nrh, long ncl, long nch);
void dilatation3SSE_bin(ulong64 ** Et, ulong64 **EtD, long nrl, long nrh, long ncl, long nch);
void ouverture3SSE_bin(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture3SSE_bin(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch);


void erosion3SSE_line_bin(ulong64 ** Et, ulong64 **EtE, int i, long nrl, long nrh, long ncl, long nch);
void dilatation3SSE_line_bin(ulong64 ** Et, ulong64 **EtD, int i, long nrl, long nrh, long ncl, long nch);
void ouverture3SSE_pipe_bin(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture3SSE_pipe_bin(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch);

void erosion3SSE_line_binOMP(ulong64 ** Et, ulong64 **EtE, int i, long nrl, long nrh, long ncl, long nch);
void dilatation3SSE_line_binOMP(ulong64 ** Et, ulong64 **EtD, int i, long nrl, long nrh, long ncl, long nch);
void ouverture3SSE_pipe_binOMP(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch);
void fermeture3SSE_pipe_binOMP(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch);

#endif // __MORPHO_SSE_H__

