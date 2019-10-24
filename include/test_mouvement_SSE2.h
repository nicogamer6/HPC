#ifndef __TEST_MOUVEMENT_SSE_H__
#define __TEST_MOUVEMENT_SSE_H__

#include "nrutil.h"
#include "nrdef.h"
#include "vnrdef.h"

//#define NBIMAGES 199
//#define SEUILFD 20

void uint_to_vuint(uint8 ** scalaire, vuint ** vecteur, int vi0, int vi1, int vj0, int vj1);

void vuint_to_uint(uint8 ** scalaire, vuint ** vecteur, int vi0, int vi1, int vj0, int vj1);

void uint_to_vuint64(ulong64 ** scalaire, vuint ** vecteur, int vi0, int vi1, int vj0, int vj1);

void vuint_to_uint64(ulong64 ** scalaire, vuint ** vecteur, int vi0, int vi1, int vj0, int vj1);

void test_routineFD_SSE(int seuil);
void test_routineFD_SSE_OMP(int seuil);
void test_unitaire_FD_SSE();

void test_routineSD_SSE();
void test_routineSD_SSE_OMP();
void test_unitaire_SD_SSE();


#endif

