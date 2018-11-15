#ifndef __TEST_MOUVEMENT_SSE_H__
#define __TEST_MOUVEMENT_SSE_H__

#include "nrutil.h"
#include "nrdef.h"
#include "vnrdef.h"

//#define NBIMAGES 299
//#define SEUILFD 20

void uint_to_vuint(uint8 ** scal, vuint ** vect, int vi0, int vi1, int vj0, int vj1);

void vuint_to_uint(uint8 ** scal, vuint ** vect, int vi0, int vi1, int vj0, int vj1);

void test_routineFD_SSE(int seuil);

void test_unitaire_fd_SSE();

void test_routineSD_SSE();
void test_unitaire_SD_SSE();


#endif

