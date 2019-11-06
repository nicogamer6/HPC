#ifndef __TEST_MORPHO_SSE_H__
#define __TEST_MORPHO_SSE_H__

#include "nrutil.h"
#include "nrdef.h"
#include "vnrdef.h"
#include "vnrutil.h"


void test_EtapemorphoSSE(void);

void test_routineFD_SSEmorpho3xOuv(int seuil);
void test_routineFD_SSEmorpho3xFerm(int seuil);
void test_routineFD_SSEmorpho3xOuvFerm(int seuil);
void test_routineFD_SSEmorpho3xFermOuv(int seuil);

void test_routineSD_SSEmorpho3xOuv();
void test_routineSD_SSEmorpho3xFerm();
void test_routineSD_SSEmorpho3xOuvFerm();
void test_routineSD_SSEmorpho3xFermOuv();

void test_routineSD_SSEmorpho3xOuv_bin();
void test_routineSD_SSEmorpho3xFerm_bin();
void test_routineSD_SSEmorpho3xOuvFerm_bin();
void test_routineSD_SSEmorpho3xFermOuv_bin();

void test_routineSD_SSEmorpho3xOuv_pipebin();
void test_routineSD_SSEmorpho3xFerm_pipebin();
void test_routineSD_SSEmorpho3xOuvFerm_pipebin();
void test_routineSD_SSEmorpho3xFermOuv_pipebin();

void test_unitaire_erosion3SSE();
void test_unitaire_dilatation3SSE();
void test_unitaire_erosion3SSE_bin();
void test_unitaire_dilatation3SSE_bin();

#endif // __TEST_MORPHO_SSE_H__

