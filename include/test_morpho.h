#ifndef __TEST_MORPHO_H__
#define __TEST_MORPHO_H__

#include "nrutil.h"
#include "nrdef.h"

void convBinToChar(ulong64 ** Et, uint8 ** Etout, long nrl, long nrh, long ncl, long nch);
void convCharToBin(uint8 ** Et, ulong64 ** Etout, long nrl, long nrh, long ncl, long nch);

void test_routineFDmorpho3xOuv(int seuil);
void test_routineFDmorpho3xFerm(int seuil);
void test_routineFDmorpho3xOuvFerm(int seuil);
void test_routineFDmorpho3xFermOuv(int seuil);

void test_routineFDmorpho3xOuv_opti(int seuil);
void test_routineFDmorpho3xFerm_opti(int seuil);
void test_routineFDmorpho3xOuvFerm_opti(int seuil);
void test_routineFDmorpho3xFermOuv_opti(int seuil);

void test_routineFDmorpho3xOuv_pipe(int seuil);
void test_routineFDmorpho3xFerm_pipe(int seuil);

void test_routineFDmorpho3xOuv_bin(int seuil);
void test_routineFDmorpho3xFerm_bin(int seuil);
    
void test_routineSDmorpho3xOuv();
void test_routineSDmorpho3xFerm();
void test_routineSDmorpho3xOuvFerm();
void test_routineSDmorpho3xFermOuv();

void test_routineSDmorpho3xOuv_opti();
void test_routineSDmorpho3xFerm_opti();
void test_routineSDmorpho3xOuvFerm_opti();
void test_routineSDmorpho3xFermOuv_opti();

void test_routineSDmorpho3xOuv_pipe();
void test_routineSDmorpho3xFerm_pipe();

void test_routineSDmorpho3xOuv_bin();
void test_routineSDmorpho3xFerm_bin();

void test_Etapemorpho(void);

    
#endif // __TEST_MORPHO_H__

