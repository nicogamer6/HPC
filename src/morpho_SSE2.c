#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "vnrdef.h"
#include "nrdef.h"
#include "morpho_SSE2.h"
#include "vnrutil.h"

#define BORD 2


void erosion3SSE(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch){
	int i, j;
    int m;
    vuint8 l_1, l0, l1;
    vuint8 tmp_1,tmp0,tmp1;
    vuint8 res;
    vuint8 xl, xr;
    vuint8 y;
    vuint8 ero;

    // Parcours de l'image
    for(i = nrl; i <= nrh; i++)
    {
        for(j = ncl; j < nch; j+=14)
        {

                l_1 = _mm_loadu_si128((__m128i*)&Et[i-1][j-1]);
                l0  = _mm_loadu_si128((__m128i*)&Et[i][j-1]);
                l1  = _mm_loadu_si128((__m128i*)&Et[i+1][j-1]);


                res=_mm_and_si128(_mm_and_si128(l_1,l0),l1); //On met dans res l1 & l0 & l_1
                
                tmp0=_mm_srli_si128(res,1);
                tmp1=_mm_slli_si128(res,1);
                res=_mm_and_si128(_mm_and_si128(tmp0,tmp1),res);

                res = _mm_srli_si128(res,1);
                /*
                

                l_1 = _mm_loadu_si128((__m128i*)&Et[i-1][j-1]);
				l0  = _mm_loadu_si128((__m128i*)&Et[i][j-1]);
				l1  = _mm_loadu_si128((__m128i*)&Et[i+1][j-1]);
				res=_mm_and_si128(_mm_and_si128(l_1,l0),l1); //On met dans res l1 & l0 & l_1

				tmp_1=_mm_srli_si128(l_1,1);
				tmp0=_mm_srli_si128(l0,1);
				tmp1=_mm_srli_si128(l1,1);
				res=_mm_and_si128(_mm_and_si128(_mm_and_si128(tmp_1,tmp0),tmp1),res);
				tmp_1=_mm_slli_si128(l_1,1);
				tmp0=_mm_slli_si128(l0,1);
				tmp1=_mm_slli_si128(l1,1);
				res=_mm_and_si128(_mm_and_si128(_mm_and_si128(tmp_1,tmp0),tmp1),res);
				res=_mm_srli_si128(res,1);
*/

            _mm_storeu_si128((__m128i*)&EtE[i][j], res); //Store dans Et la valeur de l'erosion
        }
        //j=nch;
        //EtE[i][j]=(255&Et[i-1][j-1]&Et[i-1][j]&Et[i-1][j+1]&Et[i][j-1]&Et[i][j]&Et[i][j+1]&Et[i+1][j-1]&Et[i+1][j]&Et[i+1][j+1]);
    }

}

void dilatation3SSE(uint8 ** Et, uint8 **EtD, long nrl, long nrh, long ncl, long nch){
	int i, j;
    int m;
    vuint8 l_1, l0, l1;
    vuint8 tmp_1,tmp0,tmp1;
    vuint8 res;
    //vuint8 xl, xr;
    //vuint8 y;
    //vuint8 ero;


    // Parcours de l'image
    for(i = nrl; i <= nrh; i++)
    {
        for(j = ncl; j < nch; j+=14)
        {


        		l_1 = _mm_loadu_si128((__m128i*)&Et[i-1][j-1]);
				l0  = _mm_loadu_si128((__m128i*)&Et[i][j-1]);
				l1  = _mm_loadu_si128((__m128i*)&Et[i+1][j-1]);


				res=_mm_or_si128(_mm_or_si128(l_1,l0),l1); //On met dans res l1 & l0 & l_1

				tmp0=_mm_srli_si128(res,1);
				tmp1=_mm_slli_si128(res,1);
				res=_mm_or_si128(_mm_or_si128(tmp0,tmp1),res);

				res = _mm_srli_si128(res,1);


        	/*l_1 = _mm_loadu_si128((__m128i*)&Et[i-1][j]);
			l0  = _mm_loadu_si128((__m128i*)&Et[i][j]);
			l1  = _mm_loadu_si128((__m128i*)&Et[i+1][j]);
			res=_mm_or_si128(_mm_or_si128(l_1,l0),l1); //On met dans res l1 & l0 & l_1

			tmp_1=_mm_srli_si128(l_1,1);
			tmp0=_mm_srli_si128(l0,1);
			tmp1=_mm_srli_si128(l1,1);
			res=_mm_or_si128(_mm_or_si128(_mm_or_si128(tmp_1,tmp0),tmp1),res);
			tmp_1=_mm_slli_si128(l_1,1);
			tmp0=_mm_slli_si128(l0,1);
			tmp1=_mm_slli_si128(l1,1);
			res=_mm_or_si128(_mm_or_si128(_mm_or_si128(tmp_1,tmp0),tmp1),res);
*/




            _mm_storeu_si128((__m128i*)&EtD[i][j], res); //Store dans Et la valeur de l'erosion
        }
        //j=nch;
        //EtD[i][j]=(Et[i-1][j-1]|Et[i-1][j]|Et[i-1][j+1]|Et[i][j-1]|Et[i][j]|Et[i][j+1]|Et[i+1][j-1]|Et[i+1][j]|Et[i+1][j+1]);
    }

}
void ouverture3SSE(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch){
	uint8 ** tmp = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    erosion3SSE(Et, tmp, nrl, nrh, ncl, nch);
    dilatation3SSE(tmp, Etout, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}
void fermeture3SSE(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch){
	uint8 ** tmp = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    dilatation3SSE(Et, tmp, nrl, nrh, ncl, nch);
    erosion3SSE(tmp, Etout, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}


void erosion5SSE(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch){
	uint8 ** tmp = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    erosion3SSE(Et, tmp, nrl, nrh, ncl, nch);
    erosion3SSE(tmp, EtE, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}
void dilatation5SSE(uint8 ** Et, uint8 **EtD, long nrl, long nrh, long ncl, long nch){
	uint8 ** tmp = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    dilatation3SSE(Et, tmp, nrl, nrh, ncl, nch);
    dilatation3SSE(tmp, EtD, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}
void ouverture5SSE(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch){
	uint8 ** tmp = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    erosion5SSE(Et, tmp, nrl, nrh, ncl, nch);
    dilatation5SSE(tmp, Etout, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}
void fermeture5SSE(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch){
	uint8 ** tmp = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    dilatation5SSE(Et, tmp, nrl, nrh, ncl, nch);
    erosion5SSE(tmp, Etout, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}


void erosion3SSE_bin(ulong64 ** Et, ulong64 **EtE, long nrl, long nrh, long ncl, long nch){

}

void dilatation3SSE_bin(ulong64 ** Et, ulong64 **EtD, long nrl, long nrh, long ncl, long nch){

}

void ouverture3SSE_bin(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch){

}

void fermeture3SSE_bin(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch){

}



