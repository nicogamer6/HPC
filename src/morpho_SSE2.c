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

    uint8 un = 255;

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
    uint8 zero = 0;

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
/*
	int i,j;
	ulong64 res;
	ulong64 binleft, binright; // Colonne gauche et droite i.e ulong64 avant et après

	for(i=nrl;i<=nrh;i++){
		for(j=ncl;j<=nch;j++){
			res = 0;
			res = ~res; //Passer tout à 1
			res &= (Et[i-1][j] & Et[i][j] & Et[i+1][j]);	// 3 premiers ulong64 pour les 3 1ères lignes

			binleft = (Et[i-1][j-1] & Et[i][j-1] & Et[i+1][j-1]);	// 3 ulong64 de gauche pour faire la morpho du dernier pixel du ulong64 res (bit poid fort)
			binleft = (res << 1ULL) | (binleft >> (NBBITS-1) & 1ULL);	// pour ajouter les 3 pixels du ulong de gauche pour faire l'erosion

			binright = (Et[i-1][j+1] & Et[i][j+1] & Et[i+1][j+1]);	// 3 ulong64 de droite pour faire la morpho du premier pixel du ulong64 res (bit poid faible)
			binright = ((res >> 1ULL) & ~(1ULL<<(NBBITS-1))) | (binright & 1ULL) << (NBBITS-1); // pour ajouter les 3 pixels du ulong de droite pour faire l'erosion


			EtE[i][j] = (res & binright & binleft);
		}
	}
	*/

	/*
	 * Il faut récupérer pour chaque ulong64(donc 1 truc sur 16 de chaque vuint) le pixel suivant et le pixel d'avant dans le ulong64 suivant ou précédent
	 * Il faut donc récup le poids faible du ulong64 suivant et le poids fort du ulong64 d'avant (si 1er du vuint, ça sera le dernier ulong64 du vuint d'avant)
	 */

	int i, j;
	vulong64 l_1, l0, l1;

	vulong64 res;
	vulong64 binleft, binright;

	vulong64 zero = init_vulong64(0); // Problème c'est que ce n'est pas des unsigned, faut faire attention pour les 64 bits à les ramener en positif
	//display_vulong64(zero," %d ","\ntest\n");
	vulong64 un = init_vulong64(1);
	vulong64 test = init_vulong64(~(1ULL<<(NBBITS-1)));


	// Parcours de l'image
	for(i = nrl; i <= nrh; i++)
	{
		for(j = ncl; j <= nch; j+=2) // pcqu'il y a 2 ulong64 par vulong64
		{
			res = ~zero; // tout à 1

			l_1 = _mm_loadu_si128((__m128i*)&Et[i-1][j]);
			l0  = _mm_loadu_si128((__m128i*)&Et[i][j]);
			l1  = _mm_loadu_si128((__m128i*)&Et[i+1][j]);

			res = _mm_and_si128(_mm_and_si128(_mm_and_si128(res,l_1),l0),l1);

			l_1 = _mm_loadu_si128((__m128i*)&Et[i-1][j-1]);
			l0  = _mm_loadu_si128((__m128i*)&Et[i][j-1]);
			l1  = _mm_loadu_si128((__m128i*)&Et[i+1][j-1]);

			binleft = _mm_and_si128(_mm_and_si128(l_1,l0),l1);
			binleft = _mm_or_si128((res<< 1ULL),_mm_and_si128((binleft>>(NBBITS-1)),un));


			binright = _mm_srli_si128(res,8); //2ème ulong64 de res qui est le vecteur de droite
			binright = _mm_or_si128(_mm_and_si128((res>> 1ULL),test),((_mm_and_si128(binright,un))<<(NBBITS-1)));

			res = _mm_and_si128(_mm_and_si128(res,binright),binleft);

			_mm_storeu_si128((__m128i*)&EtE[i][j], res); //Store dans Et la valeur de l'erosion

		}

	}

}

void dilatation3SSE_bin(ulong64 ** Et, ulong64 **EtD, long nrl, long nrh, long ncl, long nch){
/*
	int i,j;
	ulong64 res;
	ulong64 binleft, binright; // Colonne gauche et droite i.e ulong64 avant et après

	for(i=nrl;i<=nrh;i++){
		for(j=ncl;j<=nch;j++){
			res = 0;
			res |= (Et[i-1][j] | Et[i][j] | Et[i+1][j]); // 3 premiers ulong64 pour les 3 1ères lignes

			binleft = (Et[i-1][j-1] | Et[i][j-1] | Et[i+1][j-1]); // 3 ulong64 de gauche pour faire la morpho du dernier pixel du ulong64 res (bit poid fort)
			binleft = (res << 1ULL) | (binleft >> (NBBITS - 1) & 1ULL);	// pour ajouter les 3 pixels du ulong de gauche pour faire la dilatation

			binright = (Et[i-1][j+1] | Et[i][j+1] | Et[i+1][j+1]); // 3 ulong64 de droite pour faire la morpho du premier pixel du ulong64 res (bit poid faible)
			binright = ((res >> 1ULL) & ~(1ULL << (NBBITS - 1))) | (binright & 1ULL) << (NBBITS-1);	// pour ajouter les 3 pixels du ulong de droite pour faire la dilatation


			EtD[i][j] = (res | binright | binleft);
		}
	}
 */

	int i, j;
	vulong64 l_1, l0, l1;

	vulong64 res;
	vulong64 binleft, binright;

	vulong64 zero = init_vulong64(0);
	vulong64 un = init_vulong64(1);
	vulong64 test = init_vulong64(~(1ULL<<(NBBITS-1)));

	// Parcours de l'image
	for(i = nrl; i <= nrh; i++)
	{
		for(j = ncl; j <= nch; j+=2) // pcqu'il y a 2 ulong64 par vulong64
		{
			res = zero; // tout à 0

			l_1 = _mm_loadu_si128((__m128i*)&Et[i-1][j]);
			l0  = _mm_loadu_si128((__m128i*)&Et[i][j]);
			l1  = _mm_loadu_si128((__m128i*)&Et[i+1][j]);

			res = _mm_or_si128(_mm_or_si128(_mm_or_si128(res,l_1),l0),l1);

			l_1 = _mm_loadu_si128((__m128i*)&Et[i-1][j-1]);
			l0  = _mm_loadu_si128((__m128i*)&Et[i][j-1]);
			l1  = _mm_loadu_si128((__m128i*)&Et[i+1][j-1]);

			binleft = _mm_or_si128(_mm_or_si128(l_1,l0),l1);
			binleft = _mm_or_si128((res<< 1ULL),_mm_and_si128((binleft>>(NBBITS-1)),un));


			binright = _mm_srli_si128(res,8); //2ème ulong64 de res qui est le vecteur de droite
			binright = _mm_or_si128(_mm_and_si128((res>> 1ULL),test),((_mm_and_si128(binright,un))<<(NBBITS-1)));

			res = _mm_and_si128(_mm_and_si128(res,binright),binleft);

			_mm_storeu_si128((__m128i*)&EtD[i][j], res); //Store dans Et la valeur de la dilatation

		}

	}

}

void ouverture3SSE_bin(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch){
	erosion3SSE_bin(Et, tmp, nrl, nrh, ncl, nch);
	dilatation3SSE_bin(tmp, Etout, nrl, nrh, ncl, nch);
}

void fermeture3SSE_bin(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch){
	dilatation3SSE_bin(Et, tmp, nrl, nrh, ncl, nch);
	erosion3SSE_bin(tmp, Etout, nrl, nrh, ncl, nch);
}



