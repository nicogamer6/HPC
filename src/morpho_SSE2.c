#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "vnrdef.h"
#include "nrdef.h"
#include "morpho_SSE2.h"
#include "vnrutil.h"
#include <omp.h>

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
	 * Il faut récupérer pour chaque ulong64(donc 1 truc sur 16 de chaque vuint) le pixel suivant et le pixel d'avant dans le ulong64 suivant ou précédent
	 * Il faut donc récup le poids faible du ulong64 suivant et le poids fort du ulong64 d'avant (si 1er du vuint, ça sera le dernier ulong64 du vuint d'avant)
	 */

	int i, j;
	vulong64 l_1, l0, l1;

	vulong64 res;
	vulong64 binleft, binright;

	vulong64 zero = init_vulong64(0);
	//display_vulong64(zero," %d ","\ntest\n");
	vulong64 un = init_vulong64(1);

	//vulong64 bitfort = init_vulong64(~(1ULL<<(NBBITS-1)));
	vulong64 bitfort = _mm_slli_epi64(un,(NBBITS-1));
	bitfort = ~bitfort;


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
			//binleft = _mm_or_si128((res<< 1ULL),_mm_and_si128((binleft>>(NBBITS-1)),un));
			binleft = _mm_or_si128(_mm_slli_epi64(res,1),_mm_and_si128(_mm_slli_epi64(binleft,(NBBITS-1)),un));


			binright = _mm_srli_si128(res,8); //2ème ulong64 de res qui est le vecteur de droite
			//binright = _mm_or_si128(_mm_and_si128((res>> 1ULL),bitfort),((_mm_and_si128(binright,un))<<(NBBITS-1)));
			binright = _mm_or_si128(_mm_and_si128(_mm_srli_epi64(res,1),bitfort),(_mm_slli_epi64((_mm_and_si128(binright,un)),(NBBITS-1))));

			res = _mm_and_si128(_mm_and_si128(res,binright),binleft);

			_mm_storeu_si128((__m128i*)&EtE[i][j], res); //Store dans Et la valeur de l'erosion

		}

	}

}

void dilatation3SSE_bin(ulong64 ** Et, ulong64 **EtD, long nrl, long nrh, long ncl, long nch){

	int i, j;
	vulong64 l_1, l0, l1;

	vulong64 res;
	vulong64 binleft, binright;

	vulong64 zero = init_vulong64(0);
	vulong64 un = init_vulong64(1);

	//vulong64 bitfort = init_vulong64(~(1ULL<<(NBBITS-1)));
	vulong64 bitfort = (_mm_slli_epi64(un,(NBBITS-1)));
	bitfort = ~bitfort;

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
			//binleft = _mm_or_si128((res<< 1ULL),_mm_and_si128((binleft>>(NBBITS-1)),un));
			binleft = _mm_or_si128(_mm_slli_epi64(res,1),_mm_and_si128(_mm_slli_epi64(binleft,(NBBITS-1)),un));


			binright = _mm_srli_si128(res,8); //2ème ulong64 de res qui est le vecteur de droite
			//binright = _mm_or_si128(_mm_and_si128((res>> 1ULL),bitfort),((_mm_and_si128(binright,un))<<(NBBITS-1)));
			binright = _mm_or_si128(_mm_and_si128(_mm_srli_epi64(res,1),bitfort),(_mm_slli_epi64((_mm_and_si128(binright,un)),(NBBITS-1))));

			res = _mm_or_si128(_mm_or_si128(res,binright),binleft);

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





void erosion3SSE_line_bin(ulong64 ** Et, ulong64 **EtE, int i, long nrl, long nrh, long ncl, long nch){
	int j;
	vulong64 l_1, l0, l1;

	vulong64 res;
	vulong64 binleft, binright;

	vulong64 zero = init_vulong64(0);
	//display_vulong64(zero," %d ","\ntest\n");
	vulong64 un = init_vulong64(1);

	//vulong64 bitfort = init_vulong64(~(1ULL<<(NBBITS-1)));
	vulong64 bitfort = _mm_slli_epi64(un,(NBBITS-1));
	bitfort = ~bitfort;

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
		//binleft = _mm_or_si128((res<< 1ULL),_mm_and_si128((binleft>>(NBBITS-1)),un));
		binleft = _mm_or_si128(_mm_slli_epi64(res,1),_mm_and_si128(_mm_slli_epi64(binleft,(NBBITS-1)),un));


		binright = _mm_srli_si128(res,8); //2ème ulong64 de res qui est le vecteur de droite
		//binright = _mm_or_si128(_mm_and_si128((res>> 1ULL),bitfort),((_mm_and_si128(binright,un))<<(NBBITS-1)));
		binright = _mm_or_si128(_mm_and_si128(_mm_srli_epi64(res,1),bitfort),(_mm_slli_epi64((_mm_and_si128(binright,un)),(NBBITS-1))));

		res = _mm_and_si128(_mm_and_si128(res,binright),binleft);

		_mm_storeu_si128((__m128i*)&EtE[i][j], res); //Store dans Et la valeur de l'erosion

	}

}

void dilatation3SSE_line_bin(ulong64 ** Et, ulong64 **EtD, int i, long nrl, long nrh, long ncl, long nch){
	int j;
	vulong64 l_1, l0, l1;

	vulong64 res;
	vulong64 binleft, binright;

	vulong64 zero = init_vulong64(0);
	vulong64 un = init_vulong64(1);

	//vulong64 bitfort = init_vulong64(~(1ULL<<(NBBITS-1)));
	vulong64 bitfort = (_mm_slli_epi64(un,(NBBITS-1)));
	bitfort = ~bitfort;

	// Parcours de l'image
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
		//binleft = _mm_or_si128((res<< 1ULL),_mm_and_si128((binleft>>(NBBITS-1)),un));
		binleft = _mm_or_si128(_mm_slli_epi64(res,1),_mm_and_si128(_mm_slli_epi64(binleft,(NBBITS-1)),un));


		binright = _mm_srli_si128(res,8); //2ème ulong64 de res qui est le vecteur de droite
		//binright = _mm_or_si128(_mm_and_si128((res>> 1ULL),bitfort),((_mm_and_si128(binright,un))<<(NBBITS-1)));
		binright = _mm_or_si128(_mm_and_si128(_mm_srli_epi64(res,1),bitfort),(_mm_slli_epi64((_mm_and_si128(binright,un)),(NBBITS-1))));

		res = _mm_or_si128(_mm_or_si128(res,binright),binleft);

		_mm_storeu_si128((__m128i*)&EtD[i][j], res); //Store dans Et la valeur de la dilatation

	}

}


void ouverture3SSE_pipe_bin(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch){
	int i=nrl, j=ncl;

	erosion3SSE_line_bin(Et, tmp, i, nrl, nrh, ncl, nch); //1ère ligne

	for(i=nrl;i<=nrh-1;i++){
		erosion3SSE_line_bin(Et, tmp, i+1, nrl, nrh, ncl, nch);
		dilatation3SSE_line_bin(tmp, Etout, i, nrl, nrh, ncl, nch);
	}
	dilatation3SSE_line_bin(tmp, Etout, nrh, nrl, nrh, ncl, nch); //Dernière ligne
}

void fermeture3SSE_pipe_bin(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch){
	int i=nrl, j=ncl;

	dilatation3SSE_line_bin(Et, tmp, i, nrl, nrh, ncl, nch); //1ère ligne

	for(i=nrl;i<=nrh-1;i++){
		dilatation3SSE_line_bin(Et, tmp, i+1, nrl, nrh, ncl, nch);
		erosion3SSE_line_bin(tmp, Etout, i, nrl, nrh, ncl, nch);
	}
	erosion3SSE_line_bin(tmp, Etout, nrh, nrl, nrh, ncl, nch); //Dernière ligne
}



void erosion3SSE_line_binOMP(ulong64 ** Et, ulong64 **EtE, int i, long nrl, long nrh, long ncl, long nch){
	int j;
	vulong64 l_1, l0, l1;

	vulong64 res;
	vulong64 binleft, binright;

	vulong64 zero = init_vulong64(0);
	//display_vulong64(zero," %d ","\ntest\n");
	vulong64 un = init_vulong64(1);

	//vulong64 bitfort = init_vulong64(~(1ULL<<(NBBITS-1)));
	vulong64 bitfort = _mm_slli_epi64(un,(NBBITS-1));
	bitfort = ~bitfort;

	#pragma omp parallel for
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
		//binleft = _mm_or_si128((res<< 1ULL),_mm_and_si128((binleft>>(NBBITS-1)),un));
		binleft = _mm_or_si128(_mm_slli_epi64(res,1),_mm_and_si128(_mm_slli_epi64(binleft,(NBBITS-1)),un));


		binright = _mm_srli_si128(res,8); //2ème ulong64 de res qui est le vecteur de droite
		//binright = _mm_or_si128(_mm_and_si128((res>> 1ULL),bitfort),((_mm_and_si128(binright,un))<<(NBBITS-1)));
		binright = _mm_or_si128(_mm_and_si128(_mm_srli_epi64(res,1),bitfort),(_mm_slli_epi64((_mm_and_si128(binright,un)),(NBBITS-1))));

		res = _mm_and_si128(_mm_and_si128(res,binright),binleft);

		_mm_storeu_si128((__m128i*)&EtE[i][j], res); //Store dans Et la valeur de l'erosion

	}

}

void dilatation3SSE_line_binOMP(ulong64 ** Et, ulong64 **EtD, int i, long nrl, long nrh, long ncl, long nch){
	int j;
	vulong64 l_1, l0, l1;

	vulong64 res;
	vulong64 binleft, binright;

	vulong64 zero = init_vulong64(0);
	vulong64 un = init_vulong64(1);

	//vulong64 bitfort = init_vulong64(~(1ULL<<(NBBITS-1)));
	vulong64 bitfort = (_mm_slli_epi64(un,(NBBITS-1)));
	bitfort = ~bitfort;

	#pragma omp parallel for
	// Parcours de l'image
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
		//binleft = _mm_or_si128((res<< 1ULL),_mm_and_si128((binleft>>(NBBITS-1)),un));
		binleft = _mm_or_si128(_mm_slli_epi64(res,1),_mm_and_si128(_mm_slli_epi64(binleft,(NBBITS-1)),un));


		binright = _mm_srli_si128(res,8); //2ème ulong64 de res qui est le vecteur de droite
		//binright = _mm_or_si128(_mm_and_si128((res>> 1ULL),bitfort),((_mm_and_si128(binright,un))<<(NBBITS-1)));
		binright = _mm_or_si128(_mm_and_si128(_mm_srli_epi64(res,1),bitfort),(_mm_slli_epi64((_mm_and_si128(binright,un)),(NBBITS-1))));

		res = _mm_or_si128(_mm_or_si128(res,binright),binleft);

		_mm_storeu_si128((__m128i*)&EtD[i][j], res); //Store dans Et la valeur de la dilatation

	}

}


void ouverture3SSE_pipe_binOMP(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch){
	int i=nrl, j=ncl;

	erosion3SSE_line_binOMP(Et, tmp, i, nrl, nrh, ncl, nch); //1ère ligne

	#pragma omp ordered
	for(i=nrl;i<=nrh-1;i++){
		erosion3SSE_line_binOMP(Et, tmp, i+1, nrl, nrh, ncl, nch);
		dilatation3SSE_line_binOMP(tmp, Etout, i, nrl, nrh, ncl, nch);
	}
	dilatation3SSE_line_binOMP(tmp, Etout, nrh, nrl, nrh, ncl, nch); //Dernière ligne
}

void fermeture3SSE_pipe_binOMP(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch){
	int i=nrl, j=ncl;

	dilatation3SSE_line_binOMP(Et, tmp, i, nrl, nrh, ncl, nch); //1ère ligne

	#pragma omp ordered
	for(i=nrl;i<=nrh-1;i++){
		dilatation3SSE_line_binOMP(Et, tmp, i+1, nrl, nrh, ncl, nch);
		erosion3SSE_line_binOMP(tmp, Etout, i, nrl, nrh, ncl, nch);
	}
	erosion3SSE_line_binOMP(tmp, Etout, nrh, nrl, nrh, ncl, nch); //Dernière ligne
}

