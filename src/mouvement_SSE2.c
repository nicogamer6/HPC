#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "vnrdef.h"
#include "nrdef.h"
#include "mouvement_SSE2.h"
#include "vnrutil.h"


////////////////////////////////////
//      FRAME DIFFERENCE SSE      //
////////////////////////////////////



void routine_FrameDifference_SSE2(vuint8** It, vuint8** It_1, vuint8** Ot, vuint8** Et, long nrl, long nrh, long ncl, long nch)
{
	int i, j;

	vuint8 v_SEUILFD = init_vuint8(SEUILFD);
	vuint8 v_255 = init_vuint8(255);
	vuint8 v_128 = init_vuint8(128);

	vuint8 a,b;


	vuint8 v_it;
	vuint8 v_it_1;

	for (i=nrl;i<=nrh;i++) {
		for(j=ncl;j<=nch;j++) {
		
			// Ici on load les images
			v_it = _mm_load_si128(&It[i][j]);
 			v_it_1 = _mm_load_si128(&It_1[i][j]);
		
			// On fait la valeur absolue
            a = _mm_abs_epi8 (_mm_sub_epi8(v_it , v_it_1));
 			// On sauvegarde la valeur de abs dans Ot
 			_mm_store_si128(&Ot[i][j], a);
		}
	}
    
	for (i=nrl;i<=nrh;i++) {
		for(j=ncl;j<=nch;j++) {
		
            // si Ot < SEUILFD donc dépasse le seuil alors a_0 reçoit 255 sinon 0

            a = _mm_cmplt_epi8(Ot[i][j],v_SEUILFD);

            //  On sauvegarde la valeur de b dans Et
            _mm_store_si128(&Et[i][j], a);
		}
	}
	// return Et;
}



////////////////////////////////////////
//       SIGMA DELTA STEP0 SSE        //
////////////////////////////////////////

void SigmaDelta_1step0_SSE2 (vuint8** V, vuint8** M, vuint8** It,long nrl, long nrh, long ncl, long nch)
{
    int i,j;
    vuint8 v_vmin = init_vuint8(VMIN);
    
    for(i=nrl; j<=ncl; j++)
    {
        for (j=ncl; j<=nch; j++)
        {
            // On stocke VMIN dans V
            _mm_store_si128(&V[i][j], v_vmin);
            //On stocke It dans M
            _mm_store_si128(&M[i][j], It[i][j]);
            
        }
    }
}



////////////////////////////////////////
//       SIGMA DELTA STEP1 SSE        //
////////////////////////////////////////


void SigmaDelta_1step_SSE2 (vuint8** It_1, vuint8** Ot, vuint8** M, vuint8** V, vuint8** Et,long nrl, long nrh, long ncl, long nch)
{
    int i,j;
    vuint8 v_01 = _mm_set1_epi8(0x01);
    vuint a,b,c;
    vuint v_m, v_it_1;
    
    
    //Step 1 Estimation
    for (i=nrl; i<=nrh; i++)
    {
        for (j=ncl; j<=nch; j++)
        {
            //on fait les loads
            v_m = _mm_load_si128(&M[i][j]);
            v_it_1 = _mm_load_si128(&It_1[i][j]);
            
            //si M < I. Stocke 1 dans a si M < I, sinon 0
            a = _mm_cmplt_epi8 (v_m,v_it_1);
            
            // Si l'inégalité es vérifiée, on rajoute +1 donc on fait un and entre 1 et le resultat de la comparaison
            a = _mm_and_si128(a, v_01);
            
            //On fait l'addition pour v_m et le resultat du and
            v_m = _mm_add_epi8 (v_m , a);
            
              //si M > I. Stocke 1 dans a si M > I, sinon 0
            b = _mm_cmpgt_epi8(_mm_sub_epi8(v_m, v_128), _mm_sub_epi8(v_it_1, v_128));
            
            // Si l'inégalité es vérifiée, on rajoute +1 donc on fait un and entre 1 et le resultat de la comparaison
            b = _mm_and_si128(b, v_01);
            
            //On fait l'addition pour v_m et le resultat du and
            v_m = _mm_add_epi8 (v_m , b);
            
            //Enfin on fait le store de la valeur de v_m
            _mm_store_si128(&M[i][j], v_m);
        }
    }
    
    
        //Step 2 Difference Computation
    
    for (i=nrl; i<=nrh; i++)
    {
        for (j=ncl; j<=nch; j++)
        {
            //On fait les loads
            v_M = _mm_load_si128(&M[i][j]);
            v_It_1 = _mm_load_si128(&It_1[i][j]);
            
            // Valeur absolue
            c = _mm_abs_epi8 (_mm_sub_epi8(v_M, v_It_1));
            
            // on stocke dans Ot
            _mm_store_si128(&Ot[i][j], c);
        }
    }
    
        //Step 3 Update and clamping
    for (i=nrl; i<=nrh; i++)
    {
        for (j=ncl; j<=nch; j++)
        {
            //On fait les loads
            v_M = _mm_load_si128(&M[i][j]);
            v_It_1 = _mm_load_si128(&It_1[i][j]);
            
        }
        
    }
        
    
    
    
}
