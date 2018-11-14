#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "vnrdef.h"
#include "nrdef.h"
#include "mouvement_SSE2.h"
#include "vnrutil.h"

#define NB_IMAGE 299
#define SEUILFD 25
#define VMIN 1
#define VMAX 254
#define N 2
//SI a = 1 alors resultat = a sinon resultat = b
#define vec_sel(a,b,c) _mm_or_si128(_mm_and_si128(c,a),_mm_andnot_si128(c,b));


////////////////////////////////////
//      FRAME DIFFERENCE SSE      //
////////////////////////////////////



void routine_FrameDifference_SSE2(vuint8** It, vuint8** It_1, vuint8** Et, long nrl, long nrh, long ncl, long nch, int seuil)
{
	int i, j;

	vuint8 v_SEUILFD = init_vuint8(seuil);
	vuint8 v_255 = init_vuint8(255);
	vuint8 v_128 = init_vuint8(128);

	vuint8 a,b;


	vuint8 v_it;
	vuint8 v_it_1;
	vuint8 **Ot;

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


void SigmaDelta_1step_SSE2 (vuint8** V,vuint8** Vtm1, vuint8** M, vuint8** Mtm1, vuint8** It, vuint8** Et,long nrl, long nrh, long ncl, long nch)
{
    int i,j,k;
    vuint8 **Ot = vui8matrix(nrl,nrh,ncl,nch);
    
    vuint8 one = _mm_set1_epi8(0x01);
   // vuint8 v_N = _mm_set1_epi8((char)N);
    vuint a,b,c;
    vuint v_M, v_It, v_Mtm1, v_Ot, v_V, v_Vtm1;
    
    
    
    // MAX ET MIN
    vuint8 v_max = init_vuint8(VMAX);
    vuint8 v_min = init_vuint8(VMIN);
    vuint8 v_128 = init_vuint8(128);
    
    //Step 1 Estimation
    for (i=nrl; i<=nrh; i++)
    {
        for (j=ncl; j<=nch; j++)
        {
            //on fait les loads
            v_M = _mm_load_si128(&M[i][j]);
            v_Mtm1 = _mm_load_si128(&Mtm1[i][j]);
            v_It = _mm_load_si128(&It[i][j]);
            
            
            
            // égalité
            a = _mm_cmpeq_epi8 (_mm_sub_epi8(v_M,v_128) , _mm_sub_epi8(v_It,v_128));
            v_M = vec_sel (v_Mtm1 , _mm_sub_epi8 (v_Mtm1 , one) , a);
            
            //Si M < I. Stocke 1 dans a si M < I, sinon 0
        
            b = _mm_cmplt_epi8 (_mm_sub_epi8(v_M,v_128) , _mm_sub_epi8(v_It,v_128));
            v_M = vec_sel (_mm_add_epi8 (v_Mtm1 , one),_mm_sub_epi8 (v_Mtm1 , one),b);
            
            c = _mm_cmpgt_epi8 (_mm_sub_epi8(v_M,v_128) , _mm_sub_epi8(v_It,v_128));
            v_M = vec_sel (_mm_sub_epi8 (v_Mtm1 , one),_mm_add_epi8 (v_Mtm1 , one),b);
            

            //Enfin on fait le store de la valeur de v_m
            _mm_store_si128(&M[i][j], v_M);
        }
    }
    
    
        //Step 2 Difference Computation
    
    for (i=nrl; i<=nrh; i++)
    {
        for (j=ncl; j<=nch; j++)
        {
            //On fait les loads
            v_M = _mm_load_si128(&M[i][j]);
            v_It = _mm_load_si128(&It[i][j]);
            
            // Valeur absolue
            c = _mm_abs_epi8 (_mm_sub_epi8(v_M, v_It));
            
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
            v_Ot = _mm_load_si128(&Ot[i][j]);
            v_V = _mm_load_si128(&V[i][j]);
            v_Vtm1 = _mm_load_si128(&Vtm1[i][j]);
            
            //Ici ilfaut d'abord faire la multiplication de n et Ot
            a = v_Ot;
            for(k = 0; k < N; k++)    {
                v_Ot = _mm_adds_epu8(v_Ot, a);
            }
            
            // Egalité
            a = _mm_cmpeq_epi8 (_mm_sub_epi8(v_Vtm1,v_128) , _mm_sub_epi8(v_Ot,v_128));
            v_V = vec_sel (v_Vtm1 , _mm_sub_epi8 (v_Vtm1 , one) , a);

            //Si V < n *Ot. Stocke 1 dans a si M < I, sinon 0
            
            b = _mm_cmplt_epi8 (_mm_sub_epi8(v_Vtm1,v_128) , _mm_sub_epi8(v_Ot,v_128));
            v_V = vec_sel (_mm_add_epi8 (v_Vtm1 , one),_mm_sub_epi8 (v_Vtm1 , one),b);
            
            c = _mm_cmpgt_epi8 (_mm_sub_epi8(v_Vtm1,v_128) , _mm_sub_epi8(v_Ot,v_128));
            v_V = vec_sel (_mm_sub_epi8 (v_Vtm1 , one),_mm_add_epi8 (v_Vtm1 , one),b);
            
           
 /*
            //si V < n* Ot. Stocke 1 dans a si V < n* Ot, sinon 0
            a = _mm_cmplt_epi8 (v_V , v_Ot);
            a = _mm_and_si128(a , v_01);
            v_V = _mm_add_epi8(v_V , a);
            
            b = _mm_cmpgt_epi8(v_V , v_Ot);
            b = _mm_and_si128(b , v_01);
            v_V = _mm_sub_epi8(v_V , b);*/
            
            //Clamp to [VMIN,VMAX]
            
            // a reçoit (V, VMAX)
            a = _mm_min_epu8(v_V, v_max);
            // a <- max(a, VMIN)
            a = _mm_max_epu8(a, v_min);
        
            _mm_store_si128(&V[i][j] , a);
            
        }
    }
    
    vuint8 v_255 = init_vuint8(255);
    
        //Step 4 Estimation
    for (i=nrl; i<=nrh; i++)
    {
        for (j=ncl; j<=nch; j++)
        {
            //On fait les loads
            v_Ot = _mm_load_si128(&Ot[i][j]);
            v_V = _mm_load_si128(&V[i][j]);
            
            a = _mm_cmplt_epi8(v_Ot , v_V);
            
            //Le retour de cmplt est de 1 si l'inégalité est vraie
            a = _mm_andnot_si128(a , v_255);
            
            _mm_store_si128(&Et[i][j] , a);
        }
    }
    

}


