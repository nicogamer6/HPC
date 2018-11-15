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
            
            // si a < SEUILFD donc dépasse le seuil alors a_0 reçoit 255 sinon 0

            b = _mm_cmplt_epi8(v_SEUILFD,a);

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
    vuint8 sel;
    vuint8 a = init_vuint8(0);
    vuint8 b = init_vuint8(0);
    vuint8 c = init_vuint8(0);

    vuint v_M, v_It, v_Mtm1, v_Ot, v_V, v_Vtm1;
    
    
    // MAX ET MIN
    vuint8 v_max = init_vuint8(VMAX);
    vuint8 v_min = init_vuint8(VMIN);
    vuint8 v_128 = init_vuint8(128);
    vuint8 v_255 = init_vuint8(255);
    int n1, n2, n3, n4;
    
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    
  ////////////////////Step 1 Estimation/////////////////////////////
    for (i=nrl; i<=nrh; i++)
    {
        for (j=ncl; j<=nch; j++)
        {
            //on fait les loads
            // pas besoin de load le M car on store à l'intérieur
           v_M = _mm_load_si128(&M[i][j]);
            v_Mtm1 = _mm_load_si128(&Mtm1[i][j]);
            v_It = _mm_load_si128(&It[i][j]);
        
            //Si M < I. Stocke 1 dans a si M < I, sinon 0

            sel = _mm_cmplt_epi8 (_mm_sub_epi8(v_Mtm1,v_128) , _mm_sub_epi8(v_It,v_128));
            v_M = vec_sel (_mm_add_epi8 (v_Mtm1 , one),v_Mtm1, sel);
            
            // Cas le plus grand
            sel = _mm_cmpgt_epi8 (_mm_sub_epi8(v_Mtm1,v_128) , _mm_sub_epi8(v_It,v_128));
            v_M = vec_sel (_mm_sub_epi8 (v_Mtm1 , one),v_M , sel);
            
            // égalité
            sel = _mm_cmpeq_epi8 (_mm_sub_epi8(v_Mtm1,v_128) , _mm_sub_epi8(v_It,v_128));
            v_M= vec_sel (v_Mtm1 , v_M , sel);
            
    
            
            //Enfin on fait le store de la valeur de v_m
            _mm_store_si128(&M[i][j], v_M);
            

            
            
       //////////////////// //Step 2 Difference Computation//////////////////
    
 
            //On fait les loads
            v_M = _mm_load_si128(&M[i][j]);
            v_It = _mm_load_si128(&It[i][j]);
            
            // Valeur absolue
            vuint8 max = _mm_max_epu8(v_M,v_It);
            vuint8 min = _mm_min_epu8(v_M,v_It);
            c = _mm_sub_epi8(max,min);
            
            // on stocke dans Ot
            _mm_store_si128(&Ot[i][j], c);
    
        
    
     ///////////////   //Step 3 Update and clamping//////////////////////////////
  
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

            //Si V < n *Ot. Stocke 1 dans a si M < I, sinon 0
            
            sel = _mm_cmplt_epi8 (_mm_sub_epi8(v_Vtm1,v_128) , _mm_sub_epi8(v_Ot,v_128));
            v_V = vec_sel (_mm_add_epi8 (v_Vtm1 , one), v_Vtm1 ,sel);
            
            sel = _mm_cmpgt_epi8 (_mm_sub_epi8(v_Vtm1,v_128) , _mm_sub_epi8(v_Ot,v_128));
            v_V = vec_sel (_mm_sub_epi8 (v_Vtm1 , one), v_V ,sel);
            
            
            sel = _mm_cmpeq_epi8 (_mm_sub_epi8(v_Vtm1,v_128) , _mm_sub_epi8(v_Ot,v_128));
            v_V = vec_sel (v_Vtm1 , v_V, sel);
     

           
 /*
            //si V < n* Ot. Stocke 1 dans a si V < n* Ot, sinon 0
            a = _mm_cmplt_epi8 (v_V , v_Ot);
            a = _mm_and_si128(a , v_01);
            v_V = _mm_add_epi8(v_V , a);
            
            b = _mm_cmpgt_epi8(v_V , v_Ot);
            b = _mm_and_si128(b , v_01);
            v_V = _mm_sub_epi8(v_V , b);
            
            Clamp to [VMIN,VMAX]
            
             a reçoit (V, VMAX)*/

            b = _mm_min_epu8(v_V, v_max);
            // a <- max(a, VMIN)
            b = _mm_max_epu8(a, v_min);
        
            _mm_store_si128(&V[i][j] , b);
            
  
    

    
       /////////////////// //Step 4 Estimation///////////////////////////
   
            //On fait les loads
            v_Ot = _mm_load_si128(&Ot[i][j]);
            v_V = _mm_load_si128(&V[i][j]);
            
            a = _mm_cmplt_epi8(_mm_sub_epi8(v_Ot,v_128) , _mm_sub_epi8(v_V,v_128));
            
            //Le retour de cmplt est de 1 si l'inégalité est vraie
            a = _mm_andnot_si128(a , v_255);
            
            _mm_store_si128(&Et[i][j] , a);
        }
    }
    display_vui8matrix(Et, n1, n2, n3, n4," %d ","\nimage scal\n");

    

}


