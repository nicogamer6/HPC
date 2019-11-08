#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "vnrdef.h"
#include "nrdef.h"
#include "mouvement_SSE2.h"
#include "vnrutil.h"
#include <omp.h>

#define NB_IMAGE 299
#define SEUILFD 25
#define VMIN 20
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
            _mm_store_si128(&Et[i][j], b);
        }
    }
    // return Et;
}

void routine_FrameDifference_SSE2_OMP(vuint8** It, vuint8** It_1, vuint8** Et, long nrl, long nrh, long ncl, long nch, int seuil)
{
    int i, j;

    vuint8 v_SEUILFD = init_vuint8(seuil);
    vuint8 v_255 = init_vuint8(255);
    vuint8 v_128 = init_vuint8(128);

    vuint8 a,b;
    vuint8 v_it;
    vuint8 v_it_1;

    #pragma omp parallel for schedule(static)
    
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
                _mm_store_si128(&Et[i][j], b);
            }
        }
    
    // return Et;
}


////////////////////////////////////////
//       SIGMA DELTA STEP0 SSE        //
////////////////////////////////////////

void SigmaDelta_step0_SSE2 (vuint8** M, vuint8** V,  vuint8** It,long nrl, long nrh, long ncl, long nch)
{
    int i,j;
    vuint8 v_vmin = init_vuint8(VMIN);
    vuint8 I;
    
    for(i=nrl; i<=nrh; i++)
    {
        for (j=ncl; j<=nch; j++)
        {
             //On stocke It dans M
            I = _mm_load_si128((vuint8*)&It[i][j]);
            _mm_store_si128(&M[i][j], I);
            // On stocke VMIN dans V
            V[i][j] = v_vmin;
        }
    }
}


void SigmaDelta_step0_SSE2_OMP (vuint8** M, vuint8** V,  vuint8** It,long nrl, long nrh, long ncl, long nch)
{
    int i,j;
    vuint8 v_vmin = init_vuint8(VMIN);
    vuint8 I;
    
    #pragma omp parallel for
    
        for(i=nrl; i<=nrh; i++)
        {
            for (j=ncl; j<=nch; j++)
            {
                 //On stocke It dans M
                I = _mm_load_si128((vuint8*)&It[i][j]);
                _mm_store_si128(&M[i][j], I);
                // On stocke VMIN dans V
                V[i][j] = v_vmin;
            }
        }
    
}


////////////////////////////////////////
//       SIGMA DELTA STEP1 SSE        //
////////////////////////////////////////


void SigmaDelta_1step_SSE2 (vuint8** V,vuint8** Vtm1, vuint8** M, vuint8** Mtm1, vuint8** It, vuint8** Et,long nrl, long nrh, long ncl, long nch)
{
  //  vuint8 Mtmoins1, Mt, It, Vtmoins1, Vt, Et, res;
    vuint8 v_M, v_It, v_Mtm1, v_V, v_Vtm1, sel;
    vuint8 Ot;
    
    vuint8 v_min = init_vuint8(VMIN);
    vuint8 v_max = init_vuint8(VMAX);
    vuint8 un = init_vuint8(1);
    vuint8 zero = init_vuint8(0);
    vuint8 v_128 = init_vuint8(128);
    vuint8 v_255 = init_vuint8(255);

    
    for (int i = nrl ; i <= nrh ; i++){
        for (int j = ncl ; j <= nch ; j++){
          
            v_Vtm1 = _mm_load_si128((vuint8*) &Vtm1[i][j]);
            v_It = _mm_load_si128((vuint8*) &It[i][j]);
            v_Mtm1 = _mm_load_si128((vuint8*) &Mtm1[i][j]);
            
            //display_vuint8(v_Vtm1," %d ","\ndebug\n");

            
            ////////////////////Step 1 Estimation/////////////////////////////

 
            //Si M < I. Stocke 1 dans a si M < I, sinon 0

            sel = _mm_cmplt_epi8(_mm_sub_epi8(v_Mtm1,v_128), _mm_sub_epi8(v_It, v_128));
            v_M = vec_sel( _mm_add_epi8(v_Mtm1, un), v_Mtm1, sel);
            
            // Cas le plus grand

            sel = _mm_cmpgt_epi8(_mm_sub_epi8(v_Mtm1,v_128), _mm_sub_epi8(v_It, v_128));
            v_M = vec_sel( _mm_sub_epi8(v_Mtm1, un), v_M, sel);
            
            // égalité

            sel = _mm_cmpeq_epi8(_mm_sub_epi8(v_Mtm1,v_128), _mm_sub_epi8(v_It, v_128));
            v_M = vec_sel(v_Mtm1,  v_M, sel);
            
            _mm_store_si128(&M[i][j],v_M);
       
            
            //////////////////// //Step 2 Difference Computation//////////////////
            
            vuint8 max = init_vuint8(0);
            max = _mm_max_epu8(v_M,v_It);
            vuint8 min = init_vuint8(0);
            min = _mm_min_epu8(v_M,v_It);
            Ot = _mm_sub_epi8(max,min);
           

            ///////////////   //Step 3 Update and clamping//////////////////////////////


            vuint8 v_Ot = init_vuint8(0);;
            
            //Ici ilfaut d'abord faire la multiplication de n et Ot

            
            for(int k = 0; k < N; k++)    {
                v_Ot = _mm_adds_epu8(v_Ot, Ot);
            }
            

            sel = _mm_cmplt_epi8(_mm_sub_epi8(v_Vtm1,v_128), _mm_sub_epi8(v_Ot,v_128));
            v_V = vec_sel(_mm_add_epi8(v_Vtm1, un), v_Vtm1, sel);
            
            sel = _mm_cmpgt_epi8(_mm_sub_epi8(v_Vtm1,v_128), _mm_sub_epi8(v_Ot,v_128));
            v_V = vec_sel(_mm_sub_epi8(v_Vtm1, un), v_V, sel);
            
            sel = _mm_cmpeq_epi8(_mm_sub_epi8(v_Vtm1,v_128), _mm_sub_epi8(v_Ot, v_128));
            v_V = vec_sel(v_Vtm1, v_V, sel);
            
            
            //MAX ET MIN
            
            v_V = _mm_min_epu8(v_V, v_max);
            // a <- max(a, VMIN)
            v_V = _mm_max_epu8(v_V, v_min);
            
            _mm_store_si128(&V[i][j],v_V);
           
    
            
       /////////////////// //Step 4 Estimation///////////////////////////
            
            vuint8 a  = _mm_cmplt_epi8(_mm_sub_epi8(Ot,v_128), _mm_sub_epi8(v_V,v_128));
            
            //Le retour de cmplt est de 1 si l'inégalité est vraie
            a = _mm_andnot_si128(a , v_255);
            _mm_store_si128(&Et[i][j],a);
        }
    }
}







////////////////////////////////////////
//       SIGMA DELTA STEP1 SSE AoSoA       //
////////////////////////////////////////

/*
void SigmaDelta_1step_SSE2_AoSoA (SoA *Vm, vuint8** V, vuint8** M, vuint8** Et,long nrl, long nrh, long ncl, long nch)
{
  //  vuint8 Mtmoins1, Mt, It, Vtmoins1, Vt, Et, res;
    vuint8 v_M, v_It, v_Mtm1, v_V, v_Vtm1, sel;
    vuint8 Ot;

    // Créer un tableau de vuint8 et stocker au fur et à mesure V I et M
    
    vuint8 v_min = init_vuint8(VMIN);
    vuint8 v_max = init_vuint8(VMAX);
    vuint8 un = init_vuint8(1);
    vuint8 zero = init_vuint8(0);
    vuint8 v_128 = init_vuint8(128);
    vuint8 v_255 = init_vuint8(255);

    
    for (int i = nrl ; i <= nrh ; i++){
        for (int j = ncl ; j <= nch ; j++){
          
            
            v_Vtm1 = _mm_load_si128((vuint8*) &Vtm1[i][j]); Vm[0].p1[i][j]
            v_It = _mm_load_si128((vuint8*) &It[i][j]) Vm[0].p2[i][j]
            v_Mtm1 = _mm_load_si128((vuint8*) &Mtm1[i][j]);
            
            //display_vuint8(v_Vtm1," %d ","\ndebug\n");

            
            ////////////////////Step 1 Estimation/////////////////////////////

 
            //Si M < I. Stocke 1 dans a si M < I, sinon 0

            sel = _mm_cmplt_epi8(_mm_sub_epi8(v_Mtm1,v_128), _mm_sub_epi8(v_It, v_128));
            v_M = vec_sel( _mm_add_epi8(v_Mtm1, un), v_Mtm1, sel);
            
            // Cas le plus grand

            sel = _mm_cmpgt_epi8(_mm_sub_epi8(v_Mtm1,v_128), _mm_sub_epi8(v_It, v_128));
            v_M = vec_sel( _mm_sub_epi8(v_Mtm1, un), v_M, sel);
            
            // égalité

            sel = _mm_cmpeq_epi8(_mm_sub_epi8(v_Mtm1,v_128), _mm_sub_epi8(v_It, v_128));
            v_M = vec_sel(v_Mtm1,  v_M, sel);
            
            _mm_store_si128(&M[i][j],v_M);
       
            
            //////////////////// //Step 2 Difference Computation//////////////////
            
            vuint8 max = init_vuint8(0);
            max = _mm_max_epu8(v_M,v_It);
            vuint8 min = init_vuint8(0);
            min = _mm_min_epu8(v_M,v_It);
            Ot = _mm_sub_epi8(max,min);
           

            ///////////////   //Step 3 Update and clamping//////////////////////////////


            vuint8 v_Ot = init_vuint8(0);;
            
            //Ici ilfaut d'abord faire la multiplication de n et Ot

            
            for(int k = 0; k < N; k++)    {
                v_Ot = _mm_adds_epu8(v_Ot, Ot);
            }
            

            sel = _mm_cmplt_epi8(_mm_sub_epi8(v_Vtm1,v_128), _mm_sub_epi8(v_Ot,v_128));
            v_V = vec_sel(_mm_add_epi8(v_Vtm1, un), v_Vtm1, sel);
            
            sel = _mm_cmpgt_epi8(_mm_sub_epi8(v_Vtm1,v_128), _mm_sub_epi8(v_Ot,v_128));
            v_V = vec_sel(_mm_sub_epi8(v_Vtm1, un), v_V, sel);
            
            sel = _mm_cmpeq_epi8(_mm_sub_epi8(v_Vtm1,v_128), _mm_sub_epi8(v_Ot, v_128));
            v_V = vec_sel(v_Vtm1, v_V, sel);
            
            
            //MAX ET MIN
            
            v_V = _mm_min_epu8(v_V, v_max);
            // a <- max(a, VMIN)
            v_V = _mm_max_epu8(v_V, v_min);
            
            _mm_store_si128(&V[i][j],v_V);
           
    
            
       /////////////////// //Step 4 Estimation///////////////////////////
            
            vuint8 a  = _mm_cmplt_epi8(_mm_sub_epi8(Ot,v_128), _mm_sub_epi8(v_V,v_128));
            
            //Le retour de cmplt est de 1 si l'inégalité est vraie
            a = _mm_andnot_si128(a , v_255);
            _mm_store_si128(&Et[i][j],a);
        }
    }
}

*/



void SigmaDelta_1step_SSE2_OMP (vuint8** V,vuint8** Vtm1, vuint8** M, vuint8** Mtm1, vuint8** It, vuint8** Et,long nrl, long nrh, long ncl, long nch)
{
  //  vuint8 Mtmoins1, Mt, It, Vtmoins1, Vt, Et, res;
    vuint8 v_M, v_It, v_Mtm1, v_V, v_Vtm1, sel;
    vuint8 Ot;

    
    vuint8 v_min = init_vuint8(VMIN);
    vuint8 v_max = init_vuint8(VMAX);
    vuint8 un = init_vuint8(1);
    vuint8 zero = init_vuint8(0);
    vuint8 v_128 = init_vuint8(128);
    vuint8 v_255 = init_vuint8(255);

    #pragma omp parallel for
    
        for (int i = nrl ; i <= nrh ; i++){
            for (int j = ncl ; j <= nch ; j++){
              
                
                v_Vtm1 = _mm_load_si128((vuint8*) &Vtm1[i][j]);
                v_It = _mm_load_si128((vuint8*) &It[i][j]);
                v_Mtm1 = _mm_load_si128((vuint8*) &Mtm1[i][j]);
                
                //display_vuint8(v_Vtm1," %d ","\ndebug\n");

                
                ////////////////////Step 1 Estimation/////////////////////////////

     
                //Si M < I. Stocke 1 dans a si M < I, sinon 0

                sel = _mm_cmplt_epi8(_mm_sub_epi8(v_Mtm1,v_128), _mm_sub_epi8(v_It, v_128));
                v_M = vec_sel( _mm_add_epi8(v_Mtm1, un), v_Mtm1, sel);
                
                // Cas le plus grand

                sel = _mm_cmpgt_epi8(_mm_sub_epi8(v_Mtm1,v_128), _mm_sub_epi8(v_It, v_128));
                v_M = vec_sel( _mm_sub_epi8(v_Mtm1, un), v_M, sel);
                
                // égalité

                sel = _mm_cmpeq_epi8(_mm_sub_epi8(v_Mtm1,v_128), _mm_sub_epi8(v_It, v_128));
                v_M = vec_sel(v_Mtm1,  v_M, sel);
                
                _mm_store_si128(&M[i][j],v_M);
           
                
                //////////////////// //Step 2 Difference Computation//////////////////
                
                vuint8 max = init_vuint8(0);
                max = _mm_max_epu8(v_M,v_It);
                vuint8 min = init_vuint8(0);
                min = _mm_min_epu8(v_M,v_It);
                Ot = _mm_sub_epi8(max,min);
               

                ///////////////   //Step 3 Update and clamping//////////////////////////////


                vuint8 v_Ot = init_vuint8(0);;
                
                //Ici ilfaut d'abord faire la multiplication de n et Ot

                
                for(int k = 0; k < N; k++)    {
                    v_Ot = _mm_adds_epu8(v_Ot, Ot);
                }
                

                sel = _mm_cmplt_epi8(_mm_sub_epi8(v_Vtm1,v_128), _mm_sub_epi8(v_Ot,v_128));
                v_V = vec_sel(_mm_add_epi8(v_Vtm1, un), v_Vtm1, sel);
                
                sel = _mm_cmpgt_epi8(_mm_sub_epi8(v_Vtm1,v_128), _mm_sub_epi8(v_Ot,v_128));
                v_V = vec_sel(_mm_sub_epi8(v_Vtm1, un), v_V, sel);
                
                sel = _mm_cmpeq_epi8(_mm_sub_epi8(v_Vtm1,v_128), _mm_sub_epi8(v_Ot, v_128));
                v_V = vec_sel(v_Vtm1, v_V, sel);
                
                
                //MAX ET MIN
                
                v_V = _mm_min_epu8(v_V, v_max);
                // a <- max(a, VMIN)
                v_V = _mm_max_epu8(v_V, v_min);
                
                _mm_store_si128(&V[i][j],v_V);
               
        
                
           /////////////////// //Step 4 Estimation///////////////////////////
                
                vuint8 a  = _mm_cmplt_epi8(_mm_sub_epi8(Ot,v_128), _mm_sub_epi8(v_V,v_128));
                
                //Le retour de cmplt est de 1 si l'inégalité est vraie
                a = _mm_andnot_si128(a , v_255);
                _mm_store_si128(&Et[i][j],a);
            }
        }
    
}



