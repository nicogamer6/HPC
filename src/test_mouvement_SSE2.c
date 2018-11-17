#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "vnrdef.h"
#include "nrdef.h"
#include "vnrutil.h"
#include "mymacro.h"

#include "mouvement_SSE2.h"
#include "test_mouvement_SSE2.h"

#define NBIMAGES 299
#define BORD 2
#define VMIN 1
#define VMAX 255
#define N 4
#define vec_sel(a,b,c) _mm_or_si128(_mm_and_si128(c, a), _mm_andnot_si128(c, b));


/////////////////////////////////////////////////
//      CONVERSION SCALIRE EN VECTEUR          //
/////////////////////////////////////////////////




void uint_to_vuint(uint8 ** scalaire, vuint ** vecteur, int n1, int n2, int n3, int n4)
{
    int i,j,k;
    vuint8 tmp[1];
    uint8 *p = (uint8*) tmp;
    for (i = n1 ; i <= n2 ; i++)
    {
        for (j = n3 ;j <= n4  ;j++)
        {
            for (k = 0 ;k < 16;k++)
            {
                p[k] = scalaire[i][j*16+k];
            }
            vecteur[i][j]= tmp[0];
        }
    }
}

/////////////////////////////////////////////////
//      CONVERSION VECTEUR EN SCALAIRE         //
/////////////////////////////////////////////////



void vuint_to_uint(uint8 ** scalaire, vuint ** vecteur, int n1, int n2, int n3, int n4)
{
    int i,j,k;
    vuint8 a;
    vuint8 tmp[1];
    uint8 * b = (uint8*) tmp;
    for (i = n1 ; i <= n2 ; i++)
    {
        for (j = n3 ;j <= n4 ;j++)
        {
            a = _mm_load_si128(&vecteur[i][j]);
            _mm_store_si128(tmp, a);
            for (k = 0 ;k < 16 ;k++)
            {
                scalaire[i][j*16+k]= b[k];
            }
        }
    }
}



////////////////////////////////////////
//      TEST ROUTINE FD_SSE           //
////////////////////////////////////////



void test_routineFD_SSE(int seuil)
{
    //Cycle par point//
    double cycles, totalcy=0;
    int iter, niter=2;
    int run, nrun = 5;
    double t0,t1,dt,tmin,t;

    char *format = "%6.2f\n";
    ///////////////////

    long nrl, nrh, ncl, nch;
    int n1, n2, n3, n4;
    char nameload1[100], nameload2[100], namesave[100];
    //vuint8 vseuil=init_vuint8(seuil);


    sprintf(nameload1,"hall/hall000000.pgm");
    
    
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    
    vuint8 ** m = vui8matrix(n1,n2,n3,n4);
    vuint8 ** It = vui8matrix(n1,n2,n3,n4);
    vuint8 ** I0 = vui8matrix(n1,n2,n3,n4);
    uint_to_vuint(Itm1, I0, n1, n2, n3, n4);
    uint8** a = ui8matrix(nrl, nrh, ncl, nch);
    //uint8 **m, **m1 , **m2;
    //uint8** imaget1 = ui8matrix(nrl, nrh, ncl, nch);
    uint8** out = ui8matrix(nrl, nrh, ncl, nch);
    
    
    int i;
   
    for(i=1;i<=NBIMAGES;i++){
        sprintf(nameload2,"hall/hall000%03d.pgm",i);
        
        a = LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
        //TODO convertir uint** en vuint** (m)
        

        //TODO convertir m en uint **
        uint_to_vuint(a, It, n1, n2, n3, n4);
        CHRONO(routine_FrameDifference_SSE2(I0,It,m,n1, n2,n3,n4,seuil),cycles);
        totalcy += cycles;
        vuint_to_uint(out, m, n1, n2, n3, n4);

        
        sprintf(namesave,"testFD_SSE/hall000%03d.pgm",i);
        SavePGM_ui8matrix(out,nrl,nrh,ncl,nch,namesave);
        dup_vui8matrix(It, n1, n2, n3, n4, I0);
    }
    totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
    totalcy = totalcy / ((nch+1)*(nrh+1));

    BENCH(printf("Cycles FD_SSE = "));
    BENCH(printf(format,totalcy));
    //printf("%ld",((nch+1)*(nrh+1)));

    free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
    free_ui8matrix(a,nrl,nrh,ncl,nch);
    free_ui8matrix(out,nrl,nrh,ncl,nch);
    free_vui8matrix(m,n1,n2,n3,n4);
    free_vui8matrix(It,n1,n2,n3,n4);
    free_vui8matrix(I0,n1,n2,n3,n4);
    
}



////////////////////////////////////////
//      TEST UNITAIRE FD_SSE          //
////////////////////////////////////////


void test_unitaire_FD_SSE()
{
    vuint8 It = init_vuint8_all(255, 0, 128, 75, 255, 100, 1, 254, 9, 56, 38, 124, 198, 50, 87, 220);
    
    vuint8 It_1 = init_vuint8_all(254, 1, 127, 73, 200, 47, 30, 200, 40, 58, 38, 123, 198, 51, 0, 210);
    
    vuint8 Imagefinal;
    vuint8 seuil = init_vuint8(25);
    vuint8 tempSeuil = _mm_abs_epi8(_mm_sub_epi8(It,It_1));
    Imagefinal = _mm_cmpgt_epi8(tempSeuil,seuil);
    
   // vuint8 res= _mm_cmpeq_epi8(Imagefinal, ImageReel);
    display_vuint8(Imagefinal,"%d ","res=");
    printf("\n");
    
    
}



////////////////////////////////////////
//      TEST ROUTINE SD_SSE           //
////////////////////////////////////////


void test_routineSD_SSE()
{

    //Cycle par point//
    double cycles, totalcy=0;
    int iter, niter=2;
    int run, nrun = 5;
    double t0,t1,dt,tmin,t;

    char *format = "%6.2f\n";
    ///////////////////

    int n1, n2, n3, n4;
    long nrl, nrh, ncl, nch;
    char nameload[100];     //"hall/hall000..";
    char namesave[100];     //"testSD/SD...";
    int i;
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    sprintf(nameload,"hall/hall000000.pgm");
    
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);




   
    uint8 ** a = ui8matrix(nrl, nrh, ncl, nch);
    vuint8 ** I = vui8matrix(n1,n2,n3,n4);
    vuint8 ** V  = vui8matrix(n1,n2,n3,n4);
    vuint8 ** Vtm1 = vui8matrix(n1,n2,n3,n4);
    vuint8 ** M  = vui8matrix(n1,n2,n3,n4);
    vuint8 ** Mtm1 = vui8matrix(n1,n2,n3,n4);
    vuint8 ** Et = vui8matrix(n1,n2,n3,n4);
    //uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);
    
   
    vuint8 ** Iv = vui8matrix(n1,n2,n3,n4);
    uint8 ** res = ui8matrix(nrl, nrh, ncl, nch);
    
    vuint8 ** Itm1v = vui8matrix(n1,n2,n3,n4);
    uint_to_vuint(Itm1, Itm1v, n1, n2, n3, n4);
    
    SigmaDelta_step0_SSE2 (Mtm1, Vtm1, Itm1v, n1, n2, n3, n4);
    
    
    for(i=1;i<=NBIMAGES;i++){
        
        sprintf(nameload,"hall/hall000%03d.pgm",i);
        a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
        
        
        //Conversion
        uint_to_vuint(a, Iv, n1, n2, n3, n4);
        
        SigmaDelta_1step_SSE2(V, Vtm1, M, Iv, Mtm1 , Et, n1, n2, n3, n4);
        CHRONO(SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4),cycles);
        totalcy += cycles;
        
        sprintf(namesave,"testSD_SSE/hall000%03d.pgm",i);
        
        vuint_to_uint(res, Et, n1, n2, n3, n4);
        
        SavePGM_ui8matrix(res,nrl,nrh,ncl,nch,namesave);
        //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
     
        
        dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
        dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
        dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);
        
        
        
    }

    totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
    totalcy = totalcy / ((nch+1)*(nrh+1));
    
    BENCH(printf("Cycles SD_SSE = "));
    BENCH(printf(format,totalcy));

    free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
    free_ui8matrix(a,nrl,nrh,ncl,nch);
    free_ui8matrix(res,nrl,nrh,ncl,nch);
    free_vui8matrix(I,n1,n2,n3,n4);
    free_vui8matrix(V,n1,n2,n3,n4);
    free_vui8matrix(Vtm1,n1,n2,n3,n4);
    free_vui8matrix(M,n1,n2,n3,n4);
    free_vui8matrix(Mtm1,n1,n2,n3,n4);
    free_vui8matrix(Iv,n1,n2,n3,n4);
    free_vui8matrix(Itm1v,n1,n2,n3,n4);
    free_vui8matrix(Et,n1,n2,n3,n4);
    
    
}


////////////////////////////////////////
//      TEST UNITAIRE SD SSE          //
////////////////////////////////////////

void test_unitaire_SD_SSE()
{
    //vuint8 im =
    vuint8 M = init_vuint8_all(127, 255, 0, 0, 1, 255, 254, 50,127, 255, 0, 0, 1, 255, 254, 200);
    vuint8 It = init_vuint8_all( 128, 255, 0, 1, 0, 254, 255, 60,128, 245, 0, 1, 0, 254, 255, 1);
    vuint8 Vtm1 = init_vuint8_all( 128, 245, 0, 1, 0, 254, 255, 60,128, 245, 0, 1, 0, 254, 255, 60);
    
    //vuint8 M = init_vuint8_all(127, 250, 0, 2, 1, 255, 254, 50,127, 255, 5, 0, 1, 255, 254, 50);
   // vuint8 It = init_vuint8_all( 121, 254, 0, 0, 1, 254, 253, 49,124, 252, 0, 0, 0, 251, 250, 40);
    
    
    vuint8 Vtmoins1 = init_vuint8_all( 128, 245, 0, 1, 0, 254, 255, 60,128, 245, 0, 1, 0, 254, 255, 60);
    vuint8 v_M, Vt, Et, sel, Ot;
    vuint8 one = init_vuint8(1);
    vuint8 zero = init_vuint8(0);
    vuint8 v_128 = init_vuint8(128);
    vuint8 v_255 = init_vuint8(255);
    
    vuint8 v_min = init_vuint8(VMIN);
    vuint8 v_max = init_vuint8(VMAX);
    
    
    display_vuint8(It," %.3d ","It\n");
    display_vuint8(M," %.3d ","\nmtm1\n");
    
    
   /////////////// //Estimation////////////////////////


    //Si M < I. Stocke 1 dans a si M < I, sinon 0

    //  On rajoute 128 parce qu'il n'existe pas de comparaisons signées en SSE
    
    sel = _mm_cmplt_epi8 (_mm_sub_epi8(M,v_128) , _mm_sub_epi8(It,v_128));
    v_M = vec_sel (_mm_add_epi8 (M , one), M , sel);
    
    sel = _mm_cmpgt_epi8 (_mm_sub_epi8(M,v_128),_mm_sub_epi8(It,v_128));
    v_M = vec_sel (_mm_sub_epi8 (M , one), v_M , sel);
    
   
     // égalité
    sel = _mm_cmpeq_epi8 (_mm_sub_epi8(M,v_128) , _mm_sub_epi8(It,v_128));
    v_M = vec_sel (M, v_M , sel);
     
 
    
    printf("\n/////////////////TEST ESTIMATION//////////////");
    display_vuint8(v_M," %.3d ", "\nv_M\n");
    
       printf("\n\n\n");
    
    
    //////////////// Différence computation ////////////////
    
    // Valeur absolue
    vuint8 max = _mm_max_epu8(v_M,It);
    vuint8 min = _mm_min_epu8(v_M,It);
    Ot = _mm_sub_epi8(max,min);
   
    printf("\n/////////////TEST COMPUTATION////////////////");
    display_vuint8(Ot," %.3d ", "\nOt\n");
    
    printf("\n\n\n");
    
    // UPDATE AND CLAMPING
    
    //Ici ilfaut d'abord faire la multiplication de n et Ot
    vuint8 a = init_vuint8(0);
    vuint8 v_Ot = init_vuint8(0);
    a= Ot;

    for(int k = 0; k < N; k++)    {
        v_Ot = _mm_adds_epu8(v_Ot , a);
    }
    // De base c'est non signé, on soustrait 128 pour les recentrer et faitre les comparaisons.
    //On a un problème avec le bit de signe
    
    sel = _mm_cmplt_epi8 (_mm_sub_epi8(Vtm1,v_128) , _mm_sub_epi8(v_Ot,v_128));
    Vt = vec_sel  (_mm_add_epi8 (Vtm1 , one), Vtm1 , sel);

    //Si V < n *Ot. Stocke 1 dans a si M < I, sinon 0
    
    sel = _mm_cmpgt_epi8 (_mm_sub_epi8(Vtm1,v_128) , _mm_sub_epi8(v_Ot,v_128));
    Vt = vec_sel (_mm_sub_epi8 (Vtm1 , one),Vt ,sel);

    // Egalité
    sel = _mm_cmpeq_epi8 (_mm_sub_epi8(Vtm1,v_128) , _mm_sub_epi8(v_Ot,v_128));
    Vt = vec_sel (Vtm1 , Vt , sel);
    
    vuint8 b;
    b = _mm_min_epu8(Vt, v_max);
    b = _mm_max_epu8(b, v_min);
    
    
  
    printf("\n////////////////////TEST UPDATE AND CLAMPING////////////////////");
    display_vuint8(v_Ot," %.3d ","\nOtavant\n");

    
   display_vuint8(Vtm1," %.3d ","\nVtm1\n");
    display_vuint8(v_Ot," %.3d ","\nOt\n");


    display_vuint8(Vt," %.3d ", "\nVt\n");
    display_vuint8(b," %.3d ", "\nB\n");
    printf("\n");
    
    
    
    // ESTIMATION

    Et = _mm_cmplt_epi8(_mm_sub_epi8(v_Ot,v_128) , _mm_sub_epi8(Vt,v_128));
    
    //Le retour de cmplt est de 1 si l'inégalité est vraie
    Et = _mm_andnot_si128(Et , v_255);
    display_vuint8(Et," %.3d ", "\nEt\n");
    printf("\n");
    
}
