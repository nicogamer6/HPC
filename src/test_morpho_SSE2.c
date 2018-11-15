#include <stdio.h>
#include "nrutil.h"
#include "nrdef.h"
#include "string.h"
#include "test_morpho.h"
#include "morpho.h"
#include "mouvement.h"
#include "mouvement_SSE2.h"
#include "test_mouvement_SSE2.h"
#include "test_morpho_SSE2.h"

#include "vnrdef.h"
#include "morpho_SSE2.h"
#include "vnrutil.h"

//#define SEUILFD 20
#define NBIMAGES 299
#define BORD 2


void test_EtapemorphoSSE(void){
    long nrl, nrh, ncl, nch;
    int vi0, vi1, vj0, vj1;

    vuint8 tmp;
    
    uint8 **Et;

    int i=0,j=0;
    
    Et= LoadPGM_ui8matrix("smile.pgm",&nrl,&nrh,&ncl,&nch);
    //s2v(nrl, nrh, ncl, nch, card_vuint8(), &vi0, &vi1, &vj0, &vj1);
    //vuint8 ** It = vui8matrix(vi0,vi1,vj0,vj1);
    uint8 **Etbord = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    uint8 **Etout = ui8matrix(nrl, nrh, ncl, nch);

    //uint_to_vuint(Et, It, vi0, vi1, vj0, vj1);
    

    /*for(i=vi0;i<=vi1;i++){
        for(j=vj0;j<=vj1;j++){
            //_mm_store_si128(&Itbord[i][j],It[i][j]);
            tmp = _mm_load_si128(&It[i][j]);
            _mm_store_si128(&Itbord[i][j], tmp);
        }
    }*/
    //dup_vui8matrix(It, vi0, vi1, vj0, vj1, Itbord);
    //**Etbord=Et[1][1];
    copy_ui8matrix_ui8matrix(Et,nrl,nrh,ncl,nch,Etbord);
    
    erosion3SSE(Etbord,Etout,nrl,nrh,ncl,nch);
    //display_vui8matrix(Itbord,vi0-BORD,vi1+BORD,vj0-BORD,vj1+BORD,"%3d","\nItbord\n");
    //vuint_to_uint(Etout, It, vi0, vi1, vj0, vj1);
    //display_ui8matrix(Etout ,nrl,nrh,ncl,nch,"%3d","Etout");
    SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,"testmorphoSSE/testeroSSE.pgm");
    
    dilatation3SSE(Etbord,Etout,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,"testmorphoSSE/testdilSSE.pgm");
    
    
    ouverture3(Etbord,Etout,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,"testmorphoSSE/testouvSSE.pgm");
    
    fermeture3(Etbord,Etout,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,"testmorphoSSE/testfermSSE.pgm");
    
}


///////////////////////////////////////////////////
//     TEST FRAME DIFFERENCE SSE O,OF,F,FO       //
///////////////////////////////////////////////////

void test_routineFD_SSEmorpho3xOuv(int seuil){
    
    long nrl, nrh, ncl, nch;
    int vi0, vi1, vj0, vj1;
    char nameload1[100], nameload2[100], namesave[100];


    sprintf(nameload1,"hall/hall000000.pgm");
    
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &vi0, &vi1, &vj0, &vj1);

    
    vuint8 ** m = vui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    vuint8 ** It = vui8matrix(nrl,nrh,ncl,nch);
    vuint8 ** It0 = vui8matrix(nrl,nrh,ncl,nch);
    uint_to_vuint(Itm1, It0, vi0, vi1, vj0, vj1);
    uint8** a = ui8matrix(nrl, nrh, ncl, nch);
    uint8** out = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout = ui8matrix(nrl, nrh, ncl, nch);
    
    
    int i;
   
    for(i=1;i<=NBIMAGES;i++){
        sprintf(nameload2,"hall/hall000%03d.pgm",i);
        
        a = LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
        
        uint_to_vuint(a, It, vi0, vi1, vj0, vj1);
        routine_FrameDifference_SSE2(It0,It,m,vi0, vi1,vj0,vj1,seuil);
        vuint_to_uint(out, m, vi0, vi1, vj0, vj1);


        ouverture3SSE(out,Etout,nrl,nrh,ncl,nch);

        
        sprintf(namesave,"testFD_SSEmorphoO/hall000%03d.pgm",i);
        SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
        dup_vui8matrix(It, vi0, vi1, vj0, vj1, It0);
    }
}

void test_routineFD_SSEmorpho3xFerm(int seuil){
    long nrl, nrh, ncl, nch;
    int vi0, vi1, vj0, vj1;
    char nameload1[100], nameload2[100], namesave[100];


    sprintf(nameload1,"hall/hall000000.pgm");
    
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &vi0, &vi1, &vj0, &vj1);

    
    vuint8 ** m = vui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    vuint8 ** It = vui8matrix(nrl,nrh,ncl,nch);
    vuint8 ** It0 = vui8matrix(nrl,nrh,ncl,nch);
    uint_to_vuint(Itm1, It0, vi0, vi1, vj0, vj1);
    uint8** a = ui8matrix(nrl, nrh, ncl, nch);
    uint8** out = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout = ui8matrix(nrl, nrh, ncl, nch);
    
    
    int i;
   
    for(i=1;i<=NBIMAGES;i++){
        sprintf(nameload2,"hall/hall000%03d.pgm",i);
        
        a = LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
        
        uint_to_vuint(a, It, vi0, vi1, vj0, vj1);
        routine_FrameDifference_SSE2(It0,It,m,vi0, vi1,vj0,vj1,seuil);
        vuint_to_uint(out, m, vi0, vi1, vj0, vj1);


        fermeture3SSE(out,Etout,nrl,nrh,ncl,nch);

        
        sprintf(namesave,"testFD_SSEmorphoF/hall000%03d.pgm",i);
        SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
        dup_vui8matrix(It, vi0, vi1, vj0, vj1, It0);
    }
}

void test_routineFD_SSEmorpho3xOuvFerm(int seuil){
    long nrl, nrh, ncl, nch;
    int vi0, vi1, vj0, vj1;
    char nameload1[100], nameload2[100], namesave[100];


    sprintf(nameload1,"hall/hall000000.pgm");
    
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &vi0, &vi1, &vj0, &vj1);

    
    vuint8 ** m = vui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    vuint8 ** It = vui8matrix(nrl,nrh,ncl,nch);
    vuint8 ** It0 = vui8matrix(nrl,nrh,ncl,nch);
    uint_to_vuint(Itm1, It0, vi0, vi1, vj0, vj1);
    uint8** a = ui8matrix(nrl, nrh, ncl, nch);
    uint8** out = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout2 = ui8matrix(nrl, nrh, ncl, nch);
    
    
    int i;
   
    for(i=1;i<=NBIMAGES;i++){
        sprintf(nameload2,"hall/hall000%03d.pgm",i);
        
        a = LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
        
        uint_to_vuint(a, It, vi0, vi1, vj0, vj1);
        routine_FrameDifference_SSE2(It0,It,m,vi0, vi1,vj0,vj1,seuil);
        vuint_to_uint(out, m, vi0, vi1, vj0, vj1);


        ouverture3SSE(out,Etout,nrl,nrh,ncl,nch);
        fermeture3SSE(Etout,Etout2,nrl,nrh,ncl,nch);
        
        sprintf(namesave,"testFD_SSEmorphoOF/hall000%03d.pgm",i);
        SavePGM_ui8matrix(Etout2,nrl,nrh,ncl,nch,namesave);
        dup_vui8matrix(It, vi0, vi1, vj0, vj1, It0);
    }
}

void test_routineFD_SSEmorpho3xFermOuv(int seuil){
    long nrl, nrh, ncl, nch;
    int vi0, vi1, vj0, vj1;
    char nameload1[100], nameload2[100], namesave[100];


    sprintf(nameload1,"hall/hall000000.pgm");
    
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &vi0, &vi1, &vj0, &vj1);

    
    vuint8 ** m = vui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    vuint8 ** It = vui8matrix(nrl,nrh,ncl,nch);
    vuint8 ** It0 = vui8matrix(nrl,nrh,ncl,nch);
    uint_to_vuint(Itm1, It0, vi0, vi1, vj0, vj1);
    uint8** a = ui8matrix(nrl, nrh, ncl, nch);
    uint8** out = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout2 = ui8matrix(nrl, nrh, ncl, nch);
    
    
    int i;
   
    for(i=1;i<=NBIMAGES;i++){
        sprintf(nameload2,"hall/hall000%03d.pgm",i);
        
        a = LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
        
        uint_to_vuint(a, It, vi0, vi1, vj0, vj1);
        routine_FrameDifference_SSE2(It0,It,m,vi0, vi1,vj0,vj1,seuil);
        vuint_to_uint(out, m, vi0, vi1, vj0, vj1);


        fermeture3SSE(out,Etout,nrl,nrh,ncl,nch);
        ouverture3SSE(Etout,Etout2,nrl,nrh,ncl,nch);
        
        sprintf(namesave,"testFD_SSEmorphoFO/hall000%03d.pgm",i);
        SavePGM_ui8matrix(Etout2,nrl,nrh,ncl,nch,namesave);
        dup_vui8matrix(It, vi0, vi1, vj0, vj1, It0);
    }
}

//////////////////////////////////////////////
//     TEST SIGMA DELTA SSE O,OF,F,FO       //
//////////////////////////////////////////////

void test_routineSD_SSEmorpho3xOuv(){
    int n1, n2, n3, n4;
    long nrl, nrh, ncl, nch;
    char nameload[100];     //"hall/hall000..";
    char namesave[100];     //"testSD/SD...";
    int i;
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    sprintf(nameload,"hall/hall000000.pgm");
    
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    vuint8 ** Itm1v = vui8matrix(n1,n2,n3,n4);
    uint_to_vuint(Itm1, Itm1v, n1, n2, n3, n4);


   
    uint8** a = ui8matrix(nrl, nrh, ncl, nch);
    vuint8 ** I = vui8matrix(n1,n2,n3,n4);
    vuint8 ** V  = vui8matrix(n1,n2,n3,n4);
    vuint8 ** Vtm1 = vui8matrix(n1,n2,n3,n4);
    vuint8** M  = vui8matrix(n1,n2,n3,n4);
     vuint8 ** Mtm1 = vui8matrix(n1,n2,n3,n4);
    vuint8 ** Et = vui8matrix(n1-BORD,n2+BORD,n3-BORD,n4+BORD);
    //uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);
    
    uint8** out = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout2 = ui8matrix(nrl, nrh, ncl, nch);
   
    vuint8 ** Iv = vui8matrix(n1,n2,n3,n4);
    uint8 ** res = ui8matrix(nrl, nrh, ncl, nch);
    
    SigmaDelta_step0_SSE2 (Vtm1, Mtm1, Itm1v, n1, n2, n3, n4);

    
    
    for(i=1;i<=NBIMAGES;i++){
        sprintf(nameload,"hall/hall000%03d.pgm",i);
        a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
        //Conversion
        uint_to_vuint(a, Iv, n1, n2, n3, n4);
        
        SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4);
        
        sprintf(namesave,"testSD_SSEmorphoO/hall000%03d.pgm",i);
        
        vuint_to_uint(out, Et, n1, n2, n3, n4);

        ouverture3SSE(out,Etout,nrl,nrh,ncl,nch);

        
        SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
        //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
        
        //copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
        //copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
        //copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
        
        dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
        dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
        dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);
        
        
        
    }
    
 /*   free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
    free_ui8matrix(I,nrl,nrh,ncl,nch);
    free_ui8matrix(V,nrl,nrh,ncl,nch);
    free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
    free_ui8matrix(M,nrl,nrh,ncl,nch);
    free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
    free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);*/
    //free_ui8matrix(Ot,nrl,nrh,ncl,nch);
}


void test_routineSD_SSEmorpho3xFerm(){
    int n1, n2, n3, n4;
    long nrl, nrh, ncl, nch;
    char nameload[100];     //"hall/hall000..";
    char namesave[100];     //"testSD/SD...";
    int i;
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    sprintf(nameload,"hall/hall000000.pgm");
    
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    vuint8 ** Itm1v = vui8matrix(n1,n2,n3,n4);
    uint_to_vuint(Itm1, Itm1v, n1, n2, n3, n4);


   
    uint8** a = ui8matrix(nrl, nrh, ncl, nch);
    vuint8 ** I = vui8matrix(n1,n2,n3,n4);
    vuint8 ** V  = vui8matrix(n1,n2,n3,n4);
    vuint8 ** Vtm1 = vui8matrix(n1,n2,n3,n4);
    vuint8** M  = vui8matrix(n1,n2,n3,n4);
     vuint8 ** Mtm1 = vui8matrix(n1,n2,n3,n4);
    vuint8 ** Et = vui8matrix(n1-BORD,n2+BORD,n3-BORD,n4+BORD);
    //uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);
    
    uint8** out = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout2 = ui8matrix(nrl, nrh, ncl, nch);
   
    vuint8 ** Iv = vui8matrix(n1,n2,n3,n4);
    uint8 ** res = ui8matrix(nrl, nrh, ncl, nch);
    
    SigmaDelta_step0_SSE2 (Vtm1, Mtm1, Itm1v, n1, n2, n3, n4);

    
    
    for(i=1;i<=NBIMAGES;i++){
        sprintf(nameload,"hall/hall000%03d.pgm",i);
        a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
        //Conversion
        uint_to_vuint(a, Iv, n1, n2, n3, n4);
        
        SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4);
        
        sprintf(namesave,"testSD_SSEmorphoF/hall000%03d.pgm",i);
        
        vuint_to_uint(out, Et, n1, n2, n3, n4);

        fermeture3SSE(out,Etout,nrl,nrh,ncl,nch);

        
        SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
        //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
        
        //copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
        //copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
        //copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
        
        dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
        dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
        dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);
        
        
        
    }
    
 /*   free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
    free_ui8matrix(I,nrl,nrh,ncl,nch);
    free_ui8matrix(V,nrl,nrh,ncl,nch);
    free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
    free_ui8matrix(M,nrl,nrh,ncl,nch);
    free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
    free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);*/
    //free_ui8matrix(Ot,nrl,nrh,ncl,nch);
}

void test_routineSD_SSEmorpho3xOuvFerm(){
    int n1, n2, n3, n4;
    long nrl, nrh, ncl, nch;
    char nameload[100];     //"hall/hall000..";
    char namesave[100];     //"testSD/SD...";
    int i;
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    sprintf(nameload,"hall/hall000000.pgm");
    
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    vuint8 ** Itm1v = vui8matrix(n1,n2,n3,n4);
    uint_to_vuint(Itm1, Itm1v, n1, n2, n3, n4);


   
    uint8** a = ui8matrix(nrl, nrh, ncl, nch);
    vuint8 ** I = vui8matrix(n1,n2,n3,n4);
    vuint8 ** V  = vui8matrix(n1,n2,n3,n4);
    vuint8 ** Vtm1 = vui8matrix(n1,n2,n3,n4);
    vuint8** M  = vui8matrix(n1,n2,n3,n4);
     vuint8 ** Mtm1 = vui8matrix(n1,n2,n3,n4);
    vuint8 ** Et = vui8matrix(n1-BORD,n2+BORD,n3-BORD,n4+BORD);
    //uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);
    
    uint8** out = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout2 = ui8matrix(nrl, nrh, ncl, nch);
   
    vuint8 ** Iv = vui8matrix(n1,n2,n3,n4);
    uint8 ** res = ui8matrix(nrl, nrh, ncl, nch);
    
    SigmaDelta_step0_SSE2 (Vtm1, Mtm1, Itm1v, n1, n2, n3, n4);

    
    
    for(i=1;i<=NBIMAGES;i++){
        sprintf(nameload,"hall/hall000%03d.pgm",i);
        a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
        //Conversion
        uint_to_vuint(a, Iv, n1, n2, n3, n4);
        
        SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4);
        
        sprintf(namesave,"testSD_SSEmorphoOF/hall000%03d.pgm",i);
        
        vuint_to_uint(out, Et, n1, n2, n3, n4);

        ouverture3SSE(out,Etout,nrl,nrh,ncl,nch);
        fermeture3SSE(Etout,Etout2,nrl,nrh,ncl,nch);

        
        SavePGM_ui8matrix(Etout2,nrl,nrh,ncl,nch,namesave);
        //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
        
        //copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
        //copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
        //copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
        
        dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
        dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
        dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);
        
        
        
    }
    
 /*   free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
    free_ui8matrix(I,nrl,nrh,ncl,nch);
    free_ui8matrix(V,nrl,nrh,ncl,nch);
    free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
    free_ui8matrix(M,nrl,nrh,ncl,nch);
    free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
    free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);*/
    //free_ui8matrix(Ot,nrl,nrh,ncl,nch);
}

void test_routineSD_SSEmorpho3xFermOuv(){
    int n1, n2, n3, n4;
    long nrl, nrh, ncl, nch;
    char nameload[100];     //"hall/hall000..";
    char namesave[100];     //"testSD/SD...";
    int i;
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    sprintf(nameload,"hall/hall000000.pgm");
    
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    vuint8 ** Itm1v = vui8matrix(n1,n2,n3,n4);
    uint_to_vuint(Itm1, Itm1v, n1, n2, n3, n4);


   
    uint8** a = ui8matrix(nrl, nrh, ncl, nch);
    vuint8 ** I = vui8matrix(n1,n2,n3,n4);
    vuint8 ** V  = vui8matrix(n1,n2,n3,n4);
    vuint8 ** Vtm1 = vui8matrix(n1,n2,n3,n4);
    vuint8** M  = vui8matrix(n1,n2,n3,n4);
     vuint8 ** Mtm1 = vui8matrix(n1,n2,n3,n4);
    vuint8 ** Et = vui8matrix(n1-BORD,n2+BORD,n3-BORD,n4+BORD);
    //uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);
    
    uint8** out = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout2 = ui8matrix(nrl, nrh, ncl, nch);
   
    vuint8 ** Iv = vui8matrix(n1,n2,n3,n4);
    uint8 ** res = ui8matrix(nrl, nrh, ncl, nch);
    
    SigmaDelta_step0_SSE2 (Vtm1, Mtm1, Itm1v, n1, n2, n3, n4);

    
    
    for(i=1;i<=NBIMAGES;i++){
        sprintf(nameload,"hall/hall000%03d.pgm",i);
        a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
        //Conversion
        uint_to_vuint(a, Iv, n1, n2, n3, n4);
        
        SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4);
        
        sprintf(namesave,"testSD_SSEmorphoFO/hall000%03d.pgm",i);
        
        vuint_to_uint(out, Et, n1, n2, n3, n4);

        fermeture3SSE(out,Etout,nrl,nrh,ncl,nch);
        ouverture3SSE(Etout,Etout2,nrl,nrh,ncl,nch);

        
        SavePGM_ui8matrix(Etout2,nrl,nrh,ncl,nch,namesave);
        //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
        
        //copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
        //copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
        //copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
        
        dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
        dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
        dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);
        
        
        
    }
    
 /*   free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
    free_ui8matrix(I,nrl,nrh,ncl,nch);
    free_ui8matrix(V,nrl,nrh,ncl,nch);
    free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
    free_ui8matrix(M,nrl,nrh,ncl,nch);
    free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
    free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);*/
    //free_ui8matrix(Ot,nrl,nrh,ncl,nch);
}

