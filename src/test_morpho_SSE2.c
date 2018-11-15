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
    
}

void test_routineSD_SSEmorpho3xFerm(){
    
}

void test_routineSD_SSEmorpho3xOuvFerm(){
    
}

void test_routineSD_SSEmorpho3xFermOuv(){
    
}

