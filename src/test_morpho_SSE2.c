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
#include "mymacro.h"

//#define SEUILFD 20
#define NBIMAGES 199
#define BORD 2


void test_EtapemorphoSSE(void){
    long nrl, nrh, ncl, nch;
    int n1, n2, n3, n4;

    //vuint8 tmp;
    
    uint8 **Et=LoadPGM_ui8matrix("smile.pgm",&nrl,&nrh,&ncl,&nch);

    int i,j;

    //s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);
    //vuint8 ** It = vui8matrix(n1,n2,n3,n4);
    uint8 **Etbord = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);

    //Pour la morpho binaire
	ulong64 **Etbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **Etoutbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **tmpbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);


    //uint_to_vuint(Et, It, n1, n2, n3, n4);
    for(i=nrl;i<=nrh;i++){
		for(j=ncl;j<=nch;j++){
			Etbord[i][j]=Et[i][j];
		}
	}


    erosion3SSE(Etbord,Etout,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,"testmorphoSSE/testeroSSE.pgm");

    dilatation3SSE(Etbord,Etout,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,"testmorphoSSE/testdilSSE.pgm");
    

    ouverture3SSE(Etbord,Etout,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,"testmorphoSSE/testouvSSE.pgm");
    
    fermeture3SSE(Etbord,Etout,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,"testmorphoSSE/testfermSSE.pgm");



    convCharToBin(Etbord,Etbin,nrl,nrh,ncl,nch);
   	erosion3SSE_bin(Etbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS));
   	convBinToChar(Etoutbin,Etout,nrl,nrh,ncl,nch);
   	SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,"testmorphoSSE/testerobin3.pgm");

   	convCharToBin(Etbord,Etbin,nrl,nrh,ncl,nch);
   	dilatation3SSE_bin(Etbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS));
   	convBinToChar(Etoutbin,Etout,nrl,nrh,ncl,nch);
   	SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,"testmorphoSSE/testdilbin3.pgm");

   	convCharToBin(Etbord,Etbin,nrl,nrh,ncl,nch);
    ouverture3SSE_bin(Etbin,tmpbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS));
    convBinToChar(Etoutbin,Etout,nrl,nrh,ncl,nch);
   	SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,"testmorphoSSE/testouvbin3.pgm");

   	convCharToBin(Etbord,Etbin,nrl,nrh,ncl,nch);
   	fermeture3SSE_bin(Etbin,tmpbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS));
   	convBinToChar(Etoutbin,Etout,nrl,nrh,ncl,nch);
   	SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,"testmorphoSSE/testfermbin3.pgm");


    free_ui8matrix(Etbord,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(Et,nrl,nrh,ncl,nch);
    
    free_long64matrix(Etbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(Etoutbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(tmpbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);

}


///////////////////////////////////////////////////
//     TEST FRAME DIFFERENCE SSE O,OF,F,FO       //
///////////////////////////////////////////////////

void test_routineFD_SSEmorpho3xOuv(int seuil){

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


    sprintf(nameload1,"car3/car_3000.pgm");
    
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    
    vuint8 ** m = vui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    vuint8 ** It = vui8matrix(nrl,nrh,ncl,nch);
    vuint8 ** It0 = vui8matrix(nrl,nrh,ncl,nch);
    uint_to_vuint(Itm1, It0, n1, n2, n3, n4);
    uint8** a = ui8matrix(nrl, nrh, ncl, nch);
    uint8** out = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout = ui8matrix(nrl, nrh, ncl, nch);
    
    
    int i;
   
    for(i=1;i<=NBIMAGES;i++){
        sprintf(nameload2,"car3/car_3%03d.pgm",i);
        
        a = LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
        
        uint_to_vuint(a, It, n1, n2, n3, n4);
        routine_FrameDifference_SSE2(It0,It,m,n1, n2,n3,n4,seuil);
        vuint_to_uint(out, m, n1, n2, n3, n4);

        CHRONO(ouverture3SSE(out,Etout,nrl,nrh,ncl,nch),cycles);
        totalcy += cycles;
        
        sprintf(namesave,"testFD_SSEmorphoO/car_3%03d.pgm",i);
        SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
        dup_vui8matrix(It, n1, n2, n3, n4, It0);
    }

    totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
    totalcy = totalcy / ((nch+1)*(nrh+1));

    BENCH(printf("Cycles FD_SSE, only OuvSSE = "));
    BENCH(printf(format,totalcy));

    free_vui8matrix(m,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_vui8matrix(It,nrl,nrh,ncl,nch);
    free_vui8matrix(It0,nrl,nrh,ncl,nch);
    free_ui8matrix(out,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(a,nrl,nrh,ncl,nch);
    free_ui8matrix(Etout,nrl,nrh,ncl,nch);
    free_ui8matrix(Itm1,nrl,nrh,ncl,nch);

}

void test_routineFD_SSEmorpho3xFerm(int seuil){

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


    sprintf(nameload1,"car3/car_3000.pgm");
    
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    
    vuint8 ** m = vui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    vuint8 ** It = vui8matrix(nrl,nrh,ncl,nch);
    vuint8 ** It0 = vui8matrix(nrl,nrh,ncl,nch);
    uint_to_vuint(Itm1, It0, n1, n2, n3, n4);
    uint8** a = ui8matrix(nrl, nrh, ncl, nch);
    uint8** out = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout = ui8matrix(nrl, nrh, ncl, nch);
    
    
    int i;
   
    for(i=1;i<=NBIMAGES;i++){
        sprintf(nameload2,"car3/car_3%03d.pgm",i);
        
        a = LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
        
        uint_to_vuint(a, It, n1, n2, n3, n4);
        routine_FrameDifference_SSE2(It0,It,m,n1, n2,n3,n4,seuil);
        vuint_to_uint(out, m, n1, n2, n3, n4);

        CHRONO(fermeture3SSE(out,Etout,nrl,nrh,ncl,nch),cycles);
        totalcy += cycles;
        
        sprintf(namesave,"testFD_SSEmorphoF/car_3%03d.pgm",i);
        SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
        dup_vui8matrix(It, n1, n2, n3, n4, It0);
    }

    totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
    totalcy = totalcy / ((nch+1)*(nrh+1));

    BENCH(printf("Cycles FD_SSE, only FermSSE = "));
    BENCH(printf(format,totalcy));

    free_vui8matrix(m,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_vui8matrix(It,nrl,nrh,ncl,nch);
    free_vui8matrix(It0,nrl,nrh,ncl,nch);
    free_ui8matrix(out,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(a,nrl,nrh,ncl,nch);
    free_ui8matrix(Etout,nrl,nrh,ncl,nch);
    free_ui8matrix(Itm1,nrl,nrh,ncl,nch);

}

void test_routineFD_SSEmorpho3xOuvFerm(int seuil){

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


    sprintf(nameload1,"car3/car_3000.pgm");
    
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    
    vuint8 ** m = vui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    vuint8 ** It = vui8matrix(nrl,nrh,ncl,nch);
    vuint8 ** It0 = vui8matrix(nrl,nrh,ncl,nch);
    uint_to_vuint(Itm1, It0, n1, n2, n3, n4);
    uint8** a = ui8matrix(nrl, nrh, ncl, nch);
    uint8** out = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout2 = ui8matrix(nrl, nrh, ncl, nch);
    
    
    int i;
   
    for(i=1;i<=NBIMAGES;i++){
        sprintf(nameload2,"car3/car_3%03d.pgm",i);
        
        a = LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
        
        uint_to_vuint(a, It, n1, n2, n3, n4);
        routine_FrameDifference_SSE2(It0,It,m,n1, n2,n3,n4,seuil);
        vuint_to_uint(out, m, n1, n2, n3, n4);


        CHRONO(ouverture3SSE(out,Etout,nrl,nrh,ncl,nch),cycles);
        totalcy += cycles;
        CHRONO(fermeture3SSE(Etout,Etout2,nrl,nrh,ncl,nch),cycles);
        totalcy += cycles;
        
        sprintf(namesave,"testFD_SSEmorphoOF/car_3%03d.pgm",i);
        SavePGM_ui8matrix(Etout2,nrl,nrh,ncl,nch,namesave);
        dup_vui8matrix(It, n1, n2, n3, n4, It0);
    }

    totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
    totalcy = totalcy / ((nch+1)*(nrh+1));

    BENCH(printf("Cycles FD_SSE, only OuvFermSSE = "));
    BENCH(printf(format,totalcy));

    free_vui8matrix(m,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_vui8matrix(It,nrl,nrh,ncl,nch);
    free_vui8matrix(It0,nrl,nrh,ncl,nch);
    free_ui8matrix(out,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(a,nrl,nrh,ncl,nch);
    free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(Etout2,nrl,nrh,ncl,nch);
    free_ui8matrix(Itm1,nrl,nrh,ncl,nch);

}

void test_routineFD_SSEmorpho3xFermOuv(int seuil){

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


    sprintf(nameload1,"car3/car_3000.pgm");
    
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    
    vuint8 ** m = vui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    vuint8 ** It = vui8matrix(nrl,nrh,ncl,nch);
    vuint8 ** It0 = vui8matrix(nrl,nrh,ncl,nch);
    uint_to_vuint(Itm1, It0, n1, n2, n3, n4);
    uint8** a = ui8matrix(nrl, nrh, ncl, nch);
    uint8** out = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8** Etout2 = ui8matrix(nrl, nrh, ncl, nch);
    
    
    int i;
   
    for(i=1;i<=NBIMAGES;i++){
        sprintf(nameload2,"car3/car_3%03d.pgm",i);
        
        a = LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
        
        uint_to_vuint(a, It, n1, n2, n3, n4);
        routine_FrameDifference_SSE2(It0,It,m,n1, n2,n3,n4,seuil);
        vuint_to_uint(out, m, n1, n2, n3, n4);


        CHRONO(fermeture3SSE(out,Etout,nrl,nrh,ncl,nch),cycles);
        totalcy += cycles;
        CHRONO(ouverture3SSE(Etout,Etout2,nrl,nrh,ncl,nch),cycles);
        totalcy += cycles;
        
        sprintf(namesave,"testFD_SSEmorphoFO/car_3%03d.pgm",i);
        SavePGM_ui8matrix(Etout2,nrl,nrh,ncl,nch,namesave);
        dup_vui8matrix(It, n1, n2, n3, n4, It0);
    }
    totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
    totalcy = totalcy / ((nch+1)*(nrh+1));
    //printf("%ld",((nch+1)*(nrh+1)));

    BENCH(printf("Cycles FD_SSE, only FermOuvSSE = "));
    BENCH(printf(format,totalcy));

    free_vui8matrix(m,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_vui8matrix(It,nrl,nrh,ncl,nch);
    free_vui8matrix(It0,nrl,nrh,ncl,nch);
    free_ui8matrix(out,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(a,nrl,nrh,ncl,nch);
    free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(Etout2,nrl,nrh,ncl,nch);
    free_ui8matrix(Itm1,nrl,nrh,ncl,nch);

}


//////////////////////////////////////////////
//     TEST SIGMA DELTA SSE O,OF,F,FO       //
//////////////////////////////////////////////

void test_routineSD_SSEmorpho3xOuv(){

    //Cycle par point//
    double cycles, totalcy=0;
    int iter, niter=2;
    int run, nrun = 5;
    double t0,t1,dt,tmin,t;

    char *format = "%6.2f\n";
    ///////////////////

    int n1, n2, n3, n4;
    long nrl, nrh, ncl, nch;
    char nameload[100];     //"car3/car_3..";
    char namesave[100];     //"testSD/SD...";
    int i;
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    sprintf(nameload,"car3/car_3000.pgm");
    
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
        sprintf(nameload,"car3/car_3%03d.pgm",i);
        a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
        //Conversion
        uint_to_vuint(a, Iv, n1, n2, n3, n4);
        
        SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4);
        
        sprintf(namesave,"testSD_SSEmorphoO/car_3%03d.pgm",i);
        
        vuint_to_uint(out, Et, n1, n2, n3, n4);

        CHRONO(ouverture3SSE(out,Etout,nrl,nrh,ncl,nch),cycles);
        totalcy += cycles;
        
        SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
        //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
        
        //copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
        //copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
        //copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
        
        dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
        dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
        dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);
        
        
        
    }
    
    totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
    totalcy = totalcy / ((nch+1)*(nrh+1));

    BENCH(printf("Cycles SD_SSE, only OuvSSE = "));
    BENCH(printf(format,totalcy));

    free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
    free_ui8matrix(a,nrl,nrh,ncl,nch);
    free_vui8matrix(I,n1,n2,n3,n4);
    free_vui8matrix(V,n1,n2,n3,n4);
    free_vui8matrix(Vtm1,n1,n2,n3,n4);
    free_vui8matrix(M,n1,n2,n3,n4);
    free_vui8matrix(Mtm1,n1,n2,n3,n4);
    free_vui8matrix(Et,n1-BORD,n2+BORD,n3-BORD,n4+BORD);
    free_ui8matrix(out,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(Etout2,nrl,nrh,ncl,nch);
    free_vui8matrix(Iv,n1,n2,n3,n4);
    free_vui8matrix(Itm1v,n1,n2,n3,n4);
    free_ui8matrix(res,nrl,nrh,ncl,nch);
    
}


void test_routineSD_SSEmorpho3xFerm(){

    //Cycle par point//
    double cycles, totalcy=0;
    int iter, niter=2;
    int run, nrun = 5;
    double t0,t1,dt,tmin,t;

    char *format = "%6.2f\n";
    ///////////////////

    int n1, n2, n3, n4;
    long nrl, nrh, ncl, nch;
    char nameload[100];     //"car3/car_3..";
    char namesave[100];     //"testSD/SD...";
    int i;
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    sprintf(nameload,"car3/car_3000.pgm");
    
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
        sprintf(nameload,"car3/car_3%03d.pgm",i);
        a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
        //Conversion
        uint_to_vuint(a, Iv, n1, n2, n3, n4);
        
        SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4);
        
        sprintf(namesave,"testSD_SSEmorphoF/car_3%03d.pgm",i);
        
        vuint_to_uint(out, Et, n1, n2, n3, n4);

        CHRONO(fermeture3SSE(out,Etout,nrl,nrh,ncl,nch),cycles);
        totalcy += cycles;

        
        SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
        //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
        
        //copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
        //copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
        //copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
        
        dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
        dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
        dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);
        
        
        
    }
    
    totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
    totalcy = totalcy / ((nch+1)*(nrh+1));

    BENCH(printf("Cycles SD_SSE, only FermSSE = "));
    BENCH(printf(format,totalcy));

    free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
    free_ui8matrix(a,nrl,nrh,ncl,nch);
    free_vui8matrix(I,n1,n2,n3,n4);
    free_vui8matrix(V,n1,n2,n3,n4);
    free_vui8matrix(Vtm1,n1,n2,n3,n4);
    free_vui8matrix(M,n1,n2,n3,n4);
    free_vui8matrix(Mtm1,n1,n2,n3,n4);
    free_vui8matrix(Et,n1-BORD,n2+BORD,n3-BORD,n4+BORD);
    free_ui8matrix(out,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(Etout2,nrl,nrh,ncl,nch);
    free_vui8matrix(Iv,n1,n2,n3,n4);
    free_vui8matrix(Itm1v,n1,n2,n3,n4);
    free_ui8matrix(res,nrl,nrh,ncl,nch);

}

void test_routineSD_SSEmorpho3xOuvFerm(){

    //Cycle par point//
    double cycles, totalcy=0;
    int iter, niter=2;
    int run, nrun = 5;
    double t0,t1,dt,tmin,t;

    char *format = "%6.2f\n";
    ///////////////////

    int n1, n2, n3, n4;
    long nrl, nrh, ncl, nch;
    char nameload[100];     //"car3/car_3..";
    char namesave[100];     //"testSD/SD...";
    int i;
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    sprintf(nameload,"car3/car_3000.pgm");
    
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
    uint8** Etout2 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
   
    vuint8 ** Iv = vui8matrix(n1,n2,n3,n4);
    uint8 ** res = ui8matrix(nrl, nrh, ncl, nch);
    
    SigmaDelta_step0_SSE2 (Vtm1, Mtm1, Itm1v, n1, n2, n3, n4);

    
    
    for(i=1;i<=NBIMAGES;i++){
        sprintf(nameload,"car3/car_3%03d.pgm",i);
        a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
        //Conversion
        uint_to_vuint(a, Iv, n1, n2, n3, n4);
        
        SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4);
        
        sprintf(namesave,"testSD_SSEmorphoOF/car_3%03d.pgm",i);
        
        vuint_to_uint(out, Et, n1, n2, n3, n4);

        CHRONO(ouverture3SSE(out,Etout,nrl,nrh,ncl,nch),cycles);
        totalcy += cycles;
        CHRONO(fermeture3SSE(Etout,Etout2,nrl,nrh,ncl,nch),cycles);
        totalcy += cycles;

        
        SavePGM_ui8matrix(Etout2,nrl,nrh,ncl,nch,namesave);
        //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
        
        //copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
        //copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
        //copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
        
        dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
        dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
        dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);
        
        
        
    }

    totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
    totalcy = totalcy / ((nch+1)*(nrh+1));

    BENCH(printf("Cycles SD_SSE, only OuvFermSSE = "));
    BENCH(printf(format,totalcy));
    
    free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
    free_ui8matrix(a,nrl,nrh,ncl,nch);
    free_vui8matrix(I,n1,n2,n3,n4);
    free_vui8matrix(V,n1,n2,n3,n4);
    free_vui8matrix(Vtm1,n1,n2,n3,n4);
    free_vui8matrix(M,n1,n2,n3,n4);
    free_vui8matrix(Mtm1,n1,n2,n3,n4);
    free_vui8matrix(Et,n1-BORD,n2+BORD,n3-BORD,n4+BORD);
    free_ui8matrix(out,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(Etout2,nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    free_vui8matrix(Iv,n1,n2,n3,n4);
    free_vui8matrix(Itm1v,n1,n2,n3,n4);
    free_ui8matrix(res,nrl,nrh,ncl,nch);
}

void test_routineSD_SSEmorpho3xFermOuv(){

    //Cycle par point//
    double cycles, totalcy=0;
    int iter, niter=2;
    int run, nrun = 5;
    double t0,t1,dt,tmin,t;

    char *format = "%6.2f\n";
    ///////////////////

    int n1, n2, n3, n4;
    long nrl, nrh, ncl, nch;
    char nameload[100];     //"car3/car_3..";
    char namesave[100];     //"testSD/SD...";
    int i;
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

    sprintf(nameload,"car3/car_3000.pgm");
    
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
    uint8** Etout2 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
   
    vuint8 ** Iv = vui8matrix(n1,n2,n3,n4);
    uint8 ** res = ui8matrix(nrl, nrh, ncl, nch);
    
    SigmaDelta_step0_SSE2 (Vtm1, Mtm1, Itm1v, n1, n2, n3, n4);

    
    
    for(i=1;i<=NBIMAGES;i++){
        sprintf(nameload,"car3/car_3%03d.pgm",i);
        a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
        //Conversion
        uint_to_vuint(a, Iv, n1, n2, n3, n4);
        
        SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4);
        
        sprintf(namesave,"testSD_SSEmorphoFO/car_3%03d.pgm",i);
        
        vuint_to_uint(out, Et, n1, n2, n3, n4);

        CHRONO(fermeture3SSE(out,Etout,nrl,nrh,ncl,nch),cycles);
        totalcy += cycles;
        CHRONO(ouverture3SSE(Etout,Etout2,nrl,nrh,ncl,nch),cycles);
        totalcy += cycles;

        SavePGM_ui8matrix(Etout2,nrl,nrh,ncl,nch,namesave);
        //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
        
        //copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
        //copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
        //copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
        
        dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
        dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
        dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);
        
        
        
    }
    
    totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
    totalcy = totalcy / ((nch+1)*(nrh+1));

    BENCH(printf("Cycles SD_SSE, only FermOuvSSE = "));
    BENCH(printf(format,totalcy));

    free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
    free_ui8matrix(a,nrl,nrh,ncl,nch);
    free_vui8matrix(I,n1,n2,n3,n4);
    free_vui8matrix(V,n1,n2,n3,n4);
    free_vui8matrix(Vtm1,n1,n2,n3,n4);
    free_vui8matrix(M,n1,n2,n3,n4);
    free_vui8matrix(Mtm1,n1,n2,n3,n4);
    free_vui8matrix(Et,n1-BORD,n2+BORD,n3-BORD,n4+BORD);
    free_ui8matrix(out,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(Etout2,nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    free_vui8matrix(Iv,n1,n2,n3,n4);
    free_vui8matrix(Itm1v,n1,n2,n3,n4);
    free_ui8matrix(res,nrl,nrh,ncl,nch);
}


void test_routineSD_SSEmorpho3xOuv_bin(){
	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

	int n1, n2, n3, n4;
	long nrl, nrh, ncl, nch;
	char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;
	s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

	sprintf(nameload,"car3/car_3000.pgm");

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

	ulong64 **Etbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **Etoutbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **tmpbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);

	SigmaDelta_step0_SSE2 (Vtm1, Mtm1, Itm1v, n1, n2, n3, n4);



	for(i=1;i<=NBIMAGES;i++){
		sprintf(nameload,"car3/car_3%03d.pgm",i);
		a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
		//Conversion
		uint_to_vuint(a, Iv, n1, n2, n3, n4);

		SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4);

		sprintf(namesave,"testSD_SSEmorphoO_bin/car_3%03d.pgm",i);

		vuint_to_uint(out, Et, n1, n2, n3, n4);

		convCharToBin(out,Etbin,nrl,nrh,ncl,nch);
		CHRONO(ouverture3SSE_bin(Etbin,tmpbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS)),cycles);
		convBinToChar(Etoutbin,Etout,nrl,nrh,ncl,nch);
		totalcy += cycles;

		SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
		//On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1

		//copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
		//copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
		//copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);

		dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
		dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
		dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);



	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD_SSE, only OuvSSE bin = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(a,nrl,nrh,ncl,nch);
	free_vui8matrix(I,n1,n2,n3,n4);
	free_vui8matrix(V,n1,n2,n3,n4);
	free_vui8matrix(Vtm1,n1,n2,n3,n4);
	free_vui8matrix(M,n1,n2,n3,n4);
	free_vui8matrix(Mtm1,n1,n2,n3,n4);
	free_vui8matrix(Et,n1-BORD,n2+BORD,n3-BORD,n4+BORD);
	free_ui8matrix(out,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout2,nrl,nrh,ncl,nch);
	free_vui8matrix(Iv,n1,n2,n3,n4);
	free_vui8matrix(Itm1v,n1,n2,n3,n4);
	free_ui8matrix(res,nrl,nrh,ncl,nch);

	free_long64matrix(Etbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(Etoutbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(tmpbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
}

void test_routineSD_SSEmorpho3xFerm_bin(){
	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

	int n1, n2, n3, n4;
	long nrl, nrh, ncl, nch;
	char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;
	s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

	sprintf(nameload,"car3/car_3000.pgm");

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

	ulong64 **Etbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **Etoutbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **tmpbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);

	SigmaDelta_step0_SSE2 (Vtm1, Mtm1, Itm1v, n1, n2, n3, n4);



	for(i=1;i<=NBIMAGES;i++){
		sprintf(nameload,"car3/car_3%03d.pgm",i);
		a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
		//Conversion
		uint_to_vuint(a, Iv, n1, n2, n3, n4);

		SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4);

		sprintf(namesave,"testSD_SSEmorphoF_bin/car_3%03d.pgm",i);

		vuint_to_uint(out, Et, n1, n2, n3, n4);

		convCharToBin(out,Etbin,nrl,nrh,ncl,nch);
		CHRONO(fermeture3SSE_bin(Etbin,tmpbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS)),cycles);
		convBinToChar(Etoutbin,Etout,nrl,nrh,ncl,nch);
		totalcy += cycles;

		SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
		//On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1

		//copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
		//copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
		//copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);

		dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
		dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
		dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);



	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD_SSE, only FermSSE bin = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(a,nrl,nrh,ncl,nch);
	free_vui8matrix(I,n1,n2,n3,n4);
	free_vui8matrix(V,n1,n2,n3,n4);
	free_vui8matrix(Vtm1,n1,n2,n3,n4);
	free_vui8matrix(M,n1,n2,n3,n4);
	free_vui8matrix(Mtm1,n1,n2,n3,n4);
	free_vui8matrix(Et,n1-BORD,n2+BORD,n3-BORD,n4+BORD);
	free_ui8matrix(out,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout2,nrl,nrh,ncl,nch);
	free_vui8matrix(Iv,n1,n2,n3,n4);
	free_vui8matrix(Itm1v,n1,n2,n3,n4);
	free_ui8matrix(res,nrl,nrh,ncl,nch);

	free_long64matrix(Etbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(Etoutbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(tmpbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
}

void test_routineSD_SSEmorpho3xOuvFerm_bin(){
	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

	int n1, n2, n3, n4;
	long nrl, nrh, ncl, nch;
	char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;
	s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

	sprintf(nameload,"car3/car_3000.pgm");

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

	ulong64 **Etbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **Etoutbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **tmpbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);

	SigmaDelta_step0_SSE2 (Vtm1, Mtm1, Itm1v, n1, n2, n3, n4);



	for(i=1;i<=NBIMAGES;i++){
		sprintf(nameload,"car3/car_3%03d.pgm",i);
		a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
		//Conversion
		uint_to_vuint(a, Iv, n1, n2, n3, n4);

		SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4);

		sprintf(namesave,"testSD_SSEmorphoOF_bin/car_3%03d.pgm",i);

		vuint_to_uint(out, Et, n1, n2, n3, n4);

		convCharToBin(out,Etbin,nrl,nrh,ncl,nch);
		CHRONO(ouverture3SSE_bin(Etbin,tmpbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS)),cycles);
		totalcy += cycles;
		CHRONO(fermeture3SSE_bin(Etoutbin,tmpbin,Etbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS)),cycles);
		convBinToChar(Etbin,Etout,nrl,nrh,ncl,nch);
		totalcy += cycles;

		SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
		//On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1

		//copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
		//copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
		//copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);

		dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
		dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
		dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);



	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD_SSE, only OuvFermSSE bin = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(a,nrl,nrh,ncl,nch);
	free_vui8matrix(I,n1,n2,n3,n4);
	free_vui8matrix(V,n1,n2,n3,n4);
	free_vui8matrix(Vtm1,n1,n2,n3,n4);
	free_vui8matrix(M,n1,n2,n3,n4);
	free_vui8matrix(Mtm1,n1,n2,n3,n4);
	free_vui8matrix(Et,n1-BORD,n2+BORD,n3-BORD,n4+BORD);
	free_ui8matrix(out,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout2,nrl,nrh,ncl,nch);
	free_vui8matrix(Iv,n1,n2,n3,n4);
	free_vui8matrix(Itm1v,n1,n2,n3,n4);
	free_ui8matrix(res,nrl,nrh,ncl,nch);

	free_long64matrix(Etbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(Etoutbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(tmpbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
}

void test_routineSD_SSEmorpho3xFermOuv_bin(){
	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

	int n1, n2, n3, n4;
	long nrl, nrh, ncl, nch;
	char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;
	s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

	sprintf(nameload,"car3/car_3000.pgm");

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

	ulong64 **Etbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **Etoutbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **tmpbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);

	SigmaDelta_step0_SSE2 (Vtm1, Mtm1, Itm1v, n1, n2, n3, n4);



	for(i=1;i<=NBIMAGES;i++){
		sprintf(nameload,"car3/car_3%03d.pgm",i);
		a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
		//Conversion
		uint_to_vuint(a, Iv, n1, n2, n3, n4);

		SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4);

		sprintf(namesave,"testSD_SSEmorphoFO_bin/car_3%03d.pgm",i);

		vuint_to_uint(out, Et, n1, n2, n3, n4);

		convCharToBin(out,Etbin,nrl,nrh,ncl,nch);
		CHRONO(fermeture3SSE_bin(Etbin,tmpbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS)),cycles);
		totalcy += cycles;
		CHRONO(ouverture3SSE_bin(Etoutbin,tmpbin,Etbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS)),cycles);
		convBinToChar(Etbin,Etout,nrl,nrh,ncl,nch);
		totalcy += cycles;

		SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
		//On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1

		//copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
		//copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
		//copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);

		dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
		dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
		dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);



	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD_SSE, only FermOuvSSE bin = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(a,nrl,nrh,ncl,nch);
	free_vui8matrix(I,n1,n2,n3,n4);
	free_vui8matrix(V,n1,n2,n3,n4);
	free_vui8matrix(Vtm1,n1,n2,n3,n4);
	free_vui8matrix(M,n1,n2,n3,n4);
	free_vui8matrix(Mtm1,n1,n2,n3,n4);
	free_vui8matrix(Et,n1-BORD,n2+BORD,n3-BORD,n4+BORD);
	free_ui8matrix(out,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout2,nrl,nrh,ncl,nch);
	free_vui8matrix(Iv,n1,n2,n3,n4);
	free_vui8matrix(Itm1v,n1,n2,n3,n4);
	free_ui8matrix(res,nrl,nrh,ncl,nch);

	free_long64matrix(Etbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(Etoutbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(tmpbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
}




void test_routineSD_SSEmorpho3xOuv_pipebin(){
	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

	int n1, n2, n3, n4;
	long nrl, nrh, ncl, nch;
	char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;
	s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

	sprintf(nameload,"car3/car_3000.pgm");

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

	ulong64 **Etbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **Etoutbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **tmpbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);

	SigmaDelta_step0_SSE2 (Vtm1, Mtm1, Itm1v, n1, n2, n3, n4);



	for(i=1;i<=NBIMAGES;i++){
		sprintf(nameload,"car3/car_3%03d.pgm",i);
		a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
		//Conversion
		uint_to_vuint(a, Iv, n1, n2, n3, n4);

		SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4);

		sprintf(namesave,"testSD_SSEmorphoO_pipebin/car_3%03d.pgm",i);

		vuint_to_uint(out, Et, n1, n2, n3, n4);

		convCharToBin(out,Etbin,nrl,nrh,ncl,nch);
		CHRONO(ouverture3SSE_pipe_bin(Etbin,tmpbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS)),cycles);
		convBinToChar(Etoutbin,Etout,nrl,nrh,ncl,nch);
		totalcy += cycles;

		SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
		//On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1

		//copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
		//copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
		//copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);

		dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
		dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
		dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);



	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD_SSE, only OuvSSE pipebin = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(a,nrl,nrh,ncl,nch);
	free_vui8matrix(I,n1,n2,n3,n4);
	free_vui8matrix(V,n1,n2,n3,n4);
	free_vui8matrix(Vtm1,n1,n2,n3,n4);
	free_vui8matrix(M,n1,n2,n3,n4);
	free_vui8matrix(Mtm1,n1,n2,n3,n4);
	free_vui8matrix(Et,n1-BORD,n2+BORD,n3-BORD,n4+BORD);
	free_ui8matrix(out,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout2,nrl,nrh,ncl,nch);
	free_vui8matrix(Iv,n1,n2,n3,n4);
	free_vui8matrix(Itm1v,n1,n2,n3,n4);
	free_ui8matrix(res,nrl,nrh,ncl,nch);

	free_long64matrix(Etbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(Etoutbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(tmpbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
}

void test_routineSD_SSEmorpho3xFerm_pipebin(){
	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

	int n1, n2, n3, n4;
	long nrl, nrh, ncl, nch;
	char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;
	s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

	sprintf(nameload,"car3/car_3000.pgm");

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

	ulong64 **Etbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **Etoutbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **tmpbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);

	SigmaDelta_step0_SSE2 (Vtm1, Mtm1, Itm1v, n1, n2, n3, n4);



	for(i=1;i<=NBIMAGES;i++){
		sprintf(nameload,"car3/car_3%03d.pgm",i);
		a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
		//Conversion
		uint_to_vuint(a, Iv, n1, n2, n3, n4);

		SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4);

		sprintf(namesave,"testSD_SSEmorphoF_pipebin/car_3%03d.pgm",i);

		vuint_to_uint(out, Et, n1, n2, n3, n4);

		convCharToBin(out,Etbin,nrl,nrh,ncl,nch);
		CHRONO(fermeture3SSE_pipe_bin(Etbin,tmpbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS)),cycles);
		convBinToChar(Etoutbin,Etout,nrl,nrh,ncl,nch);
		totalcy += cycles;

		SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
		//On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1

		//copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
		//copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
		//copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);

		dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
		dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
		dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);



	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD_SSE, only FermSSE pipebin = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(a,nrl,nrh,ncl,nch);
	free_vui8matrix(I,n1,n2,n3,n4);
	free_vui8matrix(V,n1,n2,n3,n4);
	free_vui8matrix(Vtm1,n1,n2,n3,n4);
	free_vui8matrix(M,n1,n2,n3,n4);
	free_vui8matrix(Mtm1,n1,n2,n3,n4);
	free_vui8matrix(Et,n1-BORD,n2+BORD,n3-BORD,n4+BORD);
	free_ui8matrix(out,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout2,nrl,nrh,ncl,nch);
	free_vui8matrix(Iv,n1,n2,n3,n4);
	free_vui8matrix(Itm1v,n1,n2,n3,n4);
	free_ui8matrix(res,nrl,nrh,ncl,nch);

	free_long64matrix(Etbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(Etoutbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(tmpbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
}

void test_routineSD_SSEmorpho3xOuvFerm_pipebin(){
	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

	int n1, n2, n3, n4;
	long nrl, nrh, ncl, nch;
	char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;
	s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

	sprintf(nameload,"car3/car_3000.pgm");

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

	ulong64 **Etbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **Etoutbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **tmpbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);

	SigmaDelta_step0_SSE2 (Vtm1, Mtm1, Itm1v, n1, n2, n3, n4);



	for(i=1;i<=NBIMAGES;i++){
		sprintf(nameload,"car3/car_3%03d.pgm",i);
		a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
		//Conversion
		uint_to_vuint(a, Iv, n1, n2, n3, n4);

		SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4);

		sprintf(namesave,"testSD_SSEmorphoOF_pipebin/car_3%03d.pgm",i);

		vuint_to_uint(out, Et, n1, n2, n3, n4);

		convCharToBin(out,Etbin,nrl,nrh,ncl,nch);
		CHRONO(ouverture3SSE_pipe_bin(Etbin,tmpbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS)),cycles);
		totalcy += cycles;
		CHRONO(fermeture3SSE_pipe_bin(Etoutbin,tmpbin,Etbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS)),cycles);
		convBinToChar(Etbin,Etout,nrl,nrh,ncl,nch);
		totalcy += cycles;

		SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
		//On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1

		//copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
		//copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
		//copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);

		dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
		dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
		dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);



	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD_SSE, only OuvFermSSE pipebin = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(a,nrl,nrh,ncl,nch);
	free_vui8matrix(I,n1,n2,n3,n4);
	free_vui8matrix(V,n1,n2,n3,n4);
	free_vui8matrix(Vtm1,n1,n2,n3,n4);
	free_vui8matrix(M,n1,n2,n3,n4);
	free_vui8matrix(Mtm1,n1,n2,n3,n4);
	free_vui8matrix(Et,n1-BORD,n2+BORD,n3-BORD,n4+BORD);
	free_ui8matrix(out,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout2,nrl,nrh,ncl,nch);
	free_vui8matrix(Iv,n1,n2,n3,n4);
	free_vui8matrix(Itm1v,n1,n2,n3,n4);
	free_ui8matrix(res,nrl,nrh,ncl,nch);

	free_long64matrix(Etbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(Etoutbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(tmpbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
}

void test_routineSD_SSEmorpho3xFermOuv_pipebin(){
	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

	int n1, n2, n3, n4;
	long nrl, nrh, ncl, nch;
	char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;
	s2v(nrl, nrh, ncl, nch, card_vuint8(), &n1, &n2, &n3, &n4);

	sprintf(nameload,"car3/car_3000.pgm");

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

	ulong64 **Etbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **Etoutbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **tmpbin=long64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);

	SigmaDelta_step0_SSE2 (Vtm1, Mtm1, Itm1v, n1, n2, n3, n4);



	for(i=1;i<=NBIMAGES;i++){
		sprintf(nameload,"car3/car_3%03d.pgm",i);
		a=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
		//Conversion
		uint_to_vuint(a, Iv, n1, n2, n3, n4);

		SigmaDelta_1step_SSE2(V, Vtm1, M, Mtm1, Iv, Et, n1, n2, n3, n4);

		sprintf(namesave,"testSD_SSEmorphoFO_pipebin/car_3%03d.pgm",i);

		vuint_to_uint(out, Et, n1, n2, n3, n4);

		convCharToBin(out,Etbin,nrl,nrh,ncl,nch);
		CHRONO(fermeture3SSE_pipe_bin(Etbin,tmpbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS)),cycles);
		totalcy += cycles;
		CHRONO(ouverture3SSE_pipe_bin(Etoutbin,tmpbin,Etbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS)),cycles);
		convBinToChar(Etbin,Etout,nrl,nrh,ncl,nch);
		totalcy += cycles;

		SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
		//On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1

		//copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
		//copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
		//copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);

		dup_vui8matrix(V, n1, n2, n3, n4, Vtm1);
		dup_vui8matrix(M, n1, n2, n3, n4, Mtm1);
		dup_vui8matrix(Iv, n1, n2, n3, n4, Itm1v);



	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD_SSE, only FermOuvSSE pipebin = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(a,nrl,nrh,ncl,nch);
	free_vui8matrix(I,n1,n2,n3,n4);
	free_vui8matrix(V,n1,n2,n3,n4);
	free_vui8matrix(Vtm1,n1,n2,n3,n4);
	free_vui8matrix(M,n1,n2,n3,n4);
	free_vui8matrix(Mtm1,n1,n2,n3,n4);
	free_vui8matrix(Et,n1-BORD,n2+BORD,n3-BORD,n4+BORD);
	free_ui8matrix(out,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout2,nrl,nrh,ncl,nch);
	free_vui8matrix(Iv,n1,n2,n3,n4);
	free_vui8matrix(Itm1v,n1,n2,n3,n4);
	free_ui8matrix(res,nrl,nrh,ncl,nch);

	free_long64matrix(Etbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(Etoutbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_long64matrix(tmpbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
}

//////////TESTS UNITAIRES///////////////////////////////////
void test_unitaire_erosion3SSE(){

    int m;
    vuint8 l_1 = init_vuint8_all(255,0, 255, 0, 0, 255, 0, 255,255, 0, 0, 255, 0, 255,255, 0);
    vuint8 l0 = init_vuint8_all(0,255, 255, 255, 0, 0, 0, 255,255, 255, 0, 255, 0, 0,255, 0);
    vuint8 l1 = init_vuint8_all(255,255, 255, 0, 0, 0, 0, 255,0, 255, 255, 255, 0, 0,255, 0);
    
    vuint8 tmp_1,tmp0,tmp1;
    vuint8 res;
    vuint8 xl, xr;
    vuint8 y;
    vuint8 ero;
    printf("/////////////////test_unitaire_erosion3SSE////////////////\n");
    
    display_vuint8(l_1," %.3d ","l_1\n");
    display_vuint8(l0," %.3d ","\nl0\n");
    display_vuint8(l1," %.3d ","\nl1\n");


    // Parcours de l'image

                res=_mm_and_si128(_mm_and_si128(l_1,l0),l1); //On met dans res l1 & l0 & l_1
                display_vuint8(res," %.3d ","\nres\n");

                tmp0=_mm_srli_si128(res,1);
                tmp1=_mm_slli_si128(res,1);
                res=_mm_and_si128(_mm_and_si128(tmp0,tmp1),res);

                res = _mm_srli_si128(res,1);
                display_vuint8(res," %.3d ","\nres\n");


}

void test_unitaire_dilatation3SSE(){

    int m;
    vuint8 l_1 = init_vuint8_all(255,0, 255, 0, 0, 255, 0, 255,255, 0, 0, 255, 0, 255,255, 0);
    vuint8 l0 = init_vuint8_all(0,255, 255, 255, 0, 0, 0, 255,255, 255, 0, 255, 0, 0,255, 0);
    vuint8 l1 = init_vuint8_all(255,255, 255, 0, 0, 0, 0, 255,0, 255, 255, 255, 0, 0,255, 0);
    
    vuint8 tmp_1,tmp0,tmp1;
    vuint8 res;
    vuint8 xl, xr;
    vuint8 y;
    vuint8 ero;
    printf("\n\n/////////////////test_unitaire_dilatation3SSE////////////////\n");

    display_vuint8(l_1," %.3d ","l_1\n");
    display_vuint8(l0," %.3d ","\nl0\n");
    display_vuint8(l1," %.3d ","\nl1\n");

    // Parcours de l'image

                res=_mm_or_si128(_mm_or_si128(l_1,l0),l1); //On met dans res l1 & l0 & l_1
                display_vuint8(res," %.3d ","\nres\n");

                tmp0=_mm_srli_si128(res,1);
                tmp1=_mm_slli_si128(res,1);
                res=_mm_and_si128(_mm_and_si128(tmp0,tmp1),res);

                res = _mm_srli_si128(res,1);
                display_vuint8(res," %.3d ","\nres\n");
    
    printf("\n");


}


void test_unitaire_erosion3SSE_bin(){

    printf("\n\n/////////////////test_unitaire_erosion3SSE_bin////////////////\n");
    
    vulong64 res;
    vulong64 binleft, binright;

    vulong64 zero = init_vulong64(0);
    //display_vulong64(zero," %d ","\ntest\n");
    vulong64 un = init_vulong64(1);

    //vulong64 bitfort = init_vulong64(~(1ULL<<(NBBITS-1)));
    vulong64 bitfort = _mm_slli_epi64(un,(NBBITS-1));
    bitfort = ~bitfort;


    // Parcours de l'image

    res = ~zero; // tout à 1

    vulong64 l_1 = init_vulong64_all(247,156);
    vulong64 l0 = init_vulong64_all(207,207);
    vulong64 l1 = init_vulong64_all(63,63);

    //l_1 : 1111 0111   1001 1100
    //l0  : 1100 1111   1100 1111
    //l1  : 0011 1111   0011 1111
    
    //and = 0000 0111   0000 1100
    //         7            12

    display_vulong64(l_1," %llu ","l_1\n");
    display_vulong64(l0," %llu","\nl0\n");
    display_vulong64(l1," %llu ","\nl1\n");
    
    res = _mm_and_si128(_mm_and_si128(_mm_and_si128(res,l_1),l0),l1);
    
    display_vulong64(res," %llu ","\nres\n");


    l_1 = init_vulong64_all(0,0);
    l0 = init_vulong64_all(0,0);
    l1 = init_vulong64_all(0,0);

    binleft = _mm_and_si128(_mm_and_si128(l_1,l0),l1);
    display_vulong64(binleft," %llu ","\nbinleft\n");

    //binleft = _mm_or_si128((res<< 1ULL),_mm_and_si128((binleft>>(NBBITS-1)),un));
    binleft = _mm_or_si128(_mm_slli_epi64(res,1),_mm_and_si128(_mm_slli_epi64(binleft,(NBBITS-1)),un));


    binright = _mm_srli_si128(res,8); //2ème ulong64 de res qui est le vecteur de droite
    display_vulong64(binright," %llu ","\nbinright\n");

    
    //binright = _mm_or_si128(_mm_and_si128((res>> 1ULL),bitfort),((_mm_and_si128(binright,un))<<(NBBITS-1)));
    binright = _mm_or_si128(_mm_and_si128(_mm_srli_epi64(res,1),bitfort),(_mm_slli_epi64((_mm_and_si128(binright,un)),(NBBITS-1))));

    res = _mm_and_si128(_mm_and_si128(res,binright),binleft);
    
    display_vulong64(res," %llu ","\nres\n");
    // A la fin l'érosion doit nous afficher :
    //    0000 0010  0000 0000
    //        2          0

}




void test_unitaire_dilatation3SSE_bin(){
    
    printf("\n\n/////////////////test_unitaire_dilatation3SSE_bin////////////////\n");

    
    vulong64 res;
    vulong64 binleft, binright;

    vulong64 zero = init_vulong64(0);
    vulong64 un = init_vulong64(1);

    //vulong64 bitfort = init_vulong64(~(1ULL<<(NBBITS-1)));
    vulong64 bitfort = (_mm_slli_epi64(un,(NBBITS-1)));
    bitfort = ~bitfort;

    // Parcours de l'image
    res = zero; // tout à 0

    vulong64 l_1 = init_vulong64_all(7,156);
    vulong64 l0 = init_vulong64_all(135,135);
    vulong64 l1 = init_vulong64_all(23,23);
    
    //l_1 : 0000 0111   1001 1100
    //l0  : 1000 0111   1000 0111
    //l1  : 0001 0111   0001 0111
    
    //or = 1001 0111   1001 1111
    //         151         159
    
    
    display_vulong64(l_1," %llu ","l_1\n");
    display_vulong64(l0," %llu","\nl0\n");
    display_vulong64(l1," %llu ","\nl1\n");

    res = _mm_or_si128(_mm_or_si128(_mm_or_si128(res,l_1),l0),l1);
    display_vulong64(res," %llu ","\nres\n");

    l_1 = init_vulong64_all(0,0);
    l0 = init_vulong64_all(0,0);
    l1 = init_vulong64_all(0,0);

    binleft = _mm_or_si128(_mm_or_si128(l_1,l0),l1);
    //binleft = _mm_or_si128((res<< 1ULL),_mm_and_si128((binleft>>(NBBITS-1)),un));
    binleft = _mm_or_si128(_mm_slli_epi64(res,1),_mm_and_si128(_mm_slli_epi64(binleft,(NBBITS-1)),un));
    display_vulong64(binleft," %llu ","\nbinleft\n");


    binright = _mm_srli_si128(res,8); //2ème ulong64 de res qui est le vecteur de droite
    //binright = _mm_or_si128(_mm_and_si128((res>> 1ULL),bitfort),((_mm_and_si128(binright,un))<<(NBBITS-1)));
    display_vulong64(binright," %llu ","\nbinright\n");

    binright = _mm_or_si128(_mm_and_si128(_mm_srli_epi64(res,1),bitfort),(_mm_slli_epi64((_mm_and_si128(binright,un)),(NBBITS-1))));
    // binright = 1<<63 = 9223372036854775883

    res = _mm_or_si128(_mm_or_si128(res,binright),binleft);
    display_vulong64(res," %llu ","\nres\n");
    printf("\n");

}
/*
1000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000010010111000000000000000000000000100111110
0000000000000000000000001001011100000000000000000000000100101111

1000000000000000000000011111111100000000000000000000000100111111
*/
