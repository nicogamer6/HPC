#include <stdio.h>
#include "nrutil.h"
#include "nrdef.h"
#include "string.h"
#include "test_morpho.h"
#include "morpho.h"
#include "mouvement.h"
#include "mymacro.h"

//#define SEUILFD 20
//#define NBBITS 64
#define NBIMAGES 199
#define BORD 2

void convCharToBin(uint8 ** Et, ulong64 ** Etout, long nrl, long nrh, long ncl, long nch){

	int i,j,k;
	for(i=nrl;i<=nrh;i++){
		k=0;
		Etout[i][k]=0;
		for(j=ncl;j<=nch;j++){
			if((j%NBBITS)==0 && j>0){	//j>0 pour commencer Kà0 sinon il va passer direct à 1
				k++;		// Changement de ulong64 à partir du 65ème pixel de la mat Char
				Etout[i][k]=0;
			}
			if(Et[i][j]==255){
				Etout[i][k] |= (1 << (j%NBBITS)) ;
			}
		}
	}
}

void convBinToChar(ulong64 ** Et, uint8 ** Etout, long nrl, long nrh, long ncl, long nch){

	int i,j,k;
	int bit;
	for(i=nrl;i<=nrh;i++){
		k=0;
		for(j=ncl;j<=nch;j++){
			Etout[i][j]=0; // On met tous les uint8 à 0
			bit=0;
			if((j%NBBITS)==0 && j>0){	//j>0 pour commencer Kà0 sinon il va passer direct à 1
				k++;		// Changement de ulong64 à partir du 65ème pixel
			}
			bit = (Et[i][k] >> (j%NBBITS)) & 1;
			if(bit==1){
				Etout[i][j] = 255;
			}
			if(bit==0){
				Etout[i][j] = 0;
			}
		}
	}

}

void test_Etapemorpho(void){
    long nrl, nrh, ncl, nch;
    
	uint8 **m=LoadPGM_ui8matrix("smile.pgm",&nrl,&nrh,&ncl,&nch);
    int i,j;
    uint8 **bord=ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    uint8 **tmp=ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    uint8 **m2=ui8matrix(nrl,nrh,ncl,nch);
    
    //Pour la morpho binaire
    ulong64 **Etbin=ui64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
    ulong64 **Etoutbin=ui64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
    ulong64 **tmpbin=ui64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);

    for(i=nrl;i<=nrh;i++){
        for(j=ncl;j<=nch;j++){
            bord[i][j]=m[i][j]; 
        }
    }
    
    //Affichage
    /*for(i=nrl;i<=nrh;i++){
        for(j=ncl;j<=nch;j++){
             printf("%0.3ld ",bord[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");*/
   
    erosion3(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testero3.pgm");

    erosion3_opti_lu_rr(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testeroopti3.pgm");

    erosion5(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testero5.pgm");

    erosion5_opti_lu_rr(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testeroopti5.pgm");
    

    dilatation3(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testdil3.pgm");
    
    dilatation3_opti_lu_rr(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testdilopti3.pgm");

    dilatation5(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testdil5.pgm");

    dilatation5_opti_lu_rr(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testdilopti5.pgm");


    ouverture3(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testouv3.pgm");
    
    ouverture3_opti(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testouvopti3.pgm");

    ouverture3_pipe(bord,tmp,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testouvpipe3.pgm");

    ouverture5(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testouv5.pgm");

    ouverture5_opti(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testouvopti5.pgm");


    fermeture3(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testferm3.pgm");

    fermeture3_opti(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testfermopti3.pgm");

    fermeture3_pipe(bord,tmp,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testfermpipe3.pgm");

    fermeture5(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testferm5.pgm");

    fermeture5_opti(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testfermopti5.pgm");


    convCharToBin(bord,Etbin,nrl,nrh,ncl,nch);
    //ouverture3_bin(Etbin,tmpbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS));
    convBinToChar(Etbin,m2,nrl,nrh,ncl,nch);
	SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testouvbin3.pgm");
	//printf("%lu\n",Etbin[1][0]);
	//printf("%lu",Etbin[2][0]);

	convCharToBin(bord,Etbin,nrl,nrh,ncl,nch);
	//fermeture3_bin(Etbin,tmpbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS));
	convBinToChar(Etbin,m2,nrl,nrh,ncl,nch);
	SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testmorpho/testfermbin3.pgm");



	for(i=nrl-BORD;i<=nrh+BORD;i++){
		for(j=(ncl/NBBITS-BORD);j<=(nch/NBBITS)+BORD;j++){
			printf("%lu  ",Etbin[i][0]);
		}
		printf("\n");
	}

    free_ui8matrix(bord,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    free_ui8matrix(m,nrl,nrh,ncl,nch);
    free_ui8matrix(m2,nrl,nrh,ncl,nch);

    free_ui64matrix(Etbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
    free_ui64matrix(Etoutbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
    free_ui64matrix(tmpbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
}


//////////////////////////////////////////////
//     TEST FRAME DIFFERENCE O,OF,F,FO	    //
//////////////////////////////////////////////

void test_routineFDmorpho3xOuv(int seuil){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];
    
    sprintf(nameload1,"car3/car_3000.pgm");
       
	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	
	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"car3/car_3%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
            
            routine_FrameDifference(Itm1,It,Et,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testFDmorphoO/car_3%03d.pgm",i);
            CHRONO(ouverture3(Et,Etout,nrl,nrh,ncl,nch),cycles);
            totalcy += cycles;
            //ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}   
	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles FD, only Ouv = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	
}

void test_routineFDmorpho3xOuv_opti(int seuil){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];

    sprintf(nameload1,"car3/car_3000.pgm");

	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);

	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"car3/car_3%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);

            routine_FrameDifference(Itm1,It,Et,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testOptiFD/car_3%03d.pgm",i);
            CHRONO(ouverture3_opti(Et,Etout,nrl,nrh,ncl,nch),cycles);
            totalcy += cycles;
            //ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}
	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles FD, only Ouv opti = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);

}

void test_routineFDmorpho3xOuv_pipe(int seuil){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];

    sprintf(nameload1,"car3/car_3000.pgm");

	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **tmp = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD); //Matrice tmp entre ero et dil pour le pipe

	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"car3/car_3%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);

            routine_FrameDifference(Itm1,It,Et,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testFDmorphoOpipe/car_3%03d.pgm",i);
            CHRONO(ouverture3_pipe(Et,tmp,Etout,nrl,nrh,ncl,nch),cycles);
            totalcy += cycles;
            //ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}
	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles FD, only Ouv pipe = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);

}

void test_routineFDmorpho3xOuv_bin(int seuil){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];

    sprintf(nameload1,"car3/car_3000.pgm");

	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);

	ulong64 **Etbin=ui64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **Etoutbin=ui64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **tmpbin=ui64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);

	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"car3/car_3%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);

            routine_FrameDifference(Itm1,It,Et,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testFDmorphoObin/car_3%03d.pgm",i);
            convCharToBin(Et,Etbin,nrl,nrh,ncl,nch);
            CHRONO(ouverture3_bin(Etbin,tmpbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS)),cycles);
            convBinToChar(Etoutbin,Etout,nrl,nrh,ncl,nch);
            totalcy += cycles;
            //ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}
	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles FD, only Ouv bin = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);

	free_ui64matrix(Etbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_ui64matrix(Etoutbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_ui64matrix(tmpbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);

}


void test_routineFDmorpho3xFerm(int seuil){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];
    
    sprintf(nameload1,"car3/car_3000.pgm");
       
	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	
	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"car3/car_3%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
            
            routine_FrameDifference(Itm1,It,Et,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testFDmorphoF/car_3%03d.pgm",i);
            CHRONO(fermeture3(Et,Etout,nrl,nrh,ncl,nch),cycles);
            totalcy += cycles;
            //fermeture5(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}   
	
	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles FD, only Ferm = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineFDmorpho3xFerm_opti(int seuil){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];

    sprintf(nameload1,"car3/car_3000.pgm");

	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);

	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"car3/car_3%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);

            routine_FrameDifference(Itm1,It,Et,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testOptiFD/car_3%03d.pgm",i);
            CHRONO(fermeture3_opti(Et,Etout,nrl,nrh,ncl,nch),cycles);
            totalcy += cycles;
            //fermeture5(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles FD, only Ferm opti = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineFDmorpho3xFerm_pipe(int seuil){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];

    sprintf(nameload1,"car3/car_3000.pgm");

	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **tmp = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD); //Matrice tmp entre dil et ero pour le pipe

	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"car3/car_3%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);

            routine_FrameDifference(Itm1,It,Et,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testFDmorphoFpipe/car_3%03d.pgm",i);
            CHRONO(fermeture3_pipe(Et,tmp,Etout,nrl,nrh,ncl,nch),cycles);
            totalcy += cycles;
            //ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}
	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles FD, only Ferm pipe = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);

}

void test_routineFDmorpho3xFerm_bin(int seuil){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];

    sprintf(nameload1,"car3/car_3000.pgm");

	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);

	ulong64 **Etbin=ui64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **Etoutbin=ui64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **tmpbin=ui64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);

	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"car3/car_3%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);

            routine_FrameDifference(Itm1,It,Et,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testFDmorphoFbin/car_3%03d.pgm",i);
            convCharToBin(Et,Etbin,nrl,nrh,ncl,nch);
            CHRONO(fermeture3_bin(Etbin,tmpbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS)),cycles);
            convBinToChar(Etoutbin,Etout,nrl,nrh,ncl,nch);
            totalcy += cycles;
            //ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}
	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles FD, only Ferm bin = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);

	free_ui64matrix(Etbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_ui64matrix(Etoutbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_ui64matrix(tmpbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);

}


void test_routineFDmorpho3xOuvFerm(int seuil){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];
    
    sprintf(nameload1,"car3/car_3000.pgm");
       
	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	
	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"car3/car_3%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
            
            routine_FrameDifference(Itm1,It,Et,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testFDmorphoOF/car_3%03d.pgm",i);
            CHRONO(ouverture3(Et,Etout,nrl,nrh,ncl,nch),cycles);
            totalcy += cycles;
            CHRONO(fermeture3(Etout,Et,nrl,nrh,ncl,nch),cycles);
            totalcy += cycles;
            //ouverture5(Et,Etout,nrl,nrh,ncl,nch);
            //fermeture5(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}   

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles FD, only OuvFerm = "));
	BENCH(printf(format,totalcy));
	
	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineFDmorpho3xOuvFerm_opti(int seuil){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];

    sprintf(nameload1,"car3/car_3000.pgm");

	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);

	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"car3/car_3%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);

            routine_FrameDifference(Itm1,It,Et,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testOptiFD/car_3%03d.pgm",i);
            CHRONO(ouverture3_opti(Et,Etout,nrl,nrh,ncl,nch),cycles);
            totalcy += cycles;
            CHRONO(fermeture3_opti(Etout,Et,nrl,nrh,ncl,nch),cycles);
            totalcy += cycles;
            //ouverture5(Et,Etout,nrl,nrh,ncl,nch);
            //fermeture5(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles FD, only OuvFerm opti = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineFDmorpho3xFermOuv(int seuil){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];
    
    sprintf(nameload1,"car3/car_3000.pgm");
       
	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	
	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"car3/car_3%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
            
            routine_FrameDifference(Itm1,It,Et,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testFDmorphoFO/car_3%03d.pgm",i);
            CHRONO(fermeture3(Et,Etout,nrl,nrh,ncl,nch),cycles);
            totalcy += cycles;
            CHRONO(ouverture3(Etout,Et,nrl,nrh,ncl,nch),cycles);
            totalcy += cycles;
            //fermeture5(Et,Etout,nrl,nrh,ncl,nch);
            //ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}   
	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles FD, only FermOuv = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}
    
void test_routineFDmorpho3xFermOuv_opti(int seuil){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];

    sprintf(nameload1,"car3/car_3000.pgm");

	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);

	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"car3/car_3%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);

            routine_FrameDifference(Itm1,It,Et,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testOptiFD/car_3%03d.pgm",i);
            CHRONO(fermeture3_opti(Et,Etout,nrl,nrh,ncl,nch),cycles);
            totalcy += cycles;
            CHRONO(ouverture3_opti(Etout,Et,nrl,nrh,ncl,nch),cycles);
            totalcy += cycles;
            //fermeture5(Et,Etout,nrl,nrh,ncl,nch);
            //ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}
	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles FD, only FermOuv opti = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}
    
    
//////////////////////////////////////////
//     TEST SIGMA DELTA O,OF,F,FO	    //
//////////////////////////////////////////
    
       
void test_routineSDmorpho3xOuv(){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;
	
	sprintf(nameload,"car3/car_3000.pgm");
	
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    uint8 **I = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **V  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Vtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **M  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Mtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	//uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);
	
	//printf("name=%s",nameload);
	
	routine_SigmaDelta_step0(Vtm1, Mtm1, Itm1, nrl, nrh, ncl, nch);
	
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload,"car3/car_3%03d.pgm",i);
	        I= LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        CHRONO(ouverture3(Et,Etout,nrl,nrh,ncl,nch),cycles);
	        totalcy += cycles;
	        //ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            sprintf(namesave,"testSDmorphoO/car_3%03d.pgm",i);
            SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
           
	} 

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD, only Ouv = "));
	BENCH(printf(format,totalcy));
	
	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(I,nrl,nrh,ncl,nch);
	free_ui8matrix(V,nrl,nrh,ncl,nch);
	free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(M,nrl,nrh,ncl,nch);
	free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineSDmorpho3xOuv_opti(){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;

	sprintf(nameload,"car3/car_3000.pgm");

    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    uint8 **I = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **V  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Vtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **M  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Mtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	//uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);

	//printf("name=%s",nameload);

	routine_SigmaDelta_step0(Vtm1, Mtm1, Itm1, nrl, nrh, ncl, nch);

	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload,"car3/car_3%03d.pgm",i);
	        I= LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        CHRONO(ouverture3_opti(Et,Etout,nrl,nrh,ncl,nch),cycles);
	        totalcy += cycles;
	        //ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            sprintf(namesave,"testOptiSD/car_3%03d.pgm",i);
            SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);

	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD, only Ouv opti = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(I,nrl,nrh,ncl,nch);
	free_ui8matrix(V,nrl,nrh,ncl,nch);
	free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(M,nrl,nrh,ncl,nch);
	free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineSDmorpho3xOuv_pipe(){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;

	sprintf(nameload,"car3/car_3000.pgm");

    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    uint8 **I = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **V  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Vtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **M  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Mtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **tmp = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD); //Matrice tmp entre dil et ero pour le pipe
	//uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);

	//printf("name=%s",nameload);

	routine_SigmaDelta_step0(Vtm1, Mtm1, Itm1, nrl, nrh, ncl, nch);

	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload,"car3/car_3%03d.pgm",i);
	        I= LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        CHRONO(ouverture3_pipe(Et,tmp,Etout,nrl,nrh,ncl,nch),cycles);
	        totalcy += cycles;
	        //ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            sprintf(namesave,"testSDmorphoOpipe/car_3%03d.pgm",i);
            SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);

	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD, only Ouv pipe = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(I,nrl,nrh,ncl,nch);
	free_ui8matrix(V,nrl,nrh,ncl,nch);
	free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(M,nrl,nrh,ncl,nch);
	free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineSDmorpho3xOuv_bin(){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;

	sprintf(nameload,"car3/car_3000.pgm");

    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    uint8 **I = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **V  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Vtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **M  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Mtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	//uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);

	ulong64 **Etbin=ui64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **Etoutbin=ui64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **tmpbin=ui64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);

	//printf("name=%s",nameload);

	routine_SigmaDelta_step0(Vtm1, Mtm1, Itm1, nrl, nrh, ncl, nch);

	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload,"car3/car_3%03d.pgm",i);
	        I= LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        convCharToBin(Et,Etbin,nrl,nrh,ncl,nch);
			CHRONO(ouverture3_bin(Etbin,tmpbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS)),cycles);
			convBinToChar(Etoutbin,Etout,nrl,nrh,ncl,nch);
	        totalcy += cycles;
	        //ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            sprintf(namesave,"testSDmorphoObin/car_3%03d.pgm",i);
            SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);

	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD, only Ouv bin = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(I,nrl,nrh,ncl,nch);
	free_ui8matrix(V,nrl,nrh,ncl,nch);
	free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(M,nrl,nrh,ncl,nch);
	free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);

	free_ui64matrix(Etbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_ui64matrix(Etoutbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_ui64matrix(tmpbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
}

void test_routineSDmorpho3xFerm(){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;
	
	sprintf(nameload,"car3/car_3000.pgm");
	
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    uint8 **I = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **V  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Vtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **M  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Mtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	//uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);
	
	//printf("name=%s",nameload);
	
	routine_SigmaDelta_step0(Vtm1, Mtm1, Itm1, nrl, nrh, ncl, nch);
	
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload,"car3/car_3%03d.pgm",i);
	        I= LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        CHRONO(fermeture3(Et,Etout,nrl,nrh,ncl,nch),cycles);
	        totalcy += cycles;
	        //fermeture5(Etout,Et,nrl,nrh,ncl,nch);
            sprintf(namesave,"testSDmorphoF/car_3%03d.pgm",i);
            SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
           
	} 

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD, only Ferm = "));
	BENCH(printf(format,totalcy));
	
	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(I,nrl,nrh,ncl,nch);
	free_ui8matrix(V,nrl,nrh,ncl,nch);
	free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(M,nrl,nrh,ncl,nch);
	free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineSDmorpho3xFerm_opti(){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;

	sprintf(nameload,"car3/car_3000.pgm");

    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    uint8 **I = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **V  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Vtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **M  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Mtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	//uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);

	//printf("name=%s",nameload);

	routine_SigmaDelta_step0(Vtm1, Mtm1, Itm1, nrl, nrh, ncl, nch);

	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload,"car3/car_3%03d.pgm",i);
	        I= LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        CHRONO(fermeture3_opti(Et,Etout,nrl,nrh,ncl,nch),cycles);
	        totalcy += cycles;
	        //fermeture5(Etout,Et,nrl,nrh,ncl,nch);
            sprintf(namesave,"testOptiSD/car_3%03d.pgm",i);
            SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);

	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD, only Ferm opti = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(I,nrl,nrh,ncl,nch);
	free_ui8matrix(V,nrl,nrh,ncl,nch);
	free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(M,nrl,nrh,ncl,nch);
	free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineSDmorpho3xFerm_pipe(){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;

	sprintf(nameload,"car3/car_3000.pgm");

    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    uint8 **I = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **V  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Vtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **M  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Mtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **tmp = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD); //Matrice tmp entre dil et ero pour le pipe
	//uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);

	//printf("name=%s",nameload);

	routine_SigmaDelta_step0(Vtm1, Mtm1, Itm1, nrl, nrh, ncl, nch);

	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload,"car3/car_3%03d.pgm",i);
	        I= LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        CHRONO(fermeture3_pipe(Et,tmp,Etout,nrl,nrh,ncl,nch),cycles);
	        totalcy += cycles;
	        //ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            sprintf(namesave,"testSDmorphoFpipe/car_3%03d.pgm",i);
            SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);

	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD, only Ferm pipe = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(I,nrl,nrh,ncl,nch);
	free_ui8matrix(V,nrl,nrh,ncl,nch);
	free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(M,nrl,nrh,ncl,nch);
	free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineSDmorpho3xFerm_bin(){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;

	sprintf(nameload,"car3/car_3000.pgm");

    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    uint8 **I = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **V  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Vtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **M  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Mtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	//uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);

	ulong64 **Etbin=ui64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **Etoutbin=ui64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	ulong64 **tmpbin=ui64matrix(nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);

	//printf("name=%s",nameload);

	routine_SigmaDelta_step0(Vtm1, Mtm1, Itm1, nrl, nrh, ncl, nch);

	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload,"car3/car_3%03d.pgm",i);
	        I= LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        convCharToBin(Et,Etbin,nrl,nrh,ncl,nch);
			CHRONO(fermeture3_bin(Etbin,tmpbin,Etoutbin,nrl,nrh,(ncl/NBBITS),(nch/NBBITS)),cycles);
			convBinToChar(Etoutbin,Etout,nrl,nrh,ncl,nch);
	        totalcy += cycles;
	        //ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            sprintf(namesave,"testSDmorphoFbin/car_3%03d.pgm",i);
            SavePGM_ui8matrix(Etout,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);

	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD, only Ferm bin = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(I,nrl,nrh,ncl,nch);
	free_ui8matrix(V,nrl,nrh,ncl,nch);
	free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(M,nrl,nrh,ncl,nch);
	free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);

	free_ui64matrix(Etbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_ui64matrix(Etoutbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
	free_ui64matrix(tmpbin,nrl-BORD,nrh+BORD,(ncl/NBBITS)-BORD,(nch/NBBITS)+BORD);
}

void test_routineSDmorpho3xOuvFerm(){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;
	
	sprintf(nameload,"car3/car_3000.pgm");
	
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    uint8 **I = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **V  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Vtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **M  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Mtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	//uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);
	
	//printf("name=%s",nameload);
	
	routine_SigmaDelta_step0(Vtm1, Mtm1, Itm1, nrl, nrh, ncl, nch);
	
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload,"car3/car_3%03d.pgm",i);
	        I= LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        CHRONO(ouverture3(Et,Etout,nrl,nrh,ncl,nch),cycles);
	        totalcy += cycles;
	        CHRONO(fermeture3(Etout,Et,nrl,nrh,ncl,nch),cycles);
	        totalcy += cycles;
	        //ouverture5(Et,Etout,nrl,nrh,ncl,nch);
	        //fermeture5(Etout,Et,nrl,nrh,ncl,nch);
            sprintf(namesave,"testSDmorphoOF/car_3%03d.pgm",i);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
           
	} 

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD, only OuvFerm = "));
	BENCH(printf(format,totalcy));
	
	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(I,nrl,nrh,ncl,nch);
	free_ui8matrix(V,nrl,nrh,ncl,nch);
	free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(M,nrl,nrh,ncl,nch);
	free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineSDmorpho3xOuvFerm_opti(){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;

	sprintf(nameload,"car3/car_3000.pgm");

    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    uint8 **I = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **V  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Vtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **M  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Mtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	//uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);

	//printf("name=%s",nameload);

	routine_SigmaDelta_step0(Vtm1, Mtm1, Itm1, nrl, nrh, ncl, nch);

	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload,"car3/car_3%03d.pgm",i);
	        I= LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        CHRONO(ouverture3_opti(Et,Etout,nrl,nrh,ncl,nch),cycles);
	        totalcy += cycles;
	        CHRONO(fermeture3_opti(Etout,Et,nrl,nrh,ncl,nch),cycles);
	        totalcy += cycles;
	        //ouverture5(Et,Etout,nrl,nrh,ncl,nch);
	        //fermeture5(Etout,Et,nrl,nrh,ncl,nch);
            sprintf(namesave,"testOptiSD/car_3%03d.pgm",i);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);

	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD, only OuvFerm opti = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(I,nrl,nrh,ncl,nch);
	free_ui8matrix(V,nrl,nrh,ncl,nch);
	free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(M,nrl,nrh,ncl,nch);
	free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineSDmorpho3xFermOuv(){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;
	
	sprintf(nameload,"car3/car_3000.pgm");
	
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    uint8 **I = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **V  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Vtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **M  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Mtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	//uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);
	
	//printf("name=%s",nameload);
	
	routine_SigmaDelta_step0(Vtm1, Mtm1, Itm1, nrl, nrh, ncl, nch);
	
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload,"car3/car_3%03d.pgm",i);
	        I= LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        CHRONO(fermeture3(Et,Etout,nrl,nrh,ncl,nch),cycles);
	        totalcy += cycles;
	        CHRONO(ouverture3(Etout,Et,nrl,nrh,ncl,nch),cycles);
	        totalcy += cycles;
	        //fermeture5(Et,Etout,nrl,nrh,ncl,nch);
	        //ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            sprintf(namesave,"testSDmorphoFO/car_3%03d.pgm",i);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
           
	} 
	
	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD, only FermOuv = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(I,nrl,nrh,ncl,nch);
	free_ui8matrix(V,nrl,nrh,ncl,nch);
	free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(M,nrl,nrh,ncl,nch);
	free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineSDmorpho3xFermOuv_opti(){

	//Cycle par point//
	double cycles, totalcy=0;
	int iter, niter=2;
	int run, nrun = 5;
	double t0,t1,dt,tmin,t;

	char *format = "%6.2f\n";
	///////////////////

    long nrl, nrh, ncl, nch;
    char nameload[100];     //"car3/car_3..";
	char namesave[100];     //"testSD/SD...";
	int i;

	sprintf(nameload,"car3/car_3000.pgm");

    uint8 **Itm1 = LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    uint8 **I = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **V  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Vtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **M  = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Mtm1 = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	//uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);

	//printf("name=%s",nameload);

	routine_SigmaDelta_step0(Vtm1, Mtm1, Itm1, nrl, nrh, ncl, nch);

	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload,"car3/car_3%03d.pgm",i);
	        I= LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        CHRONO(fermeture3_opti(Et,Etout,nrl,nrh,ncl,nch),cycles);
	        totalcy += cycles;
	        CHRONO(ouverture3_opti(Etout,Et,nrl,nrh,ncl,nch),cycles);
	        totalcy += cycles;
	        //fermeture5(Et,Etout,nrl,nrh,ncl,nch);
	        //ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            sprintf(namesave,"testOptiSD/car_3%03d.pgm",i);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);

	}

	totalcy = totalcy / NBIMAGES; //on doit rediviser par le nombre de points pour l'avoir par point
	totalcy = totalcy / ((nch+1)*(nrh+1));

	BENCH(printf("Cycles SD, only FermOuv opti = "));
	BENCH(printf(format,totalcy));

	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(I,nrl,nrh,ncl,nch);
	free_ui8matrix(V,nrl,nrh,ncl,nch);
	free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(M,nrl,nrh,ncl,nch);
	free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}


