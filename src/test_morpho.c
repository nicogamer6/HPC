#include <stdio.h>
#include "nrutil.h"
#include "nrdef.h"
#include "string.h"
#include "test_morpho.h"
#include "morpho.h"
#include "mouvement.h"

//#define SEUILFD 20
#define NBIMAGES 299
#define BORD 2


//////////////////////////////////////////////
//     TEST FRAME DIFFERENCE O,OF,F,FO	    //
//////////////////////////////////////////////

void test_routineFDmorpho3xOuv(int seuil){
    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];
    
    sprintf(nameload1,"hall/hall000000.pgm");
       
	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	
	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"hall/hall000%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
            
            routine_FrameDifference(Itm1,It,Et,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testFDmorphoO/hall000%03d.pgm",i);
            ouverture3(Et,Etout,nrl,nrh,ncl,nch);
            ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}   
	
	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	
}

void test_routineFDmorpho3xFerm(int seuil){
    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];
    
    sprintf(nameload1,"hall/hall000000.pgm");
       
	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	
	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"hall/hall000%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
            
            routine_FrameDifference(Itm1,It,Et,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testFDmorphoF/hall000%03d.pgm",i);
            fermeture3(Et,Etout,nrl,nrh,ncl,nch);
            fermeture5(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}   
	
	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineFDmorpho3xOuvFerm(int seuil){
    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];
    
    sprintf(nameload1,"hall/hall000000.pgm");
       
	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	
	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"hall/hall000%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
            
            routine_FrameDifference(Itm1,It,Et,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testFDmorphoOF/hall000%03d.pgm",i);
            ouverture3(Et,Etout,nrl,nrh,ncl,nch);
            fermeture3(Etout,Et,nrl,nrh,ncl,nch);
            ouverture5(Et,Etout,nrl,nrh,ncl,nch);
            fermeture5(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}   
	
	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineFDmorpho3xFermOuv(int seuil){
    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];
    
    sprintf(nameload1,"hall/hall000000.pgm");
       
	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **Et = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	uint8 **Etout = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	
	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"hall/hall000%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
            
            routine_FrameDifference(Itm1,It,Et,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testFDmorphoFO/hall000%03d.pgm",i);
            fermeture3(Et,Etout,nrl,nrh,ncl,nch);
            ouverture3(Etout,Et,nrl,nrh,ncl,nch);
            fermeture5(Et,Etout,nrl,nrh,ncl,nch);
            ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}   
	
	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}
    
    
    
//////////////////////////////////////////
//     TEST SIGMA DELTA O,OF,F,FO	    //
//////////////////////////////////////////
    
       
void test_routineSDmorpho3xOuv(){
    long nrl, nrh, ncl, nch;
    char nameload[100];     //"hall/hall000..";
	char namesave[100];     //"testSD/SD...";
	int i;
	
	sprintf(nameload,"hall/hall000000.pgm");
	
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
	        sprintf(nameload,"hall/hall000%03d.pgm",i);
	        I= LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        ouverture3(Et,Etout,nrl,nrh,ncl,nch);
	        ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            sprintf(namesave,"testSDmorphoO/hall000%03d.pgm",i);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
           
	} 
	
	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(I,nrl,nrh,ncl,nch);
	free_ui8matrix(V,nrl,nrh,ncl,nch);
	free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(M,nrl,nrh,ncl,nch);
	free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineSDmorpho3xFerm(){
    long nrl, nrh, ncl, nch;
    char nameload[100];     //"hall/hall000..";
	char namesave[100];     //"testSD/SD...";
	int i;
	
	sprintf(nameload,"hall/hall000000.pgm");
	
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
	        sprintf(nameload,"hall/hall000%03d.pgm",i);
	        I= LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        fermeture3(Et,Etout,nrl,nrh,ncl,nch);
	        fermeture5(Etout,Et,nrl,nrh,ncl,nch);
            sprintf(namesave,"testSDmorphoF/hall000%03d.pgm",i);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
           
	} 
	
	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(I,nrl,nrh,ncl,nch);
	free_ui8matrix(V,nrl,nrh,ncl,nch);
	free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(M,nrl,nrh,ncl,nch);
	free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}

void test_routineSDmorpho3xOuvFerm(){
    long nrl, nrh, ncl, nch;
    char nameload[100];     //"hall/hall000..";
	char namesave[100];     //"testSD/SD...";
	int i;
	
	sprintf(nameload,"hall/hall000000.pgm");
	
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
	        sprintf(nameload,"hall/hall000%03d.pgm",i);
	        I= LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        ouverture3(Et,Etout,nrl,nrh,ncl,nch);
	        fermeture3(Etout,Et,nrl,nrh,ncl,nch);
	        ouverture5(Et,Etout,nrl,nrh,ncl,nch);
	        fermeture5(Etout,Et,nrl,nrh,ncl,nch);
            sprintf(namesave,"testSDmorphoOF/hall000%03d.pgm",i);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
           
	} 
	
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
    long nrl, nrh, ncl, nch;
    char nameload[100];     //"hall/hall000..";
	char namesave[100];     //"testSD/SD...";
	int i;
	
	sprintf(nameload,"hall/hall000000.pgm");
	
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
	        sprintf(nameload,"hall/hall000%03d.pgm",i);
	        I= LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        fermeture3(Et,Etout,nrl,nrh,ncl,nch);
	        ouverture3(Etout,Et,nrl,nrh,ncl,nch);
	        fermeture5(Et,Etout,nrl,nrh,ncl,nch);
	        ouverture5(Etout,Et,nrl,nrh,ncl,nch);
            sprintf(namesave,"testSDmorphoFO/hall000%03d.pgm",i);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
           
	} 
	
	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(I,nrl,nrh,ncl,nch);
	free_ui8matrix(V,nrl,nrh,ncl,nch);
	free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(M,nrl,nrh,ncl,nch);
	free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}



