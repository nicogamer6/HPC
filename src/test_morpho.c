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
            sprintf(namesave,"testFDmorphoO/FDmorpho3xOUV%03d.pgm",i);
            ouverture3(Et,Etout,nrl,nrh,ncl,nch);
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
            sprintf(namesave,"testFDmorphoF/FDmorpho3xFERM%03d.pgm",i);
            fermeture3(Et,Etout,nrl,nrh,ncl,nch);
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
            sprintf(namesave,"testFDmorphoOF/FDmorpho3xOUVFERM%03d.pgm",i);
            ouverture3(Et,Etout,nrl,nrh,ncl,nch);
            fermeture3(Etout,Et,nrl,nrh,ncl,nch);
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
            sprintf(namesave,"testFDmorphoFO/FDmorpho3xFERMOUV%03d.pgm",i);
            fermeture3(Et,Etout,nrl,nrh,ncl,nch);
            ouverture3(Etout,Et,nrl,nrh,ncl,nch);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}   
	
	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	free_ui8matrix(Etout,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
}
    
void test_routineSDmorpho3xOuv(){
    
}

void test_routineSDmorpho3xFerm(){

}

void test_routineSDmorpho3xOuvFerm(){

}

void test_routineSDmorpho3xFermOuv(){

}



