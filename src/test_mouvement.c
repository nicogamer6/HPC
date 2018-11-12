#include <stdio.h>
#include "nrutil.h"
#include "nrdef.h"
#include "string.h"
#include "test_mouvement.h"
#include "mouvement.h"

#define NBIMAGES 299
#define BORD 2

//Permet de remplir le dossier /testFD et de voir les nouvelles images avec algo FD sur l'ensemble du dossier /hall
void test_routineFD(int seuil){
    long nrl, nrh, ncl, nch;
    char nameload1[100], nameload2[100], namesave[100];
    
    sprintf(nameload1,"hall/hall000000.pgm");
       
	uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	uint8 **m = ui8matrix(nrl,nrh,ncl,nch);
	uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
	//uint8 **m, **m1 , **m2;
		
	int i;
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload2,"hall/hall000%03d.pgm",i);
            It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
            routine_FrameDifference(Itm1,It,m,nrl,nrh,ncl,nch,seuil);
            sprintf(namesave,"testFD/hall000%03d.pgm",i);
            SavePGM_ui8matrix(m,nrl,nrh,ncl,nch,namesave);
            copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
	}   
	
	free_ui8matrix(m,nrl,nrh,ncl,nch);
	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
}

//Permet de remplir le dossier /testSD et de voir les nouvelles images avec algo SD sur l'ensemble du dossier /hall
void test_routineSD(void){
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
	//uint8 **Ot = ui8matrix(nrl,nrh,ncl,nch);
	
	//printf("name=%s",nameload);
	
	routine_SigmaDelta_step0(Vtm1, Mtm1, Itm1, nrl, nrh, ncl, nch);
	
	for(i=1;i<=NBIMAGES;i++){
	        sprintf(nameload,"hall/hall000%03d.pgm",i);
	        I=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	        routine_SigmaDelta_1step(V, Vtm1, M, Mtm1, I, Et, nrl, nrh, ncl, nch);
	        
            sprintf(namesave,"testSD/hall000%03d.pgm",i);
            SavePGM_ui8matrix(Et,nrl,nrh,ncl,nch,namesave);
            //On doit copier M dan Mtm1, V dans Vtm1 et I dans Itm1
            
            copy_ui8matrix_ui8matrix(V, nrl, nrh, ncl, nch, Vtm1);
            copy_ui8matrix_ui8matrix(M, nrl, nrh, ncl, nch, Mtm1);
            copy_ui8matrix_ui8matrix(I, nrl, nrh, ncl, nch, Itm1);
           
	} 
	
	free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
	free_ui8matrix(I,nrl,nrh,ncl,nch);
	free_ui8matrix(V,nrl,nrh,ncl,nch);
	free_ui8matrix(Vtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(M,nrl,nrh,ncl,nch);
	free_ui8matrix(Mtm1,nrl,nrh,ncl,nch);
	free_ui8matrix(Et,nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
	//free_ui8matrix(Ot,nrl,nrh,ncl,nch);
}
