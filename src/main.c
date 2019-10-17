#include <stdio.h>
#include "nrutil.h"
#include "nrdef.h"
#include "string.h"
#include <omp.h>

#include "mouvement.h"
#include "test_mouvement.h"
#include "test_mouvement_SSE2.h"
#include "morpho.h"
#include "test_morpho.h"
#include "matriceROC.h"

#include "morpho_SSE2.h"
#include "test_morpho_SSE2.h"

#define NB 2
#define SEUILFD 20
#define BORD 2

#define NBIMAGES 199


void diffImages(char dossier1[], char dossier2[]){
    long nrl,nrh,ncl,nch;
    int k;
    int i,j;
    int cpt=0;
    char nameload1[100], nameload2[100];

    uint8 **It0 = LoadPGM_ui8matrix("car3/car_3000.pgm",&nrl,&nrh,&ncl,&nch);
	uint8 **It1 = ui8matrix(nrl,nrh,ncl,nch);

	printf("\nDIFF entre %s et %s \n\n", dossier1, dossier2);

    for(k=1; k<=NBIMAGES;k++){
    	cpt=0;
        sprintf(nameload1,"%s/car_3%03d.pgm",dossier1, k);
        sprintf(nameload2,"%s/car_3%03d.pgm",dossier2, k);

	    It0=LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	    It1=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);

        for(i=nrl;i<=nrh;i++){
            for(j=ncl;j<=nch;j++){
            	if(It0[i][j] != It1[i][j]){
            		cpt++;
            		printf("Pixel different : i,j = [%d][%d], Image n° %d\t",i,j,k);
            	}
            }
        }
        printf("\nImage n° %d : pixel diff = %d\n\n",k,cpt);
    }

    free_ui8matrix(It0,nrl,nrh,ncl,nch);
	free_ui8matrix(It1,nrl,nrh,ncl,nch);
}



int main()
{
     
    // Crée les images FD et SD
    // AFFICHE le CPP juste pour la morpho
    //printf("SEUIL FD = %d\n",SEUILFD);
    test_routineFD(SEUILFD);                // Dossier "testFD"
    test_routineSD();                       // Dossier "testSD"
    test_routineSD_opti();                  // Dossier "testOptiSD"
    test_routineSD_soa();                   // Dossier "testSD_SoA"

    printf("\n");
    ///////////////////////////////////////////////////////////////////////////

    // Crée les images FD et SD SSE
    // AFFICHE le CPP juste pour la morpho

    //test_routineFD_SSE(SEUILFD);            // Dossier "testFD_SSE"
    //test_routineFD_SSE_OMP(SEUILFD);        // Dossier "testFD_SSE_OMP"
    //test_routineSD_SSE();                   // Dossier "testSD_SSE"
    //test_routineSD_SSE_OMP();               // Dossier "testSD_SSE_OMP"
    //printf("\n");
    ///////////////////////////////////////////////////////////////////////////
    
    // Crée les images FD avec les morpho 3x3, O, F, OF et FO. 
    // AFFICHE le CPP juste pour la morpho et morphoSSE

    test_routineFDmorpho3xOuv(SEUILFD);     // Dossier "testFDmorphoO"
    //test_routineFDmorpho3xFerm(SEUILFD);    // Dossier "testFDmorphoF"
    //test_routineFDmorpho3xOuvFerm(SEUILFD); // Dossier "testFDmorphoOF"
    //test_routineFDmorpho3xFermOuv(SEUILFD); // Dossier "testFDmorphoFO"
    //printf("\n");

	test_routineFDmorpho3xOuv_opti(SEUILFD);		// Dossier "testOptiFD"
	//test_routineFDmorpho3xFerm_opti(SEUILFD);		// Dossier "testOptiFD"
	//test_routineFDmorpho3xOuvFerm_opti(SEUILFD);	// Dossier "testOptiFD"
	//test_routineFDmorpho3xFermOuv_opti(SEUILFD);	// Dossier "testOptiFD"

	test_routineFDmorpho3xOuv_pipe(SEUILFD);	// Dossier "testFDmorphoOpipe"
	//test_routineFDmorpho3xFerm_pipe(SEUILFD);	// Dossier "testFDmorphoFpipe"

	test_routineFDmorpho3xOuv_bin(SEUILFD);	// Dossier "testFDmorphoObin"
	//test_routineFDmorpho3xFerm_bin(SEUILFD);	// Dossier "testFDmorphoFbin"


    test_routineFD_SSEmorpho3xOuv(SEUILFD);     // Dossier "testFD_SSEmorphoO"
    //test_routineFD_SSEmorpho3xFerm(SEUILFD);    // Dossier "testFD_SSEmorphoF"
    //test_routineFD_SSEmorpho3xOuvFerm(SEUILFD); // Dossier "testFD_SSEmorphoOF"
    //test_routineFD_SSEmorpho3xFermOuv(SEUILFD); // Dossier "testFD_SSEmorphoFO"
    printf("\n");



    ///////////////////////////////////////////////////////////////////////////
    // Crée les images FD avec les morpho 3x3, O, F, OF et FO.
    // Affiche le CPP juste pour la morpho et morphoSSE

    test_routineSDmorpho3xOuv();            // Dossier "testSDmorphoO"
    //test_routineSDmorpho3xFerm();           // Dossier "testSDmorphoF"
    //test_routineSDmorpho3xOuvFerm();        // Dossier "testSDmorphoOF"
    //test_routineSDmorpho3xFermOuv();        // Dossier "testSDmorphoFO"
    //printf("\n");

	test_routineSDmorpho3xOuv_opti();		// Dossier "testOptiSD"
	//test_routineSDmorpho3xFerm_opti();	// Dossier "testOptiSD"
	//test_routineSDmorpho3xOuvFerm_opti();	// Dossier "testOptiSD"
	//test_routineSDmorpho3xFermOuv_opti();	// Dossier "testOptiSD"

	test_routineSDmorpho3xOuv_pipe();	// Dossier "testSDmorphoOpipe"
	//test_routineSDmorpho3xFerm_pipe();	// Dossier "testSDmorphoFpipe"

	test_routineSDmorpho3xOuv_bin();	// Dossier "testSDmorphoObin"
	//test_routineSDmorpho3xFerm_bin();	// Dossier "testSDmorphoFbin"

    test_routineSD_SSEmorpho3xOuv();        // Dossier "testSD_SSEmorphoO"
    //test_routineSD_SSEmorpho3xFerm();       // Dossier "testSD_SSEmorphoF"
    //test_routineSD_SSEmorpho3xOuvFerm();    // Dossier "testSD_SSEmorphoOF"
    //test_routineSD_SSEmorpho3xFermOuv();    // Dossier "testSD_SSEmorphoFO"
    //printf("\n");


    ///////////////////////////////////////////////////////////////////////////

    // Crée les images morpho et morpho SSE Ero, Dil, ouv et ferm dans des sous dossiers 
    test_Etapemorpho();                     // Dossier "testmorpho"
    //test_EtapemorphoSSE();                  // Dossier "testmorphoSSE"


    diffImages("testSDmorphoO", "testSDmorphoOpipe");



    ///////////////////////////////////////////////////////////////////////////

    // Pour créer la matrice ROC à partir d'un dossier d'image comparé au dossier "verite" terrain: 
    // Utilisation : matriceROC("nomDossier")
    // Paramètre de comparaison : MCC entre -1 et 1
    // 1 = perfect prediction, -1 = Total disagreement between prediction and observation
    
    //matriceROC("testFD");
    //matriceROC("testFDmorphoO");
    //matriceROC("testFDmorphoF");
    //matriceROC("testFDmorphoOF");
    //matriceROC("testFDmorphoFO");
    //printf("\n");
    //matriceROC("testSD");
    //matriceROC("testSDmorphoO");
    //matriceROC("testSDmorphoF");
    //matriceROC("testSDmorphoOF");
    //matriceROC("testSDmorphoFO");
    //printf("\n");
    
    //matriceROC("testFD_SSE");
    //matriceROC("testFD_SSEmorphoO");
    //matriceROC("testFD_SSEmorphoF");
    //matriceROC("testFD_SSEmorphoOF");
    //matriceROC("testFD_SSEmorphoFO");
    //printf("\n");
    //matriceROC("testSD_SSE");
    //matriceROC("testSD_SSEmorphoO");
    //matriceROC("testSD_SSEmorphoF");
    //matriceROC("testSD_SSEmorphoOF");
    //matriceROC("testSD_SSEmorphoFO");
    //printf("\n");
    ///////////////////////////////////////////////////////////////////////////
    
    // Tests unitaires pour SD ET FD en SSE
    //test_unitaire_SD_SSE();
    //test_unitaire_FD_SSE();
    

    return 0;
}



