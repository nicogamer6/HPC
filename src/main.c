#include <stdio.h>
#include "nrutil.h"
#include "nrdef.h"
#include "string.h"

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


int main()
{
     
    // Crée les images FD et SD. 
    // AFFICHE le CPP juste pour la morpho
    test_routineFD(SEUILFD);                // Dossier "testFD"
    test_routineSD();                       // Dossier "testSD"
    printf("\n");   
    ///////////////////////////////////////////////////////////////////////////

    // Crée les images FD et SD SSE
    // AFFICHE le CPP juste pour la morpho
    test_routineFD_SSE(SEUILFD);            // Dossier "testFD_SSE"
    test_routineSD_SSE();                   // Dossier "testSD_SSE"
    printf("\n");
    ///////////////////////////////////////////////////////////////////////////
    
    // Crée les images FD avec les morpho 3x3, O, F, OF et FO. 
    // AFFICHE le CPP juste pour la morpho et morphoSSE
    test_routineFDmorpho3xOuv(SEUILFD);     // Dossier "testFDmorphoO"
    test_routineFDmorpho3xFerm(SEUILFD);    // Dossier "testFDmorphoF"
    test_routineFDmorpho3xOuvFerm(SEUILFD); // Dossier "testFDmorphoOF"
    test_routineFDmorpho3xFermOuv(SEUILFD); // Dossier "testFDmorphoFO"
    printf("\n");

    test_routineFD_SSEmorpho3xOuv(SEUILFD);     // Dossier "testFD_SSEmorphoO"
    test_routineFD_SSEmorpho3xFerm(SEUILFD);    // Dossier "testFD_SSEmorphoF"
    test_routineFD_SSEmorpho3xOuvFerm(SEUILFD); // Dossier "testFD_SSEmorphoOF"
    test_routineFD_SSEmorpho3xFermOuv(SEUILFD); // Dossier "testFD_SSEmorphoFO"
    printf("\n");
    ///////////////////////////////////////////////////////////////////////////

    // Crée les images FD avec les morpho 3x3, O, F, OF et FO.
    // Affiche le CPP juste pour la morpho et morphoSSE
    test_routineSDmorpho3xOuv();            // Dossier "testSDmorphoO"
    test_routineSDmorpho3xFerm();           // Dossier "testSDmorphoF"
    test_routineSDmorpho3xOuvFerm();        // Dossier "testSDmorphoOF"
    test_routineSDmorpho3xFermOuv();        // Dossier "testSDmorphoFO"
    printf("\n");

    test_routineSD_SSEmorpho3xOuv();        // Dossier "testSD_SSEmorphoO"
    test_routineSD_SSEmorpho3xFerm();       // Dossier "testSD_SSEmorphoF"
    test_routineSD_SSEmorpho3xOuvFerm();    // Dossier "testSD_SSEmorphoOF"
    test_routineSD_SSEmorpho3xFermOuv();    // Dossier "testSD_SSEmorphoFO"
    printf("\n");
    ///////////////////////////////////////////////////////////////////////////

    // Crée les images morpho et morpho SSE Ero, Dil, ouv et ferm dans des sous dossiers 
    test_Etapemorpho();                     // Dossier "testmorpho"
    test_EtapemorphoSSE();                  // Dossier "testmorphoSSE"

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
    printf("\n");
    //matriceROC("testSD");
    //matriceROC("testSDmorphoO");
    //matriceROC("testSDmorphoF");
    //matriceROC("testSDmorphoOF");
    //matriceROC("testSDmorphoFO");
    //printf("\n");
    
    matriceROC("testFD_SSE");
    matriceROC("testFD_SSEmorphoO");
    matriceROC("testFD_SSEmorphoF");
    matriceROC("testFD_SSEmorphoOF");
    matriceROC("testFD_SSEmorphoFO");
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



