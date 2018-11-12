#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "nrdef.h"
#include "string.h"
#include "test_morpho.h"
#include "morpho.h"
#include "mouvement.h"
#include "matriceROC.h"

#define NBIMAGES 299

void matriceROC(char dossier[]){
    long nrl,nrh,ncl,nch;
    int matROC[2][2]={0};
    
    char nameload[100];
    double mcc;
    double VP,FN,FP,VN;
    
    uint8 **Itverite = LoadPGM_ui8matrix("hall/hall000000.pgm",&nrl,&nrh,&ncl,&nch); //Pour récup les tailles images
    uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
    
    int k; 
    
    VP = matROC[0][0];
    FN = matROC[0][1];
    FP = matROC[1][0];
    VN = matROC[1][1];
    
    for(k=0; k<=NBIMAGES;k+10){
        sprintf(nameload,"verite/hall000%03d.pgm",k); //Image verite terrain à comparer avec image du dossier arg
	    Itverite=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
	    
	    sprintf(nameload,"%s/hall000%03d.pgm",dossier,k); //Image à tester du dossier en arg
	    It=LoadPGM_ui8matrix(nameload,&nrl,&nrh,&ncl,&nch);
    }
    
    int i,j;
    
    for(i=nrl;i<=nrh;i++){
        for(j=ncl;j<=nch;j++){
            if(It[i][j] == 255 && Itverite[i][j] == 255) //VP
                VP++;
            else if(It[i][j] == 0 && Itverite[i][j] == 255) //FN
                FN++;
            else if(It[i][j] == 255 && Itverite[i][j] == 0) //FP
                FP++;
            else if(It[i][j] == 0 && Itverite[i][j] == 0) //VN
                VN++;
        }
    }
    
    mcc = (VP * VN - FP * FN)/(sqrt((VP+FP)*(VP+FN)*(VN+FP)*(VN+FN)));
    
    printf("Matrice ROC pour %s : \n VP FN \t %f %f \n FP VN \t %f %f \n", dossier, VP, FN, FP, VN);
    
    printf("MCC = %f",mcc); //Valeur entre -1 et 1, 1 = perfect prediction, -1 = Total disagreement between prediction and observation
    
    free_ui8matrix(Itverite,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
}
