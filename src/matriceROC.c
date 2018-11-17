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

#define NBIMAGESVERITE 29
#define NBIMAGES 290

void matriceROC(char dossier[]){
    long nrl,nrh,ncl,nch;
    int matROC[2][2]={0};
    
    char nameload1[100], nameload2[100];
    double mcc,acc;
    double VP,FN,FP,VN;
    
    uint8 **Itverite = LoadPGM_ui8matrix("hall/hall000000.pgm",&nrl,&nrh,&ncl,&nch); //Pour récup les tailles images
    uint8 **It = ui8matrix(nrl,nrh,ncl,nch);
    
    int k; 
    int i,j;
    int perdu;
    
    VP = matROC[0][0];
    FN = matROC[0][1];
    FP = matROC[1][0];
    VN = matROC[1][1];
    
    for(k=0; k<=NBIMAGESVERITE;k++){
        sprintf(nameload1,"verite/hall000%03d.pgm",k*10); //Image verite terrain à comparer avec image du dossier arg
	    //printf("%s\n",nameload1);
	    Itverite=LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
	    
	    sprintf(nameload2,"%s/hall000%03d.pgm",dossier,k*10+1); //Image à tester du dossier en arg
	    //printf("%s",nameload2);
	    It=LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
    
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
                else perdu++;
            }
        }
    }
    
    mcc = (VP * VN - FP * FN)/(sqrt((VP+FP)*(VP+FN)*(VN+FP)*(VN+FN)));
    acc = (VP + VN) / (VP + VN + FP + FN);
    
    printf("Matrice ROC pour %s : \n \tVP FN \t %0.0f %0.0f \n \tFP VN \t %0.0f %0.0f \n", dossier, VP, FN, FP, VN);
    printf("\tpixel analysé : %0.0f\n",VP+FN+FP+VN);
    printf("\tpixel perdu : %d\n",perdu);
    
    printf("\tMCC = %f\n",mcc); //Valeur entre -1 et 1, 1 = perfect prediction, -1 = Total disagreement between prediction and observation
    printf("\tACCURACY = %f\n",acc); //Accuracy de ce dossier
    free_ui8matrix(Itverite,nrl,nrh,ncl,nch);
	free_ui8matrix(It,nrl,nrh,ncl,nch);
}
