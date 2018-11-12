#include <stdio.h>
#include "nrutil.h"
#include "nrdef.h"
#include "string.h"
#include "mouvement.h"
#include "test_mouvement.h"
#include "test_morpho.h"
#include "matriceROC.h"

#include "morpho.h"

#define NB 2
#define SEUILFD 20


int main()
{
    /*long nrl, nrh, ncl, nch;
	uint8 **m;
	uint8 **m1;
	uint8 **m2;
	*/

   //m= LoadPGM_ui8matrix("hall000000.pgm",&nrl,&nrh,&ncl,&nch);
	
	 //MLoadPGM_ui8matrix("hall000291.pgm",nrl,nrh,ncl,nch,m);
	
	//SavePGM_ui8matrix(m,nrl,nrh,ncl,nch,"test");
	    
	//char * filename [NB]= {"hall/hall000000.pgm","hall/hall000291.pgm"};
	//char nom[4];
	//char res[20]="test/";
	/*for (int i=0;i<2;i++)
	{
		m= LoadPGM_ui8matrix(filename[i],&nrl,&nrh,&ncl,&nch);
		sprintf(nom,"%d",i);
		SavePGM_ui8matrix(m,nrl,nrh,ncl,nch, strcat (res,nom));
		strcpy(res,"test/");
	}*/

  
    //test_routineFD(SEUILFD);
    //test_routineSD();
   
   
    //test_routineFDmorpho3xOuv(SEUILFD);
    //test_routineFDmorpho3xFerm(SEUILFD);
    //test_routineFDmorpho3xOuvFerm(SEUILFD);
    //test_routineFDmorpho3xFermOuv(SEUILFD);
    
    //test_routineSDmorpho3xOuv();
    //test_routineSDmorpho3xFerm();
    //test_routineSDmorpho3xOuvFerm();
    //test_routineSDmorpho3xFermOuv();
    
    //matriceROC(testFD);
    //matriceROC(testFDmorphoO);
    //matriceROC(testFDmorphoF);
    //matriceROC(testFDmorphoOF);
    //matriceROC(testFDmorphoFO);
    
    //matriceROC(testSD);
    //matriceROC(testSDmorphoO);
    //matriceROC(testSDmorphoF);
    //matriceROC(testSDmorphoOF);
    //matriceROC(testSDmorphoFO);
    
    long nrl, nrh, ncl, nch;
    
	uint8 **m;
    int i,j;
    
    m= LoadPGM_ui8matrix("test.pgm",&nrl,&nrh,&ncl,&nch);
    uint8 **bord=ui8matrix(nrl-2,nrh+2,ncl-2,nch+2);
    uint8 **m2=ui8matrix(nrl-2,nrh+2,ncl-2,nch+2);
    for(i=nrl;i<=nrh;i++){
        for(j=ncl;j<=nch;j++){
            bord[i][j]=m[i][j];
        }
    }
    erosion3(bord,m2,nrl,nrh,ncl,nch);
    SavePGM_ui8matrix(m2,nrl,nrh,ncl,nch,"testero.pgm");


    return 0;
}



