#include <stdio.h>
#include "nrutil.h"
#include "nrdef.h"
#include "string.h"
#include "mouvement.h"
#include "test_mouvement.h"
#include "test_morpho.h"

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

  
    test_routineFD(SEUILFD);
    test_routineSD();
   
   
    test_routineFDmorpho3xOuv(SEUILFD);
    test_routineFDmorpho3xFerm(SEUILFD);
    test_routineFDmorpho3xOuvFerm(SEUILFD);
    test_routineFDmorpho3xFermOuv(SEUILFD);
    
    test_routineSDmorpho3xOuv();
    test_routineSDmorpho3xFerm();
    test_routineSDmorpho3xOuvFerm();
    test_routineSDmorpho3xFermOuv();
    


    return 0;
}



