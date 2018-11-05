#include <stdio.h>
#include "nrutil.h"
#include "nrdef.h"
#include "string.h"
#define nb 2


int main()
{
    long nrl, nrh, ncl, nch;
	uint8 **m;

   //m= LoadPGM_ui8matrix("hall000000.pgm",&nrl,&nrh,&ncl,&nch);
	
	 //MLoadPGM_ui8matrix("hall000291.pgm",nrl,nrh,ncl,nch,m);
	
	//SavePGM_ui8matrix(m,nrl,nrh,ncl,nch,"test");
	    
	char * filename [nb]= {"hall/hall000000.pgm","hall/hall000291.pgm"};
	char nom[4];
	char res[20]="test/";
	for (int i=0;i<2;i++)
	{
		m= LoadPGM_ui8matrix(filename[i],&nrl,&nrh,&ncl,&nch);
		sprintf(nom,"%d",i);
		SavePGM_ui8matrix(m,nrl,nrh,ncl,nch, strcat (res,nom));
		strcpy(res,"test/");
	}

    return 0;
}



