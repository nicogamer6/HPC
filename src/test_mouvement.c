#include <stdio.h>
#include "nrutil.h"
#include "nrdef.h"
#include "string.h"
#include "test_mouvement.h"
#include "mouvement.h"


void test_routineFD(void)
{
    long nrl, nrh, ncl, nch;
	uint8 **m;
	uint8 **m1;
	uint8 **m2;
	
	m1= LoadPGM_ui8matrix("hall/hall000000.pgm",&nrl,&nrh,&ncl,&nch);
    m2= LoadPGM_ui8matrix("hall/hall000291.pgm",&nrl,&nrh,&ncl,&nch);
     
    m=routine_FrameDifference(m1,m2,nrl,nrh,ncl,nch,20);
    SavePGM_ui8matrix(m,nrl,nrh,ncl,nch,"test/test.pgm");

}
