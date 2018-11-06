#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nrutil.h"
#include "nrdef.h"
#include "mouvement.h"

uint8** routine_FrameDifference(uint8 **in1, uint8 **in2,  long nrl, long nrh, long ncl, long nch, int seuil){
    uint8 ** res=ui8matrix(nrl,nrh,ncl,nch);
    
    int i,j;
    
    for(i=nrl;i<=nrh;i++){
        for(j=ncl;j<=nch;j++){
            if(abs(in1[i][j]-in2[i][j]) < seuil)
                res[i][j]=255; //blanc
            else res[i][j]=0; //noir
        }
    }
    
    return res;

}

//YO 


//TEST
