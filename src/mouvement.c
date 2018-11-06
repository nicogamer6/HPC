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


uint8** SigmaDelta_step0(uint8** V, uint8 ** M	, uint8** I, long nrl, long nrh, long ncl, long nch)
{
	int i,j;
	for(i = nrl; i<=nrh;i++)
	{
		for(j=ncl;j<=nch;j++)
		{
			V[i][j] = VMIN;
			M[i][j] =  I[i][j];
		}
	
	}
	return M;
}

uint8** SigmaDelta_1step(uint8** V, uint8 **M, uint8 **I, uint8 **Et, uint8** Ot, long nrl, long nrh, long ncl, long nch)
{
    int i,j;
    
    //Step 1 Estimation
    for(i = nrl; i<=nrh;i++)
	{
		for(j=ncl;j<=nch;j++)
		{
			if(M[i][j] < I[i][j])
			    M[i][j] = M[i][j]+1;
			else if(M[i][j] > I[i][j])
			    M[i][j] = M[i][j]-1;
			else M[i][j] = M[i][j];
		}
	}
	
	//Step 2 Difference Computation
	for(i = nrl; i<=nrh;i++)
	{
		for(j=ncl;j<=nch;j++)
		{
		    Ot[i][j]=abs(M[i][j]-I[i][j]);
		}
	}
    
    //Step 3 Update and clamping
    for(i = nrl; i<=nrh;i++)
	{
		for(j=ncl;j<=nch;j++)
		{
			if(V[i][j] < N * Ot[i][j])
			    V[i][j] = V[i][j]+1;
			else if(V[i][j] > N * Ot[i][j])
			    V[i][j] = V[i][j]-1;
			else V[i][j] = V[i][j];
			//Clamp to [VMIN,VMAX]
			V[i][j]=max(min(V[i][j],VMAX),VMIN);
		}
	}
	
	//Step 4 Estimation
	for(i = nrl; i<=nrh;i++)
	{
		for(j=ncl;j<=nch;j++)
		{
			if(Ot[i][j] < V[i][j])
			    Et[i][j] = 0;
			else Et[i][j] = 1; //ou 255
		}
	}
}	 


int min(int a, int b)	{
	if(a < b)
		return a;
	else
		return b;
}
int max(int a, int b)	{
	if(a > b)
		return a;
	else
		return b;
}

