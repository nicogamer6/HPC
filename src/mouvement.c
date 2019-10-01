#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nrutil.h"
#include "nrdef.h"
#include "mouvement.h"

#define VMIN 20
#define VMAX 254
#define N 4

//////////////////////////////
//     FRAME DIFFERENCE	    //
//////////////////////////////


void routine_FrameDifference(uint8 **in1, uint8 **in2, uint8 **res,  long nrl, long nrh, long ncl, long nch, int seuil){
    //uint8 ** res=ui8matrix(nrl,nrh,ncl,nch);
    
    int i,j;
    uint8 Ot;
    
    for(i=nrl;i<=nrh;i++){
        for(j=ncl;j<=nch;j++){
            Ot = abs(in2[i][j]-in1[i][j]);
            if(Ot < seuil)
                res[i][j]=0; //noir
            else res[i][j]=255; //blanc
        }
    }
    
    //return res;

}

void routine_FrameDifference_opti(uint8 **in1, uint8 **in2, uint8 **res,  long nrl, long nrh, long ncl, long nch, int seuil){
    //uint8 ** res=ui8matrix(nrl,nrh,ncl,nch);
    
    int i,j;
    uint8 Ot;
    uint8 a = 0;
    
    for(i=nrl;i<=nrh;i++){
        for(j=ncl;j<=nch;j++){
            Ot = abs(in2[i][j]-in1[i][j]);
            if(Ot < seuil)
                a=0; //noir
            else a=255; //blanc
        }
        res[i][j]=a; //Scalarisation 
    }
    
    //return res;

}

///////////////////////////////////
//	   SIGMA DELTA STEP0	     //
///////////////////////////////////


void routine_SigmaDelta_step0(uint8 **V, uint8 **M, uint8 **I, long nrl, long nrh, long ncl, long nch)
{
	int i,j;
	uint8 vmin = (uint8) VMIN;
	for(i = nrl; i<=nrh;i++)
	{
		for(j=ncl;j<=nch;j++)
		{
			M[i][j] =  I[i][j];
			V[i][j] = vmin;
		}
	
	}
	//return M;
}


///////////////////////////////////
//	   SIGMA DELTA STEP1	     //
///////////////////////////////////


void routine_SigmaDelta_1step(uint8 **V, uint8 **Vtm1, uint8 **M, uint8 **Mtm1, uint8 **I, uint8 **Et, long nrl, long nrh, long ncl, long nch)
{
    int i,j;
    uint8 **Ot=ui8matrix(nrl,nrh,ncl,nch);
    
    uint8 n = (uint8) N;
    uint8 vmax = (uint8) VMAX;
    uint8 vmin = (uint8) VMIN;
    
    for(i = nrl; i<=nrh;i++)
	{
		for(j=ncl;j<=nch;j++)
		{
		    //Step 1 Estimation
			if(Mtm1[i][j] < I[i][j])
			    M[i][j] = Mtm1[i][j]+1;
			else if(Mtm1[i][j] > I[i][j])
			    M[i][j] = Mtm1[i][j]-1;
			else M[i][j] = Mtm1[i][j];
	
	        //Step 2 Difference Computation
		    Ot[i][j]=abs(M[i][j]-I[i][j]);
	
            //Step 3 Update and clamping    
			if(Vtm1[i][j] < (n * Ot[i][j]))
			    V[i][j] = Vtm1[i][j]+1;
			else if(Vtm1[i][j] > (n * Ot[i][j]))
			    V[i][j] = Vtm1[i][j]-1;
			else V[i][j] = Vtm1[i][j];
			//Clamp to [VMIN,VMAX]
			V[i][j]=max(min(V[i][j],vmax),vmin);
	
	        //Step 4 Estimation
			if(Ot[i][j] < V[i][j])
			    Et[i][j] = 0;
			else Et[i][j] = 255; //ou 1
			
		}
	}
	
    /*//Step 1 Estimation
    for(i = nrl; i<=nrh;i++)
	{
		for(j=ncl;j<=nch;j++)
		{
			if(Mtm1[i][j] < I[i][j])
			    M[i][j] = Mtm1[i][j]+1;
			else if(Mtm1[i][j] > I[i][j])
			    M[i][j] = Mtm1[i][j]-1;
			else M[i][j] = Mtm1[i][j];
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
			if(Vtm1[i][j] < (n * Ot[i][j]))
			    V[i][j] = Vtm1[i][j]+1;
			else if(Vtm1[i][j] > n * Ot[i][j])
			    V[i][j] = Vtm1[i][j]-1;
			else V[i][j] = Vtm1[i][j];
			//Clamp to [VMIN,VMAX]
			V[i][j]=max(min(V[i][j],vmax),vmin);
		}
	}
	
	//Step 4 Estimation
	for(i = nrl; i<=nrh;i++)
	{
		for(j=ncl;j<=nch;j++)
		{
			if(Ot[i][j] < V[i][j])
			    Et[i][j] = 0;
			else Et[i][j] = 255; //ou 1
		}
	}*/
	//return Et;
}	 


///////////////////////////////
//	   MIN ET MAX		     //
///////////////////////////////

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

