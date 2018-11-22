#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "nrdef.h"
#include "morpho.h"


#define BORD 2

//Pour chaque pixel on va comparer la valeur du voisinnage du pixel avec 255
//Faire attention aux BORD pour pas que erosion ou dilatation 3 ou 5 accède pas à un pixel voisin hors de la matrice 
void erosion3(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch){  
    int i,j; // pour tous les pixels
    int m,n; // pour le voisinnage
    uint8 val;
    for(i=nrl;i<=nrh;i++){
        for(j=ncl;j<=nch;j++){
            val = 255;
            for(m=i-1;m<=i+1;m++){
                for(n=j-1;n<=j+1;n++){
                    val &= Et[m][n];
                }
            }
            EtE[i][j]=val;
        }
    }
    //printf("nrl,nrh %d %d",ncl,nch);
}

void dilatation3(uint8 ** Et, uint8 **EtD, long nrl, long nrh, long ncl, long nch){
    int i,j; // pour tous les pixels
    int m,n; // pour le voisinnage
    uint8 val;
    for(i=nrl;i<=nrh;i++){
        for(j=ncl;j<=nch;j++){
            val = 0;
            for(m=i-1;m<=i+1;m++){
                for(n=j-1;n<=j+1;n++){
                    val |= Et[m][n];
                }
            }
            EtD[i][j]=val;
        }
    }
}

void ouverture3(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch){
    uint8 ** tmp = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    erosion3(Et, tmp, nrl, nrh, ncl, nch);
    dilatation3(tmp, Etout, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}

void fermeture3(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch){
    uint8 ** tmp = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    dilatation3(Et, tmp, nrl, nrh, ncl, nch);
    erosion3(tmp, Etout, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}


void erosion5(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch){
    int i,j; // pour tous les pixels
    int m,n; // pour le voisinnage
    uint8 val;
    for(i=nrl;i<=nrh;i++){
        for(j=ncl;j<=nch;j++){
            val = 255;
            for(m=i-2;m<=i+2;m++){
                for(n=j-2;n<=j+2;n++){
                    val &= Et[m][n];
                }
            }
            EtE[i][j]=val;
        }
    }
}

void dilatation5(uint8 ** Et, uint8 **EtD, long nrl, long nrh, long ncl, long nch){
    int i,j; // pour tous les pixels
    int m,n; // pour le voisinnage
    uint8 val;
    for(i=nrl;i<=nrh;i++){
        for(j=ncl;j<=nch;j++){
            val = 0;
            for(m=i-2;m<=i+2;m++){
                for(n=j-2;n<=j+2;n++){
                    val |= Et[m][n];
                }
            }
            EtD[i][j]=val;
        }
    }
}

void ouverture5(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch){
    uint8 ** tmp = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    erosion5(Et, tmp, nrl, nrh, ncl, nch);
    dilatation5(tmp, Etout, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}

void fermeture5(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch){
    uint8 ** tmp = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    dilatation5(Et, tmp, nrl, nrh, ncl, nch);
    erosion5(tmp, Etout, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}



