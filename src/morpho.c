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


// Suppression des boucles pour le voisinnage et déroulage de 2 avec RR
//  val0    val3    val6
//  val1    val4    val7
//  val2    val5    val8

void erosion3_opti_lu_rr(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch){  
    int i = nrl, j=ncl; // pour tous les pixels

    uint8 val, valdefaut;
    uint8 val0, val1, val2, val3, val4, val5, val6, val7, val8;
    
    int r = (nrh+1-nrl) % 2;
    
    for(j=ncl;j<=nch;j++){
    	i=nrl;
    	val0 = Et[i-1][j-1];
    	val1 = Et[i-1][j];
    	val2 = Et[i-1][j+1];

    	val3 = Et[i][j-1];
    	val4 = Et[i][j];
    	val5 = Et[i][j+1];

        for(i=nrl;i<=(nrh-r);i+=2){
            valdefaut = 255;
            
            val6 = Et[i+1][j-1];
            val7 = Et[i+1][j];
            val8 = Et[i+1][j+1];
            
            val = valdefaut & val0 & val1 & val2 & val3 & val4 & val5 & val6 & val7 & val8;
            
            EtE[i][j]=val;
            
            //Rotation de variables
            val0 = val3;
            val1 = val4;
            val2 = val5; 
            
            val3 = val6; 
            val4 = val7;
            val5 = val8;
            

            val6 = Et[i+2][j-1];
            val7 = Et[i+2][j];
            val8 = Et[i+2][j+1];
            
            val = valdefaut & val0 & val1 & val2 & val3 & val4 & val5 & val6 & val7 & val8;
            
            EtE[i+1][j]=val;
            
            // RR
            val0 = val3;
            val1 = val4;
            val2 = val5;

            val3 = val6;
            val4 = val7;
            val5 = val8;

        }
        switch(r){
                case 0 : 	break;
                case 1 : 	val6 = Et[i+1][j-1];
                			val7 = Et[i+1][j];
                			val8 = Et[i+1][j+1];

                			val = valdefaut & val0 & val1 & val2 & val3 & val4 & val5 & val6 & val7 & val8;

                			EtE[i][j]=val;
        }
    }

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


void dilatation3_opti_lu_rr(uint8 ** Et, uint8 **EtD, long nrl, long nrh, long ncl, long nch){
    int i = nrl, j=ncl; // pour tous les pixels

    uint8 val, valdefaut, val_lu2;
    uint8 val0, val1, val2, val3, val4, val5, val6, val7, val8;

    int r = (nrh+1-nrl) % 2;

    for(j=ncl;j<=nch;j++){
    	i=nrl;
    	val0 = Et[i-1][j-1];
    	val1 = Et[i-1][j];
    	val2 = Et[i-1][j+1];

    	val3 = Et[i][j-1];
    	val4 = Et[i][j];
    	val5 = Et[i][j+1];

    	val = valdefaut | val0 | val1 | val2 | val3 | val4 | val5; //réduction par colonne
    	val_lu2 = valdefaut | val3 | val4 | val5;

        for(i=nrl;i<=(nrh-r);i+=2){
        	valdefaut = 0;

            val6 = Et[i+1][j-1];
            val7 = Et[i+1][j];
            val8 = Et[i+1][j+1];

            val = val | val6 | val7 | val8;

            EtD[i][j]=val;

            //Rotation de variables
            val0 = val3;
            val1 = val4;
            val2 = val5;

            val3 = val6;
            val4 = val7;
            val5 = val8;

            val6 = Et[i+2][j-1];
            val7 = Et[i+2][j];
            val8 = Et[i+2][j+1];

            val = val_lu2 | val3 | val4 | val5 | val6 | val7 | val8;

            EtD[i+1][j]=val;

            // RR
            val0 = val3;
            val1 = val4;
            val2 = val5;

            val3 = val6;
            val4 = val7;
            val5 = val8;

        }
        switch(r){
                case 0 : 	break;
                case 1 : 	val6 = Et[i+1][j-1];
                			val7 = Et[i+1][j];
                			val8 = Et[i+1][j+1];

                			val = valdefaut | val0 | val1 | val2 | val3 | val4 | val5 | val6 | val7 | val8;
                			EtD[i][j]=val;
        }
    }


    //Fonctionne sans le déroulage de 2 !
    /*for(j=ncl;j<=nch;j++){

        	val0 = Et[i-1][j-1];
        	val1 = Et[i-1][j];
        	val2 = Et[i-1][j+1];

        	val3 = Et[i][j-1];
        	val4 = Et[i][j];
        	val5 = Et[i][j+1];

            for(i=nrl;i<=(nrh);i++){
            	valdefaut = 0;

                val6 = Et[i+1][j-1];
                val7 = Et[i+1][j];
                val8 = Et[i+1][j+1];

                val = valdefaut | val0 | val1 | val2 | val3 | val4 | val5 | val6 | val7 | val8;

                EtD[i][j]=val;

                //Rotation de variables
                val0 = val3;
                val1 = val4;
                val2 = val5;

                val3 = val6;
                val4 = val7;
                val5 = val8;
            }
        }*/

}


// Cette version est lente, c'est un "slug", la refaire plus rapide au niveau accès mémoire de ui8matrix
void ouverture3(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch){
    uint8 ** tmp = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    erosion3(Et, tmp, nrl, nrh, ncl, nch);
    dilatation3(tmp, Etout, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}

void ouverture3_opti(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch){
    uint8 ** tmp = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    erosion3_opti_lu_rr(Et, tmp, nrl, nrh, ncl, nch);
    dilatation3_opti_lu_rr(tmp, Etout, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}

void fermeture3(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch){
    uint8 ** tmp = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    dilatation3(Et, tmp, nrl, nrh, ncl, nch);
    erosion3(tmp, Etout, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}

void fermeture3_opti(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch){
    uint8 ** tmp = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    dilatation3_opti_lu_rr(Et, tmp, nrl, nrh, ncl, nch);
    erosion3_opti_lu_rr(tmp, Etout, nrl, nrh, ncl, nch);
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



