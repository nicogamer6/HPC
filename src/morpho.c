#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "nrdef.h"
#include "morpho.h"

//#define NBBITS 64
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


void erosion3_bin(ulong64 ** Et, ulong64 **EtE, long nrl, long nrh, long ncl, long nch){
	int i,j;
	ulong64 res;
	ulong64 leftc, rightc; // Colonne gauche et droite i.e ulong64 avant et après


	for(i=nrl;i<=nrh;i++){
		for(j=ncl;j<=nch;j++){
			res = 0;
			res = ~res;
			res &= (Et[i-1][j] & Et[i][j] & Et[i+1][j]);

			res &= (res>>1) & (res <<1);

			//Manque la gestion des fin de ulong64 pour faire l'ero avec les points du ulong64 suivants

			EtE[i][j] = res;
		}
	}
}


// Suppression des boucles pour le voisinnage et déroulage de 3 avec RR
//  val0    val1    val2
//  val3    val4    val5
//  val6    val7    val8

void erosion3_opti_lu_rr(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch){
    int i = nrl, j=ncl; // pour tous les pixels

    uint8 val, valdefaut = 255;
    uint8 val_lig1, val_lig2, val_lig3;
    uint8 val0, val1, val2, val3, val4, val5, val6, val7, val8;
    
    int r = (nrh+1-nrl) % 3;
    
    for(j=ncl;j<=nch;j++){
    	i=nrl;
    	val0 = Et[i-1][j-1];
    	val1 = Et[i-1][j];
    	val2 = Et[i-1][j+1];

    	val3 = Et[i][j-1];
    	val4 = Et[i][j];
    	val5 = Et[i][j+1];

    	val_lig1 = val0 & val1 & val2;
    	val_lig2 = val3 & val4 & val5;

        for(i=nrl;i<=(nrh-r);i+=3){
            
            val6 = Et[i+1][j-1];
            val7 = Et[i+1][j];
            val8 = Et[i+1][j+1];
            
            val_lig3 = val6 & val7 & val8;
            val = valdefaut & val_lig1 & val_lig2 & val_lig3;
            
            EtE[i][j]=val;
            
            //Rotation de variables
            val0 = val3;
            val1 = val4;
            val2 = val5; 
            
            val3 = val6; 
            val4 = val7;
            val5 = val8;
            
            val_lig1 = val_lig2;
            val_lig2 = val_lig3;


            val6 = Et[i+2][j-1];
            val7 = Et[i+2][j];
            val8 = Et[i+2][j+1];
            
            val_lig3 = val6 & val7 & val8;

            val = valdefaut & val_lig1 & val_lig2 & val_lig3;
            
            EtE[i+1][j]=val;
            
            // RR
            val0 = val3;
            val1 = val4;
            val2 = val5;

            val3 = val6;
            val4 = val7;
            val5 = val8;

            val_lig1 = val_lig2;
            val_lig2 = val_lig3;


            val6 = Et[i+3][j-1];
			val7 = Et[i+3][j];
			val8 = Et[i+3][j+1];

			val_lig3 = val6 & val7 & val8;

			val = valdefaut & val_lig1 & val_lig2 & val_lig3;

			EtE[i+2][j]=val;

			// RR
			val0 = val3;
			val1 = val4;
			val2 = val5;

			val3 = val6;
			val4 = val7;
			val5 = val8;

			val_lig1 = val_lig2;
			val_lig2 = val_lig3;

        }
        switch(r){
                case 0 : 	break;
                case 1 : 	val6 = Et[i+1][j-1];
							val7 = Et[i+1][j];
							val8 = Et[i+1][j+1];

							val_lig3 = val6 & val7 & val8;
							val = valdefaut & val_lig1 & val_lig2 & val_lig3;

							EtE[i][j]=val;

							//Rotation de variables
							val0 = val3;
							val1 = val4;
							val2 = val5;

							val3 = val6;
							val4 = val7;
							val5 = val8;

							val_lig1 = val_lig2;
							val_lig2 = val_lig3;
							break;

                case 2 : 	val6 = Et[i+1][j-1];
							val7 = Et[i+1][j];
							val8 = Et[i+1][j+1];

							val_lig3 = val6 & val7 & val8;
							val = valdefaut & val_lig1 & val_lig2 & val_lig3;

							EtE[i][j]=val;

							//Rotation de variables
							val0 = val3;
							val1 = val4;
							val2 = val5;

							val3 = val6;
							val4 = val7;
							val5 = val8;

							val_lig1 = val_lig2;
							val_lig2 = val_lig3;


							val6 = Et[i+2][j-1];
							val7 = Et[i+2][j];
							val8 = Et[i+2][j+1];

							val_lig3 = val6 & val7 & val8;

							val = valdefaut & val_lig1 & val_lig2 & val_lig3;

							EtE[i+1][j]=val;

							// RR
							val0 = val3;
							val1 = val4;
							val2 = val5;

							val3 = val6;
							val4 = val7;
							val5 = val8;

							val_lig1 = val_lig2;
							val_lig2 = val_lig3;
							break;
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

void dilatation3_bin(ulong64 ** Et, ulong64 **EtD, long nrl, long nrh, long ncl, long nch){
	int i,j;
	ulong64 res;
	ulong64 leftc, rightc; // Colonne gauche et droite i.e ulong64 avant et après
	for(i=nrl;i<=nrh;i++){
		for(j=ncl;j<=nch;j++){
			res = 0;
			res |= (Et[i-1][j] | Et[i][j] | Et[i+1][j]);

			res |= (res>>1) | (res <<1);

			//Manque la gestion des fin de ulong64 pour faire l'ero avec les points du ulong64 suivants

			EtD[i][j] = res;
		}
	}
}

void dilatation3_opti_lu_rr(uint8 ** Et, uint8 **EtD, long nrl, long nrh, long ncl, long nch){
    int i = nrl, j=ncl; // pour tous les pixels

    uint8 val, valdefaut = 0;
    uint8 val_lig1, val_lig2, val_lig3;
    uint8 val0, val1, val2, val3, val4, val5, val6, val7, val8;

    int r = (nrh+1-nrl) % 3;

    for(j=ncl;j<=nch;j++){
    	i=nrl;
    	val0 = Et[i-1][j-1];
    	val1 = Et[i-1][j];
    	val2 = Et[i-1][j+1];

    	val3 = Et[i][j-1];
    	val4 = Et[i][j];
    	val5 = Et[i][j+1];

    	val_lig1 = val0 | val1 | val2;
    	val_lig2 = val3 | val4 | val5;

        for(i=nrl;i<=(nrh-r);i+=3){

            val6 = Et[i+1][j-1];
            val7 = Et[i+1][j];
            val8 = Et[i+1][j+1];

            val_lig3 = val6 | val7 | val8;
            val = valdefaut | val_lig1 | val_lig2 | val_lig3;

            EtD[i][j]=val;

            //Rotation de variables
            val0 = val3;
            val1 = val4;
            val2 = val5;

            val3 = val6;
            val4 = val7;
            val5 = val8;

            val_lig1 = val_lig2;
            val_lig2 = val_lig3;


            val6 = Et[i+2][j-1];
            val7 = Et[i+2][j];
            val8 = Et[i+2][j+1];

            val_lig3 = val6 | val7 | val8;

            val = valdefaut | val_lig1 | val_lig2 | val_lig3;

            EtD[i+1][j]=val;

            // RR
            val0 = val3;
            val1 = val4;
            val2 = val5;

            val3 = val6;
            val4 = val7;
            val5 = val8;

            val_lig1 = val_lig2;
            val_lig2 = val_lig3;


            val6 = Et[i+3][j-1];
			val7 = Et[i+3][j];
			val8 = Et[i+3][j+1];

			val_lig3 = val6 | val7 | val8;

			val = valdefaut | val_lig1 | val_lig2 | val_lig3;

			EtD[i+2][j]=val;

			// RR
			val0 = val3;
			val1 = val4;
			val2 = val5;

			val3 = val6;
			val4 = val7;
			val5 = val8;

			val_lig1 = val_lig2;
			val_lig2 = val_lig3;

        }
        switch(r){
                case 0 : 	break;
                case 1 : 	val6 = Et[i+1][j-1];
							val7 = Et[i+1][j];
							val8 = Et[i+1][j+1];

							val_lig3 = val6 | val7 | val8;
							val = valdefaut | val_lig1 | val_lig2 | val_lig3;

							EtD[i][j]=val;

							//Rotation de variables
							val0 = val3;
							val1 = val4;
							val2 = val5;

							val3 = val6;
							val4 = val7;
							val5 = val8;

							val_lig1 = val_lig2;
							val_lig2 = val_lig3;
							break;

				case 2 : 	val6 = Et[i+1][j-1];
							val7 = Et[i+1][j];
							val8 = Et[i+1][j+1];

							val_lig3 = val6 | val7 | val8;
							val = valdefaut | val_lig1 | val_lig2 | val_lig3;

							EtD[i][j]=val;

							//Rotation de variables
							val0 = val3;
							val1 = val4;
							val2 = val5;

							val3 = val6;
							val4 = val7;
							val5 = val8;

							val_lig1 = val_lig2;
							val_lig2 = val_lig3;


							val6 = Et[i+2][j-1];
							val7 = Et[i+2][j];
							val8 = Et[i+2][j+1];

							val_lig3 = val6 | val7 | val8;
							val = valdefaut | val_lig1 | val_lig2 | val_lig3;

							EtD[i+1][j]=val;

							// RR
							val0 = val3;
							val1 = val4;
							val2 = val5;

							val3 = val6;
							val4 = val7;
							val5 = val8;

							val_lig1 = val_lig2;
							val_lig2 = val_lig3;
							break;
        }
    }

}


//Fonctions pour effectuer ouv et ferm pipe pour les enchaîner plus rapidement
void erosion3_column(uint8 ** Et, uint8 **EtE, int j, long nrl, long nrh, long ncl, long nch){
	int i = nrl; // pour tous les pixels

	uint8 val, valdefaut = 255;
	uint8 val_lig1, val_lig2, val_lig3;
	uint8 val0, val1, val2, val3, val4, val5, val6, val7, val8;

	int r = (nrh+1-nrl) % 3;

	val0 = Et[i-1][j-1];
	val1 = Et[i-1][j];
	val2 = Et[i-1][j+1];

	val3 = Et[i][j-1];
	val4 = Et[i][j];
	val5 = Et[i][j+1];

	val_lig1 = val0 & val1 & val2;
	val_lig2 = val3 & val4 & val5;

	for(i=nrl;i<=(nrh-r);i+=3){

		val6 = Et[i+1][j-1];
		val7 = Et[i+1][j];
		val8 = Et[i+1][j+1];

		val_lig3 = val6 & val7 & val8;
		val = valdefaut & val_lig1 & val_lig2 & val_lig3;

		EtE[i][j]=val;

		//Rotation de variables
		val0 = val3;
		val1 = val4;
		val2 = val5;

		val3 = val6;
		val4 = val7;
		val5 = val8;

		val_lig1 = val_lig2;
		val_lig2 = val_lig3;


		val6 = Et[i+2][j-1];
		val7 = Et[i+2][j];
		val8 = Et[i+2][j+1];

		val_lig3 = val6 & val7 & val8;

		val = valdefaut & val_lig1 & val_lig2 & val_lig3;

		EtE[i+1][j]=val;

		// RR
		val0 = val3;
		val1 = val4;
		val2 = val5;

		val3 = val6;
		val4 = val7;
		val5 = val8;

		val_lig1 = val_lig2;
		val_lig2 = val_lig3;


		val6 = Et[i+3][j-1];
		val7 = Et[i+3][j];
		val8 = Et[i+3][j+1];

		val_lig3 = val6 & val7 & val8;

		val = valdefaut & val_lig1 & val_lig2 & val_lig3;

		EtE[i+2][j]=val;

		// RR
		val0 = val3;
		val1 = val4;
		val2 = val5;

		val3 = val6;
		val4 = val7;
		val5 = val8;

		val_lig1 = val_lig2;
		val_lig2 = val_lig3;

	}
	switch(r){
			case 0 : 	break;
			case 1 : 	val6 = Et[i+1][j-1];
						val7 = Et[i+1][j];
						val8 = Et[i+1][j+1];

						val_lig3 = val6 & val7 & val8;
						val = valdefaut & val_lig1 & val_lig2 & val_lig3;

						EtE[i][j]=val;

						//Rotation de variables
						val0 = val3;
						val1 = val4;
						val2 = val5;

						val3 = val6;
						val4 = val7;
						val5 = val8;

						val_lig1 = val_lig2;
						val_lig2 = val_lig3;
						break;

			case 2 : 	val6 = Et[i+1][j-1];
						val7 = Et[i+1][j];
						val8 = Et[i+1][j+1];

						val_lig3 = val6 & val7 & val8;
						val = valdefaut & val_lig1 & val_lig2 & val_lig3;

						EtE[i][j]=val;

						//Rotation de variables
						val0 = val3;
						val1 = val4;
						val2 = val5;

						val3 = val6;
						val4 = val7;
						val5 = val8;

						val_lig1 = val_lig2;
						val_lig2 = val_lig3;


						val6 = Et[i+2][j-1];
						val7 = Et[i+2][j];
						val8 = Et[i+2][j+1];

						val_lig3 = val6 & val7 & val8;
						val = valdefaut & val_lig1 & val_lig2 & val_lig3;

						EtE[i+1][j]=val;

						// RR
						val0 = val3;
						val1 = val4;
						val2 = val5;

						val3 = val6;
						val4 = val7;
						val5 = val8;

						val_lig1 = val_lig2;
						val_lig2 = val_lig3;
						break;
	}


}

void dilatation3_column(uint8 ** Et, uint8 **EtD, int j, long nrl, long nrh, long ncl, long nch){
	int i = nrl; // pour tous les pixels

	uint8 val, valdefaut = 0;
	uint8 val_lig1, val_lig2, val_lig3;
	uint8 val0, val1, val2, val3, val4, val5, val6, val7, val8;

	int r = (nrh+1-nrl) % 3;

	val0 = Et[i-1][j-1];
	val1 = Et[i-1][j];
	val2 = Et[i-1][j+1];

	val3 = Et[i][j-1];
	val4 = Et[i][j];
	val5 = Et[i][j+1];

	val_lig1 = val0 | val1 | val2;
	val_lig2 = val3 | val4 | val5;

	for(i=nrl;i<=(nrh-r);i+=3){

		val6 = Et[i+1][j-1];
		val7 = Et[i+1][j];
		val8 = Et[i+1][j+1];

		val_lig3 = val6 | val7 | val8;
		val = valdefaut | val_lig1 | val_lig2 | val_lig3;

		EtD[i][j]=val;

		//Rotation de variables
		val0 = val3;
		val1 = val4;
		val2 = val5;

		val3 = val6;
		val4 = val7;
		val5 = val8;

		val_lig1 = val_lig2;
		val_lig2 = val_lig3;


		val6 = Et[i+2][j-1];
		val7 = Et[i+2][j];
		val8 = Et[i+2][j+1];

		val_lig3 = val6 | val7 | val8;

		val = valdefaut | val_lig1 | val_lig2 | val_lig3;

		EtD[i+1][j]=val;

		// RR
		val0 = val3;
		val1 = val4;
		val2 = val5;

		val3 = val6;
		val4 = val7;
		val5 = val8;

		val_lig1 = val_lig2;
		val_lig2 = val_lig3;


		val6 = Et[i+3][j-1];
		val7 = Et[i+3][j];
		val8 = Et[i+3][j+1];

		val_lig3 = val6 | val7 | val8;

		val = valdefaut | val_lig1 | val_lig2 | val_lig3;

		EtD[i+2][j]=val;

		// RR
		val0 = val3;
		val1 = val4;
		val2 = val5;

		val3 = val6;
		val4 = val7;
		val5 = val8;

		val_lig1 = val_lig2;
		val_lig2 = val_lig3;

	}
	switch(r){
			case 0 : 	break;
			case 1 : 	val6 = Et[i+1][j-1];
						val7 = Et[i+1][j];
						val8 = Et[i+1][j+1];

						val_lig3 = val6 | val7 | val8;
						val = valdefaut | val_lig1 | val_lig2 | val_lig3;

						EtD[i][j]=val;

						//Rotation de variables
						val0 = val3;
						val1 = val4;
						val2 = val5;

						val3 = val6;
						val4 = val7;
						val5 = val8;

						val_lig1 = val_lig2;
						val_lig2 = val_lig3;
						break;

			case 2 : 	val6 = Et[i+1][j-1];
						val7 = Et[i+1][j];
						val8 = Et[i+1][j+1];

						val_lig3 = val6 | val7 | val8;
						val = valdefaut | val_lig1 | val_lig2 | val_lig3;

						EtD[i][j]=val;

						//Rotation de variables
						val0 = val3;
						val1 = val4;
						val2 = val5;

						val3 = val6;
						val4 = val7;
						val5 = val8;

						val_lig1 = val_lig2;
						val_lig2 = val_lig3;


						val6 = Et[i+2][j-1];
						val7 = Et[i+2][j];
						val8 = Et[i+2][j+1];

						val_lig3 = val6 | val7 | val8;
						val = valdefaut | val_lig1 | val_lig2 | val_lig3;

						EtD[i+1][j]=val;

						// RR
						val0 = val3;
						val1 = val4;
						val2 = val5;

						val3 = val6;
						val4 = val7;
						val5 = val8;

						val_lig1 = val_lig2;
						val_lig2 = val_lig3;
						break;
	}
}

// Cette version est lente, c'est un "sligg", la refaire pligs rapide au niveau accès mémoire de ui8matrix
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

void ouverture3_pipe(uint8 ** Et, uint8 ** tmp, uint8 **Etout, long nrl, long nrh, long ncl, long nch){

	int i=nrl, j=ncl;

	erosion3_column(Et, tmp, j, nrl, nrh, ncl, nch); //1ère colonne

	for(j=ncl;j<=nch-1;j++){
		erosion3_column(Et, tmp, j+1, nrl, nrh, ncl, nch);
		dilatation3_column(tmp, Etout, j, nrl, nrh, ncl, nch);
	}
	dilatation3_column(tmp, Etout, j+1, nrl, nrh, ncl, nch); //Dernière colonne
}

void ouverture3_bin(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch){
	erosion3_bin(Et, tmp, nrl, nrh, ncl, nch);
	dilatation3_bin(tmp, Etout, nrl, nrh, ncl, nch);
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

void fermeture3_pipe(uint8 ** Et, uint8 ** tmp, uint8 **Etout, long nrl, long nrh, long ncl, long nch){

	int i=nrl, j=ncl;

	dilatation3_column(Et, tmp, j, nrl, nrh, ncl, nch); //1ère colonne

	for(j=ncl;j<=nch-1;j++){
		dilatation3_column(Et, tmp, j+1, nrl, nrh, ncl, nch);
		erosion3_column(tmp, Etout, j, nrl, nrh, ncl, nch);
	}
	erosion3_column(tmp, Etout, j+1, nrl, nrh, ncl, nch); //Dernière colonne

}

void fermeture3_bin(ulong64 ** Et, ulong64 ** tmp, ulong64 **Etout, long nrl, long nrh, long ncl, long nch){
	dilatation3_bin(Et, tmp, nrl, nrh, ncl, nch);
	erosion3_bin(tmp, Etout, nrl, nrh, ncl, nch);
}





/////////////////////////////////////////////////////////////////////////////
/////////////////////////////PARTIE MORPHO x5////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


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


// Suppression des boucles pour le voisinnage et déroulage de 3 avec RR
//	j-2		j-1		j		j+1		j+2
//  val0    val1    val2	val3	val4	i-2
//  val5	val6	val7	val8	val9	i-1
//  val10	val11	val12	val13	val14	i
//	val15	val16	val17	val18	val19	i+1
//	val20	val21	val22	val23	val24	i+2


void erosion5_opti_lu_rr(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch){
    int i = nrl, j=ncl; // pour tous les pixels

    uint8 val, valdefaut = 255;
    uint8 val_lig1, val_lig2, val_lig3, val_lig4, val_lig5;

    uint8 val0, val1, val2, val3, val4, val5, val6, val7, val8, val9;
    uint8 val10, val11, val12, val13, val14, val15, val16, val17, val18, val19;
    uint8 val20, val21, val22, val23, val24;

    int r = (nrh+1-nrl) % 3;

    for(j=ncl;j<=nch;j++){
    	i=nrl;
    	val0 = Et[i-2][j-2];
    	val1 = Et[i-2][j-1];
    	val2 = Et[i-2][j];
    	val3 = Et[i-2][j+1];
		val4 = Et[i-2][j+2];

		val5 = Et[i-1][j-2];
		val6 = Et[i-1][j-1];
		val7 = Et[i-1][j];
		val8 = Et[i-1][j+1];
		val9 = Et[i-1][j+2];

		val10 = Et[i][j-2];
		val11 = Et[i][j-1];
		val12 = Et[i][j];
		val13 = Et[i][j+1];
		val14 = Et[i][j+2];

		val15 = Et[i+1][j-2];
		val16 = Et[i+1][j-1];
		val17 = Et[i+1][j];
		val18 = Et[i+1][j+1];
		val19 = Et[i+1][j+2];

    	val_lig1 = val0 & val1 & val2 & val3 & val4;
    	val_lig2 = val5 & val6 & val7 & val8 & val9;
    	val_lig3 = val10 & val11 & val12 & val13 & val14;
    	val_lig4 = val15 & val16 & val17 & val18 & val19;

        for(i=nrl;i<=(nrh-r);i+=3){

        	val20 = Et[i+2][j-2];
			val21 = Et[i+2][j-1];
			val22 = Et[i+2][j];
			val23 = Et[i+2][j+1];
			val24 = Et[i+2][j+2];

            val_lig5 = val20 & val21 & val22 & val23 & val24;
            val = valdefaut & val_lig1 & val_lig2 & val_lig3 & val_lig4 & val_lig5;

            EtE[i][j]=val;

            //Rotation de variables
            val0 = val5;
			val1 = val6;
			val2 = val7;
			val3 = val8;
			val4 = val9;

			val5 = val10;
			val6 = val11;
			val7 = val12;
			val8 = val13;
			val9 = val14;

			val10 = val15;
			val11 = val16;
			val12 = val17;
			val13 = val18;
			val14 = val19;

            val_lig1 = val_lig2;
            val_lig2 = val_lig3;
            val_lig3 = val_lig4;
            val_lig4 = val_lig5;


            val20 = Et[i+3][j-2];
			val21 = Et[i+3][j-1];
			val22 = Et[i+3][j];
			val23 = Et[i+3][j+1];
			val24 = Et[i+3][j+2];

            val_lig5 = val20 & val21 & val22 & val23 & val24;
            val = valdefaut & val_lig1 & val_lig2 & val_lig3 & val_lig4 & val_lig5;

            EtE[i+1][j]=val;

            // RR
            val0 = val5;
			val1 = val6;
			val2 = val7;
			val3 = val8;
			val4 = val9;

			val5 = val10;
			val6 = val11;
			val7 = val12;
			val8 = val13;
			val9 = val14;

			val10 = val15;
			val11 = val16;
			val12 = val17;
			val13 = val18;
			val14 = val19;

            val_lig1 = val_lig2;
			val_lig2 = val_lig3;
			val_lig3 = val_lig4;
			val_lig4 = val_lig5;


			val20 = Et[i+4][j-2];
			val21 = Et[i+4][j-1];
			val22 = Et[i+4][j];
			val23 = Et[i+4][j+1];
			val24 = Et[i+4][j+2];

			val_lig5 = val20 & val21 & val22 & val23 & val24;
			val = valdefaut & val_lig1 & val_lig2 & val_lig3 & val_lig4 & val_lig5;

			EtE[i+2][j]=val;

			// RR
			val0 = val5;
			val1 = val6;
			val2 = val7;
			val3 = val8;
			val4 = val9;

			val5 = val10;
			val6 = val11;
			val7 = val12;
			val8 = val13;
			val9 = val14;

			val10 = val15;
			val11 = val16;
			val12 = val17;
			val13 = val18;
			val14 = val19;

			val_lig1 = val_lig2;
			val_lig2 = val_lig3;
			val_lig3 = val_lig4;
			val_lig4 = val_lig5;

        }
        switch(r){
                case 0 : 	break;
                case 1 : 	val20 = Et[i+2][j-2];
							val21 = Et[i+2][j-1];
							val22 = Et[i+2][j];
							val23 = Et[i+2][j+1];
							val24 = Et[i+2][j+2];

							val_lig5 = val20 & val21 & val22 & val23 & val24;
							val = valdefaut & val_lig1 & val_lig2 & val_lig3 & val_lig4 & val_lig5;

							EtE[i][j]=val;

							//Rotation de variables
							val0 = val5;
							val1 = val6;
							val2 = val7;
							val3 = val8;
							val4 = val9;

							val5 = val10;
							val6 = val11;
							val7 = val12;
							val8 = val13;
							val9 = val14;

							val10 = val15;
							val11 = val16;
							val12 = val17;
							val13 = val18;
							val14 = val19;

							val_lig1 = val_lig2;
							val_lig2 = val_lig3;
							val_lig3 = val_lig4;
							val_lig4 = val_lig5;
							break;

                case 2 : 	val20 = Et[i+2][j-2];
							val21 = Et[i+2][j-1];
							val22 = Et[i+2][j];
							val23 = Et[i+2][j+1];
							val24 = Et[i+2][j+2];

							val_lig5 = val20 & val21 & val22 & val23 & val24;
							val = valdefaut & val_lig1 & val_lig2 & val_lig3 & val_lig4 & val_lig5;

							EtE[i][j]=val;

							//Rotation de variables
							val0 = val5;
							val1 = val6;
							val2 = val7;
							val3 = val8;
							val4 = val9;

							val5 = val10;
							val6 = val11;
							val7 = val12;
							val8 = val13;
							val9 = val14;

							val10 = val15;
							val11 = val16;
							val12 = val17;
							val13 = val18;
							val14 = val19;

							val_lig1 = val_lig2;
							val_lig2 = val_lig3;
							val_lig3 = val_lig4;
							val_lig4 = val_lig5;


							val20 = Et[i+3][j-2];
							val21 = Et[i+3][j-1];
							val22 = Et[i+3][j];
							val23 = Et[i+3][j+1];
							val24 = Et[i+3][j+2];

							val_lig5 = val20 & val21 & val22 & val23 & val24;
							val = valdefaut & val_lig1 & val_lig2 & val_lig3 & val_lig4 & val_lig5;

							EtE[i+1][j]=val;

							// RR
							val0 = val5;
							val1 = val6;
							val2 = val7;
							val3 = val8;
							val4 = val9;

							val5 = val10;
							val6 = val11;
							val7 = val12;
							val8 = val13;
							val9 = val14;

							val10 = val15;
							val11 = val16;
							val12 = val17;
							val13 = val18;
							val14 = val19;

							val_lig1 = val_lig2;
							val_lig2 = val_lig3;
							val_lig3 = val_lig4;
							val_lig4 = val_lig5;
							break;
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

void dilatation5_opti_lu_rr(uint8 ** Et, uint8 **EtE, long nrl, long nrh, long ncl, long nch){
    int i = nrl, j=ncl; // pour tous les pixels

    uint8 val, valdefaut = 0;
    uint8 val_lig1, val_lig2, val_lig3, val_lig4, val_lig5;

    uint8 val0, val1, val2, val3, val4, val5, val6, val7, val8, val9;
    uint8 val10, val11, val12, val13, val14, val15, val16, val17, val18, val19;
    uint8 val20, val21, val22, val23, val24;

    int r = (nrh+1-nrl) % 3;

    for(j=ncl;j<=nch;j++){
    	i=nrl;
    	val0 = Et[i-2][j-2];
    	val1 = Et[i-2][j-1];
    	val2 = Et[i-2][j];
    	val3 = Et[i-2][j+1];
		val4 = Et[i-2][j+2];

		val5 = Et[i-1][j-2];
		val6 = Et[i-1][j-1];
		val7 = Et[i-1][j];
		val8 = Et[i-1][j+1];
		val9 = Et[i-1][j+2];

		val10 = Et[i][j-2];
		val11 = Et[i][j-1];
		val12 = Et[i][j];
		val13 = Et[i][j+1];
		val14 = Et[i][j+2];

		val15 = Et[i+1][j-2];
		val16 = Et[i+1][j-1];
		val17 = Et[i+1][j];
		val18 = Et[i+1][j+1];
		val19 = Et[i+1][j+2];

    	val_lig1 = val0 | val1 | val2 | val3 | val4;
    	val_lig2 = val5 | val6 | val7 | val8 | val9;
    	val_lig3 = val10 | val11 | val12 | val13 | val14;
    	val_lig4 = val15 | val16 | val17 | val18 | val19;

        for(i=nrl;i<=(nrh-r);i+=3){

        	val20 = Et[i+2][j-2];
			val21 = Et[i+2][j-1];
			val22 = Et[i+2][j];
			val23 = Et[i+2][j+1];
			val24 = Et[i+2][j+2];

            val_lig5 = val20 | val21 | val22 | val23 | val24;
            val = valdefaut | val_lig1 | val_lig2 | val_lig3 | val_lig4 | val_lig5;

            EtE[i][j]=val;

            //Rotation de variables
            val0 = val5;
			val1 = val6;
			val2 = val7;
			val3 = val8;
			val4 = val9;

			val5 = val10;
			val6 = val11;
			val7 = val12;
			val8 = val13;
			val9 = val14;

			val10 = val15;
			val11 = val16;
			val12 = val17;
			val13 = val18;
			val14 = val19;

            val_lig1 = val_lig2;
            val_lig2 = val_lig3;
            val_lig3 = val_lig4;
            val_lig4 = val_lig5;


            val20 = Et[i+3][j-2];
			val21 = Et[i+3][j-1];
			val22 = Et[i+3][j];
			val23 = Et[i+3][j+1];
			val24 = Et[i+3][j+2];

            val_lig5 = val20 | val21 | val22 | val23 | val24;
            val = valdefaut | val_lig1 | val_lig2 | val_lig3 | val_lig4 | val_lig5;

            EtE[i+1][j]=val;

            // RR
            val0 = val5;
			val1 = val6;
			val2 = val7;
			val3 = val8;
			val4 = val9;

			val5 = val10;
			val6 = val11;
			val7 = val12;
			val8 = val13;
			val9 = val14;

			val10 = val15;
			val11 = val16;
			val12 = val17;
			val13 = val18;
			val14 = val19;

            val_lig1 = val_lig2;
			val_lig2 = val_lig3;
			val_lig3 = val_lig4;
			val_lig4 = val_lig5;


			val20 = Et[i+4][j-2];
			val21 = Et[i+4][j-1];
			val22 = Et[i+4][j];
			val23 = Et[i+4][j+1];
			val24 = Et[i+4][j+2];

			val_lig5 = val20 | val21 | val22 | val23 | val24;
			val = valdefaut | val_lig1 | val_lig2 | val_lig3 | val_lig4 | val_lig5;

			EtE[i+2][j]=val;

			// RR
			val0 = val5;
			val1 = val6;
			val2 = val7;
			val3 = val8;
			val4 = val9;

			val5 = val10;
			val6 = val11;
			val7 = val12;
			val8 = val13;
			val9 = val14;

			val10 = val15;
			val11 = val16;
			val12 = val17;
			val13 = val18;
			val14 = val19;

			val_lig1 = val_lig2;
			val_lig2 = val_lig3;
			val_lig3 = val_lig4;
			val_lig4 = val_lig5;

        }
        switch(r){
                case 0 : 	break;
                case 1 : 	val20 = Et[i+2][j-2];
							val21 = Et[i+2][j-1];
							val22 = Et[i+2][j];
							val23 = Et[i+2][j+1];
							val24 = Et[i+2][j+2];

							val_lig5 = val20 | val21 | val22 | val23 | val24;
							val = valdefaut | val_lig1 | val_lig2 | val_lig3 | val_lig4 | val_lig5;

							EtE[i][j]=val;

							//Rotation de variables
							val0 = val5;
							val1 = val6;
							val2 = val7;
							val3 = val8;
							val4 = val9;

							val5 = val10;
							val6 = val11;
							val7 = val12;
							val8 = val13;
							val9 = val14;

							val10 = val15;
							val11 = val16;
							val12 = val17;
							val13 = val18;
							val14 = val19;

							val_lig1 = val_lig2;
							val_lig2 = val_lig3;
							val_lig3 = val_lig4;
							val_lig4 = val_lig5;
							break;

                case 2 : 	val20 = Et[i+2][j-2];
							val21 = Et[i+2][j-1];
							val22 = Et[i+2][j];
							val23 = Et[i+2][j+1];
							val24 = Et[i+2][j+2];

							val_lig5 = val20 | val21 | val22 | val23 | val24;
							val = valdefaut | val_lig1 | val_lig2 | val_lig3 | val_lig4 | val_lig5;

							EtE[i][j]=val;

							//Rotation de variables
							val0 = val5;
							val1 = val6;
							val2 = val7;
							val3 = val8;
							val4 = val9;

							val5 = val10;
							val6 = val11;
							val7 = val12;
							val8 = val13;
							val9 = val14;

							val10 = val15;
							val11 = val16;
							val12 = val17;
							val13 = val18;
							val14 = val19;

							val_lig1 = val_lig2;
							val_lig2 = val_lig3;
							val_lig3 = val_lig4;
							val_lig4 = val_lig5;


							val20 = Et[i+3][j-2];
							val21 = Et[i+3][j-1];
							val22 = Et[i+3][j];
							val23 = Et[i+3][j+1];
							val24 = Et[i+3][j+2];

							val_lig5 = val20 | val21 | val22 | val23 | val24;
							val = valdefaut | val_lig1 | val_lig2 | val_lig3 | val_lig4 | val_lig5;

							EtE[i+1][j]=val;

							// RR
							val0 = val5;
							val1 = val6;
							val2 = val7;
							val3 = val8;
							val4 = val9;

							val5 = val10;
							val6 = val11;
							val7 = val12;
							val8 = val13;
							val9 = val14;

							val10 = val15;
							val11 = val16;
							val12 = val17;
							val13 = val18;
							val14 = val19;

							val_lig1 = val_lig2;
							val_lig2 = val_lig3;
							val_lig3 = val_lig4;
							val_lig4 = val_lig5;
							break;
        }
    }

}


void ouverture5(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch){
    uint8 ** tmp = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    erosion5(Et, tmp, nrl, nrh, ncl, nch);
    dilatation5(tmp, Etout, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}

void ouverture5_opti(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch){
    uint8 ** tmp = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    erosion5_opti_lu_rr(Et, tmp, nrl, nrh, ncl, nch);
    dilatation5_opti_lu_rr(tmp, Etout, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}

void fermeture5(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch){
    uint8 ** tmp = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    dilatation5(Et, tmp, nrl, nrh, ncl, nch);
    erosion5(tmp, Etout, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}

void fermeture5_opti(uint8 ** Et, uint8 **Etout, long nrl, long nrh, long ncl, long nch){
    uint8 ** tmp = ui8matrix(nrl-BORD,nrh+BORD,ncl-BORD,nch+BORD);
    dilatation5_opti_lu_rr(Et, tmp, nrl, nrh, ncl, nch);
    erosion5_opti_lu_rr(tmp, Etout, nrl, nrh, ncl, nch);
    free_ui8matrix(tmp, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}



