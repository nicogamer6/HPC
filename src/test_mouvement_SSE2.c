#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "vnrdef.h"
#include "nrdef.h"
#include "vnrutil.h"

#include "mouvement_SSE2.h"
#include "test_mouvement_SSE2.h"

#define NBIMAGES 299


void uint_to_vuint(uint8 ** scal, vuint ** vect, int vi0, int vi1, int vj0, int vj1){
    vuint8 tmp[1];
    uint8 *p = (uint8*) tmp;
    for (int i = vi0 ; i < vi1 + 1 ; i++){
        for (int j = vj0 ; j < vj1 + 1 ; j++){
            for (int k = 0 ; k < 16; k++){
                p[k] = scal[i][j*16+k];
            }
            vect[i][j]= tmp[0];
        }
    }
}

void vuint_to_uint(uint8 ** scal, vuint ** vect, int vi0, int vi1, int vj0, int vj1){
    vuint8 tmp[1];
    vuint8 x;
    uint8 * p = (uint8*) tmp;
    for (int i = vi0 ; i < vi1 + 1 ; i++){
        for (int j = vj0 ; j < vj1 + 1 ; j++){
            x = _mm_load_si128((vuint*)&vect[i][j]);
            _mm_store_si128(tmp, x);
            for (int k = 0 ; k < 16; k++){
                scal[i][j*16+k]= p[k];
            }
        }
    }
}





//Permet de remplir le dossier /testFD et de voir les nouvelles images avec algo FD sur l'ensemble du dossier /hall
void test_routineFD_SSE(int seuil)
{
    
    //printf("TEST1");
    long nrl, nrh, ncl, nch;
    int vi0, vi1, vj0, vj1;
    char nameload1[100], nameload2[100], namesave[100];


    sprintf(nameload1,"hall/hall000000.pgm");
    
    
    
    uint8 **Itm1 = LoadPGM_ui8matrix(nameload1,&nrl,&nrh,&ncl,&nch);
    s2v(nrl, nrh, ncl, nch, card_vuint8(), &vi0, &vi1, &vj0, &vj1);

    
    vuint8 ** m = vui8matrix(nrl,nrh,ncl,nch);
    vuint8 ** It = vui8matrix(nrl,nrh,ncl,nch);
    vuint8 ** I0 = vui8matrix(nrl,nrh,ncl,nch);
    uint_to_vuint(Itm1, I0, vi0, vi1, vj0, vj1);
    uint8** a = ui8matrix(nrl, nrh, ncl, nch);
    //uint8 **m, **m1 , **m2;
    //uint8** imaget1 = ui8matrix(nrl, nrh, ncl, nch);
    uint8** out = ui8matrix(nrl, nrh, ncl, nch);
    
    
    int i;
   
    for(i=1;i<=NBIMAGES;i++){
        sprintf(nameload2,"hall/hall000%03d.pgm",i);
        
        a = LoadPGM_ui8matrix(nameload2,&nrl,&nrh,&ncl,&nch);
        //todo convertir uint** en vuint** (m)
        

        //todo convertir m en uint **
        uint_to_vuint(a, It, vi0, vi1, vj0, vj1);
        routine_FrameDifference_SSE2(I0,It,m,vi0, vi1,vj0,vj1,seuil);
        vuint_to_uint(out, m, vi0, vi1, vj0, vj1);

        
        sprintf(namesave,"testFD_SSE/hall000%03d.pgm",i);
        SavePGM_ui8matrix(out,nrl,nrh,ncl,nch,namesave);
        dup_vui8matrix(It, vi0, vi1, vj0, vj1, I0);
    }
    
   // free_ui8matrix(Itm1,nrl,nrh,ncl,nch);
    //free_ui8matrix(imaget1,nrl,nrh,ncl,nch);
    //free_ui8matrix(out,nrl,nrh,ncl,nch);
}
