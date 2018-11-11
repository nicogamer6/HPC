#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "vnrdef.h"
#include "nrdef.h"
#include "mouvement_SSE2.h"
#include "vnrutil.h"



vuint8** routine_FrameDifference_SSE2(vuint8** It, vuint8** It_1, vuint8** Ot, vuint8** Et, long nrl, long nrh, long ncl, long nch)
{
	int i, j;

	vuint8 v_theta = init_vuint8(THETA);
	vuint8 v_255 = init_vuint8(255);
	vuint8 v_128 = init_vuint8(128);

	vuint8 a_0;
	vuint8 b_0;

	vuint8 v_it;
	vuint8 v_it_1;


	for (i=nrl;i<=nrh;i++) {
		for(j=ncl;j=nch;j++) {
		
			// Ici on load les images
			v_it = _mm_load_si128(&It[i][j]);
 			v_it_1 = _mm_load_si128(&It_1[i][j]);
		
			// On fait la valeur absolue
 			a_0 = _mm_min_epu8(v_it_1, v_it);
 			b_0 = _mm_max_epu8(v_it_1, v_it);
 			a_0 = _mm_sub_epi8(b_0, a_0);
 			// On sauvegarde la valeur de abs dans Ot
 			_mm_store_si128(&Ot[i][j], a_0);
		}		

	}

	for (i=nrl;i<=nrh;i++) {
		for(j=ncl;j=nch;j++) {
		
            // si Ot < THETA donc dépasse le seuil alors a_0 reçoit 255 sinon 0
            a_0 = _mm_sub_epi8(Ot[i][j], v_128);
            b_0 = _mm_sub_epi8(v_theta, v_128);
            a_0 = _mm_cmplt_epi8(a_0, b_0);
            // Si a_0 vaut 0 alors b_0 reçoit 255
            b_0 = _mm_andnot_si128(a_0, v_255);
            //  On sauvegarde la valeur de b dans Et
            _mm_store_si128(&Et[i][j], b_0);
		}
	}
		
	return Et;
}

