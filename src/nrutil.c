
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <ctype.h> 
#include <string.h>
#include <math.h> /* fabs */

#include "nrdef.h"
#include "nrutil.h"

//NR_END est maintenant defini dans nrutil.h

//#define NR_END 1
//#define FREE_ARG char*

long nr_end = NR_END;

/* ------------------------ */
/* -- PGM IO for bmatrix -- */
/* ------------------------ */

char *readitem   (FILE *file, char *buffer);
void  ReadPGMrow (FILE *file, long width, uint8  *line);
void  WritePGMrow(uint8 *line, long width, FILE  *file);

/* --------------------------------- */
char *readitem(FILE *file,char *buffer)
/* --------------------------------- */
/* lecture d'un mot */
{
  char *aux;
  int k;

  k=0;
  aux=buffer;
  while (!feof(file))
    {
      *aux=fgetc(file);
      switch(k)
        {
        case 0:
          if (*aux=='#') k=1;
          if (isalnum(*aux)) k=2,aux++;
          break;
        case 1:
          if (*aux==0xA) k=0;
          break;
        case 2:
          if (!isalnum(*aux))
            {
              *aux=0;
              return buffer;
            }
          aux++;
          break;
        }
    }
  *aux=0;
  return buffer;
}
/* ---------------------------------------------- */
void ReadPGMrow(FILE *file, long width, uint8  *line)
/* ---------------------------------------------- */
{
    /* Le fichier est ouvert (en lecture) et ne sera pas ferme a la fin */
     fread(&(line[0]), sizeof(uint8), width, file);
}
/* ----------------------------------------------- */
void WritePGMrow(uint8 *line, long width, FILE  *file)
/* ----------------------------------------------- */
{
/* Le fichier est deja ouvert et ne sera pas ferme a la fin */

   fwrite(&(line[0]), sizeof(uint8), width, file);
}
/* ------------------------------------------------------------------------------ */
uint8** LoadPGM_ui8matrix(char *filename, long *nrl, long *nrh, long *ncl, long *nch)
/* ------------------------------------------------------------------------------ */
{
  /* cette version ne lit plus que le type P5 */

  long height, width, gris;
  uint8 **m;
  FILE *file;
  /*int   format;/**/

  char *buffer;
  /*char  c;/**/
  int i;
  
  buffer = (char*) calloc(80, sizeof(char));
  /* ouverture du fichier */
  file = fopen(filename,"rb");
  if (file==NULL)
    nrerror("ouverture du fichier impossible\n");
    //nrerror("ouverture du fichier %s impossible\n", filename);

  /* lecture de l'entete du fichier pgm */
  readitem(file, buffer);
  /*fscanf(fichier, "%s", buffer);*/
  if(strcmp(buffer, "P5") != 0)
    nrerror("entete du fichier %s invalide\n");
    //nrerror("entete du fichier %s invalide\n", filename);

  width  = atoi(readitem(file, buffer));
  height = atoi(readitem(file, buffer));
  gris   = atoi(readitem(file, buffer));

  *nrl = 0;
  *nrh = height - 1;
  *ncl = 0;
  *nch = width - 1;
  m = ui8matrix(*nrl, *nrh, *ncl, *nch);
  
  for(i=0; i<height; i++) {
    ReadPGMrow(file, width, m[i]);
  }

  fclose(file);
  free(buffer);

  return m;
}
/* ------------------------------------------------------------------------------- */
void MLoadPGM_ui8matrix(char *filename, int nrl, int nrh, int ncl, int nch, uint8 **m)
/* ------------------------------------------------------------------------------- */
{
    /* cette version ne lit plus que le type P5 */
    
    int height, width, gris;
    FILE *file;
    
    char *buffer;
    int i;
    
    buffer = (char*) calloc(80, sizeof(char));
    /* ouverture du fichier */
    file = fopen(filename,"rb");
    if (file==NULL)
        nrerror("ouverture du fichier impossible\n");
    //nrerror("ouverture du fichier %s impossible\n", filename);
    
    /* lecture de l'entete du fichier pgm */
    readitem(file, buffer);
    /*fscanf(fichier, "%s", buffer);*/
    if(strcmp(buffer, "P5") != 0)
        nrerror("entete du fichier %s invalide\n");
    //nrerror("entete du fichier %s invalide\n", filename);
    
    width  = atoi(readitem(file, buffer));
    height = atoi(readitem(file, buffer));
    gris   = atoi(readitem(file, buffer));
    
    for(i=0; i<height; i++) {
        ReadPGMrow(file, width, m[i]);
    }
    
    fclose(file);
    free(buffer);
}
/* ----------------------------------------------------------------------------------- */
void SavePGM_ui8matrix(uint8 **m, long nrl, long nrh, long ncl, long nch, char *filename)
/* ----------------------------------------------------------------------------------- */
{
  long nrow = nrh-nrl+1;
  long ncol = nch-ncl+1;

  char buffer[80];
  
  FILE *file;
  int  i;

  file = fopen(filename, "wb");
  if (file == NULL)
    //nrerror("ouverture du fichier %s impossible dans SavePGM_bmatrix\n", filename);
    nrerror("ouverture du fichier %s impossible dans SavePGM_ui8matrix\n");

  /* enregistrement de l'image au format rpgm */

  sprintf(buffer,"P5\n%d %d\n255\n",ncol, nrow);
  fwrite(buffer,strlen(buffer),1,file);
  for(i=nrl; i<=nrh; i++)
    WritePGMrow(m[i], ncol, file);

  /* fermeture du fichier */
  fclose(file);
}
/* --------------------------- */
/* -- PNM IO for rgb8matrix -- */
/* --------------------------- */

/* ----------------------------------------------- */
void ReadPNMrow(FILE  *file, long width, uint8  *line)
/* ----------------------------------------------- */
{
    /* Le fichier est ouvert (en lecture) et ne sera pas ferme a la fin */
     fread(&(line[0]), sizeof(uint8), 3*sizeof(uint8)*width, file);
}
/* ------------------------------------------------ */
void WritePNMrow(uint8  *line, long width, FILE  *file)
/* ------------------------------------------------ */
{
/* Le fichier est deja ouvert et ne sera pas ferme a la fin */

   fwrite(&(line[0]), sizeof(uint8), 3*sizeof(uint8)*width, file);
}
/* ------------------------------------------------------------------------------- */
rgb8** LoadPPM_rgb8matrix(char *filename, long *nrl, long *nrh, long *ncl, long *nch)
/* ------------------------------------------------------------------------------- */
{
  /* cette version ne lit plus que le type P6 */

  long height, width, gris;
  rgb8 **m;
  FILE *file;
  /*int   format;/**/

  char *buffer;
  /*char  c;/**/
  int i;
  
  buffer = (char*) calloc(80, sizeof(char));
  /* ouverture du fichier */
  file = fopen(filename,"rb");
  if (file==NULL)
    nrerror("ouverture du fichier impossible\n");
    //nrerror("ouverture du fichier %s impossible\n", filename);

  /* lecture de l'entete du fichier pgm */
  readitem(file, buffer);
  /*fscanf(fichier, "%s", buffer);*/
  if(strcmp(buffer, "P6") != 0)
    nrerror("entete du fichier %s invalide\n");
    //nrerror("entete du fichier %s invalide\n", filename);

  width  = atoi(readitem(file, buffer));
  height = atoi(readitem(file, buffer));
  gris   = atoi(readitem(file, buffer));

  *nrl = 0;
  *nrh = height - 1;
  *ncl = 0;
  *nch = width - 1;
  m = rgb8matrix(*nrl, *nrh, *ncl, *nch);
  
  for(i=0; i<height; i++) {
    ReadPNMrow(file, width, (uint8*)m[i]);
  }

  fclose(file);
  free(buffer);

  return m;
}
/* ----------------------------------------------------------------------------------- */
void SavePPM_rgb8matrix(rgb8 **m, long nrl, long nrh, long ncl, long nch, char *filename)
/* ----------------------------------------------------------------------------------- */
{
  long nrow = nrh-nrl+1;
  long ncol = nch-ncl+1;

  char buffer[80];
  
  FILE *file;
  int  i;

  file = fopen(filename, "wb");
  if (file == NULL)
    //nrerror("ouverture du fichier %s impossible dans SavePGM_bmatrix\n", filename);
    nrerror("ouverture du fichier %s impossible dans SavePPM_bmatrix\n");

  /* enregistrement de l'image au format rpgm */

  sprintf(buffer,"P6\n%d %d\n255\n",ncol, nrow);
  fwrite(buffer,strlen(buffer),1,file);
  for(i=nrl; i<=nrh; i++)
    WritePNMrow((uint8*)m[i], ncol, file);

  /* fermeture du fichier */
  fclose(file);
}
