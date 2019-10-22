/* --------------- */
/* --- nrdef.h --- */
/* --------------- */


#ifndef __NR_DEF_H__
#define __NR_DEF_H__

#include <stdint.h>

#define NBBITS (sizeof(ulong) * 8)

typedef uint8_t byte;
typedef uint32_t ulong32;
typedef uint64_t ulong64; // Penser à la fin à remettre tous les 32 ou 64 adéquats, pour les fonctions free et long64matrix ou 32

typedef uint64_t ulong;

typedef          char  int8;
typedef unsigned char uint8;
typedef   signed char sint8;

typedef          short  int16;
typedef unsigned short uint16;
typedef   signed short sint16;

typedef          int  int32;
typedef unsigned int uint32;
typedef   signed int sint32;

typedef float         float32;
typedef double        float64;
typedef struct { byte r; byte g; byte b;} rgb8;

#endif // __NR_DEF_H__
