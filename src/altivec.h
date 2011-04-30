#ifndef __ALTIVEC_INCLUDE__
#define __ALTIVEC_INCLUDE__
#ifdef ALTIVEC
#include "migration.h"

// ALTIVEC specific utility functions
// tips from Doug
// from http://home.san.rr.com/altivec/Pages/TipArchive.html
//First we have the two missing logical operators:
#define vec_not(A) vec_nor(A,A)
#define vec_nand(A,B) vec_not(vec_and(A,B))

//Next we have the missing arithmetic operators:
#define vec_abs(A) vec_sel(A,vec_sub( vec_sub(A,A),A),vec_cmpgt(vec_sub(A,A),A))
#define vec_2sComp(x) vec_sub(vec_sub(x,x), x)
//Finally here a two of the missing multiplies:
#define vec_halfmul_16(A,B) vec_mladd(A,B,vec_sub(A,A))
#define vec_halfmul_8(A,B) vec_mergel(vec_pack(vec_mule(A,B),vec_mule(A,B)),vec_pack(vec_mulo(A,B),vec_mulo(A,B)))

typedef union {
    float f[4];
    vector float v;
} FloatVec;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//functions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* load a single float and spread it to all vector elements */
extern MYINLINE vector float vec_float_splat(float *thisptr);
extern MYINLINE vector float vec_zero(void);
extern MYINLINE vector float vec_float_one(void);
extern MYINLINE void load_double_floatvec(FloatVec *out, double *in, long size);
extern MYINLINE void load_float_floatvec(FloatVec *out, float *in, long size);
extern MYINLINE vector float vector_square_root( vector float v );
extern MYINLINE MYREAL vector_multiply(vector float *v1, vector float *v2);
extern MYREAL   vdot_product (FloatVec *v1, FloatVec *v2, int size);

extern MYINLINE FloatVec load_float_splat(float *thisptr);
extern MYINLINE FloatVec load_double_splat(float *thisptr);
extern MYINLINE FloatVec vector_add(FloatVec a, FloatVec b );
#endif /*ALTIVEC*/
#endif /*ALTIVEC_INCLUDE*/
