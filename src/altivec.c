/*! \file altivec.c */
//
// ALTIVEC specific utility functions
// tips from Doug
// from http://home.san.rr.com/altivec/Pages/TipArchive.html
//First we have the two missing logical operators:
#ifdef ALTIVEC
#include "migration.h"
#include "altivec.h"
#include "sighandler.h"

MYINLINE vector float vector_square_root( vector float v );
MYINLINE vector float vector_divide( vector float a, vector float b );
MYINLINE vector float vector_reciprocal( vector float v );
MYINLINE vector float vector_reciprocal_square_root( vector float v );

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//functions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* load a single float and spread it to all vector elements */
MYINLINE FloatVec load_float_splat(float *thisptr)
{
    vector unsigned char tempc;
    FloatVec temp;
    temp.v = vec_lde(0,thisptr);
    tempc  = vec_lvsl(0,thisptr);
    temp.v  = vec_perm(temp.v,temp.v,tempc);
    temp.v = vec_splat(temp.v,0);
    return temp;
}

MYINLINE vector float vec_float_splat(float *thisptr)
{
    vector unsigned char tempc;
    vector float temp;
    temp = vec_lde(0,thisptr);
    tempc  = vec_lvsl(0,thisptr);
    temp  = vec_perm(temp,temp,tempc);
    temp = vec_splat(temp,0);
    return temp;
}


MYINLINE vector float vec_zero(void)
{
    return vec_ctf(vec_splat_u32(0),0);
}

MYINLINE vector float vec_float_one(void)
{
    return vec_ctf(vec_splat_u32(1),0);
}

MYINLINE void load_double_floatvec(FloatVec *out, double *in, long size)
{
    long i;
    float *tmp = (float *) mycalloc(size, sizeof(float));
    for(i=0; i < size; i++)
        tmp[i] = (float) in[i];
    memcpy(out,tmp,sizeof(float)*size);
    myfree(tmp);
}

MYINLINE void load_float_floatvec(FloatVec *out, float *in, long size)
{
    long i;
    float *tmp = (float *) mycalloc(size, sizeof(float));
    for(i=0; i < size; i++)
        tmp[i] = (float) in[i];
    memcpy(out,tmp,sizeof(float)*size);
    myfree(tmp);
}

///
/// Vector multiply operation, for vectors 4 floats long,
/// this is not very efficient because the pipelines are not full, try to use
/// the vdot_product() or then use standard FPU and for-loops.
MYINLINE MYREAL vector_multiply(vector float *v1, vector float *v2)
{
    float result;
    vector float temp;
    vector float zero = vec_zero();
    temp = vec_madd(*v1, *v2, zero);
    temp = vec_add(temp, vec_sld(temp,temp,4));
    temp = vec_add(temp, vec_sld(temp,temp,8));
    vec_ste(temp,0,&result);
    return (MYREAL) result;
}

///
/// Calculates the dot product of two vectors
/// of size of 4 floats each
MYREAL
vdot_product (FloatVec *v1, FloatVec *v2, int size)
{
    // size needs to be a multiple of 4 !!!!!!
    vector float temp1 = (vector float) vec_splat_s8(0);
    vector float temp2 = temp1;
    vector float temp3 = temp1;
    vector float temp4 = temp1;
    float result;
    int i;
    for (i = 0; i < size; i += 4)
    {
        temp1 = vec_madd(v1[i].v, v2[i].v, temp1);
        temp2 = vec_madd(v1[i+1].v, v2[i+1].v, temp2);
        temp3 = vec_madd(v1[i+2].v, v2[i+2].v, temp3);
        temp4 = vec_madd(v1[i+3].v, v2[i+3].v, temp4);
    }
    // sum the vectors
    temp1 = vec_add(temp1,temp2);
    temp3 = vec_add(temp3,temp4);
    temp1 = vec_add(temp1,temp3);
    
    // sum across vector
    temp1 = vec_add(temp1, vec_sld(temp1,temp1,4));
    temp1 = vec_add(temp1, vec_sld(temp1,temp1,8));
    
    // copy result to the stach to return to FPU
    vec_ste(temp1,0,&result);
    return (MYREAL) result;
}

//result = a + b
MYINLINE FloatVec vector_add(FloatVec a, FloatVec b )
{
    FloatVec temp;
    temp.v = vec_add( a.v, b.v);
    return temp;    
}

//result = a/b
MYINLINE vector float vector_divide( vector float a, vector float b )
{
    return vec_madd( a, vector_reciprocal( b ), (vector float)(0) );
}


//result = v^0.5
MYINLINE vector float vector_square_root( vector float v )
{
    
    return vec_madd( v, vector_reciprocal_square_root( v ), (vector float)(0) );
    
    
}

///
///result = b^-1
MYINLINE vector float vector_reciprocal( vector float v )
{
    
    //Get the reciprocal estimate
    vector float estimate = vec_re( v );
    
    //One round of Newton-Raphson refinement
    return vec_madd( vec_nmsub( estimate, v, (vector float) (1.0) ), estimate, estimate );
    
    
    
}

///
///result = v^-0.5
MYINLINE vector float vector_reciprocal_square_root( vector float v )
{
    
    //Get the square root reciprocal estimate
    vector float zero = (vector float)(0);
    vector float oneHalf = (vector float)(0.5);
    vector float one = (vector float)(1.0);
    vector float estimate = vec_rsqrte( v );
    
    //One round of Newton-Raphson refinement
    vector float estimateSquared = vec_madd( estimate, estimate, zero );
    vector float halfEstimate = vec_madd( estimate, oneHalf, zero );
    return vec_madd( vec_nmsub( v, estimateSquared, one ), halfEstimate, estimate );
}

#endif /*ALTIVEC*/
