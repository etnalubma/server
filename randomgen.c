#include <stdlib.h>
#include "randomgen.h"

#define UNI32 2.328306437e-10
#define UNI64 5.42101086243e-20

typedef unsigned long int unlong;

float ran2(long *idum);
unlong mzran13();
void ran13set(unlong xx, unlong yy, unlong zz, unlong nn);

struct srandgen{
    genmethod method;
    long seed;
    int bits;
};

randgen create_rg(genmethod method, long seed, int bits){
    randgen rg;
    rg = (randgen)malloc(1*sizeof(struct srandgen));
    
    rg->method = method;
    rg->seed = seed;
    rg->bits = bits;
    
    return rg;
}

randgen destroy_rg(randgen rg){
    free(rg);
    rg = NULL;
    return rg;
}

float random_gen(randgen rg){

    if(rg->method == RAN2){
        return ran2((long *)&(rg->seed));
    }
    else{ /*MZRAN13*/
        return rg->bits == 32 ? mzran13()*UNI32 : mzran13()*UNI64;
    }
          
}

/**
* RAN2
*
*/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


/**
* MZRAN13
*
*/
unlong x=521288629, y=362436069, z=16163801, c=1, n=1131199209;
unlong mzran13()
{ 
	long int s;
  	if (y>x+c) {
		s=y-(x+c); 
		c=0; 
	}
	else {
		s=y-(x+c)-18; 
		c=1; 
	}
	x=y; 
	y=z; 
	return (z=s) + (n=60069*n+1013904243);
}

void ran13set(unlong xx, unlong yy, unlong zz, unlong nn)
{ 
	x = xx; 
	y = yy; 
	z = zz; 
	c = y > z;
}

