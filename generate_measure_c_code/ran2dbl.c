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

#define NULL 0//see if necessary for time.
#include<time.h>

double ran2dbl(void)
{
	//zpart start
	static long z_spaceofalongint=0;//const seed
	static long *idum=&z_spaceofalongint;
	static long ran2count=0;//count how many times this function is used.
	if(ran2count==0)z_spaceofalongint=(long)time(NULL);//if first time, seed it with time//this sentence can turn on time seed.
	ran2count++;
	//zpart end
	
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;
	
	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum <0) *idum += IM1;
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

int warmup_ran2dbl()
{
	//warm up ran2
	int i_warmup,n_warmup=1e3;
	for(i_warmup=0;i_warmup<n_warmup;i_warmup++)ran2dbl();
	return(n_warmup);
}