#include <math.h>
#define NRANSI
#include "nrutil.h"
#define MAXITS 100
#define TOLF 1.0e-4
#define TOLMIN 1.0e-6
#define TOLX 1.0e-7
#define STPMX 100.0

//int nn;
//float *fvec;
//void (*nrfuncv)(int n, float v[], float f[]);
int nrITSN;

#define FREERETURN {free_vector(fvec,1,n);free_vector(xold,1,n);\
	free_vector(p,1,n);free_vector(g,1,n);free_matrix(fjac,1,n,1,n);\
	free_ivector(indx,1,n);return;}

void newt(float x[], int n, int *check,
	void (*vecfunc)(int, float [], float []))
{
	void fdjac(int n, float x[], float fvec[], float **df,
		void (*vecfunc)(int, float [], float []));
	float f_nr_min(float x[]);
	void lnsrch(int n, float xold[], float fold, float g[], float p[], float x[],
		 float *f, float stpmax, int *check, float (*func)(float []));
	void lubksb(float **a, int n, int *indx, float b[]);
	void ludcmp(float **a, int n, int *indx, float *d);
	//void qrupdt(float **r, float **qt, int n, float u[], float v[]);
	int i,its,j,*indx;
	float d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;

	indx=ivector(1,n);
	fjac=matrix(1,n,1,n);
	g=vector(1,n);
	p=vector(1,n);
	xold=vector(1,n);
	fvec=vector(1,n);
	nn=n;
	nrfuncv=vecfunc;
	f=f_nr_min(x);
	test=0.0;
	for (i=1;i<=n;i++)
		if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
	if (test<0.01*TOLF) FREERETURN
	for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
	stpmax=STPMX*FMAX(sqrt(sum),(float)n);
	for (its=1;its<=MAXITS;its++) {
		nrITSN = its;
		fdjac(n,x,fvec,fjac,vecfunc);
		for (i=1;i<=n;i++) {
			for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec[j];
			g[i]=sum;
		}
		for (i=1;i<=n;i++) xold[i]=x[i];
		fold=f;
		for (i=1;i<=n;i++) p[i] = -fvec[i];
		ludcmp(fjac,n,indx,&d);
		//qrupdt(fjac,);

		if(errLUDCMP) {
			*check = 5;
			return;
		}
		lubksb(fjac,n,indx,p);
		lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,f_nr_min);
		test=0.0;
		for (i=1;i<=n;i++)
			if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
		if (test < TOLF) {
			*check=0;
			FREERETURN
		}
		if (*check) {
			test=0.0;
			den=FMAX(f,0.5*n);
			for (i=1;i<=n;i++) {
				temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
				if (temp > test) test=temp;
			}
			*check=(test < TOLMIN ? 1 : 0);
			FREERETURN
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) FREERETURN
	}
	//nrerror("MAXITS exceeded in newt");
    *check = 2;
    //FREERETURN
}
#undef MAXITS
#undef TOLF
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
#undef NRANSI
