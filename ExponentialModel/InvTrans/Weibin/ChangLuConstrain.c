#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include"ran2dbl.h"

int main(int argc, const char *argv[])
{
	warmup_ran2dbl();

	int N=atoi(argv[1]);
	double gam=2;//atof(argv[2]);
	int SqMax=1e4;//atoi(argv[3]);
	int RANDOM=atoi(argv[4]);

	int i;
	for(i=0;i<RANDOM;i++)ran2dbl();//warmup using RANDOM


	////////// filename
	char fname[256];
	//sprintf(fname,"ChangLuConstrainTest_N%dr%.2fSqMax%d.txt",N,gam,SqMax);
	//sprintf(fname,"ChangLuConstrainTest_RG_N%dSqMax%d.txt",N,SqMax);
	//FILE *fp=fopen(fname,"w");

	double logsum, logsumtmp, logprob, logcumultmp, logcumul;
	int n,Sq,dn,dmax,good,goodsum;
	long dsum;
	double goodrate;
	for(gam=2.00;gam<4.00;gam+=0.01)
	{
		////////// logsum
		logsum = 0.0;
		for (n=2; n<N; n++)
		{
			logsumtmp = -gam*log(n);
			logsum = fmax(logsum,logsumtmp) + log1p(exp(-fabs(logsum-logsumtmp)));
		}

		////////// seqgen
		goodsum=0;
		for(Sq=0;Sq<SqMax;Sq++)
		{
			dmax=0;
			dsum=0;

			for (n=0; n<N; n++)
			{
				logprob = log1p(ran2dbl());
				logcumul = 0.0;
				i = 0;
				do
				{
					logcumultmp = -gam*log(++i)-logsum;
					logcumul = fmax(logcumul,logcumultmp) + log1p(exp(-fabs(logcumul-logcumultmp)));
				} while (logprob>logcumul);
				dn = i;
				printf("%d\n", i);
				dmax=(dmax>dn)?(dmax):(dn);
				dsum+=dn;
			}
			good=(dmax*dmax<dsum);
			goodsum+=good;
			//fprintf(fp,"%d\t%ld\t%d\n",dmax,dsum,good);
		}
		goodrate=(double)goodsum/(double)SqMax;
		//fprintf(fp,"%f\t%f\n",gam,goodrate);
		//printf("%f\t%f\n",gam,goodrate);
	}
	//fclose(fp);

	return 0;
}
