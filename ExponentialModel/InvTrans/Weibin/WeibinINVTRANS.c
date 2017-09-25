#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include"ran2dbl.h"

int main(int argc, const char *argv[])
{
	warmup_ran2dbl();

	int N=atoi(argv[1]);
	double gam= atof(argv[2]);//atof(argv[2]);
	int SqMax= atoi(argv[3]);//atoi(argv[3]);
	int RANDOM=atoi(argv[4]);

	int i;
	for(i=0;i<RANDOM;i++)ran2dbl();//warmup using RANDOM

	////////// filename
	char fname[256];
	//sprintf(fname,"ChangLuConstrainTest_N%dr%.2fSqMax%d.txt",N,gam,SqMax);
	sprintf(fname,"ChangLuConstrainTest_RG_N%dSqMax%d.txt",N,SqMax);
	FILE *fp=fopen(fname,"w");

	double logsum, logsumtmp, logprob, logcumultmp, logcumul;
	int n,Sq,dmax,good,goodsum, *dn;

  dn = (int*) calloc(N, sizeof(int));

		////////// logsum
	logsum = 0.0;
	for (n=2; n<N; n++)
	{
		logsumtmp = gam*log(n);
		logsum = fmax(logsum,logsumtmp) + log1p(exp(-fabs(logsum-logsumtmp)));
	}

  int j;
  for(j = 0; j < SqMax; j++) {
  	////////// seqgen
  	for (n=0; n<N; n++)
  	{
  		logprob = log1p(ran2dbl());
  		logcumul = 0.0;
  		i = 0;
  		do
  		{
  			logcumultmp = gam*log(++i)-logsum;
  			logcumul = fmax(logcumul,logcumultmp) + log1p(exp(-fabs(logcumul-logcumultmp)));
  		} while (logprob>logcumul);

      //printf("%d\n", i);
      dn[i]+=1.0;

  	}

  }

  for(n = 1; n < N; n++) {
    fprintf(fp, "%d\t%lf\n", n, dn[n]/(double)(N*SqMax));
    //fprintf(fp, "%d\t%d\n", n, dn[n]);

  }


	fclose(fp);

	return 0;
}
