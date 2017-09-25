#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include"ran2dbl.h"
#include"GSamp.h"

double CharoNode(double doffset, double dj, double m)
{
	return 1;
}

double CharoStub(double doffset, double dj, double m)
{
	return dj;
}

double dd4m(double doffset, double dj, double m)
{
	double myweight=dj*pow(1-(dj-1)*(doffset-1)/(4*(m-1)*(m-1)+0.01),m-1);
	return myweight;
}

double dd3m(double doffset, double dj, double m)
{
	double myweight=dj*pow(1-(dj-1)*(doffset-1)/(3*(m-1)*(m-1)+0.01),m-1);
	return myweight;
}

double dd2m(double doffset, double dj, double m)
{
	double myweight=dj*pow(1-(dj-1)*(doffset-1)/((2*m-3)*(m-1)+0.01),m-1);
	return myweight;
}

int L3(int *zdegseq, int **zat, int **zam, int N, FILE *fp_l)
{
		// zdegseq: is the expected degree sequence? in my case might be the generated degree sequence
		// zat: adjtable
		// zam: adj matrix
		// N: size of network
		// fp_l: where the number of local triangles is written
    int i,j,k,a,b;
    int L3_sum=0;
		int L3local=0;
    for(k=0;k<N;k++)
    {
				L3local=0;
				for(i=0;i<zdegseq[k];i++)
        {
            a=zat[k][i];
            for(j=i+1;j<zdegseq[k];j++)
            {
                b=zat[k][j];
								L3_sum += zam[a][b];
								L3local += zam[a][b];
            }
        }
				fprintf(fp_l,"%d\n",L3local);
    }
    L3_sum=L3_sum/3;
    //printf("L3=%d\n",L3_sum);
    return L3_sum;
}

int L4(int *zdegseq, int **zat, int **zam, int **zel, int N, int edgenum)
{
		// zdegseq: the expected or generated degree sequence
		// zat: the adj table
		// zam: adj matrix
		// zel: edge list
		// N: size of network
		// edgenum: number of edges in graph
    int i,j,k,p,q,a,b;
    int L4_sum=0;
    for(k=0;k<edgenum;k++)
    {
        p=zel[k][0];
        q=zel[k][1];

        for(i=0;i<zdegseq[p];i++)
        {
            a=zat[p][i];
            if(a==q) continue;
            else
            {
                for(j=0;j<zdegseq[q];j++)
                {
                    b=zat[q][j];
                    if(b==p)continue;
                    else L4_sum += zam[a][b];
                }
            }
        }
    }
    L4_sum=L4_sum/4;
    //printf("L4=%d\n",L4_sum);
    return L4_sum;
}

int main(int argc, const char *argv[])
{
	warmup_ran2dbl();

	int N=atoi(argv[1]);
	double gamma=atof(argv[2]);
	int Sq=atoi(argv[3]);
	//int Sd=atoi(argv[4]);
	int RANDOM=atoi(argv[4]);
	int GperSq=1e1;

	int i,j,k,g;
	for(i=0;i<RANDOM;i++)ran2dbl();//warmup using RANDOM

	int *degseq=(int *)malloc(N*sizeof(int));

	char fds[256];
	char flgwt[256];
	char f3[256];
	char f4[256];
	char fg[256];
	char fl[256];

	int L3value;
	int L4value;
	double ccgvalue;

	////scan degseq
	sprintf(fds,"degseq_N%dr%.1fSq%d.txt",N,gamma,Sq);
	FILE *fp_ds=fopen(fds,"r");
	for(i=0;i<N;i++)
	{
		fscanf(fp_ds,"%d\n",&degseq[i]);
	}
	fclose(fp_ds);

	////calculate sumd1 and sumd2
	int sumd1=0;
	int sumd2=0;
	for(i=0;i<N;i++)
	{
			sumd1 += degseq[i];
			sumd2 += degseq[i]*degseq[i];
	}
	int edgenum = sumd1/2;

	////el
	int **zel=malloc(edgenum*sizeof(int*));
	zel[0]=malloc(edgenum*2*sizeof(int));
	for (k=1; k<edgenum; k++) zel[k] = zel[k-1]+ 2;
	////am
	int **zam=malloc(N*sizeof(int*));
	zam[0]=calloc(N*N,sizeof(int));
	for(i=1;i<N;i++) zam[i]=zam[i-1]+N;


	////generate and store data
	sprintf(flgwt,"lgwt_dd3m_N%dr%.1fSq%d.txt",N,gamma,Sq);
	FILE *fp_lgwt=fopen(flgwt,"w");
	sprintf(f3,"dd3m_L3_N%dr%.1fSq%d.txt",N,gamma,Sq);
	FILE *fp_3=fopen(f3,"w");
	sprintf(f4,"dd3m_L4_N%dr%.1fSq%d.txt",N,gamma,Sq);
	FILE *fp_4=fopen(f4,"w");
	sprintf(fg,"dd3m_ccg_N%dr%.1fSq%d.txt",N,gamma,Sq);
	FILE *fp_g=fopen(fg,"w");
	sprintf(fl,"dd3m_localTri_N%dr%.1fSq%d.txt",N,gamma,Sq);
	FILE *fp_l=fopen(fl,"w");

	graph G;
	gsaminit(degseq,N);
	for(g=0;g<GperSq;g++)
	{
		G=gsam(ran2dbl,dd3m);

		////el
		k=0;
		for(i=0;i<N;i++)
		{
			for(j=0;j<degseq[i];j++)
			{
				if(G.list[i][j]>i)//write each edge only once
				{
					zel[k][0]=i;
					zel[k][1]=G.list[i][j];
					k++;
				}
			}
			//fprintf(fp,"\n");
		}
		////am clear
		for(i=0;i<N;i++)
			for(j=0;j<N;j++)
				zam[i][j]=0;
		////am
		for(k=0;k<edgenum;k++)
		{
				i=zel[k][0];
				j=zel[k][1];
				zam[i][j]=1;
				zam[j][i]=1;
		}
		//G.list: adj table
		L3value=L3(degseq,G.list,zam,N,fp_l);
		L4value=L4(degseq,G.list,zam,zel,N,edgenum);
		ccgvalue=6.0*L3value/(sumd2-sumd1); // nice quick way to calculate gcc using degree sequence

		fprintf(fp_lgwt,"%f\n",G.weight);
		fprintf(fp_3,"%d\n",L3value);
		fprintf(fp_4,"%d\n",L4value);
		fprintf(fp_g,"%f\n",ccgvalue);
	}
	gsamclean();
	fclose(fp_lgwt);
	fclose(fp_3);
	fclose(fp_4);
	fclose(fp_g);
	fclose(fp_l);

	free(zam[0]);
	free(zam);
	free(zel[0]);
	free(zel);
	free(degseq);

	return 0;
}
