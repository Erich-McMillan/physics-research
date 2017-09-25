library(Rcpp)
library(igraph)

#Erdos Gallai test written by Charo
cppFunction('
int GTest(NumericVector DS)
{
	int n=DS.size();
	int *s=(int*)malloc(n*sizeof(int));
	int *xk=(int*)malloc(n*sizeof(int));
	int iz;
	for(iz=0;iz<n;iz++)s[iz]=DS[iz];
	
	int i, what, flag=0, k;
	long int degsum=0, minsum, c;

	for (i=0;i<n;i++) degsum += s[i];
	if (degsum%2) what = -1;
	else {
		c = -1;

		for (k=n-1; k>=(*s); k--) *(xk+k) = 0;									// Hehehe...
		for (k=1;k<n;k++) for (c=(*(s+k-1))-1; c>=(*(s+k)); c--) *(xk+c)=k;		// Dehihiho
		for (k=c;k>=0;k--) *(xk+k)=n;

		degsum = s[0];
		minsum = n-1;

		k = 1;
		while (k<n-1 && degsum<=minsum) {										// Test!
	        degsum += s[k];
	        if ((*(xk+k))<k+1) flag=1;
	        if (!flag) minsum += (*(xk+k))-1;
	        else minsum += 2*k-(*(s+k));
	        k++;
	    }

	    if (degsum>minsum) what = 0;
	    else what = 1;
	}
	
	free(s);
	free(xk);
    return what;
}
')

zERSG=function(N,p)#random degseq
{
	s=sort(rbinom(N-1,N,p),decreasing=T)
	while(GTest(s)!=1)s=sort(rbinom(N-1,N,p),decreasing=T)
	return(s)
}

cppFunction('
List el_ChungLu(NumericVector ds)
{
	int N=ds.size();
	int m=sum(ds)/2;
	
	std::vector<int> startlist;
	std::vector<int> endlist;
	
	int i,j;
	
	for(i=0;i<N-1;i++)
	{
		for(j=i+1;j<N;j++)
		{
			if(R::runif(0,1)<(ds[i]*ds[j]/2.0/m))
			{
				startlist.push_back(i);
				endlist.push_back(j);
			}
		}
	}
	return List::create(startlist,endlist);
}
')


zChungLu=function(ds)
{
	el_list=el_ChungLu(ds)
	el=cbind(el_list[[1]],el_list[[2]])+1
	g=graph_from_edgelist(el,directed=F)
	return(g)
}

N=316
D=3
Sq=0
SqTotal=10

tCL=rep(-1,SqTotal)
triple_CL=rep(-1,SqTotal)
tSNM=rep(-1,SqTotal)
triple_SNM=rep(-1,SqTotal)

ds=rep(D,N)
zdraw=T
for(Sq in seq(SqTotal))
{
	#ds=zERSG(N,p)
	g=zChungLu(ds)
	tCL[Sq]=sum(count_triangles(g))/3
	triple_CL[Sq]=sum(degree(g)*(degree(g)-1))/2
	if(zdraw==T)
	{
		fname=sprintf("gcl_N%dD%dSq%d.png",N,D,Sq)
		png(fname)
		plot(g,vertex.size=1,vertex.label=NA)
		dev.off()
	}
	
	g=sample_degseq(ds,method="simple.no.multiple")
	tSNM[Sq]=sum(count_triangles(g))/3
	triple_SNM[Sq]=sum(degree(g)*(degree(g)-1))/2	
	if(zdraw==T)
	{
		fname=sprintf("gsnm_N%dD%dSq%d.png",N,D,Sq)
		png(fname)
		plot(g,vertex.size=1,vertex.label=NA)
		dev.off()
	}
}
ccg_CL=tCL/triple_CL
ccg_SNM=tSNM/triple_SNM
d=data.frame(tCL,tSNM,triple_CL,triple_SNM,ccg_CL,ccg_SNM)
fname=sprintf("triangle_compare_ChungLu_Vs_SNM_N%dD%dSqTotal%d.txt",N,D,SqTotal)
write.table(d,file=fname,row.names=F)

