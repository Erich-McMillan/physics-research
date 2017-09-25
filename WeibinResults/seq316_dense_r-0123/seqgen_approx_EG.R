library(Rcpp)

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

zPLSG=function(N,gamma)
{
	kmax=N-0.5
	kmin=0.5
	kmax1r=kmax^(1-gamma)
	kmin1r=kmin^(1-gamma)
	k1rdiff=kmax1r-kmin1r
	
	s=sort(round((kmin1r+k1rdiff*runif(N))^(1/(1-gamma))),decreasing=T)
	while(GTest(s)!=1)
	{
		s=sort(round((kmin1r+k1rdiff*runif(N))^(1/(1-gamma))),decreasing=T)
	}
	return(s)
}

N=316
gamma=1.0
Sq=0

#for(gamma in seq(0.5,1.5,0.5))
#for(gamma in seq(2.0,3.5,0.1))
for(gamma in seq(-3,0,1))
{
	#dirname=sprintf("/Users/wbzhang/Desktop/Erich_solve_beta/degseq_N%dr%.1f/",N,gamma)
	#dir.create(dirname)
	#setwd(dirname)
for(Sq in seq(1000))
{
	s=zPLSG(N,gamma)
	fname=sprintf("degseq_N%dr%.1fSq%03d.txt",N,gamma,Sq)
	write.table(s,file=fname,row.names=F,col.names=F)
}
}

