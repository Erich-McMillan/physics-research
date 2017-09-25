library(Rcpp)
cppFunction('
NumericVector diffc(NumericVector b, NumericVector deg, NumericVector freq)
{
	int n=b.size();
	NumericVector diff(n);
	
	int i,j;
	for(i=0;i<n;i++)
	{
		diff[i]=deg[i];
		for(j=0;j<n;j++)
		{
			diff[i]-=freq[j]/(1+exp(-b[i]-b[j]));
		}
	}
	return diff;
}
')

library(rootSolve)
solvebeta=function(deg,freq)
{
	bg=log(deg/sqrt(sum(deg*freq)))
	zo=multiroot(f=diffc,deg=deg,freq=freq,start=bg)
	while(zo$estim.precis>1e-3)zo=multiroot(f=diffc,deg=deg,freq=freq,start=runif(length(deg),-0.1,0.1))
	return(zo)
}

N=316
gamma=2.0
Sq=0
for(gamma in seq(2.0,3.5,0.5))
{
for(Sq in 1:1000)
{
	fin=sprintf("degseq_N%dr%.1fSq%d.txt",N,gamma,Sq)
	d=read.table(fin)$V1
	t=as.data.frame(table(d))
	deg=as.numeric(as.vector(t$d))
	freq=t$Freq
	
	s=solvebeta(deg,freq)
	b=s[[1]]
	dout=data.frame(deg,b,freq,s[[2]])
	dout=dout[order(deg,decreasing=T),]
	fout=sprintf("dbroot_N%dr%.1fSq%d.txt",N,gamma,Sq)
	write.table(dout,file=fout,row.names=F,col.names=F)
}
}
