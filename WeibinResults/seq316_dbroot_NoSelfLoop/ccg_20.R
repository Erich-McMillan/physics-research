library(igraph)
library(Rcpp)
cppFunction('
NumericMatrix pij(NumericVector b, NumericVector freq)
{
	int i,j,k;
	int N=0;
	for(i=0;i<freq.size();i++)N+=freq[i];
	NumericMatrix p(N,N);
	NumericVector blist(N);
	i=0;
	for(j=0;j<freq.size();j++)
	{
			for(k=0;k<freq[j];k++)
			{
				blist[i]=b[j];
				i++;
			}
	}
	
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			p(i,j)=1/(1+exp(-blist[i]-blist[j]));
		}
	}
	return p;
}
')

cppFunction('
List el_gen(NumericMatrix pij, int N)
{
	std::vector<int> startlist;
	std::vector<int> endlist;
	
	int i,j;
	
	for(i=0;i<N-1;i++)
	{
		for(j=i+1;j<N;j++)
		{
			if(R::runif(0,1)<pij(i,j))
			{
				startlist.push_back(i);
				endlist.push_back(j);
			}
		}
	}
	return List::create(startlist,endlist);
}
')



N=316
gamma=2.0
Sq=1
GperSq=1e4
#for(gamma in seq(2.0,3.5,0.5))
#{
for(Sq in 1:1e3)
{
	fin=sprintf("dbroot_N%dr%.1fSq%d.txt",N,gamma,Sq)
	df=read.table(fin)
	b=rev(df$V2)
	freq=rev(df$V3)
	P=pij(b,freq)
	
	edgeNumber=rep(-1,GperSq)
	ccgList=rep(-1,GperSq)
	triangleNumber=rep(-1,GperSq)
	for(Sd in 1:GperSq)
	{
		el_list=el_gen(P,N)
		el=cbind(el_list[[1]],el_list[[2]])+1
		g=graph_from_edgelist(el,directed=F)
		g=add_vertices(g,N-length(V(g)))
		#fout=paste("el_",gamma,"_",Sq,"_",Sd,".txt",sep="")
		#write.table(el,file=fout,row.names=F,col.names=F)
		#fout=paste("g_",gamma,"_",Sq,"_",Sd,".png",sep="")
		#png(fout)
		#plot(g,vertex.size=1,vertex.label=NA)
		#dev.off()
		edgeNumber[Sd]=length(E(g))
		ccgList[Sd]=transitivity(g)
		triangleNumber[Sd]=sum(count_triangles(g))/3
	}
	fout=paste("ccg_triangle_",gamma,"_",Sq,".txt",sep="")
	write.table(data.frame(edgeNumber,triangleNumber,ccgList),file=fout,row.names=F)
}

#}
