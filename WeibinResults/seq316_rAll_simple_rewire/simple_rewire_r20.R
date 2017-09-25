library(igraph)

#### Havel Hakimi

hh_addnode=function(degseq)
{
	nodenum=length(degseq)
	g=make_empty_graph(n=0,directed=F)
	for(i in seq(nodenum))
	{
		g=add_vertices(g,1,resdeg=degseq[i])
	}
	
	return(g)
}

hh_addedge=function(g)
{
	vlist=as_ids(V(g))
	rdlist=V(g)$resdeg
	d=data.frame(vlist,rdlist)
	dsorted=d[order(d$rdlist,decreasing=T),]
	
	i=dsorted$vlist[1]
	ri=dsorted$rdlist[1]
	for(k in seq(ri))
	{
		j=dsorted$vlist[k+1]
		g=add_edges(g,c(i,j))
		V(g)[i]$resdeg=V(g)[i]$resdeg-1
		V(g)[j]$resdeg=V(g)[j]$resdeg-1
	}
	
	if(max(V(g)$resdeg)>0)
	{
		hh_addedge(g)
	}
	else
	{
		return(g)
	}
	
}

hh=function(degseq)
{
	g=hh_addnode(degseq)
	g=hh_addedge(g)
	return(g)
}

N=316
gamma=2.0
Sq=1
GperSq=1e4

for(Sq in 1:1000)
{	
	print(Sq)
	fi=sprintf("degseq_N%dr%.1fSq%d.txt",N,gamma,Sq)
	a=read.table(fi)
	zdegseq=a$V1
	sweeplen=sum(zdegseq)

	ccgList=rep(-1,GperSq)
	triangleNumber=rep(-1,GperSq)

	g=hh(zdegseq)
	#print("hh done")
	for(Sd in 1:GperSq)
	{
		#print(Sd)
		g=rewire(g,with=keeping_degseq(loops=F,niter=sweeplen))
		ccgList[Sd]=transitivity(g)
		triangleNumber[Sd]=sum(count_triangles(g))/3
	}
	
	fout=sprintf("SimpleRewire_ccg_triangle_N%dr%.1fSq%d.txt",N,gamma,Sq)
	write.table(data.frame(triangleNumber,ccgList),file=fout,row.names=F)
}
