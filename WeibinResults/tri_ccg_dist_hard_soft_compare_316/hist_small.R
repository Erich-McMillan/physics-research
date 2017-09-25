N=316
D=3
Sq=0
SqTotal=1000
for(D in 2:10)
{
	finame=sprintf("triangle_compare_ChungLu_Vs_SNM_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	tSNM=read.table(finame,header=T)$tSNM
	tCL=read.table(finame,header=T)$tCL
	finame=sprintf("triangle_compare_Rewire_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	tRW=read.table(finame,header=T)$tRW
	finame=sprintf("triangle_compare_Rvl_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	tRvl=read.table(finame,header=T)$tRW

	tRange=range(tSNM,tRW,tRvl,tCL)
	zstep=(tRange[2]-tRange[1])/100
	zbreaks=seq(tRange[1],tRange[2],1)
	hSNM=hist(tSNM,breaks=zbreaks,plot=F)
	hRW=hist(tRW,breaks=zbreaks,plot=F)
	hRvl=hist(tRvl,breaks=zbreaks,plot=F)
	hCL=hist(tCL,breaks=zbreaks,plot=F)
	foname=sprintf("triangle_compare_SNM_RW_Rvl_N%dD%dSqTotal%d.png",N,D,SqTotal)
	png(foname)
	plot(hSNM$mids,hSNM$counts,col="red",type="l",main=foname)
	lines(hRW$mids,hRW$counts,col="green")
	lines(hRvl$mids,hRvl$counts,col="blue")
	lines(hCL$mids,hCL$counts,col="black")
	dev.off()
}

