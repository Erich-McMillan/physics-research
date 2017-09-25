N=316
D=3
Sq=0
SqTotal=1000
for(D in 67:79)
{
	finame=sprintf("triangle_compare_ChungLu_Vs_SNM_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	tSNM=read.table(finame,header=T)$tSNM
	tCL=read.table(finame,header=T)$tCL
	finame=sprintf("triangle_compare_Rewire_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	tRW=read.table(finame,header=T)$tRW
	finame=sprintf("triangle_compare_Rvl_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	tRvl=read.table(finame,header=T)$tRW
	finame=sprintf("triangle_compare_gnp_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	tNP=read.table(finame,header=T)$tNP
	finame=sprintf("triangle_compare_gnm_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	tNM=read.table(finame,header=T)$tNM
	finame=sprintf("triangle_compare_mygnm_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	tmyNM=read.table(finame,header=T)$tNM

	tRange=range(tSNM,tRW,tRvl,tCL,tNP,tNM,tmyNM)
	zstep=(tRange[2]-tRange[1])/50
	zbreaks=seq(tRange[1],tRange[2],zstep)
	hSNM=hist(tSNM,breaks=zbreaks,plot=F)
	hRW=hist(tRW,breaks=zbreaks,plot=F)
	hRvl=hist(tRvl,breaks=zbreaks,plot=F)
	hCL=hist(tCL,breaks=zbreaks,plot=F)
	hNP=hist(tNP,breaks=zbreaks,plot=F)
	hNM=hist(tNM,breaks=zbreaks,plot=F)
	hmyNM=hist(tmyNM,breaks=zbreaks,plot=F)
	
	foname=sprintf("triangle_compare_SNM_RW_Rvl_N%dD%dSqTotal%d.png",N,D,SqTotal)
	png(foname)
	plot(hSNM$mids,hSNM$counts,col="red",type="l",main=foname,xlab="number of triangles",ylab="frequency")
	lines(hRW$mids,hRW$counts,col="green")
	lines(hRvl$mids,hRvl$counts,col="blue")
	lines(hCL$mids,hCL$counts,col="black")
	lines(hNP$mids,hNP$counts,col="purple")
	lines(hNM$mids,hNM$counts,col="orange")
	lines(hmyNM$mids,hmyNM$counts,col="pink")
	#legend("topright",c("SNM","RW","Rvl","CL","NP","NM"),col=c("red","green","blue","black","purple","orange"),lty=1)
	legend("topright",c("Hard_SimpleNoMultiple","Hard_RewireHH","Hard_RewireVL","Soft_ChungLu","Soft_GNP","Soft_GNM","Soft_myNM"),col=c("red","green","blue","black","purple","orange","pink"),lty=1)
	dev.off()
}

