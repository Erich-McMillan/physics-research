N=316
D=3
Sq=0
SqTotal=1000
for(D in 67:79)
{
	finame=sprintf("triangle_compare_ChungLu_Vs_SNM_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	ccg_SNM=read.table(finame,header=T)$ccg_SNM
	ccg_CL=read.table(finame,header=T)$ccg_CL
	finame=sprintf("triangle_compare_Rewire_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	ccg_RW=read.table(finame,header=T)$ccg_RW
	finame=sprintf("triangle_compare_Rvl_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	ccg_Rvl=read.table(finame,header=T)$ccg_RW
	finame=sprintf("triangle_compare_gnp_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	ccg_NP=read.table(finame,header=T)$ccg_NP
	finame=sprintf("triangle_compare_gnm_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	ccg_NM=read.table(finame,header=T)$ccg_NM
	finame=sprintf("triangle_compare_mygnm_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	ccg_myNM=read.table(finame,header=T)$ccg_NM

	ccg_Range=range(ccg_SNM,ccg_RW,ccg_Rvl,ccg_CL,ccg_NP,ccg_NM,ccg_myNM)
	zstep=(ccg_Range[2]-ccg_Range[1])/50
	zbreaks=seq(ccg_Range[1],ccg_Range[2],zstep)
	hSNM=hist(ccg_SNM,breaks=zbreaks,plot=F)
	hRW=hist(ccg_RW,breaks=zbreaks,plot=F)
	hRvl=hist(ccg_Rvl,breaks=zbreaks,plot=F)
	hCL=hist(ccg_CL,breaks=zbreaks,plot=F)
	hNP=hist(ccg_NP,breaks=zbreaks,plot=F)
	hNM=hist(ccg_NM,breaks=zbreaks,plot=F)
	hmyNM=hist(ccg_myNM,breaks=zbreaks,plot=F)
	
	foname=sprintf("ccg_compare_SNM_RW_Rvl_N%dD%dSqTotal%d.png",N,D,SqTotal)
	png(foname)
	plot(hSNM$mids,hSNM$counts,col="red",type="l",main=foname,xlab="GlobalClusteringCoefficient", ylab="frequency")
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

