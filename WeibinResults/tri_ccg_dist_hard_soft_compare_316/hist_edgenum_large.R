N=316
D=3
Sq=0
SqTotal=1000
for(D in 11:66)
{
	finame=sprintf("triangle_compare_ChungLu_Vs_SNM_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	m_CL=read.table(finame,header=T)$m_CL
	finame=sprintf("triangle_compare_gnp_N%dD%dSqTotal%d.txt",N,D,SqTotal)
	m_NP=read.table(finame,header=T)$m_NP

	tRange=range(m_CL,m_NP)
	zstep=(tRange[2]-tRange[1])/50
	zbreaks=seq(tRange[1],tRange[2],zstep)
	hCL=hist(m_CL,breaks=zbreaks,plot=F)
	hNP=hist(m_NP,breaks=zbreaks,plot=F)
	
	foname=sprintf("edgenum_compare_SNM_RW_Rvl_N%dD%dSqTotal%d.png",N,D,SqTotal)
	png(foname)
	plot(hCL$mids,hCL$counts,col="black",type="l",main=foname)
	lines(hNP$mids,hNP$counts,col="purple")
	abline(v=D*N/2,col="red")
	legend("topright",c("CL","NP"),col=c("black","purple"),lty=1)
	dev.off()
}

