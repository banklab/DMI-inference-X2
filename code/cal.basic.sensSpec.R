library(DescTools)

args=commandArgs(trailingOnly = TRUE)

paras=function(fab, faB, fAb, fAB, nSample){
	gab=fab/nSample
	gaB=faB/nSample
	gAb=fAb/nSample
	gAB=fAB/nSample

	pA=gAb+gAB
	pB=gaB+gAB
	
	h1=1-pA^2-(1-pA)^2
	h2=1-pB^2-(1-pB)^2
    
	sigma2=h1+h2-h1^2-h2^2
	sk2=(2-h1-h2)*(h1+h2-1)+2*(gab^2+gaB^2+gAb^2+gAB^2)
	delta=sk2-sigma2
	X2=delta/sigma2
	D=gab*gAB-gAb*gaB


	#Keep the sign of the data
	if(D>0){
		Dmax=min((1-pA)*pB, pA*(1-pB))
	}
	else{
		Dmax=min(pA*pB, (1-pA)*(1-pB))
	}

	D_prime=D/Dmax
	cor=D/(pA*(1-pA)*pB*(1-pB))^0.5
	D2_ab=gab^2-(1-pA)^2*(1-pB)^2
	D2_aB=gaB^2-(1-pA)^2*pB^2
	D2_Ab=gAb^2-pA^2*(1-pB)^2
	D2_AB=gAB^2-pA^2*pB^2

	expProb=c(((1-pA)*(1-pB)), ((1-pA)*pB), (pA*(1-pB)), (pA*pB))
	gtest=GTest(x=c(fab, faB, fAb, fAB), p=expProb)$p.value

	return(c(gab, gAb, gaB, gAB, D, D_prime, cor, sigma2, sk2, delta, X2, D2_aB, D2_Ab, D2_AB, gtest))
}


#obs_ab obs_Ab obs_aB obs_AB
data=read.table(args[1], header=TRUE)
output=data.frame(matrix(, nrow=length(data[,1]), ncol=17))
names(output)=c("simulation", "Gen","pab", "pAb", "paB", "pAB", "D", "D_prime", "cor", "sigma2", "sk2", "delta", "X2", "D2_aB", "D2_Ab", "D2_AB", "GTest")
nSample=sum(data[1,7:10])
gabAll=data[,7]/nSample
gAbAll=data[,8]/nSample
gaBAll=data[,9]/nSample
gABAll=data[,10]/nSample

for(i in 1:length(data[,1])){
	tmp=paras(fab=data[i,7], fAb=data[i,8], faB=data[i,9], fAB=data[i,10], nSample=nSample)
	output[i,]=c(data[i,1], data[i,2], tmp)
}

### add simulation and generation

write.table(output, file=args[2], quote=FALSE, sep="\t", row.names = FALSE)

