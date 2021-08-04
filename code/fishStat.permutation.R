library(doParallel)
library(foreach)

warning()

cores=detectCores()
cl=makeCluster(cores[1]-1)
registerDoParallel(cl)

args=commandArgs(trailingOnly = TRUE)
map=read.table(args[1])
ped=read.table(args[2])

sampleName=ped[,1]
geno=as.matrix(ped[,seq(7,length(ped[1,]),2)])
geno[geno==0] <- NA
geno=geno-1
marker=as.vector(map[,2])
nMarker=length(map[,2])

chrAall <- read.table(args[1])
chrBall <- read.table(args[2])

chrA_nMarker=length(chrAall[,1])
numCol=length(chrAall[1,])
tmp=t(chrAall[,3:(numCol-3)])
chrA=as.matrix(tmp)

chrB_nMarker=length(chrBall[,1])
numCol=length(chrBall[1,])
tmp=t(chrBall[,3:(numCol-3)])
chrB=as.matrix(tmp)

haplotype=c("11", "10", "01", "00")

if(nMarker>1){
#for(i in 1:chrA_nMarker){
#	for(j in 1:chrB_nMarker){
output=foreach(i = 1:(nMarker-1), .combine='rbind') %:%
	foreach(j=(i+1):nMarker, .combine='rbind') %dopar% {
		haplo2=cbind(geno[,i], geno[,j])
		haplo2=na.omit(haplo2)
		nPop=length(haplo2[,1])
		pA=sum(haplo2[,1])/nPop
		pB=sum(haplo2[,2])/nPop

		if(pA==1 || pA==0 || pB==1 || pB==0 || nPop<50){
			return(NULL)
		}

		g=as.data.frame(table(factor(paste(haplo2[,1], haplo2[,2], sep=""), levels=c("00", "01", "10", "11"))))
		gfreq=g
		gRaw=g
		gfreq[,2]=gfreq[,2]^0.5
		nGeno=sum(gfreq[,2])
		g[,2]=gfreq[,2]/nGeno
		pA=g[3,2]+g[4,2]
		pB=g[2,2]+g[4,2]
		expProb=c(((1-pA)*(1-pB))^2, ((1-pA)*pB)^2, (pA*(1-pB))^2, (pA*pB)^2)
		expProb=expProb/sum(expProb)

		# change the frequency in to proportion.
		h1A=2*pA*(1-pA)
		h1B=2*pB*(1-pB)
		h1=h1A+h1B
		h2=h1A^2+h1B^2

		sk2=(2-h1)*(h1-1)+2*sum(g[,2]^2)
		sigma=h1-h2
		x2=sk2/sigma-1

		### 2, 10, 11, 00, 01, cal D', as the increase of the new haplotype.
		D=g[g$Var1=='01',2]*g[g$Var1=='10',2]-g[g$Var1=='11',2]*g[g$Var1=='00',2]
		

#Keep the sign of the data
		if(D<0){
			Dmax=min((1-pA)*pB, pA*(1-pB))
		} else{
			Dmax=min(pA*pB, (1-pA)*(1-pB))
		}

		## The gametes contribute to the deviation of variance.
		D2_11=g[g$Var1=='11',2]^2-(pA*pB)^2
		D2_00=g[g$Var1=='00',2]^2-((1-pA)*(1-pB))^2
		D2_10=g[g$Var1=='10',2]^2-(pA*(1-pB))^2
		D2_01=g[g$Var1=='01',2]^2-((1-pA)*pB)^2

		Dprime=D/Dmax
		out=c(marker[i], marker[j], pA, pB, g[g$Var1=='00',2], g[g$Var1=='11',2], g[g$Var1=='01',2], g[g$Var1=='10',2], D, Dprime, sk2, sigma, x2, nPop, D2_11, D2_00, D2_10, D2_01)
		print(out)
	}
}
#}
colnames(output)=c("#Pos1", "Pos2", "P1", "P2", "00", "11", "01", "10", "D", "Dprime", "sk2", "sigma2", "x2", "sampleSize", "D2_11", "D2_00", "D2_10", "D2_01")

write.table(output, file=args[3], quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

stopCluster(cl)
