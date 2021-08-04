library(doParallel)
library(foreach)

warning()

### use multiple CPUs. This should be improved in the future.
cores=detectCores()
#cl=makeCluster(cores[1])
cl=makeCluster(4)
registerDoParallel(cl)

args=commandArgs(trailingOnly = TRUE)

#### test data. ready matrix for convience.

data <- read.table(args[1])

genome=as.matrix(data)
nPop=length(genome[,1])
nMarker=length(genome[1,])

het1=rep(0, nMarker)
het2=rep(0, nMarker)
het3=rep(0, nMarker)
het4=rep(0, nMarker)
p=rep(0, nMarker)

for(i in 1:nMarker){
	p[i]=sum(genome[,i])/nPop ### pop1 is 0, pop2 is 1.
	het1[i]=2*p[i]*(1-p[i])
	het2[i]=het1[i]^2
	het3[i]=het1[i]^3
	het4[i]=het1[i]^4
}

#output=c("#Pos1", "Pos2", "P1", "P2", "h_pos1", "h_pos2", "00", "11", "01", "10", "D", "D'","D11", "D00", "D10", "D01", "sk2", "L", "x2", "new_sk2", "new_x2", "change", "sampleSize","D2_11", "D2_00", "D2_10", "D2_01")
oNames=c("#Pos1", "Pos2", "P1", "P2", "h_pos1", "h_pos2", "00", "11", "01", "10", "D", "D'","D11", "D00", "D10", "D01", "sk2", "L", "x2", "new_sk2", "new_x2", "change", "sampleSize","D2_11", "D2_00", "D2_10", "D2_01")

haplotype=c("11", "10", "01", "00")

if(nMarker>1){
#for(i in 1:(nMarker-1)){
#	for(j in (i+1):nMarker){
output=foreach(i = 1:(nMarker-1), .combine='rbind') %:%
	foreach(j=(i+1):nMarker, .combine='rbind') %dopar% {
		g=as.data.frame(table(paste(genome[,i],genome[,j], sep="")))
# change the frequency in to proportion.
		g[,2]=g[,2]/nPop
		
#	print(g)

		h1=het1[i]+het1[j]
		h2=het2[i]+het2[j]
		h3=het3[i]+het3[j]
		h4=het4[i]+het4[j]

		sk2=(2-h1)*(h1-1)+2*sum(g[,2]^2)
		sigma=h1-h2
		x2=sk2/sigma-1
		var=(h1-h2/7+12*h3-6*h4-2*(h1-h2)^2)/2
		L=sigma+2*var^0.5

		missing_haplotype=setdiff(haplotype, g$Var1)

#		print(cbind(missing_haplotype), digits=3)
		
		if(! identical(missing_haplotype, character(0)) ){
			tmp=data.frame(Var1=missing_haplotype, Freq=rep(0, length(missing_haplotype)))
#			print(tmp)
			g=rbind(g, tmp)
#			print(g)
		}

		### 2, 10, 11, 00, 01, cal D', as the increase of the new haplotype.
		D=g[g$Var1=='01',2]*g[g$Var1=='10',2]-g[g$Var1=='11',2]*g[g$Var1=='00',2]
		

#Keep the sign of the data
		if(D<0){
			Dmax=min((1-p[i])*p[j], p[i]*(1-p[j]))
		} else{
			Dmax=min(p[i]*p[j], (1-p[i])*(1-p[j]))
		}

		## pop1 is 0, pop2 is 1. p is the frequency of allele from pop2
		D11=g[g$Var1=='11',2]-p[i]*p[j]
		D00=g[g$Var1=='00',2]-(1-p[i])*(1-p[j])
		D10=g[g$Var1=='10',2]-p[i]*(1-p[j])
		D01=g[g$Var1=='01',2]-(1-p[i])*p[j]

		D2_11=g[g$Var1=='11',2]^2-(p[i]*p[j])^2
		D2_00=g[g$Var1=='00',2]^2-((1-p[i])*(1-p[j]))^2
		D2_10=g[g$Var1=='10',2]^2-(p[i]*(1-p[j]))^2
		D2_01=g[g$Var1=='01',2]^2-((1-p[i])*p[j])^2
		
		ld=2*(D11^2+D00^2+D10^2+D01^2)+4*(D11*p[i]*p[j]+D00*(1-p[i])*(1-p[j])+D10*p[i]*(1-p[j])+D01*(1-p[i])*p[j])
		new_sk2=sigma+ld
		new_x2=new_sk2/sigma-1
#		print(D)
		Dprime=D/Dmax
		out=c(i, j, p[i], p[j], het1[i], het1[j], g[g$Var1=='00',2], g[g$Var1=='11',2], g[g$Var1=='01',2], g[g$Var1=='10',2], D, Dprime, D11, D00, D10, D01, sk2, L, x2, new_sk2, new_x2, ld, nPop, D2_11, D2_00, D2_10, D2_01)


#		output=rbind(output, out)
	}
#}

if(is.vector(output)){
	output=t(output)
}

} else {
	output=t(rep(NA, 27))
}

colnames(output)=oNames

write.table(output, file=args[2], quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

stopCluster(cl)

