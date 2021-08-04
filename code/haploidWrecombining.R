
### A haploid population with recombination.

args=commandArgs(trailingOnly = TRUE)

### parameters
alpha=as.numeric(args[1]) ### direct selection on allele A
beta=as.numeric(args[2]) ### direct selection on allele B
gamma=as.numeric(args[3]) ### epistatic coefficient on AB
r=as.numeric(args[4]) ### recombination rate
generation=1000 ### total generation for simulation.
m=0 ### migration rate can be changed here.

### define G0 and fitness. gab, gAb, gaB, gAB are haplotype frequencies. w is fitness foreach haplotype.
gab=0 # x1 0
w0=1
gaB=as.numeric(args[5]) # genotype aB,  beta
P0_B=gaB
wb=1+beta
gAb=1-gaB
#gAb=1-gaB # x3 alpha
P0_A=gAb
wa=1+alpha
gAB=0 # x4 alpha+beta-gamma
wab=(1+alpha)*(1+beta)*(1-gamma)

### Define the variables.for record results.
sigma2=0
sk2=0
delta=0
title=c("#Gen", "gab", "gaB", "gAb", "gAB", "D", "D_prime", "r", "sigma2", "sk2", "delta", "X2", "OR", "logOR", "dev_ab" ,"dev_aB", "dev_Ab", "dev_AB", "log_dev_ab", "log_dev_aB", "log_dev_Ab","log_dev_AB", "recomb_diff", "anst_diff", "recomb_diff_pro", "anst_diff_pro","pA", "pB", "ds_recomb", "ds_anst", "recomb")
before=matrix(NA, ncol=31, nrow=generation)
after_s=matrix(NA, ncol=31, nrow=generation)
after_r=matrix(NA, ncol=31, nrow=generation)
i=0
dsk2=0

### function calculate parameters
paras=function(gab, gaB, gAb, gAB){
	
	pA=gAb+gAB
	pB=gaB+gAB
	
	h1=1-pA^2-(1-pA)^2
	h2=1-pB^2-(1-pB)^2
	sigma2=h1+h2-h1^2-h2^2
	sk2=(2-h1-h2)*(h1+h2-1)+2*(gab^2+gaB^2+gAb^2+gAB^2)
	delta=sk2-sigma2
	X2=delta/sigma2
	D=gab*gAB-gAb*gaB
	delta1=delta-8*D^2
	D_prime=D/min(pA*pB, (1-pA)*(1-pB))
	cor=D/(pA*(1-pA)*pB*(1-pB))^0.5
	if(gaB*gAb!=0){
	OR=(gab*gAB)/(gaB*gAb)
	LOR=log(OR)
	}
	else{
		OR=99999
		LOR=99999
	}

	dev_ab=gab^2-(1-pA)^2*(1-pB)^2
	dev_aB=gaB^2-(1-pA)^2*pB^2
	dev_Ab=gAb^2-pA^2*(1-pB)^2
	dev_AB=gAB^2-pA^2*pB^2

	l_dev_ab=10^dev_ab
	l_dev_aB=10^dev_aB
	l_dev_Ab=10^dev_Ab
	l_dev_AB=10^dev_AB

	diffr=(dev_AB-dev_ab)/(abs(dev_aB+dev_Ab)+abs(dev_AB+dev_ab))
	diffa=(dev_Ab-dev_aB)/(abs(dev_aB+dev_Ab)+abs(dev_AB+dev_ab))

	diffa_pro=abs(dev_Ab-dev_aB)/(abs(dev_aB+dev_Ab))
	diffr_pro=abs(dev_AB-dev_ab)/(abs(dev_AB+dev_ab))

	return(c(D, D_prime, cor, sigma2, sk2, delta, X2, OR, LOR, dev_ab, dev_aB, dev_Ab, dev_AB, l_dev_ab, l_dev_aB, l_dev_Ab, l_dev_AB, diffr, diffa, diffr_pro, diffa_pro, pA, pB))
	
}

#### functioni to calculate expected trajectory under direct selection. Only valid for no-DMI sceinarios.
ds=function(s, t, p0){
        pt=p0/(p0+(1-p0)*exp(-s*t))
        return(pt)
}

#### record the negative X2 intervel (negX2start, negX2end), and the minimal X2 value (minX2) and time (minX2Gen)
negX2start=NA
negX2end=NA
minX2Gen=NA
minX2=0

#### Simulation
for(i in 1:generation){
	pA=gAb+gAB
	pB=gaB+gAB

### Stop the simulation when any alleles are almost lost/fixed in the population.
	if( pA<0.005 || pB<0.005 || pA>0.995 || pB>0.995 ){break}

#### The expected allele frequecies under direct selection model.
	exp_PA=ds(s=alpha,t=i, p0=P0_A)
	exp_PB=ds(s=beta,t=i, p0=P0_B)
	exp_diff_a=exp_PA-exp_PB ## col 30
	exp_diff_r=1-exp_PA-exp_PB ## col 29

## before selection and recombination, record the values for each paramter. Or record the parents characters.
	before[i,]=c(i, gab, gaB, gAb, gAB, paras(gab=gab, gaB=gaB, gAb=gAb, gAB=gAB), exp_diff_r, exp_diff_a, r)
	###define X2 value
	if(before[i,12]<0){
		if(is.na(negX2start)){
			negX2start=i
		}
		if(negX2start>0){
			negX2end=i
		}
		if(minX2>before[i,12]){
			minX2=before[i,12]
			minX2Gen=i
		}
	}

############ Selection section
### mean population fitness
	w_bar=w0*gab+wb*gaB+wa*gAb+wab*gAB
### new haplotype frequencies
	Gab=gab*(w0/w_bar)
	GaB=gaB*(wb/w_bar)
	GAb=gAb*(wa/w_bar)
	GAB=gAB*(wab/w_bar)

### Linkage disequilibrium
	D=Gab*GAB-GaB*GAb

#### After recombination + migration
	gab=Gab-r*D-m*Gab
	gaB=GaB+r*D+m*(1-GaB)
	gAb=GAb+r*D-m*GAb
	gAB=GAB-r*D-m*GAB
}

#### output file nanme
OF_name=paste("k1965", alpha, beta, gamma, r, args[5], sep="_")
#### output the trajectories
write.table(before, file=paste(args[6], "/",  OF_name, ".data", sep=""), quote=FALSE, col.names=title, row.names=FALSE)

#### output negative X2 intervals.
stats=cbind(r, args[5], negX2start, negX2end, minX2, minX2Gen)
colnames(stats)=c("#recombinationRate", "admixPro.", "negX2startGen", "negX2endGen", "minX2", "minX2Gen")
write.table(stats, file=paste(args[6], "/", OF_name, ".stat", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE)

