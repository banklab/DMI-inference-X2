#!/usr/local/bin/python3.7

import numpy as np
import pandas as pd
import argparse
from collections import Counter

class birthDeathProcess():
	"""
	define a populatiion in birth-death process, which to simulate the fish hybrid population.
	diploid
	Record generations
	DMI is the post-zygotic isolation
	additive effect on each locus
	the dmi is recessive, or as the h1 type, at least one locus is homozygous, the negative epistasis can be seen.
	"""
	def __init__(self, nPop, frqMPP, c, Gtime, alpha, beta, gamma, m1, m2):
		self.nPop=nPop
		self.frqMPP=frqMPP
		self.c=c
		self.Gtime=Gtime
		self.alpha=alpha
		self.beta=beta
		self.gamma=gamma
		self.gametes=["ab", "Ab", "aB", "AB"]
		###set all the gametes and genotypes.
		self.genoAll=[]
		for gamete1 in self.gametes:
			for gamete2 in self.gametes:
				geno=gamete1+gamete2
				self.genoAll.append(geno)
		#### migration
		self.m1=float(m1)
		self.m2=float(m2)
		### define the fitness matrix
		a=self.alpha
		b=self.beta
		e=self.gamma
		self.fitness={"abab":1, "abaB": (1+b), "abAb":(1+a), "abAB":(1+a)*(1+b),
			"aBab":(1+b), "aBaB":(1+b)**2, "aBAb":(1+a)*(1+b), "aBAB":(1+a)*(1+b)**2*(1+e)**2,
			"Abab":(1+a), "AbaB":(1+a)*(1+b), "AbAb":(1+a)**2, "AbAB":(1+a)**2*(1+b)*(1+e)**2,
			"ABab":(1+a)*(1+b), "ABaB":(1+a)*(1+b)**2*(1+e)**2, "ABAb":(1+a)**2*(1+b)*(1+e)**2,"ABAB":(1+a)**2*(1+b)**2*(1+e)**4}

###### basic function: calculate the genotype frequency from a population of genotypes.
	def calFrqGenotypes(self, genotypes):
		frqGenotypes={}
		for geno in self.genoAll:
			frqGenotypes[geno]=0
		tmpFrqGeno=np.unique(genotypes, return_counts=True)
		nInd=sum(tmpFrqGeno[1])
		for i in range(0,len(tmpFrqGeno[0])):
			frqGenotypes[tmpFrqGeno[0][i]]=tmpFrqGeno[1][i]/nInd
		
		return(frqGenotypes)

####the frequencies of each genotype after selection
	def afterSelFfrq(self, frqGenotypes):
		w_bar=0
		tmpFrq={}
		w={}
		for geno in self.genoAll:
			w_bar=w_bar+self.fitness[geno]*frqGenotypes[geno]
			w[geno]=0

		for geno in self.genoAll:
			tmpFrq[geno]=frqGenotypes[geno]*self.fitness[geno]/w_bar

		return tmpFrq

###### the frequency of each genotype after recombination
	def afterRecombFrq(self, frqGenotypes):
		frqGenoAR={} ## genotype frequency after recombination
		for loc11 in ("a", "A"):
			for loc12 in ("a", "A"):
				for loc21 in ("b", "B"):
					for loc22 in ("b", "B"):
						#target genotype
						gamete1=loc11+loc21
						gamete2=loc12+loc22
						geno=gamete1+gamete2
						#recombination from
						gamete1=loc11+loc22
						gamete2=loc12+loc21
						recomb_g1=gamete1+gamete2
						frqGenoAR[geno]=frqGenotypes[geno]+frqGenotypes[recomb_g1]*self.c-frqGenotypes[geno]*self.c
		return frqGenoAR

##### reproduction process: 1. recombination; 2. selection;
	def reproduction(self, frqGenotypes, nOffspring):
		### gametes frequencies after recombination
		afterR_frqGenotypes=self.afterRecombFrq(frqGenotypes=frqGenotypes)
		frqGametes=self.calGametesFrq(frqGenotypes=afterR_frqGenotypes)
		g1=list(np.random.choice(["ab", "Ab", "aB", "AB"], nOffspring, p=[frqGametes["ab"], frqGametes["Ab"], frqGametes["aB"], frqGametes["AB"]]))
		g2=list(np.random.choice(["ab", "Ab", "aB", "AB"], nOffspring, p=[frqGametes["ab"], frqGametes["Ab"], frqGametes["aB"], frqGametes["AB"]]))
		genomes=pd.DataFrame({"g1":g1, "g2":g2})
		genotypes=list(genomes["g1"]+genomes["g2"])
		tmpFrqGenotypes=self.calFrqGenotypes(genotypes=genotypes)
		### gamtes frequencies after selection
		afterS_frqGenotypes=self.afterSelFfrq(frqGenotypes=tmpFrqGenotypes)
		
		return(afterS_frqGenotypes)

##### after reproduction, change the ages of the population
	def viability(self, genotypes, ages, ageStructure):
		#the viability/survavil rate is controled by both the fitness of the individuals or the 
		#Here, the fitness can be decided by two components, 1. the age structure, 2. the population size.
		nInd=len(genotypes)
		survival=1-ageStructure
		totalFitness=0
		for age in ages:
			totalFitness += survival[age]
		fitness=[survival[age]*self.nPop/totalFitness for age in ages]
		random=np.random.uniform(0, 1, nInd)
		newGenotypes=[]
		newAges=[]
		for i in range(nInd):
			if(fitness[i]>random[i]):
				newGenotypes.append(genotypes[i])
				newAges.append(ages[i]+1)

		return(newGenotypes, newAges)

###basic function: returned in the population size, not the haplotype size.
	def calGametesFrq(self, frqGenotypes):
		frqGametes={"ab":0, "aB":0, "Ab":0, "AB":0}
		for gamete1 in self.gametes:
			for gamete2 in self.gametes:
				geno=gamete1+gamete2
				frqGametes[gamete1]=frqGametes[gamete1]+frqGenotypes[geno]
				frqGametes[gamete2]=frqGametes[gamete2]+frqGenotypes[geno]
	
		for gamete in self.gametes:
			frqGametes[gamete]=frqGametes[gamete]/2
		
		return frqGametes

#### sampling a group of individuals
	def sampling4genotyping(self, m, nPop, frqGenotypes):
		sampling={}
		nSampling=m
		nTotal=nPop ## not the defined population size, but the actual population size
		nRemaining=m
		tmpn=0
		tmpt=0
		for geno in self.genoAll:
			nGeno=int(frqGenotypes[geno]*nTotal+0.5)
			others=nPop-nGeno
			if(nRemaining<1 or others<0):
				sampling[geno]=0
				continue
			sampling[geno]=np.random.hypergeometric(nGeno, others, nRemaining)
			tmpn+=sampling[geno]
			tmpt+=nGeno
			nRemaining=nRemaining-sampling[geno]
			nPop=others
		return sampling
	
	def migration(self, frqGenotypes):
		newfrqGenotypes={}
		for geno in self.genoAll:
			if geno=="AbAb":
				newfrqGenotypes[geno]=frqGenotypes[geno]+(1-frqGenotypes[geno])*self.m1-frqGenotypes[geno]*self.m2
			elif geno=="aBaB":
				newfrqGenotypes[geno]=frqGenotypes[geno]+(1-frqGenotypes[geno])*self.m2-frqGenotypes[geno]*self.m1
			else:
				newfrqGenotypes[geno]=frqGenotypes[geno]*(1-self.m1-self.m2)

		return newfrqGenotypes



		return newfrqGenotypes
	

	def nonWF(self):
			#nPop, c, frqMPP, nGen, alpha, beta, gamma, h1, h2):
		## ten generations, low viability for the new offspring; The age structure can be changed to overlapping generation model.To maintain an equilbrium, the number of newIndividuals should be calculated by Moran.equilibrium.py. See Crow and Kimura 1973 Chapter 1.3, Model 3: overlapping generations, discrete time intervals. 
		ageStructure=np.array([0, 1, 1, 1, 1, 1, 1, 1, 1, 1.0])
		frqGametes={'ab':0, 'Ab':self.frqMPP,'aB':1-self.frqMPP, 'AB':0}
		
		#The output in the future; The frequency in population, as g_xx; the sampling frequency, as obs_xx, in the WF model, 200 individuals are sampled each generation. 
		# record from F1. I do not record G0.
		generation=[]
		g_ab=[]
		g_Ab=[]
		g_aB=[]
		g_AB=[]
		obs_ab=[]
		obs_Ab=[]
		obs_aB=[]
		obs_AB=[]
		samplingInd=[]
		frqGenotypes={}

	## set G0
		genotypes=list(np.random.choice(["AbAb", "aBaB"], self.nPop, p=[self.frqMPP, 1-self.frqMPP]))
		ages=list(np.random.choice(range(10), self.nPop))
	
#### from generation 2
		for t in range(1, self.Gtime+1):
			### calculate the frequencies of genotypes for generating the offspring
			frqGenotypes=self.calFrqGenotypes(genotypes=genotypes)
			#### maybe add migration here. if need, write a method here, change the frqGenotypes
			if self.m1+self.m2 > 0 :
				frqGenotypes=self.migration(frqGenotypes=frqGenotypes)

			### get the genotype frequencies in offspring, and sample the offspring randomly (alpha, beta, gamma, frqGenotypes, c, h1, h2)
			nOffspring=int(1*len(genotypes))
			frqGenotypesOffspring=self.reproduction(frqGenotypes=frqGenotypes, nOffspring=nOffspring)
			tmpFrq=[frqGenotypesOffspring[geno] for geno in self.genoAll]
			genotypesOffspring=list(np.random.choice(self.genoAll, nOffspring, p=tmpFrq))
			agesOffspring=list(np.repeat(0,nOffspring))
	
			### add those individuals to the new generation
			newGenotypes=genotypes+genotypesOffspring
			newAges=ages+agesOffspring
			
			### calculate the frequencies in both population and sampling population before the viablility fuction.
			frqAfterReprodGenotypes=self.calFrqGenotypes(genotypes=newGenotypes)
			frqGametes=self.calGametesFrq(frqGenotypes=frqAfterReprodGenotypes)
			## if one of the parental haplotypes almost eliminated, break the loop. 
			if(frqGametes['Ab']<0.005 or frqGametes['aB']<0.005):
				break
			
			### output generation
			generation.append(t)
			#### output population gametes frequencies
			g_ab.append(frqGametes['ab'])
			g_Ab.append(frqGametes['Ab'])
			g_aB.append(frqGametes['aB'])
			g_AB.append(frqGametes['AB'])
			## sampling for analyses.
			nInd=len(newGenotypes)
			tmpSample=self.sampling4genotyping(m=300, nPop=nInd, frqGenotypes=frqAfterReprodGenotypes)
			tmp=self.calGametesFrq(frqGenotypes=tmpSample)
			obs_ab.append(int(tmp['ab']*2))
			obs_Ab.append(int(tmp['Ab']*2))
			obs_aB.append(int(tmp['aB']*2))
			obs_AB.append(int(tmp['AB']*2))
			samplingInd.append(tmpSample.values())
			
			### Survival rate/ viability 
			(genotypes, ages)=self.viability(genotypes=newGenotypes, ages=newAges, ageStructure=ageStructure)
		
		### output
		gametes_frq=pd.DataFrame({'Gen':generation, 'gab':g_ab, 'gAb':g_Ab, 'gaB':g_aB, 'gAB':g_AB, 'obs_ab':obs_ab, 'obs_Ab':obs_Ab, 'obs_aB':obs_aB, 'obs_AB':obs_AB})
		samplingInd=pd.DataFrame(samplingInd, columns=self.genoAll)
		samplingInd.insert(0,"Gen", generation, True)
		return gametes_frq, samplingInd
