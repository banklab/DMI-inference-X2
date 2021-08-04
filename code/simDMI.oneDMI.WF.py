#!/usr/local/bin/python3.7
#### I do not know why I am using python3.7

import argparse
from WF_oneDMI import birthDeathProcess ### I changed the age structure manually here. Because I do not want to debug a huge script.
import time
start_time = time.time()

parser = argparse.ArgumentParser(description='A one DMI model')
parser.add_argument("--nPop", default=5000, type=int, help="The population size (default=5000)")
parser.add_argument("--frqMPP", default=0.5, type=float, help="The starting frequency of parental population 1 (AbAb) or the minor parental population. Pop1 is always the minor parental population (MPP) (default=0.5)")
parser.add_argument("--c", default=0.5, type=float, help="The recombination (default=0.5)")
parser.add_argument("--Gtime", default=150, type=int, help="The total generation time (default=150)")
parser.add_argument("--alpha", default=0.001, type=float, help="direct selecion on locus 1 (default=0.001)")
parser.add_argument("--beta", default=0.002, type=float, help="direct selection on locus 2 (default=0.002)")
parser.add_argument("--gamma", default=-0.5, type=float, help="epistasis (default=-0.5)")
parser.add_argument("--m1", default=0.0, type=str, help="migration from pop1 (AbAb) (default=0.0)")
parser.add_argument("--m2", default=0.0, type=str, help="migration from pop2 (aBaB) (default=0.0)")
parser.add_argument("--output", default="tmp.dmi.nonWF", type=str, help="output (default=tmp.dmi.nonWF)")

args = parser.parse_args()

for i in range(1,101):
	data=birthDeathProcess(nPop=args.nPop, c=args.c, frqMPP=args.frqMPP, Gtime=args.Gtime, alpha=args.alpha, beta=args.beta, gamma=args.gamma, m1=args.m1, m2=args.m2)
	(out1, out2)=data.nonWF()
	out1.insert(0, "simulation", [i]*len(out1), True)
	out2.insert(0, "simulation", [i]*len(out2), True)
	if i > 1 :
		out1.to_csv(args.output, sep='\t', index=False, header=False, mode="a")
		out2.to_csv(args.output+".geno", sep='\t', index=False, header=False, mode="a")
	else:
		out1.to_csv(args.output, sep='\t', index=False)
		out2.to_csv(args.output+".geno", sep='\t', index=False)

print("--- %s seconds---" % (time.time() - start_time))
