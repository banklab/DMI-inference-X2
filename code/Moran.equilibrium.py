import numpy as np

## WF_oneDMI.py is an overlapping generation model. I changed the age structure in this script to a Wright-Fisher process. We need it to simulate some scinarios like the natural hybrid population.
## Therefore, we have this script, here. 
## With a age structure, calculate the equalibrium. How many age 0 individuals should be produced at one generation? See Crow and Kimura 1973 Chapter 1.3, Model 3: overlapping generations, discrete time intervals.

def moran(birth, death):
	birthDeath=[birth]
	survival=1-death
	print(survival)
	gen=len(death)
	for i in range(gen-1):
		tmp=np.repeat(0,i)
		tmp=np.append(tmp,survival[i])
		tmp=np.append(tmp,np.repeat(0,gen-i-1))
		tmp=[tmp]
		birthDeath=birthDeath+tmp
	
	print(birthDeath)
	value,vector=np.linalg.eig(birthDeath)
	lTotal=1
	old=1
	for i in range(gen-1):
		old=survival[i]*value[0]*old
		lTotal=lTotal+old
	
	print(1/lTotal)
	return(value[0])

#AS=np.array([0.2,0.2,0,0,0,0,0.25,0.5,0.75,1])
AS=np.array([0,1,1,1,1,1,1,1,1,1])
birth=np.array(np.repeat(0.2,10))
x=moran(birth=birth, death=AS)
print(x)
