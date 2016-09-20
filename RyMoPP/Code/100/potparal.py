#!/usr/bin/python
from decimal import *
from numpy import *
import glob as gb
getcontext().prec = 6
files = gb.glob('EigenValue*.txt')
start = 48.600000
names = arange(start,51.9000000,0.3)
data = []
for i in names:
	eigen = loadtxt('EigenValue'+ str('{0:.6f}'.format(i)) + '.txt')
	if isinstance(eigen,float):
		eigenabs = [abs(eigen)]
		eigen = [eigen]
	else:
		eigenabs = abs(eigen).tolist()
		eigen = eigen.tolist()
	eigen = [i] + [eigen[eigenabs.index(min(eigenabs))]]
	data.append(eigen)
savetxt('data.txt', array(data))
#print names,str(names)
#for filet in files:
#	print filet
	#eigen = loadtxt(filet)
#	data.append([filet,eigen])
	
