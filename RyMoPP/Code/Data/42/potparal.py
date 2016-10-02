#!/usr/bin/python
from decimal import *
from numpy import *
import glob as gb
getcontext().prec = 6
files = gb.glob('EigenValue*.txt')
start = 46.000000
names = arange(start,52.600000,0.3)
data = []
for i in names:
	eigen = loadtxt('EigenValue'+ str('{0:.6f}'.format(i)) + '.txt')
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
	
