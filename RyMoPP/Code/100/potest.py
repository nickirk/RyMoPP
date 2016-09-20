#!/usr/bin/python
from numpy import *
import glob as gb
files = gb.glob('EigenValue*.txt')
start = 46.001000
names = arange(46.001000,53.001000,0.2)
data = []
for i in names:
	eigen = [loadtxt('EigenValue'+ str(i) + '000.txt')]
	print eigen
	eigenabs = abs(eigen).tolist()
	eigen = eigen.tolist()
	eigen = [i] + [eigen[eigenabs.index(min(eigenabs))]]
	data.append(eigen)
savetxt('datatest.txt', array(data))
#print names,str(names)
#for filet in files:
#	print filet
	#eigen = loadtxt(filet)
#	data.append([filet,eigen])
	
