#!/usr/bin/python
from numpy import *
import matplotlib.pyplot as plt
import glob as gb
files = gb.glob('EigenValue*.txt')
start = 20.001000
names = arange(20.001000,36.001000,0.5)
data = []
for i in names:
	eigen = loadtxt('EigenValue'+ str(i) + '000.txt')
	eigen = eigen.tolist()
	eigen = [i] + eigen
	data.append(eigen)
savetxt('data.txt', array(data))
#print names,str(names)
#for filet in files:
#	print filet
	#eigen = loadtxt(filet)
#	data.append([filet,eigen])
	
