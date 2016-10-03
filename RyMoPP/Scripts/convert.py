#!/usr/bin/python
import scipy.io as sio
import pickle
import numpy as np
import glob as gb
files=gb.glob('../MatWave/wfk*.mat')
for filet in files:
#	print filet
	mat=sio.loadmat(filet)
	wave=mat['wavefunction'][0]
	J=mat['J']
	L=mat['L']
	n=mat['n']
	xmax=mat['xmax'][0][0]
	xmin=mat['xmin'][0][0]
	xstep=mat['xstep'][0][0]
	energy=mat['energy'][0][0]
	name=filet.replace(".mat",".txt")
	name=name.replace("Matwave/","Wave/")
	x=np.arange(xmin,xmax,xstep)
	np.savetxt(name,np.transpose([x,wave]))
	with open(name, "a") as f:
    		f.write(str(energy)+" "+str(xstep)+"\n")
		f.write(str(xmin)+" "+str(xmax))

	

