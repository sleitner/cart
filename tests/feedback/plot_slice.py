#!/usr/bin/python2.4
#!/usr/bin/python2.7 -Wignore

from numpy import *
import numpy as np
from pylab import *
import sys


dpi=1024
# load data
for slice in ( sys.argv[1:] ) :
	print slice

	input = open( slice, "r" )
	endian= np.fromfile(file=input,dtype='i4',count=1,sep='') #<i little endian >i big endian
	(nx,nz)= np.fromfile(file=input,dtype='i4',count=2,sep='') #<i little endian >i big endian
	print nx,nz
	middle=nx/2
	density = np.fromfile(file=input,dtype='f4', count=nx*nz, ).reshape(nz,nx).transpose()
	int_en = np.fromfile(file=input,dtype='f4', count=nx*nz, ).reshape(nz,nx).transpose()
	momx = np.fromfile(file=input,dtype='f4', count=nx*nz, ).reshape(nz,nx).transpose()
	momz = np.fromfile(file=input,dtype='f4', count=nx*nz, ).reshape(nz,nx).transpose()
	pressure = np.fromfile(file=input,dtype='f4', count=nx*nz, ).reshape(nz,nx).transpose()
	level = np.fromfile(file=input,dtype='f4', count=nx*nz, ).reshape(nz,nx).transpose()
	input.close()

	figure( figsize=(float(nz)/dpi,float(nx)/dpi), dpi=dpi )
	a = axes([0,0,1,1])
	a.axis('off')
	a.imshow(density,cmap=cm.gist_yarg,origin="lower",interpolation="Nearest")
	a.set_autoscale_on(True)
	filename="plots/density"+slice.replace("dat","png").replace("dumps","")
	print filename
	savefig(filename,dpi=dpi)
	clf()
	
	figure( figsize=(float(nz)/dpi,float(nx)/dpi), dpi=dpi )
	a = axes([0,0,1,1])
	a.axis('off')
	a.imshow(int_en,cmap=cm.gist_yarg,origin="lower",interpolation="Nearest")
	a.set_autoscale_on(True)
	filename="plots/temp"+slice.replace("dat","png").replace("dumps","")
	print filename
	savefig(filename,dpi=dpi)
	clf()

	figure( figsize=(float(nz)/dpi,float(nx)/dpi), dpi=dpi )
	a = axes([0,0,1,1])
	a.axis('off')
	a.imshow(momx,cmap=cm.gist_yarg,origin="lower",interpolation="Nearest")
	a.set_autoscale_on(True)
	filename="plots/momx"+slice.replace("dat","png").replace("dumps","")
	print filename
	savefig(filename,dpi=dpi)
	clf()

	figure( figsize=(float(nz)/dpi,float(nx)/dpi), dpi=dpi )
	a = axes([0,0,1,1])
	a.axis('off')
	a.imshow(momz,cmap=cm.gist_yarg,origin="lower",interpolation="Nearest") #,vmin=0.001,vmax=3.31986)
	a.set_autoscale_on(True)
	filename="plots/momz"+slice.replace("dat","png").replace("dumps","")
	print filename
	savefig(filename,dpi=dpi)
	clf()

	figure( figsize=(float(nz)/dpi,float(nx)/dpi), dpi=dpi )
	a = axes([0,0,1,1])
	a.axis('off')
	a.imshow(pressure,cmap=cm.gist_yarg,origin="lower",interpolation="Nearest")
	a.set_autoscale_on(True)
	filename="plots/pressure"+slice.replace("dat","png").replace("dumps","")
	print filename
	savefig(filename,dpi=dpi)
	clf()

	figure( figsize=(float(nz)/dpi,float(nx)/dpi), dpi=dpi )
	a = axes([0,0,1,1])
	a.axis('off')
	a.imshow(level,cmap=cm.gist_yarg,origin="lower",interpolation="Nearest", vmin=0,vmax=4)
	print level[middle:(middle+20),middle]
	a.set_autoscale_on(True)
	filename="plots/level"+slice.replace("dat","png").replace("dumps","")
	print filename
	savefig(filename,dpi=dpi)
	clf()
		
