#!/usr/bin/python2.4
#!/usr/bin/python2.7 -Wignore

from numpy import *
import numpy as np
from pylab import *
import sys

# load data
for slice in ( sys.argv[1:] ) :
	print slice

	input = open( slice, "r" )
#	endian= np.fromfile(file=input,dtype='i4',count=1,sep='') #<i little endian >i big endian
	(nx,nz)= np.fromfile(file=input,dtype='i4',count=2,sep='') #<i little endian >i big endian
	density = np.fromfile(file=input,dtype='f4', count=nx*nz, ).reshape(nz,nx).transpose()
	int_en = np.fromfile(file=input,dtype='f4', count=nx*nz, ).reshape(nz,nx).transpose()
	momx = np.fromfile(file=input,dtype='f4', count=nx*nz, ).reshape(nz,nx).transpose()
	momz = np.fromfile(file=input,dtype='f4', count=nx*nz, ).reshape(nz,nx).transpose()
	input.close()

	figure( figsize=(float(nz)/72,float(nx)/72), dpi=72 )
	a = axes([0,0,1,1])
	a.axis('off')
	a.imshow(density,cmap=cm.jet,origin="lower",interpolation="Nearest",vmin=0.001,vmax=3.31986)
	a.set_autoscale_on(False)
	savefig("density_"+slice.replace("dat","png"),dpi=288)
	clf()
	
	figure( figsize=(float(nz)/72,float(nx)/72), dpi=72 )
	a = axes([0,0,1,1])
	a.axis('off')
	a.imshow(int_en,cmap=cm.jet,origin="lower",interpolation="Nearest",vmin=0.001,vmax=3.31986)
	a.set_autoscale_on(True)
	savefig("temp_"+slice.replace("dat","png"),dpi=288)
	clf()

	figure( figsize=(float(nz)/72,float(nx)/72), dpi=72 )
	a = axes([0,0,1,1])
	a.axis('off')
	a.imshow(momx,cmap=cm.jet,origin="lower",interpolation="Nearest",vmin=0.001,vmax=3.31986)
	a.set_autoscale_on(True)
	savefig("momx_"+slice.replace("dat","png"),dpi=288)
	clf()

	figure( figsize=(float(nz)/72,float(nx)/72), dpi=72 )
	a = axes([0,0,1,1])
	a.axis('off')
	a.imshow(momz,cmap=cm.jet,origin="lower",interpolation="Nearest",vmin=0.001,vmax=3.31986)
	a.set_autoscale_on(True)
	savefig("momz_"+slice.replace("dat","png"),dpi=288)
	clf()

