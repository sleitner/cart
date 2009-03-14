#!/usr/local/bin/python2.5

from scipy.io.numpyio import fread, fwrite
from pylab import *
import sys

def safe_log10(a, minimum=0.0 ) :
        return log10(a) if ( a ) else minimum

vlog10 = vectorize(safe_log10)

# load data
for slice in ( sys.argv[1:] ) :
	print slice

	input = open( slice, "r" )

	(nx,nz) = fread( input, 2, 'i' )
	density = fread( input, nx*nz, 'f' ).reshape(nz,nx).transpose()
	print min(density.ravel()), max(density.ravel())
	#temperature = fread( input, nx*nz, 'f' ).reshape(nx,nz)
	#velocity_x = fread( input, nx*nz, 'f' ).reshape(nx,nz)
	#velocity_z = fread( input, nx*nz, 'f' ).reshape(nx,nz)

	input.close()

	figure( figsize=(float(nz)/72,float(nx)/72), dpi=72 )
	a = axes([0,0,1,1])
	a.axis('off')
	a.imshow(density,cmap=cm.jet,origin="lower",interpolation="Nearest")
	a.set_autoscale_on(False)
	savefig("density_"+slice.replace("dat","png"),dpi=72)
	clf()
