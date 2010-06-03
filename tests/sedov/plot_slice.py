#!/sw/bin/python2.6 -Wignore

from scipy.io.numpyio import fread, fwrite
from pylab import *
import sys

# load data
for slice in ( sys.argv[1:] ) :
	print slice

	input = open( slice, "r" )

	(nx,nz) = fread( input, 2, 'i' )
	density = fread( input, nx*nz, 'f' ).reshape(nz,nx).transpose()
	input.close()

	figure( figsize=(float(nz)/72,float(nx)/72), dpi=72 )
	a = axes([0,0,1,1])
	a.axis('off')
	a.imshow(density,cmap=cm.jet,origin="lower",interpolation="Nearest",vmin=0.001,vmax=3.31986)
	a.set_autoscale_on(False)
	savefig("density_"+slice.replace("dat","png"),dpi=72)
	clf()