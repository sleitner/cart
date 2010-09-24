#!/usr/bin/python -Wignore

import matplotlib
matplotlib.use('Agg')

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

from numpy import sqrt, exp, log10, max, min, outer, vectorize, ones, array, arange, var, mean
from scipy.io.numpyio import fread, fwrite
from pylab import *
import sys
from matplotlib.patches import Circle
from re import search

def safe_log10(a, minimum=0.0 ) :
	if ( a ) :
		return log10(a)
	else :
		return minimum

vlog10 = vectorize(safe_log10)

min_value = 1.e50
max_value = 0.0

for file in sys.argv[1:] :
	input = open( file,"rb")
        (num_pixels) = fread( input, 1, 'i' )[0]
        (Lbox) = fread( input, 1, 'd' )[0]
        pointing = fread( input, 3, 'd' )
        (width,depth) = fread( input, 2, 'd' )
        data = fread( input, num_pixels**2, 'd' )
        input.close()

	local_min_value = min(data.take(data.nonzero()))
	local_max_value = max(data.take(data.nonzero()))

	print file, local_min_value, local_max_value 

	if ( local_min_value < min_value ) :
		min_value = local_min_value

	if ( local_max_value > max_value ) :
		max_value = local_max_value

	del data

print min_value, max_value

for file in sys.argv[1:] :
	print file

	dimension = search("_(x|y|z)_",file).group(1)
	aexpn = search("_(\d\.\d\d\d\d)_",file).group(1)

	input = open( file,"rb")
	(num_pixels) = fread( input, 1, 'i' )[0]
	(Lbox) = fread( input, 1, 'd' )[0]
	pointing = fread( input, 3, 'd' )
	(width,depth) = fread( input, 2, 'd' )
	data = fread( input, num_pixels**2, 'd' )
	input.close()

	data = data.reshape(num_pixels,num_pixels).transpose()
	#local_min_value = min(data.take(data.ravel().nonzero()))
	#local_max_value = max(data.take(data.ravel().nonzero()))
	#print local_min_value, local_max_value
	#min_value = 10.**min_value
	#min_value = 10.**(log10(max_value) - 5.0)

	#data = vlog10( data, log10(min_value)-2. )
	data = vlog10( data, log10(min_value) )

	figure( figsize=(num_pixels/150,num_pixels/150) )
	a = axes([0,0,1,1],frameon=False)
	a.axis('off')

	if ( file.find("dm") != -1 ) :
		a.imshow(data,cmap=cm.bone,origin="lower",vmin=log10(min_value)-2,vmax=log10(max_value),interpolation="Nearest")
	elif ( file.find("star") != -1 ) :
		imshow(data,cmap=cm.copper,origin="lower",vmin=log10(min_value)-1,vmax=log10(max_value),interpolation="Nearest")
	elif ( file.find("gas") != -1 ) :
		imshow(data,cmap=cm.jet,origin="lower",vmin=log10(min_value),vmax=log10(max_value),interpolation="Nearest")
	elif ( file.find("S") != -1 ) :
		imshow(data,cmap=cm.jet,origin="lower",vmin=log10(min_value),vmax=log10(max_value),interpolation="Nearest")
	else :
		imshow(data,cmap=cm.hsv,origin="lower",vmin=log10(min_value),vmax=log10(max_value),interpolation="Nearest")

	a.set_autoscale_on(False)

	if dimension == "x" :
		figtext(0.05, 0.9, aexpn, color='w',size=24)

	savefig( file.replace("dat","png"),dpi=150 )
	del data
	clf()
	close()
