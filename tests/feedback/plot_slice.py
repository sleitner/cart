#!/usr/local/Python-2.7.2/bin/python
#!/usr/bin/python -Wignore
#!/usr/bin/python2.4

def plot_field(field,clabel,filename,time):
	""" purpose: generate plot from square field; args: field,colorbar label, filename"""
	plotfig = plt.figure()
	plotwindow = plotfig.add_subplot(111)
# 	im=plt.imshow(field,cmap=py.cm.gist_yarg,origin="lower",interpolation="Nearest")
 	im=plt.imshow(field,cmap=py.cm.jet,origin="lower",interpolation="Nearest")
	title='t-t0={:.2g}Myr'.format(time)
	py.title(title)
	plt.xticks(tick_locs,tick_lbls)
	plt.yticks(tick_locs,tick_lbls)
	py.xlabel('x')
	py.ylabel('y')
	cax = py.axes([0.85, 0.1, 0.05, 0.8]) #xleft, ybot, xsz, ysz
	py.colorbar(im,cax=cax)
	py.ylabel(clabel)
	plt.savefig(filename,dpi=100)
	plt.clf()


#quit()

import matplotlib.pyplot as plt
import numpy as np
import pylab as py
import sys


dpi=1024
# load data
for slice in ( sys.argv[1:] ) :
	print slice

	input = open( slice, "r" )
	(endian)= np.fromfile(file=input,dtype='i4',count=1,sep='') #<i little endian >i big endian
	(nx,ny)= np.fromfile(file=input,dtype='i4',count=2,sep='')
	box_kpc= np.fromfile(file=input,dtype='f4',count=1)
	box_kpc=box_kpc[0]
	time = np.fromfile(file=input,dtype='f4',count=1)
	time = np.fromfile(file=input,dtype='f4',count=1)
	time=time[0]*1.0
	print 'nx=',nx,'ny=',ny,'box[kpc]=',box_kpc

	middle=nx/2
	nticks= 5
	dticks=1.0/(nticks-1.0)
	tick_locs = range(0,nticks,1) 
	tick_lbls = range(0,nticks,1) 
	for i in range(len(tick_lbls)):
#		tick_lbls[i] = int(round(tick_lbls[i]*dticks*nx,0))
		tick_locs[i] = tick_locs[i]*dticks*nx
		tick_lbls[i] = '{:.2}'.format(tick_lbls[i]*dticks*box_kpc)

	density = np.log10(np.fromfile(file=input,dtype='f4', count=nx*ny, ).reshape(ny,nx).transpose())
	temp = np.fromfile(file=input,dtype='f4', count=nx*ny, ).reshape(ny,nx).transpose()
	temp=np.log10(temp+1)
	vx = np.fromfile(file=input,dtype='f4', count=nx*ny, ).reshape(ny,nx).transpose()
	vz = np.fromfile(file=input,dtype='f4', count=nx*ny, ).reshape(ny,nx).transpose()
	pressure = np.log10(np.fromfile(file=input,dtype='f4', count=nx*ny, ).reshape(ny,nx).transpose())
	level = np.fromfile(file=input,dtype='f4', count=nx*ny, ).reshape(ny,nx).transpose()
	input.close()
	


#colorbar whose height (or width) in sync with the master axes
#
#import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import numpy as np
#
#ax = plt.subplot(111)
#im = ax.imshow(np.arange(100).reshape((10,10)))
#
## create an axes on the right side of ax. The width of cax will be 5%
## of ax and the padding between cax and ax will be fixed at 0.05 inch.
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)
#
#plt.colorbar(im, cax=cax)


###############################@@@@@@@@@@

	#fig = figure( figsize=(float(ny)/dpi,float(nx)/dpi), dpi=dpi )
	filename="plots/density"+slice.replace("dat","png").replace("dumps","")
	print filename
	clabel='log(n [1/cc])'
	plot_field(density,clabel,filename, time)

	filename="plots/temp"+slice.replace("dat","png").replace("dumps","")
	print filename
	clabel='log(T [K])'
	plot_field(temp,clabel,filename, time)
	
	filename="plots/level"+slice.replace("dat","png").replace("dumps","")
	print filename
	clabel='level'
	plot_field(level,clabel,filename, time)

