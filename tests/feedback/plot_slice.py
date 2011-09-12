#!/usr/local/Python-2.7.2/bin/python -Wignore
#!/usr/bin/python

import sys
import matplotlib.pyplot as plt
import numpy as np
import pylab as py
import sys

def mycmap_rgba(x):
    #scale for colormap
    lo=vmin
    hi=vmax
    if(lo == hi):
        hi=(x).max()
        lo=(x).min()
        
    dx=hi-lo
    x=x-lo
    x=x/dx
    transparency_function=1
    tmp=x.reshape( (x.shape[0]*x.shape[1]) )
#THIS COLORMAP IS IGNORED -- RESET TO BLACK BELOW!
    tmp = py.cm.jet(tmp)

    #scale for transparency
    hi=(x).max()
    lo=(x).min()
    dx=hi-lo
    x=x-lo
    x=x/dx
    for i in xrange(x.shape[0]):
        for j in xrange(x.shape[1]):
            k=j+(x.shape[0]*i)
            tmp[k][3] = x[i,j]
#HERE! make black
            tmp[k][0] = 0
            tmp[k][1] = 0
            tmp[k][2] = 0

    return tmp

#input 2 mesh output 2+4 mesh+rgba
def mycmap(x):
    lo=vmin
    hi=vmax
    if(lo == hi):
        hi=(x).max()
        lo=(x).min()
        
    dx=hi-lo
    x=x-lo
    x=x/dx
    transparency_function=1
    tmp = py.cm.jet(x)
#    for i in xrange(tmp.shape[0]):
#        for j in xrange(tmp.shape[1]):
#            if(x[i,j]>2):
#                tmp[i,j][3] = 1.
#            else:
#                tmp[i,j][3] = 1.
    return tmp

def findmax(a):
    if len(a) == 0:
        return 0
    curr_max = a[0]
    for i in a:
        if i > curr_max:
            curr_max = i
    return curr_max

def add_arrows(fieldx,fieldy):
	norm=abs(np.concatenate((fieldx,fieldy)).max())/1. #1d norm
	freq=1
	M = py.sqrt(pow(fieldx, 2) + pow(fieldy, 2))
	X,Y = py.meshgrid( py.arange(0,nx,1),py.arange(0,ny,1) )

        if(1==1):
            Q = py.quiver( 
		X[::freq, ::freq], Y[::freq, ::freq], 
		fieldx[::freq, ::freq], fieldy[::freq, ::freq], 
                #		M[::freq, ::freq],cmap=py.cm.jet,
                #		,color='yellow'
		#M[::freq, ::freq],cmap=py.cm.gray, alpha=0.2,
#		facecolors=mycmap(M[::freq, ::freq]),
                #		width=0.4,linewidths=(.5), edgecolors=('k'), 
#                color=colors,
		width=0.2,linewidths=(.2),
		pivot='mid', 
		units='x',
		headaxislength=5 ,
		scale=norm)
            colors=mycmap_rgba(M)
            Q.set_facecolor(colors)
            Q.set_edgecolor('k')
            Q.set_edgecolor(colors)
        else:
            Q = py.quiver( 
		X[::freq, ::freq], Y[::freq, ::freq], 
		fieldx[::freq, ::freq], fieldy[::freq, ::freq], 
                #		M[::freq, ::freq],cmap=py.cm.gray, alpha=1,
		width=0.2,
		units='x',
		headaxislength=5 ,
		scale=norm)

        legend='{:.2g}km/s'.format(norm)
        qk = py.quiverkey(Q, 0.9, 1.05, 1, legend,
                          labelpos='E',
                          fontproperties={'weight': 'bold'})

#
#		im=plt.imshow(field,cmap=py.cm.jet,origin="lower",interpolation="Nearest",vmin=vmin,vmax=vmax)
	
def plot_field(field,clabel,filename,time,vmin,vmax):
	""" purpose: generate plot from square field; args: field,colorbar label, filename"""
	plotfig = plt.figure()
	plotwindow = plotfig.add_subplot(111)
	
	if(arrows==1):
		add_arrows(vx,vy)
	
        if(vmin==vmax):
            im=plt.imshow(field,cmap=py.cm.jet,origin="lower",interpolation="Nearest")        
        else:
            im=plt.imshow(field,cmap=py.cm.jet,origin="lower",interpolation="Nearest",vmin=vmin,vmax=vmax) 

	title='t-t0={:.2g}Myr,dt={:.2g}Myr'.format(time,dtl)
	py.title(title)
	plt.xticks(tick_locs,tick_lbls)
	plt.yticks(tick_locs,tick_lbls)
	py.xlabel('x[Kpc]')
	py.ylabel('y[Kpc]')
	cax = py.axes([0.85, 0.1, 0.05, 0.8]) #xleft, ybot, xsz, ysz
	py.colorbar(im,cax=cax)
	py.ylabel(clabel)
	plt.savefig(filename,dpi=dpi)
	plt.clf()

def read_header(input):
	(endian)= np.fromfile(file=input,dtype='i4',count=1,sep='') #<i little endian >i big endian
	if(endian!=-99):
		sys.stderr.write('read issue endian=',endian,'expected -99')
		quit()

	(nx,ny)= np.fromfile(file=input,dtype='i4',count=2,sep='')
	box_kpc= np.fromfile(file=input,dtype='f4',count=1,sep='')
	box_kpc=box_kpc[0]
	time = np.fromfile(file=input,dtype='f4',count=1,sep='')
	dtl = np.fromfile(file=input,dtype='f4',count=1,sep='')
	time=time[0]*1.0
	dtl=dtl[0]*1.0
	print 'nx=',nx,'ny=',ny,'box[kpc]=',box_kpc
	return (nx,ny,box_kpc,time,dtl)
	
def set_axis(np):
	middle=np/2
	nticks= 5
	dticks=1.0/(nticks-1.0)
	tick_locs = range(0,nticks,1) 
	tick_lbls = range(0,nticks,1) 
	for i in range(len(tick_lbls)):
		tick_locs[i] = tick_locs[i]*dticks*np
		tick_lbls[i] = tick_lbls[i]*dticks*box_kpc
		tick_lbls[i] = '%.2f' % tick_lbls[i]
	return (tick_locs, tick_lbls)

# load data
for slice in ( sys.argv[1:] ) :
	print slice

	input = open( slice, "r" )


	(nx,ny,box_kpc,time,dtl)=read_header(input)
	
	(tick_locs,tick_lbls) = set_axis(nx-1)
	density = np.log10(np.fromfile(file=input,dtype='f4',count=nx*ny,sep='' ).reshape(nx,ny))

	read_header(input)
	temp = np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny)
	read_header(input)
	vx = np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny)
	read_header(input)
	vy = np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny)
	read_header(input)
	vz = np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny)

	read_header(input)
	mach = (np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny))

	read_header(input)
	soundspeed = (np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny))

	read_header(input)
	level = np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny)

	for i in range(nx/2,nx/2+10):
		print 'rho',density[i,nx/2],'vx',vx[i,nx/2], 'vy',vy[i,nx/2], 'vz',vz[i,nx/2],'mach#:',mach[i,nx/2]
		
		

	input.close()

	filename="plots/density/"+slice.replace("dat","png").replace("out/","")
	print filename
	clabel='log(n [1/cc])'
	vmin=0
	vmax=4
	arrows=1
	dpi=400
	plot_field(density,clabel,filename, time,vmin,vmax)
	dpi=100
#	arrows=0


	filename="plots/vx/"+slice.replace("dat","png").replace("out/","")
	print filename
	clabel='vx [km/s])'
	vmin=0
	vmax=50
	plot_field(vx,clabel,filename, time,0,0)

	filename="plots/mach/"+slice.replace("dat","png").replace("out/","")
	print filename
	clabel='|v/cs|'
	plot_field(mach,clabel,filename, time,0,0)

	filename="plots/cs/"+slice.replace("dat","png").replace("out/","")
	print filename
	clabel='cs [km/s]'
	plot_field(soundspeed,clabel,filename, time,0,0)

	filename="plots/temp/"+slice.replace("dat","png").replace("out/","")
	print filename
	clabel='log(T [K])'
	vmin=2
	vmax=8
	plot_field(temp,clabel,filename, time,vmin,vmax)
	
	filename="plots/level/"+slice.replace("dat","png").replace("out/","")
	print filename
	clabel='level'
	vmin=0
	vmax=10
	plot_field(level,clabel,filename, time,vmin,vmax)


#	varr = arrow(pos=ball.pos, axis=vscale*ball.velocity, color=color.yellow)
#	X,Y = py.meshgrid( py.arange(0,2*py.pi,(2*py.pi/nx)),py.arange(0,2*py.pi,(2*py.pi/nx)) )
#	U = py.cos(X)
#	V = py.sin(Y)
#	add_arrows(U,V)
#	l,r,b,t = py.axis()
#	dx, dy = r-l, t-b
#	py.axis([l-0.05*dx, r+0.05*dx, b-0.05*dy, t+0.05*dy])
#	py.axis([-1, 7, -1, 7])


#		tick_lbls[i] = int(round(tick_lbls[i]*dticks*np,0))
	#fig = figure( figsize=(float(ny)/dpi,float(nx)/dpi), dpi=dpi )

#	title='t-t0=',time,'Myr'
#		add_arrows(vx,vy,box_kpc*1000) #arrow length km/s~pc/Myr->box/Myr
#	X,Y = py.meshgrid( py.arange(0,2*py.pi,(2*py.pi/nx)),py.arange(0,2*py.pi,(2*py.pi/nx)) )
#	vx = py.cos(X)
#	vy = py.sin(Y)
