#!/usr/local/Python-2.7.2/bin/python -Wignore
#!/usr/bin/python

import sys
import matplotlib.pyplot as plt
import numpy as np
import pylab as py

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
            tmp[k][0] = 1
            tmp[k][1] = 1
            tmp[k][2] = 1

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
	freq=4
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
		width=(0.2*py.sqrt(freq)),linewidths=(.2),
#		pivot='mid', 
#		pivot='tail', 
		units='x',
		headaxislength=5 ,
		scale=norm/freq)
            colors=mycmap_rgba(M[::freq, ::freq])
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
		scale=norm/2)

        legend='{:.2g}km/s'.format(norm)
        qk = py.quiverkey(Q, 0.9, 1.05, norm, legend,
                          labelpos='E',
                          fontproperties={'weight': 'bold'})

	
def plot_field(field,clabel,filename,time,vmin,vmax):
	""" purpose: generate plot from square field; args: field,colorbar label, filename"""
	plotfig = plt.figure()
	plotwindow = plotfig.add_subplot(111)
	
	if(arrows==1):
		add_arrows(vx,vy)
	
        if(vmin==vmax):
#            im=plt.imshow(field,cmap=py.cm.gist_yarg,origin="lower",interpolation="Nearest")        
#            im=plt.imshow(field,cmap=py.cm.jet,origin="lower",interpolation=None)        
            im=plt.imshow(field,cmap=py.cm.jet,origin="lower",interpolation="Nearest")        
        else:
#            im=plt.imshow(field,cmap=py.cm.gist_yarg,origin="lower",interpolation="Nearest",vmin=vmin,vmax=vmax) 
#            im=plt.imshow(field,cmap=py.cm.jet,origin="lower",interpolation="gaussian",vmin=vmin,vmax=vmax) 
            im=plt.imshow(field,cmap=py.cm.jet,origin="lower",interpolation="Nearest",vmin=vmin,vmax=vmax) 

#	title='t-t0={:.2g}Myr,dt={:.2g}Myr'.format(time,dtl)
#	title='t-t0={:.2g}Myr'.format(time)
	title='a={:.3g}, t-t0={:.2g}Myr'.format(auni,time)
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

def Ctuple_to_string(word,sizeword):
        new=''
        i=0
        while i < sizeword-1 and word[i]+word[i+1]!='\n' :
            new=new+word[i]
            i=i+1
        return new

def read_header(input):
	(endian)= np.fromfile(file=input,dtype='i4',count=1,sep='') #<i little endian >i big endian
	if(len(endian)==0):
            print 'read everything, next file?'
            return ('','','','','','','')
            
	if(endian!=-99):
            sys.stderr.write('read issue endian=',endian,'expected -99')
            quit()
            
	(sizestring)= np.fromfile(file=input,dtype='i4',count=1,sep='') 
        if(sizestring!=256):
            sys.stderr.write('read issue sizestring=',sizestring,'expected 256')
            quit()
        else:
            (name_slice) = np.fromfile(file=input,dtype='c',count=sizestring,sep='')
            name_slice=Ctuple_to_string(name_slice, sizestring)
        
	(nx,ny)= np.fromfile(file=input,dtype='i4',count=2,sep='')
	box_kpc= np.fromfile(file=input,dtype='f4',count=1,sep='')
	box_kpc=box_kpc[0]
        print 'rounding box size %4.6f for labels' % box_kpc
        if(box_kpc>2):
            box_kpc = round(box_kpc,0)
        elif(box_kpc>0.2):
            box_kpc = round(box_kpc,1)
        elif(box_kpc>0.02):
            box_kpc = round(box_kpc,2)
	auni = np.fromfile(file=input,dtype='f4',count=1,sep='')
	time = np.fromfile(file=input,dtype='f4',count=1,sep='')
	dtl = np.fromfile(file=input,dtype='f4',count=1,sep='')
	auni=auni[0]*1.0
	time=time[0]*1.0
	dtl=dtl[0]*1.0
	print 'name=',name_slice,'nx=',nx,'ny=',ny,'box[kpc]=',box_kpc
	return (name_slice,nx,ny,box_kpc,auni,time,dtl)
	
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



################################# load data
for slice in ( sys.argv[1:] ) :
	print slice

	input = open( slice, "r" )
        (name_slice,nx,ny,box_kpc,auni,time,dtl)=read_header(input)
                
        while (len(name_slice)!=0):
            if(name_slice=='density_numbercc'):
                density = np.log10(np.fromfile(file=input,dtype='f4',count=nx*ny,sep='' ).reshape(nx,ny))
            elif(name_slice=='temperature_kelvin'):
                temp = np.log10(np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny))
            elif(name_slice=='vx_kms'):
                vx = np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny)
            elif(name_slice=='vy_kms'):
                vy = np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny)
            elif(name_slice=='vz_kms'):
                vz = np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny)
            elif(name_slice=='mach_number'):
                mach = (np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny))
            elif(name_slice=='cs_kms'):
                soundspeed = (np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny))
            elif(name_slice=='level_number'):
                level = np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny)
            elif(name_slice=='tauUV_number'):
                tauUV = np.log10(np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny))
            elif(name_slice=='Urad_ergcc'):
                Urad =  np.log10(np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny))
            elif(name_slice=='radPoverP_number'):
                radPoP = np.log10(np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny))
            elif(name_slice=='pressure_ergcc'):
                Pergcc = np.log10(np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny))
            else:
                sys.stderr.write('bad name_slice value',name_slice)
                dummy = np.fromfile(file=input,dtype='f4',count=nx*ny,sep='').reshape(nx,ny)
                
            (name_slice,ndum,ndum,bdum,adum,tdum,dtdum)=read_header(input)
            
        print 'done reading '
        (tick_locs,tick_lbls) = set_axis(nx-1)
	for i in range(nx/2,nx/2+10):
            print 'rho',density[i,nx/2],'vx',vx[i,nx/2], 'vy',vy[i,nx/2], 'vz',vz[i,nx/2]
            #,'mach#:',mach[i,nx/2] ,'radPoP#:',radPoP[i,nx/2]

	input.close()


        
        if 'vx' in locals() and 'vy' in locals():
            arrows=1
        else:
            arrows=0

        arrows=0
        if 'density' in locals():
            #	filename="plots/density/"+slice.replace("dat","png").replace("out/","")
            filename=slice.replace("dat","png").replace("out/","plots/density/")
            print filename
            clabel='log(n [1/cc])'
            vmin=0
            vmax=4
            dpi=400
            plot_field(density,clabel,filename, time,vmin,vmax)
            dpi=100
            #	arrows=0

        if 'vz' in locals():
            filename=slice.replace("dat","png").replace("out/","plots/vz/")
            print filename
            clabel='vz [km/s])'
            vmin=0
            vmax=50
            plot_field(vz,clabel,filename, time,0,0)
        if 'vy' in locals():
            filename=slice.replace("dat","png").replace("out/","plots/vy/")
            print filename
            clabel='vy [km/s])'
            vmin=0
            vmax=50
            plot_field(vy,clabel,filename, time,0,0)
        if 'vx' in locals():
            filename=slice.replace("dat","png").replace("out/","plots/vx/")
            print filename
            clabel='vx [km/s])'
            vmin=0
            vmax=50
            plot_field(vx,clabel,filename, time,0,0)

        if 'mach' in locals():
            filename=slice.replace("dat","png").replace("out/","plots/mach/")
            print filename
            clabel='|v/cs|'
            plot_field(mach,clabel,filename, time,0,0)

        if 'soundspeed' in locals():
            filename=slice.replace("dat","png").replace("out/","plots/cs/")
            print filename
            clabel='cs [km/s]'
            plot_field(soundspeed,clabel,filename, time,0,0)

        if 'temp' in locals():
            filename=slice.replace("dat","png").replace("out/","plots/temp/")
            print filename
            clabel='log(T [K])'
            vmin=2
            vmax=8
            dpi=400
            plot_field(temp,clabel,filename, time,vmin,vmax)
            dpi=100

        if 'level' in locals():
            filename=slice.replace("dat","png").replace("out/","plots/level/")
            print filename
            clabel='level'
            vmin=0
            vmax=0
            plot_field(level,clabel,filename, time,vmin,vmax)

        if 'Pergcc' in locals():
            filename=slice.replace("dat","png").replace("out/","plots/pressure/")
            print filename
            clabel='log(P[erg/cc])'
            vmin=-11
            vmax=-6
            plot_field(Pergcc,clabel,filename, time,vmin,vmax)

        if 'Urad' in locals():
            filename=slice.replace("dat","png").replace("out/","plots/urad/")
            print filename
            clabel='log(Urad[erg/cc])'
            vmin=0
            vmax=0
            plot_field(Urad,clabel,filename, time,vmin,vmax)

        if 'tauUV' in locals():
            filename=slice.replace("dat","png").replace("out/","plots/tauUV/")
            print filename
            clabel='log(tauUV)'
            vmin=0
            vmax=0
            plot_field(tauUV,clabel,filename, time,vmin,vmax)
            
        if 'radPoP' in locals():
            filename=slice.replace("dat","png").replace("out/","plots/prad/")
            print filename
            clabel='log(Prad/P)'
            vmin=0
            vmax=0
            plot_field(radPoP,clabel,filename, time,vmin,vmax)
            

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
