#!/sw/bin/python2.6 -Wignore

from pylab import *
from scipy import loadtxt

import sys
from re import search
from glob import glob

output_dir = "DAT"

for output in glob(output_dir+'/sod_*.dat') :
	(x,rho,momentum,pressure) = loadtxt(output,unpack=True)

	figure(figsize=(9,3))
	subplots_adjust(left=0.08,bottom=0.12,top=0.95, right=0.95,wspace=0.3)

	subplot(131)
	plot(x,rho,marker="o")
	gca().set_ylim(0,1.1)
	gca().set_xlabel("x")
	gca().set_ylabel("density")

	subplot(132)
	plot(x,momentum/rho,marker="o")
	gca().set_ylim(-1.0,1.0)
	gca().set_xlabel("x")
	gca().set_ylabel("vx")

	subplot(133)
	plot(x,pressure,marker="o")
	gca().set_ylim(0,1.1)
	gca().set_xlabel("x")
	gca().set_ylabel("Pressure")

	file = output.replace(".dat",".png")
	savefig(file,dpi=100)
