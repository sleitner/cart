#!/usr/bin/python -Wignore

from re import search
from glob import glob
from scipy import *
from pylab import *

Lbox = float(sys.argv[1])
num_grid = 128.
r0 = Lbox/num_grid

input = open(sys.argv[2],"r")

# skip header lines
for i in xrange(17) :
	input.readline()

x = []
y = []
z = []
rvir = []
rtrunc = []
mvir = []
np = []

for line in input.readlines() :
	cols = line.split()

	x.append(float(cols[1]))
	y.append(float(cols[2]))
	z.append(float(cols[3]))
	rvir.append(float(cols[7]))
	rtrunc.append(float(cols[8]))
	mvir.append(float(cols[9]))
	np.append(int(cols[10]))

input.close()

output = open( sys.argv[3], "w" )
print >>output, len(x)
print >>output, "0.0 0.0 0.0 %f %f %f" % (Lbox,Lbox,Lbox) 
for i in xrange(len(x)) :
    print >>output, x[i], y[i], z[i], rvir[i], rtrunc[i], np[i]
output.close()
