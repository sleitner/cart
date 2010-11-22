PROG =	hf.x

SRCS =	hfind.f

OBJS =	hfind.o

LIBS =	

FC = pgf77
FFLAGS = -fast -mp -mcmodel=medium
LDFLAGS = -fast -mp -mcmodel=medium

# BulldogM
#FC = ifort
#FFLAGS = -O3 -xHOST -m64 -fpic -fpie -mcmodel=large -shared-intel -fp-model precise
3LDFLAGS = -O3 -xHOST -m64 -fpic -fpie -mcmodel=large -shared-intel -fp-model precise

# gerneic ifort
#FC = ifort
#FFLAGS = -O2 -w -fp-model precise -mcmodel=medium -i_dynamic 
#LDFLAGS = -O2 -w -fp-model precise -mcmodel=medium -i_dynamic 

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

hfind.o: hfind.h 
