#
# Makefile for smooth.
#
CC = pgcc
CFLAGS	= -O0 -mcmodel=medium
LIBS	= -lm

default:	smooth

clean:
	rm -f *.o smooth

smooth: main.o kd.o smooth.o
	$(CC) $(CFLAGS) -o smooth main.o kd.o smooth.o $(LIBS)

main.o: main.c kd.h smooth.h
	$(CC) $(CFLAGS) -c main.c

kd.o: kd.c kd.h tipsydefs.h
	$(CC) $(CFLAGS) -c kd.c

smooth.o: smooth.c kd.h smooth.h
	$(CC) $(CFLAGS) -c smooth.c



