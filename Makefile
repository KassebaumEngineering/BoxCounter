#
# Makefile
# 
# Makefile,v 1.6 1994/05/09 07:39:36 jak Exp
#
# History:
# Makefile,v
# Revision 1.6  1994/05/09  07:39:36  jak
# Added a new tool for linear_redundancy calculations! -jak
#
# Revision 1.5  1994/05/05  04:54:38  jak
# Made fixes and additions.  -jak
#
# Revision 1.4  1994/05/03  17:09:20  jak
# Added a line the the Makefile to remove old tarfiles. -jak
#
# Revision 1.3  1994/05/03  17:05:02  jak
# Added RCS stuff to the Makefile -jak
#
#

CC=gcc
CFLAGS=-g -O
LIBS=-lm
OCTAVEINCLUDE=-I/usr/local/include/octave
OCTAVELIBS=-loctave

SOURCES= Makefile dimensions.cc redundancy.cc BoxCounter.cc BoxCounter.h FastBoxCounter.cc FastBoxCounter.h cal_du.c linear_redundancy.cc scale.cc
OBJECTS= FastBoxCounter.o
EXES=dimensions redundancy linear_redundancy scale

all: $(EXES)

BoxCounter.o: BoxCounter.cc BoxCounter.h
	$(CC) -c $(CFLAGS) BoxCounter.cc -o BoxCounter.o

FastBoxCounter.o: FastBoxCounter.cc FastBoxCounter.h
	$(CC) -c $(CFLAGS) FastBoxCounter.cc -o FastBoxCounter.o

scale: scale.cc $(OBJECTS)
	$(CC) -o scale $(CFLAGS) scale.cc $(OBJECTS) $(LIBS)

dimensions: dimensions.cc $(OBJECTS)
	$(CC) -o dimensions $(CFLAGS) dimensions.cc $(OBJECTS) $(LIBS)
	
redundancy: redundancy.cc $(OBJECTS)
	$(CC) -o redundancy $(CFLAGS) redundancy.cc $(OBJECTS) $(LIBS)
	
linear_redundancy: linear_redundancy.cc $(OBJECTS)
	$(CC) -o linear_redundancy $(CFLAGS) linear_redundancy.cc $(LIBS)

clean:
	rm -f $(EXES) $(OBJECTS) BoxCounter.o
	
cal_du: cal_du.c
	$(CC) -o cal_du cal_du.c

tarfile: $(SOURCES)
	rm -f BoxCounter.tar.bz2
	mkdir BoxCounter
	cp $(SOURCES) BoxCounter
	tar cf - BoxCounter | bzip2 > BoxCounter.tar.bz2
	rm -rf BoxCounter
	touch tarfile

install: $(EXES)
	cp $(EXES) $(HOME)/Apps
