CC=gcc
CFLAGS=-O3 -g
LIBS=-I/usr/local/include -I/usr/include -L/usr/lib64 -L/usr/local/lib -lm -lcfitsio

jedi : jedicatalog jeditransform jedidistort jedigrid_a jedigrid_b jedipaste jediconvolve jedirescale jedinoise

jedicatalog :
	$(CC) $(CFLAGS) sources/jedicatalog.c -o jedicatalog $(LIBS)
	
jeditransform :
	$(CC) $(CFLAGS) sources/jeditransform.c -o jeditransform $(LIBS)
	
jedidistort :
	$(CC) $(CFLAGS) sources/jedidistort.c -o jedidistort $(LIBS)
	
jedigrid_a :
	$(CC) $(CFLAGS) sources/jedigrid_a.c -o jedigrid_a $(LIBS)

jedigrid_b :
	$(CC) $(CFLAGS) sources/jedigrid_b.c -o jedigrid_b $(LIBS)

jedipaste :
	$(CC) $(CFLAGS) sources/jedipaste3.c -o jedipaste $(LIBS)
	
jediconvolve :
	$(CC) $(CFLAGS) sources/jediconvolve4.c -o jediconvolve $(LIBS) -lfftw3f
	
jedirescale :
	$(CC) $(CFLAGS) sources/jedirescale2.c -o jedirescale $(LIBS)

jedinoise :
	$(CC) $(CFLAGS) sources/jedinoise.c -o jedinoise $(LIBS)

