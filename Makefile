COMPILER	=	gcc
CFLAGS		=	-O2 -fopenmp
LFLAGS		=	-lm
OBJECTS		=	sphlib.o sph.o utility.o 
INCLUDES	=	utility.h sphlib.h

sph:		$(OBJECTS)
	$(COMPILER) $(CFLAGS) -o sph $(OBJECTS) $(LFLAGS)

sph.o:	$(INCLUDES) sph.c
	$(COMPILER) $(CFLAGS) -c sph.c $(LFLAGS)

	
sphlib.o:		sphlib.c $(INCLUDES)
	$(COMPILER) $(CFLAGS) -c sphlib.c $(LFLAGS)

utility.o:		utility.c utility.h
	$(COMPILER) $(CFLAGS) -c utility.c $(LFLAGS)
	

clean:
	rm -f *.o *~

realclean: clean
	rm -f sph
