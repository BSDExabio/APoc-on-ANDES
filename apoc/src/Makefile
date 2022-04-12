##############################
#### Fortran compiler
##############################

### Intel fortran
FC     = ifort
FFLAGS = -extend_source 132 -real-size 64 -fast # equivalent to using -O3
#FFLAGS = -extend_source 132 -real-size 64 -O2
#FFLAGS = -extend_source 132 -real-size 64 -g -check all
#FFLAGS = -extend_source 132 -real-size 64 -O0

#### PGF fortran
#FC	= pgf90
#FFLAGS	= -Mextend -fastsse -O3 -Mr8
##FFLAGS	= -Mextend -Bstatic -Mr8 -C
##FFLAGS = -Mextend -fastsse -O3 -Mr8 -mcmodel=medium
##FFLAGS = -Mextend -Mr8 -mcmodel=medium -C

#############
PROGRAMS = apoc
P_OBJS	= palign_batch.o scoring_batch.o

all:	$(PROGRAMS)
 
#############

apoc:	$(P_OBJS)
	$(FC) $(FFLAGS) -o $@ $(P_OBJS)


#############
clean:
	rm -f *.o *.mod $(PROGRAMS)
#############


