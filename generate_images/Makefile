#SHELL = /bin/bash
#FFLAG = -O2 -fast -Mipa=fast
FFLAG = -O2 
#IDIR  = -I/home/maad-huangyuhan/fftw-2.1.5/include
#IDIR = -I/opt/software/fftw2-2.1.5-intel/include
#IDIR = -I/share/base/fftw/2.1.5/include
#LDIR  = -L/home/maad-huangyuhan/fftw-2.1.5/lib
#LDIR = -L/opt/software/fftw2-2.1.5-intel/lib
#LDIR = -L/share/base/fftw/2.1.5/lib
FCOMP = mpif90 -c ${FFLAG} ${IDIR}
LINK  = mpif90 -mkl
#LIBS  = -ldfftw -ldrfftw #-ldfftw_mpi -ldrfftw_mpi 
#LIBS = -ldrfftw -ldfftw -lm -r8
#LIBS = -lrfftw -lfftw -lm -r8
OBJ   = post_poiseuilleflow.o
PROJ = post_poiseuille

.SUFFIXES: .o .f90

.f90.o:
	${FCOMP} $*.f90

turb:  ${OBJ}
	${LINK}  ${OBJ} ${LDIR} ${LIBS} -o ${PROJ}

clean:
	rm -f core* output fort.* *.o   *~ ${PROJ} *.mod


