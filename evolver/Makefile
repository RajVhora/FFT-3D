# Author - MPGuruajan. Last update 1.12.04

# According to make manual, every makefile is supposed to contain this
# line

SHELL = /bin/sh

cc = /usr/bin/gcc

# Suffi(ces?)-xes of files used in this makefile 

.SUFFIXES: .c .h .o

# Compilation and linking options
#
# The following options are needed only while testing and optimising
#
#
#COMPOPS = -fprofile-arcs -ftest-coverage -pg -Wall -O3 \
#-mcpu=pentiumpro
#LINKOPS = -O3 -mcpu=pentiumpro -lgsl -lgslcblas -lfftw3 -lm -pg 
#
# The following are for debugging
#
#COMPOPS = -Wall -O3 -mcpu=pentiumpro -g
#LINKOPS = -O3 -mcpu=pentiumpro -lgsl -lgslcblas -lfftw3 -lm -g
#COMPOPS = -Wall -O3 -mcpu=x86-64 -march=x86-64 -m64 -g
#LINKOPS = -O3 -mcpu=x86-64 -march=x86-64 -m64 \
# -lgsl -lgslcblas -lfftw3 -lm -g
#
# The following are for production runs - The AMD and P4 versions
#
COMPOPS = -Wall -O3  -U APPROX_CALC -D debug
#LINKOPS = -O3 -static -LLIBDIR=/usr/local/lib -lgsl -lgslcblas  -lfftw3 -lfftw3_threads -#lpthread -lm
LINKOPS = -lgsl -lgslcblas  -lfftw3 -lm
#COMPOPS = -Wall -O3 -mcpu=x86-64 -march=x86-64 -m64
#LINKOPS = -O3 -mcpu=x86-64 -march=x86-64 -m64 -lgsl -lgslcblas -lfftw3 -lm 
 


# List of object files

objects = MUSE_EVOLVER.o nrutil.o generate_elast_tensors.o \
calculate_sigma_T.o calculate_omega.o \
evolve_inhom_constrained.o calculate_Del_sigma_T.o calculate_uzero.o \
refine_u.o calculate_hom_strain.o calculate_alpha.o calculate_beta.o \
calculate_S_exp.o calculate_sig_zero_exp.o refine_u_chen.o \
calculate_alpha_prime.o calculate_beta_prime.o write_ppfiles.o \
add_noise.o evolve_diffusional.o \
calculate_S.o calculate_aleph.o 

# List of header files

headers = stdio.h stdlib.h math.h complex.h fftw3.h \
functions.h nrutil.h time.h float.h 

# List of source files

sources = muse_evolver.c nrutil.c generate_elast_tensors.c \
calculate_sigma_T.c calculate_omega.c \
evolve_inhom_constrained.c calculate_Del_sigma_T.c calculate_uzero.c \
refine_u.c calculate_hom_strain.c calculate_alpha.c calculate_beta.c \
calculate_S_exp.c calculate_sig_zero_exp.c refine_u_chen.c \
calculate_alpha_prime.c calculate_beta_prime.c write_ppfiles.c \
add_noise.c evolve_diffusional.c calculate_S.c calculate_aleph.c 


# Directory paths for the source and header files 

vpath %.c ./src
vpath %.c ./src_test
vpath %.h ./headers
vpath %.h /usr/lib/gcc/i686-linux-gnu/4.8/include/
vpath %.h /usr/lib/gcc/x86_64-linux-gnu/9/include/
vpath %.h /usr/include/
vpath %.h /usr/local/include/
vpath %.h /usr/include/gsl/
vpath %.h /usr/local/include/gsl/
vpath %.h /usr/local/lib/
vpath %.h /usr/lib/x86_64-linux-gnu/
# To compile and link the source files - For microstructural evolution

muse_evolver: $(objects) $(headers)
	gcc -g -o muse_evolver.out $(objects) $(LINKOPS)

MUSE_EVOLVER.o: $(sources) $(headers)
	gcc -o $@ -c ./src/muse_evolver.c $(COMPOPS)
	
generate_elast_tensors.o: generate_elast_tensors.c 
	gcc -c $^ $(COMPOPS)

calculate_sigma_T.o: calculate_sigma_T.c
	gcc -c $(COMPOPS) $<

nrutil.o: nrutil.c
	gcc -c $^ $(COMPOPS)

calculate_omega.o: calculate_omega.c nrutil.c nrutil.h \
math.h float.h
	gcc -c $(COMPOPS) $<

evolve_diffusional.o: evolve_diffusional.c stdio.h stdlib.h \
math.h  float.h complex.h fftw3.h functions.h nrutil.h \
nrutil.c add_noise.c
	gcc -c $(COMPOPS) $<





evolve_inhom_constrained.o: evolve_inhom_constrained.c stdio.h \
stdlib.h math.h  float.h complex.h fftw3.h functions.h nrutil.h \
nrutil.c calculate_Del_sigma_T.c calculate_uzero.c refine_u_chen.c \
calculate_alpha.c calculate_beta.c calculate_alpha_prime.c \
calculate_beta_prime.c add_noise.c
	gcc -c $(COMPOPS) $<

calculate_Del_sigma_T.o: calculate_Del_sigma_T.c
	gcc -c $(COMPOPS) $<

calculate_uzero.o: calculate_uzero.c complex.h \
fftw3.h functions.h nrutil.h nrutil.c
	gcc -c $(COMPOPS) $<

refine_u.o: refine_u.c stdio.h stdlib.h math.h complex.h \
fftw3.h functions.h nrutil.h calculate_S_exp.c \
calculate_sig_zero_exp.c nrutil.c calculate_hom_strain.c \
calculate_omega.c
	gcc -c $(COMPOPS) $<

refine_u_chen.o: refine_u_chen.c stdio.h stdlib.h math.h complex.h \
fftw3.h functions.h nrutil.h calculate_S_exp.c \
calculate_sig_zero_exp.c
	gcc -g -c $(COMPOPS) $<

calculate_hom_strain.o: calculate_hom_strain.c stdio.h stdlib.h \
complex.h fftw3.h functions.h nrutil.h nrutil.c
	gcc -c $(COMPOPS) $<

calculate_alpha.o: calculate_alpha.c math.h
	gcc -c $(COMPOPS) $<

calculate_beta.o: calculate_beta.c math.h
	gcc -c $(COMPOPS) $<

calculate_alpha_prime.o: calculate_alpha_prime.c math.h
	gcc -c $(COMPOPS) $<

calculate_beta_prime.o: calculate_beta_prime.c math.h
	gcc -c $(COMPOPS) $<

calculate_S_exp.o: calculate_S_exp.c complex.h fftw3.h \
functions.h
	gcc -c $(COMPOPS) $<

calculate_sig_zero_exp.o: calculate_sig_zero_exp.c functions.h \
complex.h fftw3.h
	gcc -c $(COMPOPS) $<

write_ppfiles.o: write_ppfiles.c stdio.h stdlib.h complex.h fftw3.h
	gcc -c $(COMPOPS) $<

add_noise.o: add_noise.c stdlib.h gsl_rng.h gsl_math.h
	gcc -c $(COMPOPS) $<

calculate_S.o: calculate_S.c complex.h fftw3.h \
functions.h
	gcc -c $(COMPOPS) $<

calculate_aleph.o: calculate_aleph.c stdio.h stdlib.h functions.h \
math.h complex.h fftw3.h
	gcc -c $(COMPOPS) $<

# To remove the intermediate object files

clean:
	rm -rf *.o
	rm -rf *.bbg 
	rm -rf *.bb 
	rm -rf *.da 
	rm -rf *.gcov 
	rm -rf ../output/postprocessor/ps/src/psmakerc.c
	rm -rf ../output/postprocessor/ps/src/psmakere.c
	rm -rf ../output/postprocessor/data/src/datamakerc.c
	rm -rf ../output/postprocessor/data/src/datamakere.c

# To remove all object files

CLEAN:
	rm -rf *.out

# End of the makefile
