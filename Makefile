#### specify the target OS
OS = _WINDOWS_
OS = _LINUX_

### specify a block of compilers/flags, etc from Makefile.in (edit this file)
MODE = zero
include Makefile.in

### decide on debugging level(s)
COPTS += -DNO_DEBUG_
COPTS += -DNO_DEBUG2_

### -rpath arguments for finding .so objects in pre-specified locations
### WARNING: note the concatanation of -rpath's below needs a leading comma!
#MY_DIR = .
#RPATH =,-rpath=$(MY_DIR)


############################### Various ##############################
COPTS += -D$(OS)
COPTS += -I $(INMPI_PATH)
LOPTS += -L $(INMPI_PATH) -lINmpi

############################### Targets ##############################
all: libs
	$(MPICC) -c $(COPTS) code.c
	$(MPICC) $(COPTS) \
       main.c -Wl,-rpath=.,-rpath=$(INMPI_PATH)$(RPATH) \
              code.o \
              -L . -lINtest -lm $(LOPTS)

fortran: libs
	$(MPICC) -c -fPIC fcode_overlap.c -I .
	$(MPIF90) $(FOPTS) fcode.f -Wl,-rpath=. fcode_overlap.o -L. -lINtest

libs:
	$(MPICXX) -c $(COPTS) -D_USE_COMM_DUP incg_facematch.cpp
	$(MPICXX) -shared -Wl,-soname,libINtest.so,-rpath=$(INMPI_PATH)$(RPATH) \
             -o libINtest.so \
             incg_facematch.o \
             $(LOPTS)


clean:
	rm -f  *.o a.out libIN*.so

