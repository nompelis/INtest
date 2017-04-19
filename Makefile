
MODE = zero

include Makefile.in

#### specify the target OS
OS = _WINDOWS_
OS = _LINUX_
COPTS += -D$(OS)

############################### Various ##############################
### -rpath arguments for finding .so objects in pre-specified locations
MY_DIR = .
RPATH = -rpath=$(MY_DIR)
### decide on debugging level(s)
COPTS += -D  _DEBUG_
COPTS += -DNO_DEBUG2_

############################### Target ##############################
all: libs
	$(MPICC) $(COPTS) \
       main.c -Wl,-rpath=.,$(RPATH)  -L . -lINtest -lm


libs:
	$(MPICXX) -c $(COPTS) -D_USE_COMM_DUP intest.cpp
	$(MPICXX) -shared -Wl,-soname,libINtest.so   -o libINtest.so intest.o


clean:
	rm -f  *.o a.out libIN*.so

