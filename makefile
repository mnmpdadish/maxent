# source files, objects, headers and libraries
SOURCES      = src/main.cpp   src/loop_run.cpp  src/utilities.cpp  src/maxEnt_data.cpp  src/maxEnt_methods.cpp  src/minimize.cpp  src/preproc.cpp src/gnuplot_pipe.cpp

OBJECTS      = $(SOURCES:.cpp=.o)

EXENAME= ./maxEnt

OPT = -std=c++11 -O0 -DARMA_DONT_USE_WRAPPER

INCLUDE_PATH = -I/home/maxime/install/armadillo-5.600.2/include
#-I/usr/include

MKL          = /opt/intel/mkl/lib/intel64
LIBS         = ${MKL}/libmkl_intel_lp64.a -Wl,--start-group $(MKL)/libmkl_blas95_lp64.a $(MKL)/libmkl_lapack95_lp64.a $(MKL)/libmkl_sequential.a ${MKL}/libmkl_core.a -Wl,--end-group  -lgomp -lpthread -lm -ldl 
#LIBS         = ${F95ROOT}/lib/intel64/libmkl_blas95_ilp64.a ${F95ROOT}/lib/intel64/libmkl_lapack95_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
#-larmadillo

LIB_PATH     = -L/usr/lib

# compiler
COMP         = g++
COMP_FLAGS   = $(INCLUDE_PATH) $(OPT)

# linker
LINK         = g++
LINK_FLAGS   = $(LIB_PATH) 

all: $(EXENAME)

# remove objects and executable
clean: 
	rm -rf $(OBJECTS) $(EXENAME)

# executable
$(EXENAME) : $(OBJECTS)
	$(LINK) $(OBJECTS) $(LIBS) $(LINK_FLAGS) -o $(EXENAME)

# rule for source to object conversion
.cpp.o:
	$(COMP) $(COMP_FLAGS) -c $< -o $@


# tool that generates dependencies
#dependencies:
#	makedepend $(COMP_FLAGS) $(SOURCES)
