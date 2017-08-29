# source files, objects, headers and libraries
SOURCES      = src/main.cpp   src/loop_run.cpp  src/utilities.cpp  src/maxEnt_data.cpp  src/maxEnt_methods.cpp  src/minimize.cpp  src/preproc.cpp

OBJECTS      = $(SOURCES:.cpp=.o)

EXENAME= /home/maxime/bin/maxEnt

OPT = -std=c++11 -O0 -DARMA_DONT_USE_WRAPPER

INCLUDE_PATH = -Isrc/armadillo-5.600.2/include
#-I/usr/include

LIBS         = -llapack -lblas 
#-larmadillo

LIB_PATH     = -L/usr/lib

# compiler
COMP         = g++-4.9
COMP_FLAGS   = $(INCLUDE_PATH) $(OPT)

# linker
LINK         = g++-4.9
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
