# Makefile for fd2d

# Set compiler
FC := gfortran
LD := gfortran

# Set compiler flags based on make options
CFLAGS := -Wall -O3 -g -m64
FFLAGS := -fopenmp -fbounds-check -ffpe-trap=zero -O3 -Wall -g -J./obj -I./obj -m64
LDFLAGS := -fopenmp -m64 

# Identify directories
SUB_DIRS := para base
SRC_DIR  := $(addprefix source/,$(SUB_DIRS))

# identify object files
#common_parameter needs to be first as mod file is depended upon by nearly everything.
OBJ_FILES := obj/kind_parameters.o obj/common_parameter.o obj/common_vars.o
OBJ_FILES += obj/derivatives.o 
OBJ_FILES += obj/output.o obj/setup.o
OBJ_FILES += obj/weights.o obj/rhs.o
OBJ_FILES += obj/step.o
OBJ_FILES += $(foreach sdir,$(SRC_DIR),$(patsubst $(sdir)/%.F90,obj/%.o,$(wildcard $(sdir)/*.F90)))
OBJ_FILES += $(foreach sdir,$(SRC_DIR),$(patsubst $(sdir)/%.cpp,obj/%.o,$(wildcard $(sdir)/*.cpp)))

HDEPS := $(OBJ_FILES:.o=.d)

vpath %.F90 $(SRC_DIR)
vpath %.cpp $(SRC_DIR)

#-------
default: fd2d
fd2d: $(OBJ_FILES)
	$(LD) -o $@ $^ $(LDFLAGS)

obj/%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $<

#-include $(HDEPS)


clean:
	rm -vf ./obj/*.o
	rm -vf ./obj/*.mod
	rm -vf ./fd2d
	rm -rfv fort.*	
	rm -vf ./data_out/layer*
	rm -vf ./data_out/statistics/*.out
	rm -vf ./paraview_files/LAYER*

