CXX := mpicc
CXXFLAGS := -acc -Minline -Minfo=accel -ta=tesla,cuda8.0,cc60,pinned
LDFLAGS :=
LDLIBS :=

EXE_DIR :=  bin/
SRC_FILES := $(wildcard src/*.c)
EXECUTABLES := $(notdir $(SRC_FILES:.c=))

OBJ_DIR :=  obj/
SRC_DIR :=  src/
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC_FILES:.c=.o)))
SRC_DIR := $(dir $(SRC_FILES) $(PROB_FILES))
# Set vpath
VPATH := $(SRC_DIR)

# Makefile targets

.PHONY : all dirs clean

all: dirs $(EXECUTABLES)

dirs : $(EXE_DIR) $(OBJ_DIR)

$(EXE_DIR):
	mkdir -p $(EXE_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Link objects into executable(s)

$(EXECUTABLES): % : $(OBJ_DIR)%.o
# echo $(subst bin/,obj/, $@).o
	$(CXX) $(CXXFLAGS) $< -o $(EXE_DIR)$@ $(LDFLAGS) $(LDLIBS)

# Create objects from source files
$(OBJ_DIR)%.o : %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

#mandelbrot_mpi.o: mandelbrot_mpi.c
#	$(CXX) -c $<

# For running on summitdev
launch_session:
	bsub -n 20 -x -P CSC261 -W 120 -Is $$SHELL

run:
	mpirun -n 1 ./bin/mandelbrot_mpi_color 1000 1000

run_profile:
	mpirun -n 1 nvprof --print-gpu-summary ./bin/mandelbrot_mpi_color 1000 1000

# Cleanup

clean:
	rm -rf $(OBJ_DIR)*
	rm -rf $(addprefix $(EXE_DIR),$(notdir $(SRC_FILES:.c=)))
#rm *o mandelbrot


# Debug makefile
print-%  : ; @echo $* = $($*)
