CXX := mpicc
CXXFLAGS := -acc -Minline -Minfo=accel -ta=tesla,cuda8.0,cc60,pinned
LDFLAGS :=
LDLIBS :=

EXE_DIR :=  bin/
EXECUTABLE := $(EXE_DIR)mandelbrot_mpi
SRC_FILES := $(wildcard src/*.c)

OBJ_DIR :=  obj/
SRC_DIR :=  src/
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC_FILES:.c=.o)))
SRC_DIR := $(dir $(SRC_FILES) $(PROB_FILES))
# Set vpath
VPATH := $(SRC_DIR)

# Makefile targets

.PHONY : all dirs clean

all: dirs $(EXECUTABLE)

dirs : $(EXE_DIR) $(OBJ_DIR)

$(EXE_DIR):
	mkdir -p $(EXE_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Link objects into executable(s)

mandelbrot: $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) mandelbrot_mpi.o -o $@ $(LDFLAGS) $(LDLIBS)

# Create objects from source files
$(OBJ_DIR)%.o : %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

#mandelbrot_mpi.o: mandelbrot_mpi.c
#	$(CXX) -c $<

# For running on summitdev
launch_session:
	bsub -n 20 -x -P CSC261 -W 120 -Is $$SHELL

run:
	mpirun -n 1 ./mandelbrot_mpi 1000 1000

run_profile:
	mpirun -n 1 nvprof --print-gpu-summary ./mandelbrot_mpi 1000 1000

# Cleanup

clean:
	rm -rf $(OBJ_DIR)*
	rm -rf $(EXECUTABLE)
#rm *o mandelbrot
