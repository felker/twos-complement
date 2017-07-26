CXX := mpicc
CXXFLAGS := -acc -Minline -Minfo=accel -ta=tesla,cuda8.0,cc60,pinned

all: mandelbrot

mandelbrot: mandelbrot_mpi.o
	$(CXX) $(CXXFLAGS) -o $@

mandelbrot_mpi.o: mandelbrot_mpi.c
	$(CXX) -c $<

launch_session:
	bsub -n 20 -x -P CSC261 -W 120 -Is $$SHELL

run:
	mpirun -n 1 ./mandelbrot_mpi 1000 1000

run_profile:
	mpirun -n 1 nvprof --print-gpu-summary ./mandelbrot_mpi 1000 1000

clean:
	rm *o mandelbrot
